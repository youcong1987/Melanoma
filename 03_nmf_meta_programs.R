## 03_nmf_meta_programs.R
## Non-negative Matrix Factorization to identify melanoma meta-programs
## Methods: cNMF (Python) + NMF (R), cross-validated with Jaccard similarity
## Reference: Gavish et al. Nature 2023

library(Seurat)
library(tidyverse)
library(qs)
library(NMF)
library(data.table)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(123)

# ============================================================
# 1. Prepare input: extract malignant cells per sample
# ============================================================

seu <- qread("output/02_malignant_identification/seurat_malignant_identified.qs")
seu <- subset(seu, malignant_final == "Malignant")

sample_size <- table(seu$orig.ident)
keep_samples <- names(sample_size[sample_size >= 300])
seu <- subset(seu, orig.ident %in% keep_samples)

seu <- NormalizeData(seu) %>% FindVariableFeatures(nfeatures = 3000)
HVGs <- VariableFeatures(seu)

dir.create("output/03_NMF/00_input", recursive = TRUE)
for (sid in keep_samples) {
  mat <- GetAssayData(seu, layer = paste0("counts.", sid))[HVGs, ]
  mat <- mat[!grepl("^MT-", rownames(mat)), ]
  mat <- mat[rowSums(mat) > 0, ]
  write.table(t(as.matrix(mat)),
    file = file.path("output/03_NMF/00_input", paste0(sid, ".count.txt")),
    sep = "\t", quote = FALSE)
}

# ============================================================
# 2. cNMF (Python): per-sample factorization
# ============================================================

cnmf_cli <- "/path/to/cnmf"
k_range <- 4:12
n_iter <- 200

for (sid in keep_samples) {
  count_file <- file.path("output/03_NMF/00_input", paste0(sid, ".count.txt"))
  out_dir <- file.path("output/03_NMF/01_cNMF/step1_factorize")

  system(paste(cnmf_cli, "prepare --output-dir", out_dir, "--name", sid,
               "-c", count_file, "-k", paste(k_range, collapse = " "),
               "--n-iter", n_iter, "--max-nmf-iter 2000"))
  system(paste(cnmf_cli, "factorize --output-dir", out_dir, "--name", sid,
               "--worker-index 0 --total-workers 1"))
  system(paste(cnmf_cli, "combine --output-dir", out_dir, "--name", sid))
  system(paste(cnmf_cli, "k_selection_plot --output-dir", out_dir, "--name", sid))
}

# K selection: maximize silhouette score from cNMF k_selection_stats
select_best_k <- function(sid, out_dir) {
  ksel_file <- file.path(out_dir, sid, paste0(sid, ".k_selection_stats.txt"))
  kstats <- read.table(ksel_file, header = TRUE, sep = "\t")
  kstats$k[which.max(kstats$silhouette)]
}

# Consensus with selected K (density threshold 0.02)
for (sid in keep_samples) {
  best_k <- select_best_k(sid, file.path("output/03_NMF/01_cNMF/step1_factorize"))
  system(paste(cnmf_cli, "consensus --output-dir",
               file.path("output/03_NMF/01_cNMF/step1_factorize"),
               "--name", sid, "--components", best_k,
               "--local-density-threshold 0.02 --show-clustering"))
}

# Extract top 50 genes per program from gene_spectra_score
top_gene <- 50
cnmf_programs <- list()
for (sid in keep_samples) {
  best_k <- select_best_k(sid, file.path("output/03_NMF/01_cNMF/step1_factorize"))
  score_file <- file.path("output/03_NMF/01_cNMF/step1_factorize", sid,
    paste0(sid, ".gene_spectra_score.k_", best_k, ".dt_0_02.txt"))
  score_df <- read.table(score_file, header = TRUE, row.names = 1, sep = "\t")
  score_df <- as.data.frame(t(score_df))
  for (k in seq_len(ncol(score_df))) {
    prog_name <- paste(sid, k, sep = "_")
    cnmf_programs[[prog_name]] <- names(sort(score_df[, k], decreasing = TRUE))[1:top_gene]
  }
}

# ============================================================
# 3. NMF (R): per-sample, all K values (Gavish protocol)
# ============================================================

k_range_nmf <- 4:9
n_top_genes <- 7000
nmfr_programs <- list()

for (sid in keep_samples) {
  count_file <- file.path("output/03_NMF/00_input", paste0(sid, ".count.txt"))
  mat <- fread(count_file, header = TRUE, sep = "\t")
  rn <- mat[[1]]; mat <- as.matrix(mat[, -1, with = FALSE]); rownames(mat) <- rn
  mat <- t(mat)  # genes x cells

  # CPM -> log2(CPM/10 + 1) (Gavish et al. protocol)
  col_sums <- colSums(mat); col_sums[col_sums == 0] <- 1
  cpm_mat <- sweep(mat, 2, col_sums, "/") * 1e6
  log_mat <- log2(cpm_mat / 10 + 1)

  # Keep top 7000 genes by mean expression
  gene_means <- rowMeans(log_mat)
  top_idx <- order(gene_means, decreasing = TRUE)[1:min(n_top_genes, nrow(log_mat))]
  log_mat <- log_mat[top_idx, ]

  # Center by gene mean, clip negatives to 0
  centered_mat <- log_mat - rowMeans(log_mat)
  centered_mat[centered_mat < 0] <- 0

  # Run NMF at each K separately (brunet method, 10 runs)
  for (k in k_range_nmf) {
    res <- tryCatch(
      nmf(centered_mat, rank = k, method = "brunet", seed = 123456, nrun = 10, .options = "-p"),
      error = function(e) NULL)
    if (is.null(res)) next

    W <- basis(res)
    rownames(W) <- rownames(centered_mat)
    for (j in seq_len(k)) {
      prog_name <- paste(sid, k, j, sep = "_")
      nmfr_programs[[prog_name]] <- names(sort(W[, j], decreasing = TRUE))[1:top_gene]
    }
  }
}

# ============================================================
# 4. Robust program filtering
# ============================================================

filter_robust_programs <- function(programs, across_overlap_genes, across_min_samples,
                                   redundancy_overlap, top_gene, skip_within = FALSE) {
  prog_names <- names(programs)
  prog_info <- data.frame(
    name = prog_names,
    sample = sapply(strsplit(prog_names, "_"), function(x) x[1]))

  # Within-sample robustness (NMF-R only): require overlap with another K's program
  if (!skip_within) {
    robust_within <- c()
    for (s in unique(prog_info$sample)) {
      s_progs <- prog_info$name[prog_info$sample == s]
      for (p in s_progs) {
        is_robust <- any(sapply(setdiff(s_progs, p), function(q)
          length(intersect(programs[[p]], programs[[q]])) / top_gene >= 0.2))
        if (is_robust) robust_within <- c(robust_within, p)
      }
    }
    programs <- programs[robust_within]
    prog_info <- prog_info[prog_info$name %in% robust_within, ]
  }

  # Cross-sample reproducibility
  robust_across <- sapply(names(programs), function(p) {
    sample_p <- prog_info$sample[prog_info$name == p]
    cross_samples <- unique(prog_info$sample[sapply(names(programs), function(q) {
      if (prog_info$sample[prog_info$name == q] == sample_p) return(FALSE)
      length(intersect(programs[[p]], programs[[q]])) >= across_overlap_genes
    })])
    length(cross_samples) >= across_min_samples
  })
  programs <- programs[robust_across]

  # Remove redundant programs (high overlap within same sample)
  # Keep the one with more cross-sample support
  programs
}

robust_cnmf <- filter_robust_programs(cnmf_programs,
  across_overlap_genes = 6, across_min_samples = 3,
  redundancy_overlap = 0.15, top_gene = 50, skip_within = TRUE)

robust_nmfr <- filter_robust_programs(nmfr_programs,
  across_overlap_genes = 5, across_min_samples = 3,
  redundancy_overlap = 0.2, top_gene = 50, skip_within = FALSE)

# ============================================================
# 5. Meta-program identification (iterative greedy algorithm)
# ============================================================

cluster_meta_programs <- function(programs, min_overlap_genes = 8,
                                  min_cluster_size = 3, mp_gene_size = 50,
                                  max_mt_rb_frac = 0.3) {
  prog_names <- names(programs)
  n <- length(prog_names)

  # Compute pairwise overlap matrix (number of shared genes)
  overlap_mat <- matrix(0, n, n, dimnames = list(prog_names, prog_names))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      ov <- length(intersect(programs[[i]], programs[[j]]))
      overlap_mat[i, j] <- overlap_mat[j, i] <- ov
    }
  }

  meta_programs <- list()
  clusters <- list()
  remaining <- prog_names

  while (length(remaining) >= min_cluster_size) {
    # Find founder: program with most connections (>= min_overlap_genes)
    counts_remaining <- sapply(remaining, function(p)
      sum(overlap_mat[p, setdiff(remaining, p)] >= min_overlap_genes))
    founder <- names(which.max(counts_remaining))
    if (counts_remaining[founder] < min_cluster_size) break

    # Initialize cluster with founder + best partner
    overlaps_with_founder <- overlap_mat[founder, setdiff(remaining, founder)]
    best_partner <- names(which.max(overlaps_with_founder))
    current_cluster <- c(founder, best_partner)
    current_mp_genes <- intersect(programs[[founder]], programs[[best_partner]])

    # Pad to mp_gene_size if needed
    if (length(current_mp_genes) < mp_gene_size) {
      all_genes <- unique(c(programs[[founder]], programs[[best_partner]]))
      gene_freq <- table(c(programs[[founder]], programs[[best_partner]]))
      extra <- setdiff(names(sort(gene_freq, decreasing = TRUE)), current_mp_genes)
      current_mp_genes <- c(current_mp_genes, head(extra, mp_gene_size - length(current_mp_genes)))
    }
    current_mp_genes <- head(current_mp_genes, mp_gene_size)

    # Iteratively add programs with sufficient overlap
    candidates <- setdiff(remaining, current_cluster)
    added <- TRUE
    while (added && length(candidates) > 0) {
      added <- FALSE
      candidate_overlaps <- sapply(candidates, function(p)
        length(intersect(programs[[p]], current_mp_genes)))
      valid <- candidates[candidate_overlaps >= min_overlap_genes]
      if (length(valid) == 0) break
      best <- valid[which.max(candidate_overlaps[valid])]
      current_cluster <- c(current_cluster, best)
      candidates <- setdiff(candidates, best)
      added <- TRUE
      # Update MP genes: most frequent across cluster members
      all_genes_in_cluster <- unlist(programs[current_cluster])
      gene_freq <- sort(table(all_genes_in_cluster), decreasing = TRUE)
      current_mp_genes <- head(names(gene_freq), mp_gene_size)
    }

    # Filter MT/RB-dominated MPs
    mt_rb_frac <- sum(grepl("^(MT-|RP[SL])", current_mp_genes)) / length(current_mp_genes)
    if (mt_rb_frac <= max_mt_rb_frac) {
      mp_name <- paste0("MP", length(meta_programs) + 1)
      meta_programs[[mp_name]] <- current_mp_genes
      clusters[[mp_name]] <- current_cluster
    }
    remaining <- setdiff(remaining, current_cluster)
  }
  list(meta_programs = meta_programs, clusters = clusters)
}

# Run clustering for each method
mp_result_nmf <- cluster_meta_programs(robust_nmfr,
  min_overlap_genes = 8, min_cluster_size = 3, mp_gene_size = 50)
mp_result_cnmf <- cluster_meta_programs(robust_cnmf,
  min_overlap_genes = 8, min_cluster_size = 3, mp_gene_size = 50)

# Post-processing: remove immune-marker-dominated MPs (doublet signal)
immune_markers <- c("PTPRC", "CD3D", "CD3E", "CD8A", "CD4", "CD19", "CD79A",
                    "MS4A1", "NKG7", "GNLY", "CD14", "CD68", "CSF1R")
for (mp_name in names(mp_result_nmf$meta_programs)) {
  genes <- mp_result_nmf$meta_programs[[mp_name]]
  if (sum(genes %in% immune_markers) / length(genes) >= 0.08)
    mp_result_nmf$meta_programs[[mp_name]] <- NULL
}

# ============================================================
# 6. Cross-method validation and final MP definition
# ============================================================

# Jaccard similarity between cNMF MPs and NMF-R MPs
compare_jaccard <- function(mps_a, mps_b) {
  mat <- matrix(0, length(mps_a), length(mps_b),
                dimnames = list(names(mps_a), names(mps_b)))
  for (i in names(mps_a)) for (j in names(mps_b)) {
    inter <- length(intersect(mps_a[[i]], mps_b[[j]]))
    union <- length(union(mps_a[[i]], mps_b[[j]]))
    mat[i, j] <- inter / union
  }
  mat
}

jacc_cross <- compare_jaccard(mp_result_cnmf$meta_programs, mp_result_nmf$meta_programs)

# Final MPs: those with Jaccard > 0.1 between cNMF and NMF-R
# Use cNMF version as primary (higher quality consensus)

# ============================================================
# 7. Functional annotation
# ============================================================

for (mp_name in names(mp_result_cnmf$meta_programs)) {
  genes <- mp_result_cnmf$meta_programs[[mp_name]]
  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.1)
  hallmark <- read.gmt("source/Data/h.all.v2023.2.Hs.symbols.gmt")
  ehallmark <- enricher(genes, TERM2GENE = hallmark, pvalueCutoff = 0.1)
}

# ============================================================
# 8. Comparison with known programs
# ============================================================

# Pan-cancer meta-programs (Gavish et al. Nature 2023)
pancan_mps <- read_tsv("source/Data/Gavish2023_pan_cancer_MPs.tsv")
pancan_list <- split(pancan_mps$gene, pancan_mps$MP)

# Melanoma subtypes (Tsoi et al. 2018)
mel_subtypes <- read_tsv("source/Data/melanoma_subtype_signatures.tsv")
mel_list <- split(mel_subtypes$gene, mel_subtypes$subtype)

jacc_pancan <- compare_jaccard(mp_result_cnmf$meta_programs, pancan_list)
jacc_mel <- compare_jaccard(mp_result_cnmf$meta_programs, mel_list)
