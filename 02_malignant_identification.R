## 02_malignant_identification.R
## Malignant cell identification using CopyKAT + inferCNV + driver gene validation
## Input: Integrated Seurat object with cell type annotations
## Classification rules:
##   - Both CopyKAT & inferCNV = Malignant -> Malignant
##   - Only one method = Malignant + driver gene support -> Malignant
##   - Only one method = Malignant, no driver support -> Malignant_like
##   - Both = Non-malignant -> Non-malignant

library(Seurat)
library(qs)
library(tidyverse)
library(data.table)
library(copykat)
library(infercnv)
library(parallel)
library(Matrix)

set.seed(20260331)

# ============================================================
# 1. Data preparation
# ============================================================

seurat_obj <- qread("output/01_scRNAseq_preprocessing/integrated/seurat_cell_annotation_main.qs")

tumor_label <- "Melanoma cells"
ref_labels <- c("T/NK cells", "Endothelial cells")
sample_ids <- unique(seurat_obj$orig.ident)

# Per-sample: extract tumor + reference cells, downsample reference to 2000
# T/NK and Endothelial allocated proportionally; if Endothelial < 50, use T/NK only

# ============================================================
# 2. CopyKAT analysis (per sample, parallel)
# ============================================================

run_copykat_single <- function(input_file, result_dir, ref_labels) {
  seu <- qread(input_file)
  sample_meta <- seu@meta.data
  norm_cells <- rownames(sample_meta)[sample_meta$celltype_main %in% ref_labels]
  counts <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "counts"))
  sid <- sub("_input\\.qs$", "", basename(input_file))

  copykat_res <- copykat(
    rawmat = counts,
    id.type = "S",
    ngene.chr = 5,
    win.size = 25,
    KS.cut = 0.1,
    sam.name = sid,
    distance = "euclidean",
    genome = "hg20",
    n.cores = 4,
    norm.cell.names = norm_cells,
    output.seg = TRUE,
    plot.genes = FALSE
  )

  # Extract prediction: aneuploid = Malignant, diploid = Non-malignant
  pred <- copykat_res$prediction
  pred$copykat_malignant <- ifelse(pred$copykat.pred == "aneuploid", "Malignant",
                            ifelse(pred$copykat.pred == "diploid", "Non-malignant", "Uncertain"))
  pred
}

# Run in parallel (5 workers, 4 cores each)
input_files <- list.files("output/02_malignant_identification/01_copykat/01_single_sample_input",
                          pattern = "\\.qs$", full.names = TRUE)
copykat_results <- mclapply(input_files, run_copykat_single,
                            result_dir = "output/02_malignant_identification/01_copykat/02_single_sample_result",
                            ref_labels = ref_labels, mc.cores = 5)

# ============================================================
# 3. inferCNV analysis (per sample)
# ============================================================

gene_order_file <- "source/Data/hg38_gencode_v44_gene_pos.txt"

run_infercnv_single <- function(input_file, gene_order_file, out_base_dir) {
  seu <- qread(input_file)
  sid <- sub("_input\\.qs$", "", basename(input_file))
  out_dir <- file.path(out_base_dir, sid)

  counts <- as.matrix(GetAssayData(seu, assay = "RNA", layer = "counts"))
  annotations <- data.frame(cell = colnames(seu), group = seu$celltype_main)
  write.table(annotations, file = file.path(out_dir, "annotations.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file = file.path(out_dir, "annotations.txt"),
    gene_order_file = gene_order_file,
    ref_group_names = ref_labels
  )

  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = out_dir,
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = FALSE,
    num_threads = 4
  )

  # Score: mean absolute CNV deviation from baseline
  cnv_mat <- infercnv_obj@expr.data
  ref_cells <- colnames(cnv_mat)[annotations$group %in% ref_labels]
  obs_cells <- setdiff(colnames(cnv_mat), ref_cells)
  ref_baseline <- mean(abs(cnv_mat[, ref_cells] - 1))
  obs_scores <- colMeans(abs(cnv_mat[, obs_cells] - 1))
  threshold <- ref_baseline + 2 * sd(abs(cnv_mat[, ref_cells] - 1))

  data.frame(cell = obs_cells,
             infercnv_score = obs_scores,
             infercnv_malignant = ifelse(obs_scores > threshold, "Malignant", "Non-malignant"))
}

# ============================================================
# 4. Driver gene validation (Bailey et al. 2018, 204 genes)
# ============================================================

driver_genes <- fread("source/Data/Bailey_204_driver_genes.txt")$Gene

# For each sample: check CNV status of driver genes from CopyKAT CNV matrix
# Amplification: mean CNV > 0.2; Deletion: mean CNV < -0.2
# Sample has "driver_support" = "Yes" if >= 2 driver genes show significant CNV

check_driver_genes <- function(copykat_cnv_file, driver_genes, amp_thresh = 0.2, del_thresh = -0.2) {
  cnv_mat <- fread(copykat_cnv_file)
  gene_col <- colnames(cnv_mat)[1]
  cnv_mat <- cnv_mat %>% column_to_rownames(gene_col)

  available_drivers <- intersect(driver_genes, rownames(cnv_mat))
  if (length(available_drivers) < 5) return("Insufficient")

  driver_cnv <- cnv_mat[available_drivers, , drop = FALSE]
  mean_cnv <- rowMeans(as.matrix(driver_cnv))
  n_amp <- sum(mean_cnv > amp_thresh)
  n_del <- sum(mean_cnv < del_thresh)

  ifelse((n_amp + n_del) >= 2, "Yes", "No")
}

# ============================================================
# 5. Integrated classification
# ============================================================

# Merge CopyKAT + inferCNV predictions per cell
copykat_pred <- rbindlist(copykat_results)
infercnv_pred <- rbindlist(infercnv_results)

merged_pred <- merge(copykat_pred[, .(cell, copykat_malignant)],
                     infercnv_pred[, .(cell, infercnv_malignant)],
                     by = "cell", all = TRUE)
merged_pred[is.na(copykat_malignant), copykat_malignant := "Uncertain"]
merged_pred[is.na(infercnv_malignant), infercnv_malignant := "Uncertain"]

# Classification rules
merged_pred$malignant_final <- "Undetermined"

# Rule 1: Both methods = Malignant
both_mal <- merged_pred$copykat_malignant == "Malignant" & merged_pred$infercnv_malignant == "Malignant"
merged_pred$malignant_final[both_mal] <- "Malignant"

# Rule 2: Both methods = Non-malignant
both_nonmal <- merged_pred$copykat_malignant == "Non-malignant" & merged_pred$infercnv_malignant == "Non-malignant"
merged_pred$malignant_final[both_nonmal] <- "Non-malignant"

# Rule 3: Only one method = Malignant
one_mal <- xor(merged_pred$copykat_malignant == "Malignant", merged_pred$infercnv_malignant == "Malignant")
# 3a: with driver gene support -> Malignant
merged_pred$malignant_final[one_mal & merged_pred$driver_support == "Yes"] <- "Malignant"
# 3b: without driver support -> Malignant_like
merged_pred$malignant_final[one_mal & merged_pred$driver_support != "Yes" &
                             merged_pred$malignant_final == "Undetermined"] <- "Malignant_like"

# Rule 4: One Non-malignant + one Uncertain -> Non-malignant
one_nonmal_uncertain <- (merged_pred$copykat_malignant == "Non-malignant" & merged_pred$infercnv_malignant == "Uncertain") |
                        (merged_pred$copykat_malignant == "Uncertain" & merged_pred$infercnv_malignant == "Non-malignant")
merged_pred$malignant_final[one_nonmal_uncertain & merged_pred$malignant_final == "Undetermined"] <- "Non-malignant"

# Write back to Seurat object
seurat_obj$malignant_final <- as.character(seurat_obj$celltype_main)
matched_cells <- intersect(colnames(seurat_obj), merged_pred$cell)
matched_idx <- match(matched_cells, merged_pred$cell)
seurat_obj$malignant_final[matched_cells] <- merged_pred$malignant_final[matched_idx]

# Non-melanoma cells keep original annotation
non_melanoma <- seurat_obj$celltype_main != tumor_label
seurat_obj$malignant_final[non_melanoma] <- as.character(seurat_obj$celltype_main[non_melanoma])

qsave(seurat_obj, "output/02_malignant_identification/seurat_malignant_identified.qs")
