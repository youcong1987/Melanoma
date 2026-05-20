## 08_drug_screening.R
## Multi-strategy drug screening for MP1/MP7 core genes
## Reference: Yang et al. Cell Discovery 2022

library(tidyverse)
library(data.table)
library(GSVA)
library(limma)
library(oncoPredict)
library(httr)
library(jsonlite)
library(ggpubr)
library(patchwork)
library(beyondcell)
library(qs)

set.seed(5201314)

# ============================================================
# Part 0: Data loading and ssGSEA scoring
# ============================================================

load("output/04_Bulk_preprocessing/01_TCGA_SKCM/counts_tpm_metadata_Tumor_survival.Rdata")
tcga_expr <- Tumor_tpm_surv
tcga_meta <- metadata_surv

mp1_genes <- readLines("output/07_core_gene_prognosis/MP1_final_core_genes.txt")
mp7_genes <- readLines("output/07_core_gene_prognosis/MP7_final_core_genes.txt")

# Auto log2 transform
log_expr <- if (max(tcga_expr, na.rm = TRUE) > 30) log2(as.matrix(tcga_expr) + 1) else as.matrix(tcga_expr)

# ssGSEA scoring (not GSVA)
gs_filt <- list(MP1 = intersect(mp1_genes, rownames(log_expr)),
                MP7 = intersect(mp7_genes, rownames(log_expr)))
param <- ssgseaParam(exprData = log_expr, geneSets = gs_filt)
gsva_scores <- as.data.frame(t(gsva(param, verbose = FALSE)))
gsva_scores$MP1_group <- ifelse(gsva_scores$MP1 > median(gsva_scores$MP1), "MP1-high", "MP1-low")
gsva_scores$MP7_group <- ifelse(gsva_scores$MP7 > median(gsva_scores$MP7), "MP7-high", "MP7-low")

# ============================================================
# Part 1: Drug sensitivity prediction (GDSC1 + GDSC2 + CTRP2)
# ============================================================

# Training data from oncoPredict
# GDSC1/GDSC2: RMA log-transformed; CTRP2: raw TPM scale
db_files <- list(
  GDSC1 = list(expr = readRDS("source/oncoPredict/GDSC1_Expr.rds"),
               res = readRDS("source/oncoPredict/GDSC1_Res.rds")),
  GDSC2 = list(expr = readRDS("source/oncoPredict/GDSC2_Expr.rds"),
               res = readRDS("source/oncoPredict/GDSC2_Res.rds")),
  CTRP2 = list(expr = readRDS("source/oncoPredict/CTRP2_Expr.rds"),
               res = readRDS("source/oncoPredict/CTRP2_Res.rds"))
)

# TCGA data is already log2 scale; for CTRP2 need raw TPM
test_expr_log <- as.matrix(log_expr)
test_expr_tpm <- 2^test_expr_log - 1
test_expr_tpm[test_expr_tpm < 0] <- 0

all_ic50 <- list()
for (db_name in names(db_files)) {
  cur_test <- if (db_name == "CTRP2") test_expr_tpm else test_expr_log
  calcPhenotype(
    trainingExprData = db_files[[db_name]]$expr,
    trainingPtype = db_files[[db_name]]$res,
    testExprData = cur_test,
    batchCorrect = "eb",
    powerTransformPhenotype = TRUE,
    removeLowVaryingGenes = 0.2,
    minNumSamples = 10,
    selection = 1
  )
  all_ic50[[db_name]] <- read.csv("DrugPredictions.csv", row.names = 1, check.names = FALSE)
}

# Differential analysis: Wilcoxon + Spearman correlation
run_drug_diff <- function(ic50, scores, group_col, group_hi, group_lo, score_col) {
  common <- intersect(rownames(ic50), rownames(scores))
  ic50 <- ic50[common, , drop = FALSE]
  sc <- scores[common, ]
  hi_idx <- which(sc[[group_col]] == group_hi)
  lo_idx <- which(sc[[group_col]] == group_lo)

  diff_res <- data.frame(drug = colnames(ic50), stringsAsFactors = FALSE)
  diff_res$mean_hi <- colMeans(ic50[hi_idx, ], na.rm = TRUE)
  diff_res$mean_lo <- colMeans(ic50[lo_idx, ], na.rm = TRUE)
  diff_res$delta <- diff_res$mean_hi - diff_res$mean_lo
  diff_res$p_value <- sapply(seq_len(ncol(ic50)), function(j) {
    tryCatch(wilcox.test(ic50[hi_idx, j], ic50[lo_idx, j])$p.value, error = function(e) NA)
  })
  diff_res$fdr <- p.adjust(diff_res$p_value, method = "BH")
  diff_res$cor_rho <- sapply(seq_len(ncol(ic50)), function(j) {
    cor(sc[[score_col]], ic50[, j], method = "spearman", use = "complete.obs")
  })
  diff_res$cor_p <- sapply(seq_len(ncol(ic50)), function(j) {
    cor.test(sc[[score_col]], ic50[, j], method = "spearman")$p.value
  })
  diff_res
}

# MP7: drugs with lower IC50 in MP7-high (target MP7-high tumors)
# MP1 sensitization: lower IC50 in MP1-high
# MP1 cold-to-hot: lower IC50 in MP1-low (activate cold tumors)
for (db_name in names(all_ic50)) {
  ic50 <- all_ic50[[db_name]]
  diff_mp7 <- run_drug_diff(ic50, gsva_scores, "MP7_group", "MP7-high", "MP7-low", "MP7")
  diff_mp1_sens <- run_drug_diff(ic50, gsva_scores, "MP1_group", "MP1-high", "MP1-low", "MP1")
  diff_mp1_c2h <- run_drug_diff(ic50, gsva_scores, "MP1_group", "MP1-low", "MP1-high", "MP1")
}

# ============================================================
# Part 2: CMap / L1000CDS2 reverse gene expression matching
# ============================================================

# L1000CDS2 API (searchMethod = "geneSet")
query_l1000cds2 <- function(up_genes, dn_genes, aggravate = FALSE) {
  url <- "https://maayanlab.cloud/L1000CDS2/query"
  payload <- list(
    data = list(upGenes = as.list(up_genes), dnGenes = as.list(dn_genes)),
    config = list(aggravate = aggravate, searchMethod = "geneSet",
                  share = FALSE, combination = FALSE, `db-version` = "latest")
  )
  resp <- POST(url, body = toJSON(payload, auto_unbox = TRUE),
               content_type_json(), timeout(120))
  if (status_code(resp) != 200) return(NULL)
  result <- content(resp, as = "parsed")
  if (is.null(result$topMeta)) return(NULL)
  do.call(rbind, lapply(result$topMeta, function(e) {
    data.frame(pert_name = e$pert_desc %||% NA, pert_id = e$pert_id %||% NA,
               cell_line = e$cell_id %||% NA, score = e$score %||% NA,
               stringsAsFactors = FALSE)
  }))
}

# MP7 reversal: up=MP7 genes, dn=MP1 genes (aggravate=FALSE -> find reversals)
cmap_mp7 <- query_l1000cds2(up_genes = mp7_genes, dn_genes = head(mp1_genes, 20),
                             aggravate = FALSE)

# MP1 cold-to-hot: up=MP7, dn=MP1 (find compounds that reverse cold state)
cmap_mp1_c2h <- query_l1000cds2(up_genes = mp7_genes, dn_genes = mp1_genes,
                                 aggravate = FALSE)

# MP1 sensitization: up=MP1, dn=MP7 (aggravate=TRUE -> enhance MP1-high state)
cmap_mp1_sens <- query_l1000cds2(up_genes = mp1_genes, dn_genes = mp7_genes,
                                  aggravate = TRUE)

# ============================================================
# Part 3: DGIdb + CTD + DSigDB drug-target networks
# ============================================================

# DGIdb GraphQL API
query_dgidb <- function(genes) {
  url <- "https://dgidb.org/api/graphql"
  gene_str <- paste0('"', genes, '"', collapse = ", ")
  gql <- sprintf('{
    genes(names: [%s]) {
      nodes {
        name
        interactions {
          drug { name conceptId }
          interactionScore
          interactionTypes { type directionality }
          publications { pmid }
        }
      }
    }
  }', gene_str)

  resp <- POST(url, body = toJSON(list(query = gql), auto_unbox = TRUE),
               content_type_json(), timeout(60))
  raw <- fromJSON(content(resp, as = "text", encoding = "UTF-8"))
  nodes <- raw$data$genes$nodes

  dgidb_df <- data.frame()
  for (i in seq_len(nrow(nodes))) {
    interactions <- nodes$interactions[[i]]
    if (is.null(interactions) || nrow(interactions) == 0) next
    for (j in seq_len(nrow(interactions))) {
      int_types <- interactions$interactionTypes[[j]]
      dgidb_df <- rbind(dgidb_df, data.frame(
        gene = nodes$name[i], drug = interactions$drug$name[j],
        interaction_type = if (!is.null(int_types) && nrow(int_types) > 0) int_types$type[1] else "unknown",
        score = interactions$interactionScore[j], stringsAsFactors = FALSE))
    }
  }
  dgidb_df
}

dgidb_mp1 <- query_dgidb(mp1_genes)
dgidb_mp7 <- query_dgidb(mp7_genes)

# CTD: batch query gene-chemical interactions
query_ctd <- function(genes) {
  base_url <- "http://ctdbase.org/tools/batchQuery.go"
  gene_input <- paste(genes, collapse = "\n")
  resp <- POST(base_url, body = list(
    inputType = "gene", inputTerms = gene_input,
    report = "genes_chemicals", format = "tsv", actionTypes = "ANY"
  ), encode = "form", timeout(120))
  read.table(text = content(resp, as = "text"), sep = "\t", header = TRUE, quote = "")
}

ctd_mp1 <- query_ctd(mp1_genes)
ctd_mp7 <- query_ctd(mp7_genes)

# DSigDB enrichment (via Enrichr API)
query_dsigdb <- function(genes) {
  add_url <- "https://maayanlab.cloud/Enrichr/addList"
  resp1 <- POST(add_url, body = list(list = paste(genes, collapse = "\n")), encode = "form")
  user_list_id <- fromJSON(content(resp1, as = "text"))$userListId
  enrich_url <- paste0("https://maayanlab.cloud/Enrichr/enrich?userListId=",
                       user_list_id, "&backgroundType=DSigDB")
  resp2 <- GET(enrich_url)
  fromJSON(content(resp2, as = "text"))$DSigDB
}

dsig_mp1 <- query_dsigdb(mp1_genes)
dsig_mp7 <- query_dsigdb(mp7_genes)

# ============================================================
# Part 4: Beyondcell single-cell drug sensitivity
# ============================================================

seu_mal <- qread("output/02_malignant_identification/seurat_malignant_identified.qs")
seu_mal <- subset(seu_mal, malignant_final %in% c("Malignant", "Malignant_like"))

# Add MP group labels
seu_mal <- AddModuleScore(seu_mal, features = list(mp1_genes), name = "MP1_sc")
seu_mal <- AddModuleScore(seu_mal, features = list(mp7_genes), name = "MP7_sc")
seu_mal$MP1_group <- ifelse(seu_mal$MP1_sc1 > median(seu_mal$MP1_sc1), "MP1-high", "MP1-low")
seu_mal$MP7_group <- ifelse(seu_mal$MP7_sc1 > median(seu_mal$MP7_sc1), "MP7-high", "MP7-low")

# Beyondcell scoring with PSc
data("PSc", package = "beyondcell")
bc <- bcScore(sc = seu_mal, gs = PSc, expr.thres = 0.1)
bc <- bcUMAP(bc, pc = 10, k.neighbors = 20)

# Differential drug sensitivity between groups
bc_mp7 <- bcRanks(bc, idents = "MP7_group", ident.1 = "MP7-high", ident.2 = "MP7-low")
bc_mp1 <- bcRanks(bc, idents = "MP1_group", ident.1 = "MP1-high", ident.2 = "MP1-low")

bc_mp7_sig <- bc_mp7 %>% filter(p_val_adj < 0.05) %>% arrange(avg_diff)

# ============================================================
# Part 5: Multi-method integration and prioritization
# ============================================================

# Collect candidate drugs from all methods
drug_sources <- list(
  oncoPredict = unique(c(
    diff_mp7$drug[diff_mp7$fdr < 0.05 & diff_mp7$delta < 0],
    diff_mp1_sens$drug[diff_mp1_sens$fdr < 0.05 & diff_mp1_sens$delta < 0])),
  L1000CDS2 = unique(c(
    cmap_mp7$pert_name, cmap_mp1_c2h$pert_name, cmap_mp1_sens$pert_name)),
  DGIdb = unique(c(dgidb_mp1$drug, dgidb_mp7$drug)),
  CTD = unique(c(ctd_mp1$ChemicalName, ctd_mp7$ChemicalName)),
  Beyondcell = unique(head(bc_mp7_sig$drug, 50))
)

# Count method support per drug
all_drugs <- unique(tolower(unlist(drug_sources)))
drug_support <- data.frame(drug = all_drugs, stringsAsFactors = FALSE)
for (method in names(drug_sources)) {
  drug_support[[method]] <- tolower(drug_support$drug) %in% tolower(drug_sources[[method]])
}
drug_support$n_methods <- rowSums(drug_support[, -1])
drug_support <- drug_support %>% arrange(desc(n_methods))

# Final prioritized list: drugs supported by >= 2 methods
final_drugs <- drug_support %>% filter(n_methods >= 2)
