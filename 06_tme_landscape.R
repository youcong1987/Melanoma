## 06_tme_landscape.R
## MP1/MP7 and tumor microenvironment (TME) landscape analysis
## Explores why MP1 associates with ICB response and MP7 with resistance

library(Seurat)
library(tidyverse)
library(qs)
library(GSVA)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)

set.seed(20260331)

# ============================================================
# Part 1: Single-cell — MP scores vs TME cell composition
# ============================================================

seurat_obj <- qread("output/01_scRNAseq_preprocessing/integrated/seurat_cell_annotation_main.qs")
seurat_mal <- qread("output/02_malignant_identification/seurat_malignant_identified.qs")

mp_gene_list <- qread("output/03_NMF/meta_programs_gene_list.qs")

# Score malignant cells (including Malignant_like)
mal_cells <- subset(seurat_mal, malignant_final %in% c("Malignant", "Malignant_like"))
DefaultAssay(mal_cells) <- "RNA"
mal_cells <- JoinLayers(mal_cells)

mp1_genes_avail <- intersect(mp_gene_list[["MP1"]], rownames(mal_cells))
mp7_genes_avail <- intersect(mp_gene_list[["MP7"]], rownames(mal_cells))
mal_cells <- AddModuleScore(mal_cells, features = list(mp1_genes_avail), name = "MP1_score")
mal_cells <- AddModuleScore(mal_cells, features = list(mp7_genes_avail), name = "MP7_score")

# Per-sample mean MP scores
mp_per_sample <- mal_cells@meta.data %>%
  group_by(orig.ident) %>%
  summarise(MP1_mean = mean(MP1_score1, na.rm = TRUE),
            MP7_mean = mean(MP7_score1, na.rm = TRUE),
            n_malignant = n(), .groups = "drop")

# Per-sample cell type proportions (all cells)
cell_prop <- seurat_obj@meta.data %>%
  group_by(orig.ident, celltype_main) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n)) %>%
  select(orig.ident, celltype_main, prop) %>%
  pivot_wider(names_from = celltype_main, values_from = prop, values_fill = 0)

# Spearman correlation: MP scores vs cell type proportions
sc_cor_df <- mp_per_sample %>% inner_join(cell_prop, by = "orig.ident")
cell_types <- setdiff(colnames(sc_cor_df), c("orig.ident", "MP1_mean", "MP7_mean", "n_malignant"))

cor_results <- data.frame()
for (ct in cell_types) {
  for (mp in c("MP1_mean", "MP7_mean")) {
    ct_test <- cor.test(sc_cor_df[[mp]], sc_cor_df[[ct]], method = "spearman")
    cor_results <- rbind(cor_results, data.frame(
      MP = gsub("_mean", "", mp), CellType = ct,
      rho = ct_test$estimate, pval = ct_test$p.value))
  }
}
cor_results$padj <- p.adjust(cor_results$pval, method = "BH")

# ============================================================
# Part 2: Bulk — TCGA-SKCM immune infiltration
# ============================================================

load("output/04_Bulk_preprocessing/01_TCGA_SKCM/counts_tpm_metadata_Tumor_survival.Rdata")
tcga_tpm <- Tumor_tpm_surv
tcga_meta <- metadata_surv

# Auto log-transform (if max > 20, apply log2(x+1))
log_mat <- if (max(tcga_tpm, na.rm = TRUE) > 20) log2(as.matrix(tcga_tpm) + 1) else as.matrix(tcga_tpm)

# GSVA scoring of MP1/MP7
gs_mp <- list(
  MP1 = intersect(mp_gene_list[["MP1"]], rownames(log_mat)),
  MP7 = intersect(mp_gene_list[["MP7"]], rownames(log_mat)))
param <- gsvaParam(exprData = log_mat, geneSets = gs_mp, kcdf = "Gaussian")
mp_scores <- gsva(param, verbose = FALSE)
tcga_meta$MP1_score <- as.numeric(mp_scores["MP1", tcga_meta$sample])
tcga_meta$MP7_score <- as.numeric(mp_scores["MP7", tcga_meta$sample])

# 2.1 ESTIMATE scores
library(estimate)
tpm_for_est <- tcga_tpm[, tcga_meta$sample]
filterCommonGenes(input.f = "tcga_for_estimate.gct", output.f = "estimate_filtered.gct")
estimateScore("estimate_filtered.gct", "estimate_scores.gct")

# 2.2 MCPcounter (uses raw TPM, not log-transformed)
library(MCPcounter)
mcp_scores <- MCPcounter.estimate(as.matrix(tcga_tpm[, tcga_meta$sample]),
                                   featuresType = "HUGO_symbols")
mcp_df <- as.data.frame(t(mcp_scores))
mcp_df$sample <- rownames(mcp_df)

# 2.3 ssGSEA immune cell signatures (12 cell types)
immune_cell_markers <- list(
  CD8_Tcells = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "EOMES", "TBX21"),
  CD4_Tcells = c("CD4", "IL7R", "TCF7", "LEF1", "CCR7", "SELL"),
  Tregs = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"),
  NK_cells = c("NKG7", "KLRD1", "KLRB1", "KLRC1", "NCR1", "GNLY"),
  B_cells = c("CD19", "MS4A1", "CD79A", "CD79B", "BANK1", "PAX5"),
  Macrophages_M1 = c("CD68", "NOS2", "IL1B", "TNF", "CXCL10", "IDO1"),
  Macrophages_M2 = c("CD68", "CD163", "MRC1", "MSR1", "CD200R1", "IL10"),
  Dendritic_cells = c("ITGAX", "CD1C", "CLEC9A", "BATF3", "IRF8", "FLT3"),
  Mast_cells = c("TPSAB1", "TPSB2", "CPA3", "KIT", "HDC"),
  Neutrophils = c("FCGR3B", "CXCR2", "CSF3R", "S100A8", "S100A9"),
  Fibroblasts = c("DCN", "COL1A1", "COL3A1", "FAP", "PDGFRB", "ACTA2"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "KDR", "FLT1", "ENG")
)

gs_immune <- lapply(immune_cell_markers, function(g) intersect(g, rownames(log_mat)))
gs_immune <- gs_immune[sapply(gs_immune, length) >= 3]
param_immune <- gsvaParam(exprData = log_mat[, tcga_meta$sample], geneSets = gs_immune, kcdf = "Gaussian")
immune_scores <- gsva(param_immune, verbose = FALSE)
immune_df <- as.data.frame(t(immune_scores))
immune_df$sample <- rownames(immune_df)

# 2.4 Correlation: MP scores vs all immune scores
bulk_merged <- tcga_meta %>%
  left_join(mcp_df, by = "sample") %>%
  left_join(immune_df, by = "sample")

# ============================================================
# Part 3: Immune-related scores and checkpoints
# ============================================================

tpm_mat <- tcga_tpm[, tcga_meta$sample]

# CYT (Cytolytic Activity): sqrt(GZMA * PRF1) on raw TPM
tcga_meta$CYT <- sqrt(as.numeric(tpm_mat["GZMA", tcga_meta$sample]) *
                       as.numeric(tpm_mat["PRF1", tcga_meta$sample]))

# TIS (Tumor Inflammation Signature, 18-gene) via GSVA
tis_genes <- c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
               "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
               "PSMB10", "STAT1", "TIGIT")
param_tis <- gsvaParam(exprData = log_mat, geneSets = list(TIS = intersect(tis_genes, rownames(log_mat))),
                       kcdf = "Gaussian")
tcga_meta$TIS <- as.numeric(gsva(param_tis, verbose = FALSE)["TIS", tcga_meta$sample])

# IFN-gamma signature (Ayers et al. 6-gene) via GSVA
ifng_genes <- c("IFNG", "STAT1", "IDO1", "CXCL10", "CXCL9", "HLA-DRA")
param_ifng <- gsvaParam(exprData = log_mat, geneSets = list(IFNg = intersect(ifng_genes, rownames(log_mat))),
                        kcdf = "Gaussian")
tcga_meta$IFNg <- as.numeric(gsva(param_ifng, verbose = FALSE)["IFNg", tcga_meta$sample])

# T cell inflamed GEP (16-gene expanded) via GSVA
gep_genes <- c("CD8A", "CD8B", "CCL5", "CXCL9", "CXCL10", "CXCL11", "GZMA", "GZMB",
               "PRF1", "IFNG", "IDO1", "LAG3", "PDCD1", "CTLA4", "TIGIT", "HAVCR2")
param_gep <- gsvaParam(exprData = log_mat, geneSets = list(GEP = intersect(gep_genes, rownames(log_mat))),
                       kcdf = "Gaussian")
tcga_meta$GEP <- as.numeric(gsva(param_gep, verbose = FALSE)["GEP", tcga_meta$sample])

# Immune checkpoint gene expression (16 genes)
checkpoint_genes <- c("CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT",
                      "TNFRSF9", "ICOS", "CD28", "CD80", "CD86", "IDO1",
                      "BTLA", "VISTA", "TNFRSF4", "TNFRSF18")
checkpoint_avail <- intersect(checkpoint_genes, rownames(log_mat))
checkpoint_expr <- as.data.frame(t(log_mat[checkpoint_avail, tcga_meta$sample]))

# Correlation: MP1/MP7 vs immune scores and checkpoints
ckpt_cor <- data.frame()
for (mp in c("MP1_score", "MP7_score")) {
  for (gene in checkpoint_avail) {
    ct <- cor.test(tcga_meta[[mp]], checkpoint_expr[[gene]], method = "spearman")
    ckpt_cor <- rbind(ckpt_cor, data.frame(
      MP = gsub("_score", "", mp), Gene = gene,
      rho = ct$estimate, pval = ct$p.value))
  }
}
ckpt_cor$padj <- p.adjust(ckpt_cor$pval, method = "BH")

# ============================================================
# Part 4: Hallmark pathway differential activation
# ============================================================

hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
gs_hallmark <- lapply(hallmark_list, function(g) intersect(g, rownames(log_mat)))
gs_hallmark <- gs_hallmark[sapply(gs_hallmark, length) >= 10]

param_hm <- gsvaParam(exprData = log_mat[, tcga_meta$sample], geneSets = gs_hallmark, kcdf = "Gaussian")
hallmark_scores <- as.data.frame(t(gsva(param_hm, verbose = FALSE)))

# Stratify by MP1/MP7 median
tcga_meta$MP1_group <- ifelse(tcga_meta$MP1_score > median(tcga_meta$MP1_score, na.rm = TRUE), "High", "Low")
tcga_meta$MP7_group <- ifelse(tcga_meta$MP7_score > median(tcga_meta$MP7_score, na.rm = TRUE), "High", "Low")

# Wilcoxon test: High vs Low for each pathway
pathway_diff <- data.frame()
for (pw in colnames(hallmark_scores)) {
  for (mp_grp in c("MP1_group", "MP7_group")) {
    wt <- wilcox.test(hallmark_scores[[pw]] ~ tcga_meta[[mp_grp]])
    delta <- mean(hallmark_scores[[pw]][tcga_meta[[mp_grp]] == "High"]) -
             mean(hallmark_scores[[pw]][tcga_meta[[mp_grp]] == "Low"])
    pathway_diff <- rbind(pathway_diff, data.frame(
      MP = gsub("_group", "", mp_grp), pathway = pw,
      delta_mean = delta, p_value = wt$p.value))
  }
}
pathway_diff$padj <- p.adjust(pathway_diff$p_value, method = "BH")
