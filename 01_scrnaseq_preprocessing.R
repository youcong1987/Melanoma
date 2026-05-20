## 01_scRNAseq_preprocessing.R
## Single-cell RNA-seq preprocessing pipeline for melanoma datasets
## Datasets: GSE277165, GSE222446

library(Seurat)
library(tidyverse)
library(qs)
library(DoubletFinder)
library(parallel)

# ============================================================
# 1. Data loading and QC filtering
# ============================================================

filenames <- list.files("../rawdata/scRNAseq/GSE277165/02_cellranger/")
path <- paste0(getwd(), "/../rawdata/scRNAseq/GSE277165/02_cellranger/")

seurat.list <- mclapply(filenames, function(file) {
  seurat.data <- Read10X(paste0(path, file, "/outs/filtered_feature_bc_matrix/"))
  obj <- CreateSeuratObject(counts = seurat.data, min.cells = 5,
                            min.features = 200, project = file)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")

  # QC filtering
  obj <- subset(obj,
    subset = nCount_RNA > 500 & nCount_RNA < 40000 &
             nFeature_RNA > 200 & nFeature_RNA < 8000 &
             percent.mt < 20)

  # Standard workflow for DoubletFinder

  obj <- NormalizeData(obj) %>%
    FindVariableFeatures(nfeatures = 2000) %>%
    ScaleData(features = VariableFeatures(.)) %>%
    RunPCA()
  PCs <- 30
  obj <- RunUMAP(obj, dims = 1:PCs) %>%
    FindNeighbors(reduction = "pca", dims = 1:PCs) %>%
    FindClusters(resolution = 0.5)

  # DoubletFinder
  sweep.data <- paramSweep(obj, PCs = 1:PCs)
  sweep.stats <- summarizeSweep(sweep.data, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  homotypic.prop <- modelHomotypic(obj$RNA_snn_res.0.5)
  nExp_poi <- round(0.05 * ncol(obj))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  obj <- doubletFinder(obj, PCs = 1:PCs, pN = 0.25,
    pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])),
    nExp = nExp_poi.adj)
  colnames(obj@meta.data)[grepl("DF.class", colnames(obj@meta.data))] <- "DF_classifications"
  return(obj)
}, mc.cores = 10)

seurat <- merge(x = seurat.list[[1]], y = seurat.list[-1]) %>% JoinLayers()

# ============================================================
# 2. Standard Seurat workflow
# ============================================================

seurat <- subset(seurat, DF_classifications == "Singlet")

# Cell cycle scoring
data("cc.genes.updated.2019")
seurat <- NormalizeData(seurat) %>%
  CellCycleScoring(s.features = cc.genes.updated.2019$s.genes,
                   g2m.features = cc.genes.updated.2019$g2m.genes)
seurat$CC.Difference <- seurat$S.Score - seurat$G2M.Score

# Normalization, scaling with regression, PCA
seurat <- seurat %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt", "percent.rb",
                                 "nCount_RNA", "CC.Difference")) %>%
  RunPCA()

# Clustering
PCs <- 30
seurat <- RunUMAP(seurat, reduction = "pca", dims = 1:PCs) %>%
  FindNeighbors(reduction = "pca", dims = 1:PCs) %>%
  FindClusters(resolution = 2)

# ============================================================
# 3. Cell type annotation (major types)
# ============================================================

# Marker-based manual annotation assisted by COSG, scType, and SingleR
# Major cell types identified:
#   Melanoma cells, T/NK cells, B cells, Myeloid cells,
#   Fibroblasts, Endothelial cells, SMCs, Keratinocytes, Mast cells

library(COSG)
marker_cosg <- cosg(seurat, groups = "all", assay = "RNA", slot = "data",
                    mu = 1, expressed_pct = 0.1, n_genes_user = 200)

# Manual annotation based on canonical markers and automated tools
# [Cluster-to-celltype mapping applied via RenameIdents]
seurat$celltype_main <- Idents(seurat)

# ============================================================
# 4. Sub-clustering for each major cell type
# ============================================================

# Each major cell type was re-clustered independently to remove
# residual doublets and refine annotations.
# Sub-types annotated: Arterial/Venous/Capillary/Lymphatic ECs,
# CAFs, myofibroblasts, CD4/CD8/Treg/NK subtypes, etc.

# ============================================================
# 5. Dataset integration (GSE277165 + GSE222446)
# ============================================================

# Merge two datasets (no batch correction applied after evaluation)
seurat_GSE277165 <- qread("output/01_scRNAseq_preprocessing/GSE277165/seurat_annotated.qs")
seurat_GSE222446 <- qread("output/01_scRNAseq_preprocessing/GSE222446/seurat_annotated.qs")

seurat1 <- CreateSeuratObject(counts = GetAssayData(seurat_GSE277165, assay = "RNA", layer = "counts"),
                              meta.data = seurat_GSE277165@meta.data, project = "GSE277165")
seurat2 <- CreateSeuratObject(counts = GetAssayData(seurat_GSE222446, assay = "RNA", layer = "counts"),
                              meta.data = seurat_GSE222446@meta.data, project = "GSE222446")

seurat_merged <- merge(seurat1, seurat2) %>% JoinLayers()
seurat_merged$batch <- paste0(seurat_merged$orig.ident, "_", seurat_merged$project)

# Standard workflow on merged object (batch effect evaluated as minimal)
seurat_merged <- seurat_merged %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt", "percent.rb", "nCount_RNA", "CC.Difference")) %>%
  RunPCA()

PCs <- 30
seurat_merged <- RunUMAP(seurat_merged, reduction = "pca", dims = 1:PCs) %>%
  FindNeighbors(reduction = "pca", dims = 1:PCs) %>%
  FindClusters(resolution = 2)

# Final cell type annotation on merged object
# [Same marker-based strategy as above]
