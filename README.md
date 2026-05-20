# Code Availability

R scripts for the analysis pipeline described in the manuscript. Scripts are organized sequentially following the analytical workflow.

## Requirements

- R >= 4.3.0
- Python >= 3.9 (for cNMF)
- Key R packages: Seurat (v5), tidyverse, GSVA, survival, survminer, NMF, copykat, infercnv, DoubletFinder, harmony, oncoPredict, Mime1, beyondcell, clusterProfiler

## Scripts

| Script | Description |
|--------|-------------|
| `01_scrnaseq_preprocessing.R` | scRNA-seq data loading, QC, doublet removal, normalization, clustering, and cell type annotation (GSE277165 + GSE222446) |
| `02_malignant_identification.R` | Malignant cell identification using CopyKAT, inferCNV, and cancer driver gene validation |
| `03_nmf_meta_programs.R` | Non-negative matrix factorization (cNMF + NMF-R) to identify melanoma meta-programs, with Jaccard similarity analysis |
| `04_bulk_preprocessing.R` | Bulk RNA-seq preprocessing for TCGA-SKCM and multiple ICB cohorts |
| `05_mp_icb_response.R` | Meta-program GSVA scoring and ICB response analysis with meta-analysis across cohorts |
| `06_tme_landscape.R` | MP1/MP7 association with tumor microenvironment: immune infiltration, checkpoint expression, and pathway activation |
| `07_core_gene_prognosis.R` | Multi-method core gene screening, ICB response modeling (7 classifiers), and survival analysis |
| `08_drug_screening.R` | Multi-strategy drug screening: oncoPredict, L1000CDS2, DGIdb, Beyondcell single-cell drug sensitivity |

## Data Availability

- Single-cell RNA-seq: GSE277165, GSE222446
- Bulk RNA-seq: TCGA-SKCM (UCSC Xena), GSE78220, GSE91061, GSE100797, GSE35640
- ICB cohorts: Gide et al. 2019, Hugo et al. 2016 (PMID: 31792460)

## Notes

- Scripts are simplified to show core analytical steps. Intermediate file paths are illustrative.
- Figure generation code is available upon request.
