## 04_bulk_preprocessing.R
## Bulk RNA-seq data preprocessing for TCGA-SKCM and ICB cohorts
## Outputs: normalized expression matrices + clinical metadata per cohort

library(tidyverse)
library(data.table)
library(GEOquery)

# ============================================================
# 1. TCGA-SKCM
# ============================================================

exp_counts <- fread("../rawdata/Bulk/TCGA-SKCM/TCGA-SKCM.star_counts.tsv.gz") %>%
  column_to_rownames("Ensembl_ID")
exp_tpm <- fread("../rawdata/Bulk/TCGA-SKCM/TCGA-SKCM.star_tpm.tsv.gz") %>%
  column_to_rownames("Ensembl_ID")
clinical <- fread("../rawdata/Bulk/TCGA-SKCM/TCGA-SKCM.clinical.tsv.gz")
survival_data <- fread("../rawdata/Bulk/TCGA-SKCM/TCGA-SKCM.survival.tsv.gz")
geneids <- fread("../rawdata/Bulk/TCGA-SKCM/g2s_vm25_gencode.txt",
                 col.names = c("Ensembl_ID", "gene_name", "gene_type"))

# Restore raw counts from log2(counts+1)
counts_raw <- round(2^as.matrix(exp_counts) - 1)

# Map Ensembl to gene symbol, keep protein-coding, deduplicate by max mean
counts <- as.data.frame(counts_raw) %>%
  rownames_to_column("Ensembl_ID") %>%
  inner_join(geneids %>% filter(gene_type == "protein_coding"), by = "Ensembl_ID") %>%
  mutate(rowmean = rowMeans(select(., -c(Ensembl_ID, gene_name, gene_type)))) %>%
  arrange(desc(rowmean)) %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  select(-rowmean, -Ensembl_ID, -gene_type) %>%
  column_to_rownames("gene_name")

# Same for TPM
tpm <- as.data.frame(exp_tpm) %>%
  rownames_to_column("Ensembl_ID") %>%
  inner_join(geneids %>% filter(gene_type == "protein_coding"), by = "Ensembl_ID") %>%
  mutate(rowmean = rowMeans(select(., -c(Ensembl_ID, gene_name, gene_type)))) %>%
  arrange(desc(rowmean)) %>%
  distinct(gene_name, .keep_all = TRUE) %>%
  select(-rowmean, -Ensembl_ID, -gene_type) %>%
  column_to_rownames("gene_name")

# Clinical metadata: tumor samples only
meta_all <- clinical %>%
  mutate(sample = submitter_id.samples) %>%
  filter(tissue_type.samples == "Tumor")

# Age (days -> years), sex, stage
meta_all$Age <- round(meta_all$age_at_diagnosis.diagnoses / 365)
meta_all$Sex <- meta_all$gender.demographic
meta_all$Stage_raw <- meta_all$ajcc_pathologic_stage.diagnoses
meta_all$Stage_raw[meta_all$Stage_raw %in% c("", "Not Reported", "Stage 0")] <- NA
meta_all$Stage_group <- ifelse(
  grepl("^Stage I[AB]?$|^Stage II[ABC]?$", meta_all$Stage_raw), "I-II",
  ifelse(is.na(meta_all$Stage_raw), NA, "III-IV"))

# Merge survival (OS, DSS, PFI)
metadata <- meta_all %>%
  mutate(sample_short = substr(sample, 1, 15)) %>%
  inner_join(survival_data, by = c("sample_short" = "sample"))

# Subset expression to tumor samples
common_samples <- intersect(colnames(counts), metadata$sample)
counts <- counts[, common_samples]
tpm <- tpm[, common_samples]
metadata <- metadata %>% filter(sample %in% common_samples)

save(counts, tpm, metadata,
     file = "output/04_Bulk_preprocessing/01_TCGA_SKCM/counts_tpm_metadata_Tumor_survival.Rdata")

# ============================================================
# 2. GSE78220 (Hugo et al. 2016, anti-PD1, melanoma, n=28)
# ============================================================

gse78220 <- getGEO("GSE78220", destdir = "../rawdata/Bulk/GSE78220/")
pdata <- pData(gse78220[[1]])

# Load TPM from supplementary
tpm_78220 <- fread("../rawdata/Bulk/GSE78220/GSE78220_PatientFPKM.txt.gz") %>%
  column_to_rownames("Gene")

metadata_78220 <- data.frame(
  sample = colnames(tpm_78220),
  Response = ifelse(pdata$`response:ch1` %in% c("Complete Response", "Partial Response"),
                    "Response",
             ifelse(pdata$`response:ch1` == "Progressive Disease",
                    "Non-response", NA)),  # SD excluded
  OS_days = as.numeric(pdata$`os (days):ch1`),
  OS_status = ifelse(pdata$`vital status:ch1` == "Dead", 1, 0),
  stringsAsFactors = FALSE
)

save(tpm_78220, metadata_78220,
     file = "output/04_Bulk_preprocessing/02_GSE78220/tpm_metadata.Rdata")

# ============================================================
# 3. PMID_31792460 (Hugo et al. 2016, anti-PD1, n=121)
# ============================================================

# Directly loads TPM (not raw counts)
tpm_31792460 <- read.table("../rawdata/Bulk/PMID_31792460/RNAseq_TPM.txt",
                           header = TRUE, row.names = 1, check.names = FALSE)
meta_31792460 <- fread("../rawdata/Bulk/PMID_31792460/clinical.tsv")

meta_31792460$Response <- ifelse(meta_31792460$BOR %in% c("CR", "PR"), "Response",
                          ifelse(meta_31792460$BOR %in% c("PD"), "Non-response", NA))
meta_31792460$OS_status <- ifelse(meta_31792460$vital_status == "Dead", 1, 0)

save(tpm_31792460, meta_31792460,
     file = "output/04_Bulk_preprocessing/03_PMID_31792460/tpm_metadata.Rdata")

# ============================================================
# 4. GSE91061 (Riaz et al. 2017, anti-PD1/CTLA4, n=109)
# ============================================================

gse91061 <- getGEO("GSE91061", destdir = "../rawdata/Bulk/GSE91061/")
pdata_91061 <- pData(gse91061[[1]])

tpm_91061 <- fread("../rawdata/Bulk/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz") %>%
  column_to_rownames("V1")

metadata_91061 <- data.frame(
  sample = colnames(tpm_91061),
  Response = ifelse(pdata_91061$`response:ch1` %in% c("CR", "PR"), "Response",
             ifelse(pdata_91061$`response:ch1` %in% c("PD", "SD"), "Non-response", NA)),
  Timepoint = pdata_91061$`timepoint:ch1`,
  stringsAsFactors = FALSE
)

save(tpm_91061, metadata_91061,
     file = "output/04_Bulk_preprocessing/04_GSE91061/tpm_metadata.Rdata")

# ============================================================
# 5. Gide2019 (anti-PD1+CTLA4, melanoma, n=72)
# ============================================================

tpm_gide <- fread("../rawdata/Bulk/Gide2019/Gide2019_TPM.txt") %>%
  column_to_rownames("Gene")
meta_gide <- fread("../rawdata/Bulk/Gide2019/Gide2019_clinical.txt")
meta_gide$Response <- ifelse(meta_gide$BOR %in% c("CR", "PR"), "Response", "Non-response")

save(tpm_gide, meta_gide,
     file = "output/04_Bulk_preprocessing/05_Gide2019/tpm_metadata.Rdata")

# ============================================================
# 6. GSE100797 (Lauss 2017, TIL-ACT, n=21)
# ============================================================

gse100797 <- getGEO("GSE100797", destdir = "../rawdata/Bulk/GSE100797/")
pdata_100797 <- pData(gse100797[[1]])

tpm_100797 <- fread("../rawdata/Bulk/GSE100797/GSE100797_expression.txt") %>%
  column_to_rownames("Gene")
meta_100797 <- data.frame(
  sample = colnames(tpm_100797),
  Response = ifelse(pdata_100797$`response:ch1` %in% c("CR", "PR"), "Response", "Non-response"),
  OS_days = as.numeric(pdata_100797$`os_days:ch1`),
  OS_status = as.numeric(pdata_100797$`os_event:ch1`),
  stringsAsFactors = FALSE
)

save(tpm_100797, meta_100797,
     file = "output/04_Bulk_preprocessing/06_GSE100797/tpm_metadata.Rdata")

# ============================================================
# 7. Additional ICI cohorts (from ICI_data-symbol.rda)
# ============================================================

# 19 ICI cohorts loaded from curated collection
# Filtered by: gene count >= 10000, min group size >= 6, ratio >= 0.15
load("../rawdata/ICI_cohort/ICI_data-symbol.rda")
