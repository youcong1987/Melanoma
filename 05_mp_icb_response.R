## 05_mp_icb_response.R
## Meta-program GSVA scoring and ICB response analysis across multiple cohorts
## Scoring: GSVA (gsvaParam, kcdf="Gaussian")
## Statistics: Wilcoxon + Z-based effect size + AUC
## Meta-analysis: metagen (fixed + random effects) by cancer type

library(tidyverse)
library(data.table)
library(GSVA)
library(pROC)
library(meta)
library(ggpubr)
library(patchwork)

set.seed(123)

# ============================================================
# 1. Load MP gene lists and ICB cohort data
# ============================================================

mp_genes_tbl <- read_tsv("output/03_NMF/meta_programs_genes.tsv")
mp_gene_list <- split(mp_genes_tbl$gene, mp_genes_tbl$metaprogram)

# Cancer type mapping
cancer_type_map <- c(
  GSE100797 = "SKCM", GSE35640 = "SKCM", GSE78220 = "SKCM",
  GSE91061 = "SKCM", Nathanson = "SKCM", VanAllen = "SKCM",
  PMID_31792460 = "SKCM", Gide2019 = "SKCM",
  IMvigor210 = "UC", GSE173839 = "BRCA")

# Load curated ICI cohorts
load("../rawdata/ICI_cohort/ICI_data-symbol.rda")

# Filter cohorts: gene count >= 10000, min group size >= 6, ratio >= 0.15
MIN_GENES <- 10000; MIN_GROUP_SIZE <- 6; MIN_RATIO <- 0.15

# ============================================================
# 2. GSVA scoring for each cohort
# ============================================================

auto_log_transform <- function(mat) {
  mat <- as.matrix(mat)
  if (min(mat, na.rm = TRUE) >= 0 && max(mat, na.rm = TRUE) > 20) log2(mat + 1) else mat
}

calc_gsva_score <- function(expr_mat, gene_sets) {
  log_mat <- auto_log_transform(expr_mat)
  log_mat <- log_mat[complete.cases(log_mat), , drop = FALSE]
  gene_sets_filtered <- lapply(gene_sets, function(g) intersect(g, rownames(log_mat)))
  gene_sets_filtered <- gene_sets_filtered[sapply(gene_sets_filtered, length) >= 5]
  if (length(gene_sets_filtered) == 0) return(NULL)
  param <- gsvaParam(exprData = log_mat, geneSets = gene_sets_filtered, kcdf = "Gaussian")
  as.data.frame(t(gsva(param, verbose = FALSE)))
}

# Score all cohorts
scores_list <- list()
for (ds_name in names(datasets_expr)) {
  sc <- calc_gsva_score(datasets_expr[[ds_name]], mp_gene_list)
  if (is.null(sc)) next
  meta <- datasets_meta[[ds_name]]
  common_samples <- intersect(meta$sample, rownames(sc))
  meta_aligned <- meta[match(common_samples, meta$sample), ]
  scores_list[[ds_name]] <- cbind(meta_aligned, sc[common_samples, , drop = FALSE])
}

# ============================================================
# 3. ICB response analysis per cohort per MP
# ============================================================

run_icb_analysis <- function(data, mp_col, response_col, dataset_name) {
  df <- data[, c(mp_col, response_col), drop = FALSE]
  colnames(df) <- c("score", "response")
  df <- df %>% filter(!is.na(score) & response %in% c("Response", "Non-response"))
  if (nrow(df) < 6 || min(table(df$response)) < 3) return(NULL)

  df$response <- factor(df$response, levels = c("Response", "Non-response"))

  # Wilcoxon test
  wt <- wilcox.test(score ~ response, data = df)

  # Effect size: Z-based r = |Z| / sqrt(N)
  z_val <- qnorm(wt$p.value / 2)
  r_effect <- abs(z_val) / sqrt(nrow(df))

  # Direction sign
  direction_sign <- ifelse(mean(df$score[df$response == "Response"]) >
                           mean(df$score[df$response == "Non-response"]), 1, -1)

  # ROC-AUC
  df$resp_binary <- ifelse(df$response == "Response", 1, 0)
  roc_obj <- tryCatch(roc(df$resp_binary, df$score, direction = "<", quiet = TRUE),
                      error = function(e) NULL)
  auc_val <- auc_lo <- auc_hi <- NA
  if (!is.null(roc_obj)) {
    auc_val <- as.numeric(auc(roc_obj))
    ci_obj <- tryCatch(ci.auc(roc_obj, method = "delong"), error = function(e) NULL)
    if (!is.null(ci_obj)) { auc_lo <- ci_obj[1]; auc_hi <- ci_obj[3] }
  }

  data.frame(
    dataset = dataset_name, MP = mp_col,
    cancer_type = cancer_type_map[dataset_name],
    n_response = sum(df$response == "Response"),
    n_nonresponse = sum(df$response == "Non-response"),
    mean_response = mean(df$score[df$response == "Response"]),
    mean_nonresponse = mean(df$score[df$response == "Non-response"]),
    wilcox_p = wt$p.value,
    effect_size = r_effect,
    direction_sign = direction_sign,
    auc = auc_val, auc_lower = auc_lo, auc_upper = auc_hi)
}

# Run for all cohorts x all MPs
all_icb <- data.frame()
for (ds_name in names(scores_list)) {
  for (mp in names(mp_gene_list)) {
    if (!mp %in% colnames(scores_list[[ds_name]])) next
    res <- run_icb_analysis(scores_list[[ds_name]], mp, "Response", ds_name)
    if (!is.null(res)) all_icb <- rbind(all_icb, res)
  }
}

# Signed -log10(p)
all_icb$signed_logp <- -log10(all_icb$wilcox_p) * all_icb$direction_sign

# ============================================================
# 4. Meta-analysis by cancer type
# ============================================================

meta_results <- data.frame()
cancer_types_for_meta <- unique(all_icb$cancer_type[!is.na(all_icb$cancer_type)])

for (ct in cancer_types_for_meta) {
  ct_data <- all_icb %>%
    filter(cancer_type == ct & !is.na(wilcox_p) & !is.na(effect_size) &
           n_response >= 3 & n_nonresponse >= 3)
  if (length(unique(ct_data$dataset)) < 2) next

  for (mp in names(mp_gene_list)) {
    mp_data <- ct_data %>% filter(MP == mp)
    if (nrow(mp_data) < 2) next

    # SE approximation: 1/sqrt(N)
    mp_data$se <- 1 / sqrt(mp_data$n_response + mp_data$n_nonresponse)
    mp_data$signed_effect <- mp_data$effect_size * mp_data$direction_sign

    m <- tryCatch(
      metagen(TE = mp_data$signed_effect, seTE = mp_data$se,
              studlab = mp_data$dataset, sm = "MD",
              fixed = TRUE, random = TRUE),
      error = function(e) NULL)
    if (is.null(m)) next

    meta_results <- rbind(meta_results, data.frame(
      cancer_type = ct, MP = mp,
      n_studies = nrow(mp_data),
      fixed_effect = m$TE.fixed, fixed_p = m$pval.fixed,
      random_effect = m$TE.random, random_p = m$pval.random,
      I2 = m$I2 * 100, tau2 = m$tau2))
  }
}

# ============================================================
# 5. Visualization
# ============================================================

# Signed -log10(p) heatmap
heatmap_df <- all_icb %>%
  filter(!is.na(wilcox_p)) %>%
  select(dataset, MP, signed_logp) %>%
  pivot_wider(names_from = MP, values_from = signed_logp) %>%
  column_to_rownames("dataset")

# AUC heatmap
auc_df <- all_icb %>%
  filter(!is.na(auc)) %>%
  select(dataset, MP, auc) %>%
  pivot_wider(names_from = MP, values_from = auc) %>%
  column_to_rownames("dataset")

# Forest plots for significant MPs
for (mp in names(mp_gene_list)) {
  for (ct in cancer_types_for_meta) {
    mp_data <- all_icb %>%
      filter(MP == mp, cancer_type == ct, !is.na(effect_size), n_response >= 3, n_nonresponse >= 3)
    if (nrow(mp_data) < 2) next
    mp_data$se <- 1 / sqrt(mp_data$n_response + mp_data$n_nonresponse)
    mp_data$signed_effect <- mp_data$effect_size * mp_data$direction_sign
    m <- metagen(TE = mp_data$signed_effect, seTE = mp_data$se,
                 studlab = mp_data$dataset, sm = "MD", fixed = TRUE, random = TRUE)
    forest(m, leftcols = "studlab", smlab = paste(mp, ct, "(GSVA)"))
  }
}
