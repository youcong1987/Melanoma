## 07_core_gene_prognosis.R
## Core gene screening -> ICB response modeling -> Prognosis validation
##
## Step 1: Multi-method core gene selection (5 methods, >= 5 required)
## Step 2: ICB response modeling (7 classifiers via Mime1)
## Step 3: Single-gene Cox regression filtering
## Step 4: GSVA-based survival analysis

library(tidyverse)
library(data.table)
library(survival)
library(survminer)
library(Mime1)
library(glmnet)
library(randomForest)
library(Boruta)
library(xgboost)
library(GSVA)
library(timeROC)
library(patchwork)

set.seed(5201314)

# ============================================================
# Step 1: Core gene screening (5 methods, >= 3 methods required)
# ============================================================

# Training data: merged ICB cohorts
# Load all ICB datasets
load("output/04_Bulk_preprocessing/01_TCGA_SKCM/counts_tpm_metadata_Tumor_survival.Rdata")
tcga_expr <- Tumor_tpm_surv; tcga_meta <- metadata_surv

# MP gene lists from NMF
mp1_genes <- readLines("output/03_NMF/MP1_genes.txt")
mp7_genes <- readLines("output/03_NMF/MP7_genes.txt")

# Build screening data: merge ICB cohorts, log2 transform, binary response
build_screen_data <- function(expr_list, meta_list, gene_list) {
  all_parts <- list()
  for (nm in names(expr_list)) {
    expr <- expr_list[[nm]]; meta <- meta_list[[nm]]
    common_genes <- intersect(gene_list, rownames(expr))
    common_s <- intersect(meta$sample, colnames(expr))
    ms <- meta[meta$sample %in% common_s & meta$Response %in% c("Response", "Non-response"), ]
    if (nrow(ms) == 0) next
    mat <- as.matrix(expr[common_genes, ms$sample, drop = FALSE])
    if (max(mat, na.rm = TRUE) > 30) mat <- log2(mat + 1)
    df <- data.frame(Response = ifelse(ms$Response == "Response", 1, 0),
                     t(mat), check.names = TRUE)
    all_parts <- c(all_parts, list(df))
  }
  cols <- Reduce(intersect, lapply(all_parts, colnames))
  do.call(rbind, lapply(all_parts, function(d) d[, cols]))
}

# 5-method screening function
screen_core_genes <- function(screen_df, n_lasso = 1000, min_methods = 3) {
  X <- as.matrix(screen_df[, -1])
  y <- screen_df$Response

  # Method 1: Wilcoxon rank-sum test (FDR < 0.05)
  wilcox_genes <- sapply(colnames(X), function(g) {
    wt <- wilcox.test(X[y == 1, g], X[y == 0, g])
    wt$p.value
  })
  wilcox_genes <- names(wilcox_genes[p.adjust(wilcox_genes, "BH") < 0.05])

  # Method 2: LASSO (repeated 1000x, frequency >= 80%)
  lasso_freq <- rep(0, ncol(X)); names(lasso_freq) <- colnames(X)
  for (i in seq_len(n_lasso)) {
    cv_fit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 5)
    coefs <- coef(cv_fit, s = "lambda.min")[-1, ]
    lasso_freq[names(coefs[coefs != 0])] <- lasso_freq[names(coefs[coefs != 0])] + 1
  }
  lasso_genes <- names(lasso_freq[lasso_freq / n_lasso >= 0.8])

  # Method 3: Random Forest (importance > 0)
  rf_fit <- randomForest(x = X, y = factor(y), ntree = 1000, importance = TRUE)
  imp <- importance(rf_fit, type = 1)
  rf_genes <- names(sort(imp[, 1], decreasing = TRUE))[imp[, 1] > 0]

  # Method 4: Boruta
  bor <- Boruta(x = X, y = factor(y), maxRuns = 300, doTrace = 0)
  bor_final <- TentativeRoughFix(bor)
  boruta_genes <- names(bor_final$finalDecision[bor_final$finalDecision == "Confirmed"])

  # Method 5: XGBoost gain importance (> median gain)
  dtrain <- xgb.DMatrix(data = X, label = y)
  params <- list(objective = "binary:logistic", eval_metric = "auc",
                 max_depth = 4, eta = 0.1, subsample = 0.8, colsample_bytree = 0.8)
  xgb_fit <- xgb.train(params, dtrain, nrounds = 200, verbose = 0)
  imp_xgb <- xgb.importance(model = xgb_fit)
  xgb_genes <- imp_xgb$Feature[imp_xgb$Gain > median(imp_xgb$Gain)]

  # Consensus: genes selected by >= min_methods
  all_selected <- list(Wilcoxon = wilcox_genes, Lasso = lasso_genes,
                       RF = rf_genes, Boruta = boruta_genes, XGBoost = xgb_genes)
  gene_counts <- table(unlist(all_selected))
  core <- names(gene_counts[gene_counts >= min_methods])
  gsub("\\.", "-", core)  # restore hyphens
}

core_mp1 <- screen_core_genes(build_screen_data(expr_list, meta_list, mp1_genes))
core_mp7 <- screen_core_genes(build_screen_data(expr_list, meta_list, mp7_genes))

# ============================================================
# Step 2: ICB response modeling (Mime1 package, 7 classifiers)
# ============================================================

# Mime1 data format: columns = ID, Var (Y/N), gene1, gene2, ...
make_icb_data <- function(expr_mat, meta_df, core_genes) {
  common_genes <- intersect(make.names(core_genes), make.names(rownames(expr_mat)))
  common_s <- intersect(meta_df$sample, colnames(expr_mat))
  ms <- meta_df[meta_df$sample %in% common_s & meta_df$Response %in% c("Response", "Non-response"), ]
  mat <- as.matrix(expr_mat[core_genes[make.names(core_genes) %in% common_genes], ms$sample])
  if (max(mat, na.rm = TRUE) > 30) mat <- log2(mat + 1)
  df <- data.frame(ID = ms$sample,
                   Var = ifelse(ms$Response == "Response", "Y", "N"),
                   t(mat), check.names = TRUE, stringsAsFactors = FALSE)
  rownames(df) <- NULL
  df
}

methods_icb <- c("nb", "svmRadialWeights", "rf", "kknn",
                 "adaboost", "LogitBoost", "cancerclass")

for (mp in c("MP1", "MP7")) {
  core_genes <- if (mp == "MP1") core_mp1 else core_mp7

  # Training: merged cohorts (Gide2019 + PMID_31792460 + GSE91061)
  train_parts <- lapply(train_specs, function(s) make_icb_data(s$expr, s$meta, core_genes))
  cols <- Reduce(intersect, lapply(train_parts, colnames))
  train_df <- do.call(rbind, lapply(train_parts, function(d) d[, cols]))

  # Validation datasets
  list_icb <- list(training = train_df)
  for (v in vali_specs) {
    vdf <- make_icb_data(v$expr, v$meta, core_genes)
    if (!is.null(vdf) && nrow(vdf) >= 10) list_icb[[v$name]] <- vdf
  }

  res_icb <- ML.Dev.Pred.Category.Sig(
    train_data = train_df,
    list_train_vali_Data = list_icb,
    candidate_genes = make.names(core_genes),
    methods = methods_icb,
    seed = 5201314, cores_for_parallel = 10)
}

# ============================================================
# Step 3: Single-gene Cox regression filtering
# ============================================================

# For each core gene, univariate Cox in multiple survival datasets
# MP1 (protective): retain HR < 1, p < 0.05 in >= 1 dataset
# MP7 (risk): retain HR > 1, p < 0.05 in >= 1 dataset

single_gene_cox <- function(expr_mat, meta_df, gene_list, time_col, status_col) {
  log_expr <- if (max(expr_mat, na.rm = TRUE) > 30) log2(as.matrix(expr_mat) + 1) else as.matrix(expr_mat)
  common_genes <- intersect(gene_list, rownames(log_expr))
  common_s <- intersect(colnames(log_expr), meta_df$sample)
  meta_sub <- meta_df[meta_df$sample %in% common_s, ]
  meta_sub <- meta_sub[!is.na(meta_sub[[time_col]]) & meta_sub[[time_col]] > 0, ]

  results <- lapply(common_genes, function(g) {
    sdf <- data.frame(
      time = meta_sub[[time_col]],
      status = meta_sub[[status_col]],
      gene_expr = as.numeric(log_expr[g, meta_sub$sample]))
    cox <- tryCatch(coxph(Surv(time, status) ~ gene_expr, data = sdf), error = function(e) NULL)
    if (is.null(cox)) return(NULL)
    s <- summary(cox)
    data.frame(gene = g, HR = s$conf.int[1, 1],
               HR_lower = s$conf.int[1, 3], HR_upper = s$conf.int[1, 4],
               p = s$coefficients[1, 5])
  })
  do.call(rbind, results)
}

# MP1: TCGA-SKCM OS + PMID_31792460 OS (retain HR < 1 in any)
cox_mp1_tcga <- single_gene_cox(tcga_expr, tcga_meta, core_mp1, "OS.time", "OS")
cox_mp1_pmid <- single_gene_cox(pmid_expr, pmid_meta, core_mp1, "OS_days", "OS_status")
mp1_final <- unique(c(
  cox_mp1_tcga$gene[cox_mp1_tcga$p < 0.05 & cox_mp1_tcga$HR < 1],
  cox_mp1_pmid$gene[cox_mp1_pmid$p < 0.05 & cox_mp1_pmid$HR < 1]))

# MP7: PMID_31792460 OS + VanAllen OS (retain HR > 1 in any)
cox_mp7_pmid <- single_gene_cox(pmid_expr, pmid_meta, core_mp7, "OS_days", "OS_status")
cox_mp7_van <- single_gene_cox(vanallen_expr, vanallen_meta, core_mp7, "OS_days", "OS_status")
mp7_final <- unique(c(
  cox_mp7_pmid$gene[cox_mp7_pmid$p < 0.05 & cox_mp7_pmid$HR > 1],
  cox_mp7_van$gene[cox_mp7_van$p < 0.05 & cox_mp7_van$HR > 1]))

# ============================================================
# Step 4: GSVA-based survival analysis
# ============================================================

run_survival_analysis <- function(expr_mat, meta_df, gene_set, time_col, status_col) {
  log_mat <- if (max(expr_mat, na.rm = TRUE) > 30) log2(as.matrix(expr_mat) + 1) else as.matrix(expr_mat)
  common_genes <- intersect(gene_set, rownames(log_mat))
  common_s <- intersect(colnames(log_mat), meta_df$sample)
  log_mat <- log_mat[, common_s]
  meta_sub <- meta_df[meta_df$sample %in% common_s, ]
  meta_sub <- meta_sub[!is.na(meta_sub[[time_col]]) & meta_sub[[time_col]] > 0, ]
  log_mat <- log_mat[, meta_sub$sample]

  # ssGSEA scoring
  param <- ssgseaParam(exprData = log_mat, geneSets = list(signature = common_genes))
  scores <- as.numeric(gsva(param, verbose = FALSE))
  names(scores) <- colnames(log_mat)

  surv_df <- data.frame(
    sample = meta_sub$sample,
    time = meta_sub[[time_col]],
    status = meta_sub[[status_col]],
    score = scores[meta_sub$sample]) %>%
    filter(!is.na(time) & time > 0 & !is.na(score))

  surv_df$Group <- ifelse(surv_df$score > median(surv_df$score), "High", "Low")

  # KM + log-rank
  fit <- survfit(Surv(time, status) ~ Group, data = surv_df)
  logrank_p <- survdiff(Surv(time, status) ~ Group, data = surv_df)$pvalue

  # Cox
  cox <- coxph(Surv(time, status) ~ score, data = surv_df)

  # Time-dependent ROC (1, 3, 5 years)
  troc <- timeROC(T = surv_df$time, delta = surv_df$status,
                  marker = surv_df$score, cause = 1,
                  times = c(12, 36, 60), iid = TRUE)

  list(fit = fit, cox = cox, logrank_p = logrank_p, troc = troc)
}

# Run across multiple datasets for MP1 and MP7
# MP1: TCGA-SKCM(OS), PMID_31792460(OS/PFS), Gide2019(PFS), VanAllen(OS/PFS)
# MP7: PMID_31792460(OS/PFS), Gide2019(PFS), VanAllen(OS/PFS)

surv_mp1_tcga <- run_survival_analysis(tcga_expr, tcga_meta, mp1_final, "OS.time", "OS")
surv_mp1_pmid <- run_survival_analysis(pmid_expr, pmid_meta, mp1_final, "OS_days", "OS_status")
surv_mp7_pmid <- run_survival_analysis(pmid_expr, pmid_meta, mp7_final, "OS_days", "OS_status")
