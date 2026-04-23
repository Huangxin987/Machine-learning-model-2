# ==============================================================================
# PMI Multi-Omics Joint Modeling v7
# ==============================================================================

desktop_path <- ifelse(Sys.info()["sysname"] == "Windows",
                       file.path(Sys.getenv("USERPROFILE"), "Desktop"),
                       file.path(Sys.getenv("HOME"), "Desktop"))
output_dir <- file.path(desktop_path,
                        paste0("PMI_MultiOmics_v7_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(tidyverse); library(glmnet); library(patchwork)
  library(ggrepel);   library(randomForest)
})

# ── Global parameters ──────────────────────────────────────────────────────────
TISSUES    <- c("Heart", "Liver", "Spleen", "Lung", "Muscle")
META_COLS  <- c("SampleID", "AnimalID", "Organ", "PMI_day")
COL_TISSUE <- "Organ"; COL_ANIMAL <- "AnimalID"
COL_PMI    <- "PMI_day"; COL_ID <- "SampleID"
PREV_CUT   <- 0.20; ABUND_CUT <- 0.001; PSEUDO <- 1e-6
ALPHA      <- 0.5;  N_REP <- 50; TRAIN_FRAC <- 2/3; SEED_BASE <- 1000
TOP_K_METAB_LIST <- c(10L, 15L, 20L, 30L)
TOP_K_MICRO_LIST <- c(5L, 10L, 15L, 20L)
KM_SELECTED <- 0L; KC_SELECTED <- 0L
BEST_M <- setNames(rep("ElasticNet", length(TISSUES)), TISSUES)
BEST_C <- setNames(rep("ElasticNet", length(TISSUES)), TISSUES)

MODEL_COLORS_2     <- c("ElasticNet" = "#6BAED6", "RandomForest" = "#FD8D3C")
MODEL_COLORS_MULTI <- c("Metab_Best" = "#C6DBEF", "Micro_Best" = "#FDD0A2",
                        "EF_All(ElasticNet)" = "#BDBDBD", "Stacking(Best)" = "#D62728")
MODEL_COLORS_TOPK  <- c("Metab_TopK" = "#FDD0A2", "Micro_TopK" = "#ECEAEF",
                        "EF_TopK" = "#2C7BB6", "Stacking_TopK" = "#D7191C")
MODEL_COLORS_BASE  <- c("Metab" = "#F4A460", "Micro" = "#9370DB",
                        "EF" = "#2C7BB6", "Stacking" = "#D7191C")
OMICS_BAR_COLORS   <- c("Metabolomics" = "#2980B9", "Microbiome" = "#E67E22")
METAB_TREND_COLOR  <- "#1A5276"; MICRO_TREND_COLOR <- "#784212"

# ── Utility functions ──────────────────────────────────────────────────────────
calc_r2   <- function(o, p) 1 - sum((o-p)^2) / sum((o-mean(o))^2)
calc_rmse <- function(o, p) sqrt(mean((o-p)^2))
calc_mae  <- function(o, p) mean(abs(o-p))

metric_row <- function(model, o_tr, p_tr, o_te, p_te, nfeat = NA_integer_)
  data.frame(Model = model, R2_train = calc_r2(o_tr, p_tr),
             RMSE_train = calc_rmse(o_tr, p_tr), R2_test = calc_r2(o_te, p_te),
             RMSE_test = calc_rmse(o_te, p_te), MAE_test = calc_mae(o_te, p_te),
             n_features = as.integer(nfeat))

split_animals <- function(meta, seed) {
  set.seed(seed)
  ap <- meta %>% distinct(.data[[COL_ANIMAL]], .data[[COL_PMI]])
  tr_ids <- c()
  for (p in sort(unique(ap[[COL_PMI]]))) {
    ids    <- ap[[COL_ANIMAL]][ap[[COL_PMI]] == p]
    k      <- max(1L, min(round(length(ids) * TRAIN_FRAC), length(ids) - 1L))
    tr_ids <- c(tr_ids, sample(ids, k, replace = FALSE))
  }
  tr_ids <- unique(tr_ids)
  list(train_idx = meta[[COL_ANIMAL]] %in% tr_ids,
       test_idx  = !(meta[[COL_ANIMAL]] %in% tr_ids))
}

# ── Preprocessing ──────────────────────────────────────────────────────────────
preprocess_metab <- function(X_tr) {
  vars <- apply(X_tr, 2, var, na.rm = TRUE)
  k1   <- names(vars)[!is.na(vars) & vars > 0]
  Xt   <- X_tr[, k1, drop = FALSE]
  mu   <- colMeans(Xt, na.rm = TRUE)
  sdv  <- apply(Xt, 2, sd, na.rm = TRUE); sdv[sdv == 0] <- 1
  list(feats = k1, mu = mu, sdv = sdv,
       X_train = as.data.frame(sweep(sweep(as.matrix(Xt), 2, mu, "-"), 2, sdv, "/")))
}
apply_metab <- function(prep, X_te) {
  Xt <- X_te[, prep$feats, drop = FALSE]
  as.data.frame(sweep(sweep(as.matrix(Xt), 2, prep$mu, "-"), 2, prep$sdv, "/"))
}
clr_mat <- function(X) {
  mat <- as.matrix(X) + PSEUDO
  sweep(log(mat), 1, rowMeans(log(mat)), "-")
}
preprocess_micro <- function(X_tr) {
  prev  <- colMeans(X_tr > 0, na.rm = TRUE)
  abund <- colMeans(X_tr, na.rm = TRUE)
  k     <- names(prev)[prev >= PREV_CUT & abund > ABUND_CUT]
  Xt    <- X_tr[, k, drop = FALSE]
  Xc    <- as.data.frame(clr_mat(Xt))
  mu    <- colMeans(Xc, na.rm = TRUE)
  sdv   <- apply(Xc, 2, sd, na.rm = TRUE); sdv[sdv == 0] <- 1
  list(feats = k, mu = mu, sdv = sdv,
       X_train = as.data.frame(sweep(sweep(as.matrix(Xc), 2, mu, "-"), 2, sdv, "/")))
}
apply_micro <- function(prep, X_te) {
  Xt <- X_te[, prep$feats, drop = FALSE]
  Xc <- as.data.frame(clr_mat(Xt))
  as.data.frame(sweep(sweep(as.matrix(Xc), 2, prep$mu, "-"), 2, prep$sdv, "/"))
}

# ── Model functions ────────────────────────────────────────────────────────────
fit_enet <- function(Xtr, y_tr, Xte, alpha = ALPHA) {
  nf    <- min(5L, max(3L, nrow(Xtr)))
  cvfit <- cv.glmnet(as.matrix(Xtr), y_tr, alpha = alpha, nfolds = nf)
  p_tr  <- drop(predict(cvfit, as.matrix(Xtr), s = "lambda.1se"))
  p_te  <- drop(predict(cvfit, as.matrix(Xte), s = "lambda.1se"))
  cm    <- coef(cvfit, s = "lambda.1se")
  cv    <- setNames(as.numeric(cm)[-1], rownames(cm)[-1])
  list(pred_train = p_tr, pred_test = p_te,
       n_features = as.integer(sum(cv != 0)), importance = abs(cv))
}
fit_rf <- function(Xtr, y_tr, Xte) {
  p      <- ncol(Xtr)
  mtry_c <- unique(pmax(1L, c(floor(sqrt(p)), floor(p/3L), floor(p/2L))))
  best_oob <- Inf; best_mtry <- mtry_c[1]
  for (m in mtry_c) {
    rf_t <- randomForest(x = as.matrix(Xtr), y = y_tr, ntree = 200L, mtry = m)
    oob  <- tail(rf_t$mse, 1L)
    if (oob < best_oob) { best_oob <- oob; best_mtry <- m }
  }
  rf  <- randomForest(x = as.matrix(Xtr), y = y_tr,
                      ntree = 500L, mtry = best_mtry, importance = TRUE)
  imp <- importance(rf)[, "%IncMSE"]; imp[imp < 0] <- 0
  list(pred_train = predict(rf, as.matrix(Xtr)),
       pred_test  = predict(rf, as.matrix(Xte)),
       n_features = p, importance = imp)
}
call_model <- function(mt, Xtr, y_tr, Xte)
  switch(mt, "ElasticNet" = fit_enet(Xtr, y_tr, Xte),
         "RandomForest" = fit_rf(Xtr, y_tr, Xte))

fit_topk_model <- function(Xtr, y_tr, Xte, k, model_type = "ElasticNet") {
  k         <- min(as.integer(k), ncol(Xtr))
  fit_full  <- call_model(model_type, Xtr, y_tr, Xte)
  top_feats <- names(sort(fit_full$importance, decreasing = TRUE))[seq_len(k)]
  fit_k     <- call_model(model_type, Xtr[, top_feats, drop=FALSE],
                          y_tr, Xte[, top_feats, drop=FALSE])
  list(pred_train = fit_k$pred_train, pred_test = fit_k$pred_test,
       selected_feats = top_feats, k_actual = k)
}

fit_stacking_v2 <- function(Xm_tr, Xc_tr, y_tr, Xm_te, Xc_te,
                            model_type_m = "ElasticNet",
                            model_type_c = "ElasticNet", kfold = 5) {
  n     <- length(y_tr)
  kfold <- min(kfold, max(3L, floor(n/3L)))
  set.seed(42)
  fid   <- sample(rep(seq_len(kfold), length.out = n))
  oof_m <- oof_c <- rep(NA_real_, n)
  for (f in seq_len(kfold)) {
    val <- fid == f; trn <- !val
    r <- call_model(model_type_m, Xm_tr[trn,,drop=FALSE], y_tr[trn], Xm_tr[val,,drop=FALSE])
    oof_m[val] <- r$pred_test
    r <- call_model(model_type_c, Xc_tr[trn,,drop=FALSE], y_tr[trn], Xc_tr[val,,drop=FALSE])
    oof_c[val] <- r$pred_test
  }
  meta_lm <- lm(y ~ m + c, data = data.frame(y = y_tr, m = oof_m, c = oof_c))
  fit_m   <- call_model(model_type_m, Xm_tr, y_tr, Xm_te)
  fit_c   <- call_model(model_type_c, Xc_tr, y_tr, Xc_te)
  p_tr <- as.numeric(predict(meta_lm, data.frame(m = fit_m$pred_train, c = fit_c$pred_train)))
  p_te <- as.numeric(predict(meta_lm, data.frame(m = fit_m$pred_test,  c = fit_c$pred_test)))
  list(pred_train = p_tr, pred_test = p_te)
}

fit_ef_topk <- function(Xm_tr, Xc_tr, y_tr, Xm_te, Xc_te,
                        model_type_m, model_type_c, km, kc,
                        fit_m_full = NULL, fit_c_full = NULL) {
  km <- min(as.integer(km), ncol(Xm_tr))
  kc <- min(as.integer(kc), ncol(Xc_tr))
  top_m <- names(sort(fit_m_full$importance, decreasing = TRUE))[seq_len(km)]
  top_c <- names(sort(fit_c_full$importance, decreasing = TRUE))[seq_len(kc)]
  Xmk_tr <- Xm_tr[, top_m, drop=FALSE]; Xmk_te <- Xm_te[, top_m, drop=FALSE]
  Xck_tr <- Xc_tr[, top_c, drop=FALSE]; Xck_te <- Xc_te[, top_c, drop=FALSE]
  colnames(Xmk_tr) <- colnames(Xmk_te) <- paste0("M_", top_m)
  colnames(Xck_tr) <- colnames(Xck_te) <- paste0("C_", top_c)
  fit_ef <- fit_enet(cbind(Xmk_tr, Xck_tr), y_tr, cbind(Xmk_te, Xck_te))
  list(pred_train = fit_ef$pred_train, pred_test = fit_ef$pred_test,
       top_m = top_m, top_c = top_c, n_features = as.integer(km + kc))
}

fit_stacking_topk <- function(Xm_tr, Xc_tr, y_tr, Xm_te, Xc_te,
                              model_type_m = "ElasticNet",
                              model_type_c = "ElasticNet", km, kc, kfold = 5) {
  n  <- length(y_tr)
  km <- min(as.integer(km), ncol(Xm_tr))
  kc <- min(as.integer(kc), ncol(Xc_tr))
  fit_m_full <- call_model(model_type_m, Xm_tr, y_tr, Xm_te)
  fit_c_full <- call_model(model_type_c, Xc_tr, y_tr, Xc_te)
  top_m <- names(sort(fit_m_full$importance, decreasing = TRUE))[seq_len(km)]
  top_c <- names(sort(fit_c_full$importance, decreasing = TRUE))[seq_len(kc)]
  Xm_tr_k <- Xm_tr[, top_m, drop=FALSE]; Xm_te_k <- Xm_te[, top_m, drop=FALSE]
  Xc_tr_k <- Xc_tr[, top_c, drop=FALSE]; Xc_te_k <- Xc_te[, top_c, drop=FALSE]
  kfold <- min(kfold, max(3L, floor(n/3L)))
  set.seed(42)
  fid   <- sample(rep(seq_len(kfold), length.out = n))
  oof_m <- oof_c <- rep(NA_real_, n)
  for (f in seq_len(kfold)) {
    val <- fid == f; trn <- !val
    r <- call_model(model_type_m, Xm_tr_k[trn,,drop=FALSE], y_tr[trn], Xm_tr_k[val,,drop=FALSE])
    oof_m[val] <- r$pred_test
    r <- call_model(model_type_c, Xc_tr_k[trn,,drop=FALSE], y_tr[trn], Xc_tr_k[val,,drop=FALSE])
    oof_c[val] <- r$pred_test
  }
  meta_lm <- lm(y ~ m + c, data = data.frame(y = y_tr, m = oof_m, c = oof_c))
  fit_m_k <- call_model(model_type_m, Xm_tr_k, y_tr, Xm_te_k)
  fit_c_k <- call_model(model_type_c, Xc_tr_k, y_tr, Xc_te_k)
  p_tr <- as.numeric(predict(meta_lm, data.frame(m = fit_m_k$pred_train, c = fit_c_k$pred_train)))
  p_te <- as.numeric(predict(meta_lm, data.frame(m = fit_m_k$pred_test,  c = fit_c_k$pred_test)))
  list(pred_train = p_tr, pred_test = p_te,
       top_m = top_m, top_c = top_c,
       fit_m_full = fit_m_full, fit_c_full = fit_c_full)
}

# ── Run engines ────────────────────────────────────────────────────────────────
run_one_rep_compare <- function(meta, X_data, seed, omics_type = "metab") {
  sp   <- split_animals(meta, seed)
  tr   <- sp$train_idx; te <- sp$test_idx
  y_tr <- meta[[COL_PMI]][tr]; y_te <- meta[[COL_PMI]][te]
  prep <- if (omics_type == "metab") preprocess_metab(X_data[tr,]) else preprocess_micro(X_data[tr,])
  Xtr  <- prep$X_train
  Xte  <- if (omics_type == "metab") apply_metab(prep, X_data[te,]) else apply_micro(prep, X_data[te,])
  bind_rows(lapply(c("ElasticNet", "RandomForest"), function(nm) {
    fit <- call_model(nm, Xtr, y_tr, Xte)
    metric_row(nm, y_tr, fit$pred_train, y_te, fit$pred_test, fit$n_features)
  })) %>% mutate(Repeat = seed - SEED_BASE)
}
run_cv_compare <- function(meta, X_data, label = "", omics_type = "metab", n_rep = N_REP) {
  cat(sprintf("  %-25s n=%3d  ", label, nrow(meta)))
  res <- lapply(seq_len(n_rep), function(r) { if (r%%10==0) cat("."); run_one_rep_compare(meta, X_data, SEED_BASE+r, omics_type) })
  cat(" done\n")
  bind_rows(res) %>% mutate(Label = label)
}

run_one_rep_topk <- function(meta, X_data, seed, k_list, model_type = "ElasticNet", omics_type = "metab") {
  sp   <- split_animals(meta, seed)
  tr   <- sp$train_idx; te <- sp$test_idx
  y_tr <- meta[[COL_PMI]][tr]; y_te <- meta[[COL_PMI]][te]
  prep <- if (omics_type == "metab") preprocess_metab(X_data[tr,]) else preprocess_micro(X_data[tr,])
  Xtr  <- prep$X_train
  Xte  <- if (omics_type == "metab") apply_metab(prep, X_data[te,]) else apply_micro(prep, X_data[te,])
  fit_all <- call_model(model_type, Xtr, y_tr, Xte)
  rows <- list(All = metric_row("All", y_tr, fit_all$pred_train, y_te, fit_all$pred_test, fit_all$n_features))
  for (k in k_list) {
    fit_k <- fit_topk_model(Xtr, y_tr, Xte, k, model_type)
    rows[[paste0("Top", k)]] <- metric_row(paste0("Top",k), y_tr, fit_k$pred_train, y_te, fit_k$pred_test, fit_k$k_actual)
  }
  bind_rows(rows) %>% mutate(Repeat = seed - SEED_BASE)
}
run_cv_topk <- function(meta, X_data, label = "", k_list, model_type = "ElasticNet", omics_type = "metab", n_rep = N_REP) {
  cat(sprintf("  %-25s n=%3d  [%s]  ", label, nrow(meta), model_type))
  res <- lapply(seq_len(n_rep), function(r) { if (r%%10==0) cat("."); run_one_rep_topk(meta, X_data, SEED_BASE+r, k_list, model_type, omics_type) })
  cat(" done\n")
  bind_rows(res) %>% mutate(Label = label)
}

run_one_rep_multiomics <- function(meta, X_metab, X_micro, seed, model_type_m = "ElasticNet", model_type_c = "ElasticNet") {
  sp   <- split_animals(meta, seed)
  tr   <- sp$train_idx; te <- sp$test_idx
  y_tr <- meta[[COL_PMI]][tr]; y_te <- meta[[COL_PMI]][te]
  pm   <- preprocess_metab(X_metab[tr,]); pc <- preprocess_micro(X_micro[tr,])
  Xm_tr <- pm$X_train; Xm_te <- apply_metab(pm, X_metab[te,])
  Xc_tr <- pc$X_train; Xc_te <- apply_micro(pc, X_micro[te,])
  fit_m <- call_model(model_type_m, Xm_tr, y_tr, Xm_te)
  fit_c <- call_model(model_type_c, Xc_tr, y_tr, Xc_te)
  Xmf   <- Xm_tr; colnames(Xmf)    <- paste0("M_", colnames(Xmf))
  Xcf   <- Xc_tr; colnames(Xcf)    <- paste0("C_", colnames(Xcf))
  Xmf_te <- Xm_te; colnames(Xmf_te) <- paste0("M_", colnames(Xmf_te))
  Xcf_te <- Xc_te; colnames(Xcf_te) <- paste0("C_", colnames(Xcf_te))
  fit_ef <- fit_enet(cbind(Xmf, Xcf), y_tr, cbind(Xmf_te, Xcf_te))
  fit_st <- fit_stacking_v2(Xm_tr, Xc_tr, y_tr, Xm_te, Xc_te, model_type_m, model_type_c)
  bind_rows(
    metric_row("Metab_Best",         y_tr, fit_m$pred_train,  y_te, fit_m$pred_test,  fit_m$n_features),
    metric_row("Micro_Best",         y_tr, fit_c$pred_train,  y_te, fit_c$pred_test,  fit_c$n_features),
    metric_row("EF_All(ElasticNet)", y_tr, fit_ef$pred_train, y_te, fit_ef$pred_test, fit_ef$n_features),
    metric_row("Stacking(Best)",     y_tr, fit_st$pred_train, y_te, fit_st$pred_test)
  ) %>% mutate(Repeat = seed - SEED_BASE)
}
run_cv_multiomics <- function(meta, X_metab, X_micro, label = "", model_type_m = "ElasticNet", model_type_c = "ElasticNet", n_rep = N_REP) {
  cat(sprintf("  %-20s n=%3d  metab=%-12s micro=%-12s", label, nrow(meta), model_type_m, model_type_c))
  res <- lapply(seq_len(n_rep), function(r) { if (r%%10==0) cat("."); run_one_rep_multiomics(meta, X_metab, X_micro, SEED_BASE+r, model_type_m, model_type_c) })
  cat(" done\n")
  bind_rows(res) %>% mutate(Label = label)
}

run_one_rep_stacking_topk <- function(meta, X_metab, X_micro, seed, model_type_m, model_type_c, km, kc) {
  sp   <- split_animals(meta, seed)
  tr   <- sp$train_idx; te <- sp$test_idx
  y_tr <- meta[[COL_PMI]][tr]; y_te <- meta[[COL_PMI]][te]
  pm   <- preprocess_metab(X_metab[tr,]); pc <- preprocess_micro(X_micro[tr,])
  Xm_tr <- pm$X_train; Xm_te <- apply_metab(pm, X_metab[te,])
  Xc_tr <- pc$X_train; Xc_te <- apply_micro(pc, X_micro[te,])
  fit_stk <- fit_stacking_topk(Xm_tr, Xc_tr, y_tr, Xm_te, Xc_te, model_type_m, model_type_c, km, kc)
  metric_row("Stacking_TopK", y_tr, fit_stk$pred_train, y_te, fit_stk$pred_test, as.integer(km+kc)) %>%
    mutate(Repeat = seed - SEED_BASE)
}
run_cv_stacking_topk <- function(meta, X_metab, X_micro, label, model_type_m, model_type_c, km, kc, n_rep = N_REP) {
  cat(sprintf("  %-25s n=%3d  [Top%dm+%dc]  ", label, nrow(meta), km, kc))
  res <- lapply(seq_len(n_rep), function(r) { if (r%%10==0) cat("."); run_one_rep_stacking_topk(meta, X_metab, X_micro, SEED_BASE+r, model_type_m, model_type_c, km, kc) })
  cat(" done\n")
  bind_rows(res) %>% mutate(Label = label)
}

# ── Load data ──────────────────────────────────────────────────────────────────
cat("Select: [1/3] metabolomics.csv\n")
metab_raw <- read.csv(file.choose(), check.names = FALSE)
cat("Select: [2/3] microbiome.csv\n")
micro_raw <- read.csv(file.choose(), check.names = FALSE)
cat("Select: [3/3] metabolite_mapping.csv\n")
metab_map_raw <- read.csv(file.choose(), check.names = FALSE)

metab_raw[[COL_PMI]] <- as.numeric(metab_raw[[COL_PMI]])
micro_raw[[COL_PMI]] <- as.numeric(micro_raw[[COL_PMI]])
meta_M_all <- metab_raw %>% select(all_of(META_COLS))
feat_M_all <- metab_raw %>% select(-all_of(META_COLS)) %>% select(where(is.numeric))
meta_C_all <- micro_raw %>% select(all_of(META_COLS))
feat_C_all <- micro_raw %>% select(-all_of(META_COLS)) %>% select(where(is.numeric))
metab_map  <- setNames(as.character(metab_map_raw[[2]]), as.character(metab_map_raw[[1]]))
map_metab_names <- function(ids) { m <- metab_map[ids]; ifelse(!is.na(m) & nchar(trimws(m)) > 0, m, ids) }

joint_ids <- intersect(meta_M_all[[COL_ID]], meta_C_all[[COL_ID]])
meta_J    <- meta_M_all %>% filter(.data[[COL_ID]] %in% joint_ids) %>% arrange(.data[[COL_ID]])
feat_MJ   <- feat_M_all[match(meta_J[[COL_ID]], meta_M_all[[COL_ID]]), ]
feat_CJ   <- feat_C_all[match(meta_J[[COL_ID]], meta_C_all[[COL_ID]]), ]

# ── Step 1: Model selection ────────────────────────────────────────────────────
cat("\n---- Step 1A: Metabolomics ElasticNet vs RF ----\n")
res_Acmp <- bind_rows(lapply(setNames(TISSUES, TISSUES), function(tis) {
  idx <- meta_J[[COL_TISSUE]] == tis
  run_cv_compare(meta_J[idx,], feat_MJ[idx,], tis, "metab")
}))
summ_Acmp <- res_Acmp %>% group_by(Label, Model) %>%
  summarise(RMSE_mean = round(mean(RMSE_test, na.rm=TRUE), 4),
            RMSE_sd   = round(sd(RMSE_test,   na.rm=TRUE), 4),
            R2_mean   = round(mean(R2_test,   na.rm=TRUE), 4),
            R2_sd     = round(sd(R2_test,     na.rm=TRUE), 4), .groups="drop") %>%
  mutate(RMSE_str = sprintf("%.3f \u00b1 %.3f", RMSE_mean, RMSE_sd),
         R2_str   = sprintf("%.3f \u00b1 %.3f", R2_mean,   R2_sd)) %>%
  arrange(Label, RMSE_mean)
write.csv(res_Acmp,  file.path(output_dir, "S1_Acmp_metab_2models_all.csv"),  row.names=FALSE)
write.csv(summ_Acmp, file.path(output_dir, "S1_Acmp_metab_2models_summ.csv"), row.names=FALSE)

cat("\n---- Step 1B: Microbiome ElasticNet vs RF ----\n")
res_Bcmp <- bind_rows(lapply(setNames(TISSUES, TISSUES), function(tis) {
  idx <- meta_J[[COL_TISSUE]] == tis
  run_cv_compare(meta_J[idx,], feat_CJ[idx,], tis, "micro")
}))
summ_Bcmp <- res_Bcmp %>% group_by(Label, Model) %>%
  summarise(RMSE_mean = round(mean(RMSE_test, na.rm=TRUE), 4),
            RMSE_sd   = round(sd(RMSE_test,   na.rm=TRUE), 4),
            R2_mean   = round(mean(R2_test,   na.rm=TRUE), 4), .groups="drop") %>%
  mutate(RMSE_str = sprintf("%.3f \u00b1 %.3f", RMSE_mean, RMSE_sd)) %>%
  arrange(Label, RMSE_mean)
write.csv(res_Bcmp,  file.path(output_dir, "S2_Bcmp_micro_2models_all.csv"),  row.names=FALSE)
write.csv(summ_Bcmp, file.path(output_dir, "S2_Bcmp_micro_2models_summ.csv"), row.names=FALSE)

BEST_M <- summ_Acmp %>% group_by(Label) %>% slice_min(RMSE_mean, n=1, with_ties=FALSE) %>% select(Label, Model) %>% deframe()
BEST_C <- summ_Bcmp %>% group_by(Label) %>% slice_min(RMSE_mean, n=1, with_ties=FALSE) %>% select(Label, Model) %>% deframe()
modal_m <- names(sort(table(BEST_M), decreasing=TRUE))[1]
modal_c <- names(sort(table(BEST_C), decreasing=TRUE))[1]
write.csv(data.frame(Tissue=TISSUES, Metab_Best=BEST_M[TISSUES], Micro_Best=BEST_C[TISSUES]),
          file.path(output_dir, "S3_ModelSelection.csv"), row.names=FALSE)

cat("\n---- Step 1C: Full-feature four-strategy (-> Supplementary S4) ----\n")
res_Cbest <- bind_rows(lapply(setNames(TISSUES, TISSUES), function(tis) {
  idx <- meta_J[[COL_TISSUE]] == tis
  cat(sprintf("  %-10s [m=%-12s c=%-12s]", tis, BEST_M[tis], BEST_C[tis]))
  res <- lapply(seq_len(N_REP), function(r) { if (r%%10==0) cat(".")
    run_one_rep_multiomics(meta_J[idx,], feat_MJ[idx,], feat_CJ[idx,], SEED_BASE+r, BEST_M[tis], BEST_C[tis]) })
  cat(" done\n")
  bind_rows(res) %>% mutate(Label = tis)
}))
summ_Cbest <- res_Cbest %>% group_by(Label, Model) %>%
  summarise(RMSE_mean=round(mean(RMSE_test,na.rm=TRUE),4), RMSE_sd=round(sd(RMSE_test,na.rm=TRUE),4),
            R2_mean=round(mean(R2_test,na.rm=TRUE),4), R2_sd=round(sd(R2_test,na.rm=TRUE),4), .groups="drop") %>%
  mutate(RMSE_str=sprintf("%.3f \u00b1 %.3f",RMSE_mean,RMSE_sd), R2_str=sprintf("%.3f \u00b1 %.3f",R2_mean,R2_sd)) %>%
  arrange(Label, RMSE_mean)
write.csv(res_Cbest,  file.path(output_dir, "S4_Cbest_fullfeature_all.csv"),  row.names=FALSE)
write.csv(summ_Cbest, file.path(output_dir, "S4_Cbest_fullfeature_summ.csv"), row.names=FALSE)

# ── Step 2: Top-K determination ────────────────────────────────────────────────
cat(sprintf("\n---- Step 2A: Metabolomics Top-K (%s) ----\n", paste(TOP_K_METAB_LIST, collapse="/")))
res_TK_metab <- bind_rows(lapply(setNames(TISSUES,TISSUES), function(tis) {
  idx <- meta_J[[COL_TISSUE]] == tis
  run_cv_topk(meta_J[idx,], feat_MJ[idx,], tis, TOP_K_METAB_LIST, BEST_M[tis], "metab")
}))
summ_TK_metab <- res_TK_metab %>% group_by(Label, Model) %>%
  summarise(RMSE_mean=round(mean(RMSE_test,na.rm=TRUE),4), RMSE_sd=round(sd(RMSE_test,na.rm=TRUE),4),
            R2_mean=round(mean(R2_test,na.rm=TRUE),4), .groups="drop") %>%
  mutate(RMSE_str=sprintf("%.3f \u00b1 %.3f",RMSE_mean,RMSE_sd)) %>% arrange(Label, RMSE_mean)
write.csv(res_TK_metab,  file.path(output_dir, "S5_TopK_metab_all.csv"),  row.names=FALSE)
write.csv(summ_TK_metab, file.path(output_dir, "S5_TopK_metab_summ.csv"), row.names=FALSE)

cat(sprintf("\n---- Step 2B: Microbiome Top-K (%s) ----\n", paste(TOP_K_MICRO_LIST, collapse="/")))
res_TK_micro <- bind_rows(lapply(setNames(TISSUES,TISSUES), function(tis) {
  idx <- meta_J[[COL_TISSUE]] == tis
  run_cv_topk(meta_J[idx,], feat_CJ[idx,], tis, TOP_K_MICRO_LIST, BEST_C[tis], "micro")
}))
summ_TK_micro <- res_TK_micro %>% group_by(Label, Model) %>%
  summarise(RMSE_mean=round(mean(RMSE_test,na.rm=TRUE),4), RMSE_sd=round(sd(RMSE_test,na.rm=TRUE),4),
            R2_mean=round(mean(R2_test,na.rm=TRUE),4), .groups="drop") %>%
  mutate(RMSE_str=sprintf("%.3f \u00b1 %.3f",RMSE_mean,RMSE_sd)) %>% arrange(Label, RMSE_mean)
write.csv(res_TK_micro,  file.path(output_dir, "S5_TopK_micro_all.csv"),  row.names=FALSE)
write.csv(summ_TK_micro, file.path(output_dir, "S5_TopK_micro_summ.csv"), row.names=FALSE)

KM_SELECTED <- summ_TK_metab %>% filter(Model!="All") %>% group_by(Model) %>%
  summarise(avg=mean(RMSE_mean,na.rm=TRUE),.groups="drop") %>% slice_min(avg,n=1,with_ties=FALSE) %>%
  mutate(K=as.integer(str_extract(Model,"\\d+"))) %>% pull(K)
KC_SELECTED <- summ_TK_micro %>% filter(Model!="All") %>% group_by(Model) %>%
  summarise(avg=mean(RMSE_mean,na.rm=TRUE),.groups="drop") %>% slice_min(avg,n=1,with_ties=FALSE) %>%
  mutate(K=as.integer(str_extract(Model,"\\d+"))) %>% pull(K)
cat(sprintf("  Consensus K (metabolomics) = %d,  Consensus K (microbiome) = %d\n", KM_SELECTED, KC_SELECTED))
# Manual override: KM_SELECTED <- 20L; KC_SELECTED <- 10L

ef_topk_label  <- sprintf("EF_Top%dm+%dc",      KM_SELECTED, KC_SELECTED)
stk_topk_label <- sprintf("Stacking_Top%dm+%dc", KM_SELECTED, KC_SELECTED)

# ── Step 3: Parsimonious four-arm comparison ───────────────────────────────────
cat(sprintf("\n---- Step 3: Four-arm comparison (Top%dm + Top%dc) ----\n", KM_SELECTED, KC_SELECTED))
res_EF_topk <- list(); res_ST_topk <- list(); feat_ST_topk_log <- list()

for (tis in TISSUES) {
  idx <- meta_J[[COL_TISSUE]] == tis
  mt_m <- BEST_M[tis]; mt_c <- BEST_C[tis]
  cat(sprintf("  %-10s [m=%-12s c=%-12s]", tis, mt_m, mt_c))
  ef_tmp <- list(); stk_tmp <- list(); feat_tmp <- list()
  for (r in seq_len(N_REP)) {
    if (r%%10==0) cat(".")
    sp <- split_animals(meta_J[idx,], SEED_BASE+r)
    tr <- sp$train_idx; te <- sp$test_idx
    y_tr <- meta_J[[COL_PMI]][idx][tr]; y_te <- meta_J[[COL_PMI]][idx][te]
    pm   <- preprocess_metab(feat_MJ[idx,][tr,]); pc <- preprocess_micro(feat_CJ[idx,][tr,])
    Xm_tr <- pm$X_train; Xm_te <- apply_metab(pm, feat_MJ[idx,][te,])
    Xc_tr <- pc$X_train; Xc_te <- apply_micro(pc, feat_CJ[idx,][te,])
    fit_stk  <- fit_stacking_topk(Xm_tr, Xc_tr, y_tr, Xm_te, Xc_te, mt_m, mt_c, KM_SELECTED, KC_SELECTED)
    stk_tmp[[r]] <- metric_row(stk_topk_label, y_tr, fit_stk$pred_train, y_te, fit_stk$pred_test,
                               as.integer(KM_SELECTED+KC_SELECTED)) %>% mutate(Repeat=r)
    fit_ef <- fit_ef_topk(Xm_tr, Xc_tr, y_tr, Xm_te, Xc_te, mt_m, mt_c, KM_SELECTED, KC_SELECTED,
                          fit_m_full=fit_stk$fit_m_full, fit_c_full=fit_stk$fit_c_full)
    ef_tmp[[r]]   <- metric_row(ef_topk_label, y_tr, fit_ef$pred_train, y_te, fit_ef$pred_test,
                                fit_ef$n_features) %>% mutate(Repeat=r)
    feat_tmp[[r]] <- bind_rows(
      data.frame(Feature=fit_stk$top_m, Omics="Metabolomics", Rank=seq_along(fit_stk$top_m), Repeat=r, Tissue=tis),
      data.frame(Feature=fit_stk$top_c, Omics="Microbiome",   Rank=seq_along(fit_stk$top_c), Repeat=r, Tissue=tis))
  }
  cat(" done\n")
  res_EF_topk[[tis]]      <- bind_rows(ef_tmp)  %>% mutate(Label=tis)
  res_ST_topk[[tis]]      <- bind_rows(stk_tmp) %>% mutate(Label=tis)
  feat_ST_topk_log[[tis]] <- bind_rows(feat_tmp)
}
res_EF_topk      <- bind_rows(res_EF_topk)
res_ST_topk      <- bind_rows(res_ST_topk)
feat_ST_topk_log <- bind_rows(feat_ST_topk_log)

topk_feat_summary <- feat_ST_topk_log %>%
  group_by(Tissue, Omics, Feature) %>%
  summarise(select_freq=n()/N_REP, mean_rank=mean(Rank), .groups="drop") %>%
  mutate(Display_Name=case_when(Omics=="Metabolomics"~map_metab_names(Feature), TRUE~Feature)) %>%
  arrange(Tissue, Omics, mean_rank)

write.csv(res_EF_topk,       file.path(output_dir, "Step3_EF_TopK_all.csv"),         row.names=FALSE)
write.csv(res_ST_topk,       file.path(output_dir, "Step3_Stacking_TopK_all.csv"),   row.names=FALSE)
write.csv(feat_ST_topk_log,  file.path(output_dir, "Step3_feat_selection_log.csv"),  row.names=FALSE)
write.csv(topk_feat_summary, file.path(output_dir, "Step3_feat_selection_summ.csv"), row.names=FALSE)

summ_EF_topk <- res_EF_topk %>% group_by(Label, Model) %>%
  summarise(RMSE_mean=round(mean(RMSE_test,na.rm=TRUE),4), RMSE_sd=round(sd(RMSE_test,na.rm=TRUE),4),
            R2_mean=round(mean(R2_test,na.rm=TRUE),4), .groups="drop") %>%
  mutate(RMSE_str=sprintf("%.3f \u00b1 %.3f",RMSE_mean,RMSE_sd))
summ_ST_topk <- res_ST_topk %>% group_by(Label, Model) %>%
  summarise(RMSE_mean=round(mean(RMSE_test,na.rm=TRUE),4), RMSE_sd=round(sd(RMSE_test,na.rm=TRUE),4),
            R2_mean=round(mean(R2_test,na.rm=TRUE),4), .groups="drop") %>%
  mutate(RMSE_str=sprintf("%.3f \u00b1 %.3f",RMSE_mean,RMSE_sd))

compute_wilcox_pair <- function(res_a, model_a, res_b, model_b) {
  da <- res_a %>% filter(Model==model_a, !is.na(RMSE_test)) %>% select(Label, Repeat, RMSE_test) %>% rename(va=RMSE_test)
  db <- res_b %>% filter(Model==model_b, !is.na(RMSE_test)) %>% select(Label, Repeat, RMSE_test) %>% rename(vb=RMSE_test)
  inner_join(da, db, by=c("Label","Repeat")) %>% group_by(Label) %>%
    summarise(delta=round(mean(va-vb,na.rm=TRUE),3),
              p_val=wilcox.test(va, vb, paired=TRUE, exact=FALSE)$p.value, .groups="drop") %>%
    mutate(sig=case_when(p_val<0.001~"***",p_val<0.01~"**",p_val<0.05~"*",TRUE~"n.s."),
           pair=sprintf("%s vs %s", model_a, model_b))
}

metab_topk_label <- paste0("Top", KM_SELECTED); micro_topk_label <- paste0("Top", KC_SELECTED)
res_metab_tk <- res_TK_metab %>% filter(Model==metab_topk_label) %>% mutate(Model="Metab_TopK")
res_micro_tk <- res_TK_micro %>% filter(Model==micro_topk_label) %>% mutate(Model="Micro_TopK")
res_ef_tk    <- res_EF_topk  %>% mutate(Model="EF_TopK")
res_stk_tk   <- res_ST_topk  %>% mutate(Model="Stacking_TopK")

wlx_metab_vs_ef  <- compute_wilcox_pair(res_metab_tk,"Metab_TopK", res_ef_tk, "EF_TopK")
wlx_metab_vs_stk <- compute_wilcox_pair(res_metab_tk,"Metab_TopK", res_stk_tk,"Stacking_TopK")
wlx_ef_vs_stk    <- compute_wilcox_pair(res_ef_tk,   "EF_TopK",    res_stk_tk,"Stacking_TopK")
write.csv(wlx_metab_vs_ef,  file.path(output_dir,"Wilcoxon_Metab_vs_EF.csv"),       row.names=FALSE)
write.csv(wlx_metab_vs_stk, file.path(output_dir,"Wilcoxon_Metab_vs_Stacking.csv"), row.names=FALSE)
write.csv(wlx_ef_vs_stk,    file.path(output_dir,"Wilcoxon_EF_vs_Stacking.csv"),    row.names=FALSE)

# ── Step 3.5: Feature compression effect ──────────────────────────────────────
stage_full <- res_Cbest %>% filter(!is.na(RMSE_test)) %>%
  mutate(Stage="Full", ModelBase=case_when(
    Model=="Metab_Best"~"Metab", Model=="Micro_Best"~"Micro",
    Model=="EF_All(ElasticNet)"~"EF", Model=="Stacking(Best)"~"Stacking")) %>%
  filter(!is.na(ModelBase))
stage_topk <- bind_rows(
  res_metab_tk%>%mutate(ModelBase="Metab"), res_micro_tk%>%mutate(ModelBase="Micro"),
  res_EF_topk %>%mutate(ModelBase="EF"),    res_ST_topk %>%mutate(ModelBase="Stacking")) %>%
  filter(!is.na(RMSE_test)) %>% mutate(Stage="TopK")

compress_4arm_all <- bind_rows(
  stage_full%>%select(Label,ModelBase,Stage,Repeat,RMSE_test),
  stage_topk%>%select(Label,ModelBase,Stage,Repeat,RMSE_test)) %>%
  mutate(Label=factor(Label,levels=TISSUES),
         ModelBase=factor(ModelBase,levels=c("Metab","Micro","EF","Stacking")),
         Stage=factor(Stage,levels=c("Full","TopK")))

compress_4arm_wilcox <- compress_4arm_all %>%
  select(Label,ModelBase,Stage,Repeat,RMSE_test) %>%
  pivot_wider(names_from=Stage, values_from=RMSE_test) %>%
  group_by(Label,ModelBase) %>%
  summarise(delta=round(mean(TopK-Full,na.rm=TRUE),3),
            p_val=wilcox.test(TopK, Full, paired=TRUE, exact=FALSE)$p.value, .groups="drop") %>%
  mutate(sig=case_when(p_val<0.001~"***",p_val<0.01~"**",p_val<0.05~"*",TRUE~"n.s."))
write.csv(compress_4arm_wilcox, file.path(output_dir,"Step3.5_compress_4arm_wilcox.csv"), row.names=FALSE)

compress_4arm_summ_se <- compress_4arm_all %>%
  group_by(Label,ModelBase,Stage) %>%
  summarise(mean_rmse=mean(RMSE_test,na.rm=TRUE),
            se_rmse=sd(RMSE_test,na.rm=TRUE)/sqrt(sum(!is.na(RMSE_test))), .groups="drop") %>%
  left_join(compress_4arm_wilcox, by=c("Label","ModelBase")) %>%
  mutate(Label=factor(Label,levels=TISSUES),
         ModelBase=factor(ModelBase,levels=c("Metab","Micro","EF","Stacking")),
         Stage=factor(Stage,levels=c("Full","TopK")),
         line_color =MODEL_COLORS_BASE[as.character(ModelBase)],
         line_type  =ifelse(sig=="n.s.","dashed","solid"),
         line_width =ifelse(sig=="n.s.",0.55,1.0),
         label_color=case_when(delta<0&sig!="n.s."~"#1A6B4A",delta>0&sig!="n.s."~"#943126",TRUE~"#999999"),
         nudge_y=case_when(ModelBase=="Micro"~0.15,ModelBase=="Metab"~0.08,ModelBase=="EF"~0.03,TRUE~0.00))
label_4arm_se <- compress_4arm_summ_se %>% filter(Stage=="TopK") %>%
  mutate(label_y=mean_rmse+se_rmse+nudge_y, label_text=paste0(sig,"  \u0394=",sprintf("%.3f",delta)))
y_lower <- min(compress_4arm_summ_se$mean_rmse-compress_4arm_summ_se$se_rmse,na.rm=TRUE)*0.92
y_upper <- max(compress_4arm_summ_se$mean_rmse+compress_4arm_summ_se$se_rmse,na.rm=TRUE)

p_compress4 <- ggplot() +
  geom_line(data=compress_4arm_summ_se, aes(x=Stage,y=mean_rmse,group=ModelBase,color=line_color,linetype=line_type,linewidth=line_width),alpha=0.75) +
  scale_color_identity() + scale_linetype_identity() + scale_linewidth_identity() +
  geom_errorbar(data=compress_4arm_summ_se, aes(x=Stage,ymin=mean_rmse-se_rmse,ymax=mean_rmse+se_rmse,group=ModelBase), width=0.07,linewidth=0.4,color="grey60",alpha=0.70) +
  geom_point(data=compress_4arm_summ_se%>%filter(Stage=="Full"), aes(x=Stage,y=mean_rmse,fill=ModelBase), shape=21,size=3.8,color="white",stroke=0.6,alpha=0.45) +
  geom_point(data=compress_4arm_summ_se%>%filter(Stage=="TopK"), aes(x=Stage,y=mean_rmse,fill=ModelBase), shape=23,size=4.8,color="white",stroke=0.9) +
  scale_fill_manual(values=MODEL_COLORS_BASE,name="Strategy") +
  geom_text(data=label_4arm_se, aes(x=Stage,y=label_y,label=label_text,color=label_color), size=2.6,lineheight=0.80,hjust=0.5,fontface="bold",family="Arial") +
  facet_wrap(~Label,nrow=1) +
  scale_x_discrete(labels=c("Full"="Full features\n(All)","TopK"=sprintf("Top-K features\n(Top%d+%d)",KM_SELECTED,KC_SELECTED))) +
  coord_cartesian(ylim=c(y_lower,y_upper*1.30)) +
  scale_y_continuous(breaks=seq(0,2,by=0.2),expand=expansion(mult=c(0.02,0.05))) +
  labs(title=sprintf("Feature compression effect: Full vs Top-K (Top%d+%d)",KM_SELECTED,KC_SELECTED),
       subtitle="Diamond=Top-K  |  Circle=Full  |  Error bar=\u00b11 SE  |  Solid=significant (p<0.05)  |  Dashed=n.s.\n\u0394=Top-K RMSE\u2212Full RMSE (negative=improved)",
       x=NULL,y="Mean Test RMSE (days)") +
  theme_bw(base_size=12,base_family="Arial") +
  theme(legend.position="bottom",legend.direction="horizontal",strip.text=element_text(face="bold",size=11),
        strip.background=element_rect(fill="grey96"),panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),panel.spacing=unit(0.8,"lines"),
        plot.title=element_text(size=12,face="bold"),plot.subtitle=element_text(size=8,color="grey35",lineheight=1.3),
        axis.text.x=element_text(size=8.5,lineheight=0.85))
ggsave(file.path(output_dir,"SuppFig_compression_4arm.tiff"), p_compress4, width=16,height=6,device="tiff",dpi=300,compression="lzw")

# ── Step 4: Multi-tissue pairwise ──────────────────────────────────────────────
cat(sprintf("\n---- Step 4A: Pooled five-tissue Stacking_TopK ----\n"))
res_D  <- run_cv_stacking_topk(meta_J,feat_MJ,feat_CJ,"Pooled_5Tissues",modal_m,modal_c,KM_SELECTED,KC_SELECTED)
summ_D <- res_D %>% group_by(Label,Model) %>%
  summarise(RMSE_mean=round(mean(RMSE_test,na.rm=TRUE),4), RMSE_sd=round(sd(RMSE_test,na.rm=TRUE),4),
            R2_mean=round(mean(R2_test,na.rm=TRUE),4), R2_sd=round(sd(R2_test,na.rm=TRUE),4),
            MAE_mean=round(mean(MAE_test,na.rm=TRUE),4), .groups="drop") %>%
  mutate(RMSE_str=sprintf("%.3f \u00b1 %.3f",RMSE_mean,RMSE_sd), R2_str=sprintf("%.3f \u00b1 %.3f",R2_mean,R2_sd))
write.csv(res_D,  file.path(output_dir,"Step4_D_pooled_all.csv"),  row.names=FALSE)
write.csv(summ_D, file.path(output_dir,"Step4_D_pooled_summ.csv"), row.names=FALSE)

cat(sprintf("\n---- Step 4B: Two-tissue pairwise (10 pairs) ----\n"))
tis_pairs_all <- combn(TISSUES, 2, simplify=FALSE)
pair_names    <- sapply(tis_pairs_all, paste, collapse="+")
res_E_list    <- list()
for (pi in seq_along(tis_pairs_all)) {
  pair  <- tis_pairs_all[[pi]]; pname <- pair_names[pi]
  idx_p <- meta_J[[COL_TISSUE]] %in% pair
  mt_m_pair <- names(sort(table(BEST_M[pair]),decreasing=TRUE))[1]
  mt_c_pair <- names(sort(table(BEST_C[pair]),decreasing=TRUE))[1]
  cat(sprintf("  [%2d/10] %s [m=%s c=%s]", pi, pname, mt_m_pair, mt_c_pair))
  res_E_list[[pname]] <- run_cv_stacking_topk(meta_J[idx_p,],feat_MJ[idx_p,],feat_CJ[idx_p,],
                                              pname,mt_m_pair,mt_c_pair,KM_SELECTED,KC_SELECTED)
}
res_E <- bind_rows(res_E_list)
summ_E <- res_E %>% group_by(Label,Model) %>%
  summarise(RMSE_mean=round(mean(RMSE_test,na.rm=TRUE),4), RMSE_sd=round(sd(RMSE_test,na.rm=TRUE),4),
            R2_mean=round(mean(R2_test,na.rm=TRUE),4), .groups="drop") %>%
  mutate(RMSE_str=sprintf("%.3f \u00b1 %.3f",RMSE_mean,RMSE_sd)) %>% arrange(RMSE_mean)
write.csv(res_E,  file.path(output_dir,"Step4_E_pairwise_all.csv"),  row.names=FALSE)
write.csv(summ_E, file.path(output_dir,"Step4_E_pairwise_summ.csv"), row.names=FALSE)

# ── Step 5: Visualization ──────────────────────────────────────────────────────
wilcox_brackets <- function(res_df, pairs, x_order, facet_col="Label", y_col="RMSE_test") {
  bind_rows(lapply(seq_along(pairs), function(i) {
    ma <- pairs[[i]][1]; mb <- pairs[[i]][2]
    xa <- which(x_order==ma); xb <- which(x_order==mb)
    da <- res_df %>% filter(Model==ma,!is.na(.data[[y_col]])) %>% select(all_of(c(facet_col,"Repeat",y_col))) %>% rename(va=all_of(y_col))
    db <- res_df %>% filter(Model==mb,!is.na(.data[[y_col]])) %>% select(all_of(c(facet_col,"Repeat",y_col))) %>% rename(vb=all_of(y_col))
    inner_join(da,db,by=c(facet_col,"Repeat")) %>% group_by(.data[[facet_col]]) %>%
      summarise(p_val=wilcox.test(va,vb,paired=TRUE,exact=FALSE)$p.value, y_top=max(c(va,vb),na.rm=TRUE), .groups="drop") %>%
      mutate(sig=case_when(p_val<0.001~"***",p_val<0.01~"**",p_val<0.05~"*",TRUE~"n.s."),
             xa=xa,xb=xb,xmid=(xa+xb)/2,comp_i=i, y_seg=y_top*1.08+(i-1)*0.22, y_txt=y_top*1.08+(i-1)*0.22+0.09)
  }))
}

make_violin_wilcox <- function(res_df, x_order, colors, title, y_label, facet_col="Label", y_col="RMSE_test", compare_pairs=NULL, base_size=12) {
  df_p <- res_df %>% filter(!is.na(.data[[y_col]])) %>%
    mutate(!!facet_col:=factor(.data[[facet_col]],levels=TISSUES), Model=factor(Model,levels=x_order))
  p <- ggplot(df_p,aes(x=Model,y=.data[[y_col]],fill=Model)) +
    geom_violin(trim=FALSE,alpha=0.82,color="white",linewidth=0.3) +
    geom_boxplot(width=0.12,fill="white",color="grey40",outlier.size=0.4,linewidth=0.4) +
    stat_summary(fun=mean,geom="point",shape=18,size=2.5,color="black") +
    scale_fill_manual(values=colors,na.value="grey80") +
    facet_wrap(as.formula(paste("~",facet_col)),nrow=1) +
    labs(title=title,x=NULL,y=y_label) +
    theme_bw(base_size=base_size,base_family="Arial") +
    theme(axis.text.x=element_text(angle=35,hjust=1,size=base_size-3,family="Arial"),
          axis.text.y=element_text(family="Arial"), axis.title=element_text(family="Arial"),
          legend.position="none", strip.text=element_text(face="bold",size=base_size,family="Arial"),
          panel.grid.major.x=element_blank(), plot.title=element_text(size=base_size,face="bold",family="Arial")) +
    scale_y_continuous(expand=expansion(mult=c(0.02,0.32)))
  brk <- wilcox_brackets(res_df, compare_pairs, x_order, facet_col, y_col)
  p + geom_segment(data=brk,aes(x=xa,xend=xb,y=y_seg,yend=y_seg),inherit.aes=FALSE,color="grey20",linewidth=0.55) +
    geom_segment(data=brk,aes(x=xa,xend=xa,y=y_seg-0.07,yend=y_seg),inherit.aes=FALSE,color="grey20",linewidth=0.55) +
    geom_segment(data=brk,aes(x=xb,xend=xb,y=y_seg-0.07,yend=y_seg),inherit.aes=FALSE,color="grey20",linewidth=0.55) +
    geom_text(data=brk,aes(x=xmid,y=y_txt,label=sig),inherit.aes=FALSE,size=base_size/2.5,fontface="bold",color="grey10",family="Arial")
}

make_violin_simple <- function(df, x_order, colors, title, y_label, facet_col="Label", base_size=12) {
  df %>% filter(!is.na(RMSE_test)) %>%
    mutate(!!facet_col:=factor(.data[[facet_col]],levels=TISSUES), Model=factor(Model,levels=x_order)) %>%
    ggplot(aes(x=Model,y=RMSE_test,fill=Model)) +
    geom_violin(trim=FALSE,alpha=0.82,color="white",linewidth=0.3) +
    geom_boxplot(width=0.12,fill="white",color="grey40",outlier.size=0.4,linewidth=0.4) +
    stat_summary(fun=mean,geom="point",shape=18,size=2.5,color="black") +
    scale_fill_manual(values=colors,na.value="grey80") +
    facet_wrap(as.formula(paste("~",facet_col)),nrow=1) +
    labs(title=title,x=NULL,y=y_label) +
    theme_bw(base_size=base_size,base_family="Arial") +
    theme(axis.text.x=element_text(angle=35,hjust=1,size=base_size-3,family="Arial"),
          legend.position="none", strip.text=element_text(face="bold",family="Arial"),
          panel.grid.major.x=element_blank(), plot.title=element_text(size=base_size,face="bold",family="Arial"))
}

# Supplementary figures
ggsave(file.path(output_dir,"SuppFig1_Acmp_metab_2models.tiff"),
       make_violin_simple(res_Acmp,c("ElasticNet","RandomForest"),MODEL_COLORS_2,"Supplementary: Metabolomics model selection","Test RMSE (days)"),
       width=14,height=5,device="tiff",dpi=300,compression="lzw")

ggsave(file.path(output_dir,"SuppFig2_Bcmp_micro_2models.tiff"),
       make_violin_simple(res_Bcmp,c("ElasticNet","RandomForest"),MODEL_COLORS_2,"Supplementary: Microbiome model selection","Test RMSE (days)"),
       width=14,height=5,device="tiff",dpi=300,compression="lzw")

p_s3 <- bind_rows(data.frame(Tissue=TISSUES,Model=BEST_M[TISSUES],Omics="Metabolomics"),
                  data.frame(Tissue=TISSUES,Model=BEST_C[TISSUES],Omics="Microbiome")) %>%
  mutate(Tissue=factor(Tissue,levels=rev(TISSUES))) %>%
  ggplot(aes(x=Omics,y=Tissue,fill=Model)) +
  geom_tile(color="white",linewidth=1.2) + geom_text(aes(label=Model),size=4,fontface="bold",family="Arial") +
  scale_fill_manual(values=c("ElasticNet"="#6BAED6","RandomForest"="#FD8D3C")) +
  labs(title="Supplementary: Per-tissue optimal model selection",x=NULL,y=NULL) +
  theme_bw(base_size=13,base_family="Arial") + theme(panel.grid=element_blank(),text=element_text(family="Arial"))
ggsave(file.path(output_dir,"SuppFig3_ModelSelection_tile.tiff"),p_s3,width=7,height=5,device="tiff",dpi=300,compression="lzw")

ggsave(file.path(output_dir,"SuppFig4_Cbest_fullfeature_wilcox.tiff"),
       make_violin_wilcox(res_Cbest,c("Metab_Best","Micro_Best","EF_All(ElasticNet)","Stacking(Best)"),MODEL_COLORS_MULTI,
                          "Supplementary: Full-feature four-strategy comparison","Test RMSE (days)",
                          compare_pairs=list(c("Metab_Best","EF_All(ElasticNet)"),c("Metab_Best","Stacking(Best)"))),
       width=14,height=6,device="tiff",dpi=300,compression="lzw")

parse_k <- function(m) case_when(grepl("^All",m)~Inf, TRUE~as.numeric(str_extract(m,"\\d+(?=[^\\d]*$)")))
p_s5m <- summ_TK_metab %>% mutate(K_val=parse_k(Model),Label=factor(Label,levels=TISSUES)) %>% filter(is.finite(K_val)) %>%
  ggplot(aes(x=K_val,y=RMSE_mean,color=Label,group=Label)) + geom_point(size=2.5) + geom_line(linewidth=0.8) +
  geom_hline(data=summ_TK_metab%>%mutate(K_val=parse_k(Model))%>%filter(!is.finite(K_val))%>%mutate(Label=factor(Label,levels=TISSUES)),
             aes(yintercept=RMSE_mean,color=Label),linetype="dashed",linewidth=0.7,alpha=0.6) +
  geom_vline(xintercept=KM_SELECTED,linetype="dotted",color="grey20",linewidth=0.8) +
  annotate("text",x=KM_SELECTED,y=Inf,label=sprintf("K=%d",KM_SELECTED),vjust=1.5,hjust=-0.1,size=3.5,color="grey20",family="Arial") +
  scale_x_continuous(breaks=sort(TOP_K_METAB_LIST)) +
  labs(title="Metabolomics Top-K",subtitle="Dashed=All baseline, dotted=consensus cutoff",x="K (features)",y="Mean RMSE (days)") +
  theme_bw(base_size=12,base_family="Arial") + theme(text=element_text(family="Arial"))
p_s5c <- summ_TK_micro %>% mutate(K_val=parse_k(Model),Label=factor(Label,levels=TISSUES)) %>% filter(is.finite(K_val)) %>%
  ggplot(aes(x=K_val,y=RMSE_mean,color=Label,group=Label)) + geom_point(size=2.5) + geom_line(linewidth=0.8) +
  geom_hline(data=summ_TK_micro%>%mutate(K_val=parse_k(Model))%>%filter(!is.finite(K_val))%>%mutate(Label=factor(Label,levels=TISSUES)),
             aes(yintercept=RMSE_mean,color=Label),linetype="dashed",linewidth=0.7,alpha=0.6) +
  geom_vline(xintercept=KC_SELECTED,linetype="dotted",color="grey20",linewidth=0.8) +
  annotate("text",x=KC_SELECTED,y=Inf,label=sprintf("K=%d",KC_SELECTED),vjust=1.5,hjust=-0.1,size=3.5,color="grey20",family="Arial") +
  scale_x_continuous(breaks=sort(TOP_K_MICRO_LIST)) +
  labs(title="Microbiome Top-K",subtitle="Dashed=All baseline, dotted=consensus cutoff",x="K (features)",y="Mean RMSE (days)") +
  theme_bw(base_size=12,base_family="Arial") + theme(text=element_text(family="Arial"))
ggsave(file.path(output_dir,"SuppFig5_TopK_lineplot.tiff"),p_s5m+p_s5c+plot_layout(guides="collect"),width=14,height=5,device="tiff",dpi=300,compression="lzw")

# Figure 5(a)
x_labels_main <- setNames(c(sprintf("Metab\nTop%d",KM_SELECTED),sprintf("Micro\nTop%d",KC_SELECTED),
                            sprintf("EarlyFusion\nTop%d+%d",KM_SELECTED,KC_SELECTED),sprintf("Stacking\nTop%d+%d",KM_SELECTED,KC_SELECTED)),
                          c("Metab_TopK","Micro_TopK","EF_TopK","Stacking_TopK"))
p_fig5a <- make_violin_wilcox(bind_rows(res_metab_tk,res_micro_tk,res_ef_tk,res_stk_tk),
                              c("Metab_TopK","Micro_TopK","EF_TopK","Stacking_TopK"), MODEL_COLORS_TOPK,
                              sprintf("PMI prediction: parsimonious single-omics vs multi-omics integration\n(Top-%d metabolites + Top-%d microbial species; 50\u00d7 nested CV)",KM_SELECTED,KC_SELECTED),
                              "Test RMSE (days)",
                              compare_pairs=list(c("Metab_TopK","EF_TopK"),c("Metab_TopK","Stacking_TopK"),c("EF_TopK","Stacking_TopK")),
                              base_size=15) +
  scale_x_discrete(labels=x_labels_main) +
  theme(axis.text.x=element_text(angle=0,hjust=0.5,size=12,lineheight=0.9,family="Arial"))
ggsave(file.path(output_dir,"Fig5a_model_performance_5tissues.tiff"),p_fig5a,width=16,height=7,device="tiff",dpi=300,compression="lzw")
cat("  [OK] Fig5a saved\n")

# Figure 5(b)
pooled_rmse      <- summ_D %>% filter(Model=="Stacking_TopK") %>% pull(RMSE_mean)
best_single_rmse <- round(mean((bind_rows(
  res_metab_tk%>%group_by(Label)%>%summarise(rmse=mean(RMSE_test,na.rm=TRUE),.groups="drop")%>%mutate(src="Metab"),
  res_micro_tk%>%group_by(Label)%>%summarise(rmse=mean(RMSE_test,na.rm=TRUE),.groups="drop")%>%mutate(src="Micro"))%>%
    group_by(Label)%>%slice_min(rmse,n=1,with_ties=FALSE)%>%ungroup())$rmse),3)

pair_stacking <- res_E %>% filter(Model=="Stacking_TopK",!is.na(RMSE_test))
pair_order_e  <- pair_stacking %>% group_by(Label) %>% summarise(med=median(RMSE_test,na.rm=TRUE),.groups="drop") %>% arrange(med) %>% pull(Label)
pair_stacking <- pair_stacking %>% mutate(Label=factor(Label,levels=pair_order_e))

pair_summary_e <- pair_stacking %>% group_by(Label) %>%
  summarise(RMSE_mean=round(mean(RMSE_test,na.rm=TRUE),3), RMSE_sd=round(sd(RMSE_test,na.rm=TRUE),3),
            p_vs_pooled=wilcox.test(RMSE_test,mu=pooled_rmse,alternative="less",exact=FALSE)$p.value,
            p_vs_single=wilcox.test(RMSE_test,mu=best_single_rmse,alternative="less",exact=FALSE)$p.value, .groups="drop") %>%
  mutate(beats_pooled=RMSE_mean<pooled_rmse, beats_single=RMSE_mean<best_single_rmse,
         sig_pooled=case_when(p_vs_pooled<0.001~"***",p_vs_pooled<0.01~"**",p_vs_pooled<0.05~"*",TRUE~"n.s."),
         sig_single=case_when(p_vs_single<0.001~"***",p_vs_single<0.01~"**",p_vs_single<0.05~"*",TRUE~"n.s."),
         fill_col=case_when(beats_pooled&sig_pooled!="n.s."~"#1A6B4A",beats_single&sig_single!="n.s."~"#5B9E6F",TRUE~"#7FB3D3"))
pair_fill_vals   <- setNames(pair_summary_e$fill_col, as.character(pair_summary_e$Label))
sig_pooled_labs  <- pair_summary_e %>% filter(sig_pooled%in%c("*","**","***")) %>% mutate(Label=factor(Label,levels=pair_order_e),y_pos=RMSE_mean-RMSE_sd*0.3-0.06,lab=sig_pooled)
sig_single_labs  <- pair_summary_e %>% filter(sig_single%in%c("*","**","***"),!sig_pooled%in%c("*","**","***")) %>% mutate(Label=factor(Label,levels=pair_order_e),y_pos=RMSE_mean-RMSE_sd*0.3-0.06,lab=paste0("(",sig_single,")"))

p_fig5b <- pair_stacking %>%
  ggplot(aes(x=Label,y=RMSE_test,fill=Label)) +
  geom_violin(trim=FALSE,alpha=0.78,color="white",linewidth=0.3) +
  geom_boxplot(width=0.10,fill="white",color="grey40",outlier.size=0.5,linewidth=0.4) +
  stat_summary(fun=mean,geom="point",shape=18,size=3.5,color="black") +
  geom_hline(yintercept=pooled_rmse,linetype="solid",color="#943126",linewidth=0.9,alpha=0.85) +
  annotate("text",x=length(pair_order_e)+0.45,y=pooled_rmse+0.022,label=sprintf("5-organ pooled\n(%.3f d)",pooled_rmse),color="#943126",size=4.5,hjust=1,fontface="italic",family="Arial") +
  geom_hline(yintercept=best_single_rmse,linetype="dashed",color="#2C5F8A",linewidth=0.85,alpha=0.85) +
  annotate("text",x=length(pair_order_e)+0.45,y=best_single_rmse+0.022,label=sprintf("Best single-omics\n(%.3f d)",best_single_rmse),color="#2C5F8A",size=4.5,hjust=1,fontface="italic",family="Arial") +
  geom_text(data=sig_pooled_labs,aes(x=Label,y=y_pos,label=lab),inherit.aes=FALSE,size=6.5,fontface="bold",color="#1A6B4A",family="Arial") +
  geom_text(data=sig_single_labs,aes(x=Label,y=y_pos,label=lab),inherit.aes=FALSE,size=5.5,fontface="bold",color="#2C5F8A",family="Arial") +
  scale_fill_manual(values=pair_fill_vals) +
  scale_x_discrete(labels=function(x) gsub("\\+","\n+",x)) +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.10))) +
  labs(title=sprintf("Two-tissue pairwise Stacking_TopK vs baselines (Top%d+%d; 50\u00d7 nested CV)",KM_SELECTED,KC_SELECTED),
       subtitle=paste0("Solid red=5-organ pooled (",sprintf("%.3f",pooled_rmse)," d)  |  Dashed blue=best single-omics (",sprintf("%.3f",best_single_rmse)," d)\n",
                       "Dark green=significantly below pooled (p<0.05)  |  Mid green=below single-omics only  |  *** p<0.001  ** p<0.01  * p<0.05"),
       x="Tissue pair",y="Test RMSE (days)") +
  theme_bw(base_size=16,base_family="Arial") +
  theme(legend.position="none",axis.text.x=element_text(size=13,lineheight=0.85,family="Arial"),
        axis.text.y=element_text(size=14,family="Arial"),axis.title=element_text(size=15,family="Arial"),
        panel.grid.major.x=element_blank(),panel.grid.minor=element_blank(),
        plot.title=element_text(size=16,face="bold",family="Arial"),
        plot.subtitle=element_text(size=12,color="grey35",lineheight=1.3,family="Arial"))
ggsave(file.path(output_dir,"Fig5b_pairwise_tissue_performance.tiff"),p_fig5b,width=16,height=6,device="tiff",dpi=300,compression="lzw")
cat("  [OK] Fig5b saved\n")
write.csv(pair_summary_e, file.path(output_dir,"Step4_pair_vs_pooled_summary.csv"), row.names=FALSE)

# ── Figure 5(c): Cross-tissue stable features (Top-K selection frequency bubble chart) ──
N_SHOW_METAB <- 20L; N_SHOW_MICRO <- 10L; FREQ_CUT_SHOW <- 0.50; MIN_TISSUES_MAIN <- 2L
N_METAB_BAR  <- 10L; N_MICRO_BAR  <- 5L;  N_METAB_TREND <- 10L; N_MICRO_TREND    <- 5L
pmi_breaks   <- sort(unique(meta_J[[COL_PMI]]))

get_feat_pool <- function(omics_type, n_show, freq_cut)
  topk_feat_summary %>% filter(Omics==omics_type) %>% group_by(Tissue) %>%
  slice_max(select_freq,n=n_show,with_ties=FALSE) %>% ungroup() %>%
  group_by(Feature) %>% filter(max(select_freq)>=freq_cut) %>% ungroup() %>% pull(Feature) %>% unique()

get_feat_pool_by_ntissue <- function(omics_type, n_show, freq_cut, min_tissues)
  topk_feat_summary %>% filter(Omics==omics_type) %>% group_by(Tissue) %>%
  slice_max(select_freq,n=n_show,with_ties=FALSE) %>% ungroup() %>%
  group_by(Feature) %>% filter(max(select_freq)>=freq_cut) %>%
  mutate(n_qual=n_distinct(Tissue[select_freq>=freq_cut])) %>%
  filter(n_qual>=min_tissues) %>% ungroup() %>% pull(Feature) %>% unique()

pool_metab_main   <- get_feat_pool_by_ntissue("Metabolomics",N_SHOW_METAB,FREQ_CUT_SHOW,MIN_TISSUES_MAIN)
pool_micro_main   <- get_feat_pool_by_ntissue("Microbiome",  N_SHOW_MICRO,FREQ_CUT_SHOW,MIN_TISSUES_MAIN)
pool_metab_single <- setdiff(get_feat_pool("Metabolomics",N_SHOW_METAB,FREQ_CUT_SHOW), pool_metab_main)
pool_micro_single <- setdiff(get_feat_pool("Microbiome",  N_SHOW_MICRO,FREQ_CUT_SHOW), pool_micro_main)

full_feat_table <- topk_feat_summary %>%
  filter(Feature%in%c(get_feat_pool("Metabolomics",N_SHOW_METAB,FREQ_CUT_SHOW), get_feat_pool("Microbiome",N_SHOW_MICRO,FREQ_CUT_SHOW))) %>%
  mutate(select_freq_pct=round(select_freq*100,1), in_main_fig=Feature%in%c(pool_metab_main,pool_micro_main)) %>%
  select(Tissue,Omics,Display_Name,select_freq_pct,mean_rank,in_main_fig) %>% arrange(Omics,Tissue,desc(select_freq_pct))
write.csv(full_feat_table, file.path(output_dir,"SuppTable_full_feature_freq.csv"), row.names=FALSE)

make_bubble_data <- function(feat_pool, omics_type) {
  feat_order <- topk_feat_summary %>% filter(Omics==omics_type,Feature%in%feat_pool) %>%
    group_by(Feature,Display_Name) %>%
    summarise(avg_freq=mean(select_freq,na.rm=TRUE), n_tissues=n_distinct(Tissue[select_freq>=FREQ_CUT_SHOW]),.groups="drop") %>%
    arrange(desc(n_tissues),desc(avg_freq))
  expand.grid(Feature=feat_pool,Tissue=TISSUES,stringsAsFactors=FALSE) %>%
    left_join(topk_feat_summary%>%filter(Omics==omics_type)%>%select(Feature,Tissue,Display_Name,select_freq), by=c("Feature","Tissue")) %>%
    left_join(feat_order%>%select(Feature,Display_Name,avg_freq,n_tissues), by="Feature",suffix=c("","_pool")) %>%
    mutate(Display_Name=dplyr::coalesce(Display_Name,Display_Name_pool), select_freq=tidyr::replace_na(select_freq,0),
           Tissue=factor(Tissue,levels=TISSUES), Display_Name=factor(Display_Name,levels=rev(feat_order$Display_Name))) %>%
    select(-Display_Name_pool)
}

make_bubble_panel <- function(df, panel_title, bubble_color) {
  right_labs <- df %>% group_by(Display_Name,n_tissues) %>% summarise(.groups="drop") %>% distinct()
  ggplot(df,aes(x=Tissue,y=Display_Name)) +
    geom_tile(aes(fill=as.integer(Display_Name)%%2==0),color=NA,alpha=0.06) +
    scale_fill_manual(values=c("TRUE"="#F0F0F0","FALSE"="white"),guide="none") +
    geom_point(aes(size=ifelse(select_freq>0,select_freq*100,NA),alpha=select_freq),
               fill=bubble_color,color="white",shape=21,stroke=0.5,na.rm=TRUE) +
    geom_text(data=right_labs,aes(x=5.65,y=Display_Name,label=sprintf("%d/5",n_tissues)),
              size=3.4,hjust=0,color=bubble_color,inherit.aes=FALSE,family="Arial") +
    scale_size_continuous(name="Selection\nfrequency (%)",range=c(1.5,9),breaks=c(25,50,75,100),limits=c(0,100)) +
    scale_alpha_continuous(range=c(0.15,0.92),limits=c(0,1),guide="none") +
    scale_x_discrete(expand=expansion(add=c(0.5,1.15))) +
    coord_cartesian(clip="off") + labs(title=panel_title,x=NULL,y=NULL) +
    theme_bw(base_size=12,base_family="Arial") +
    theme(axis.text.x=element_text(face="bold",size=11,family="Arial"),
          axis.text.y=element_text(size=10,family="Arial"),
          panel.grid.major=element_line(color="grey90",linewidth=0.25), panel.grid.minor=element_blank(),
          panel.border=element_rect(color="grey70",linewidth=0.5),
          plot.title=element_text(size=13,face="bold",color=bubble_color,hjust=0.5,family="Arial"),
          legend.position="none", plot.margin=ggplot2::margin(t=8,r=24,b=8,l=8))
}

make_legend_panel <- function() {
  ggplot(data.frame(x=c(2.2,4.4,6.6,8.8),y=1,freq=c(25,50,75,100)),aes(x=x,y=y,size=freq)) +
    geom_point(fill="grey55",color="white",shape=21,stroke=0.5,alpha=0.75) +
    geom_text(aes(label=paste0(freq,"%")),vjust=-1.6,size=3.8,color="grey30",family="Arial") +
    annotate("text",x=0.25,y=1,label="Selection\nfrequency:",size=3.8,color="grey30",hjust=0,vjust=0.5,lineheight=0.9,family="Arial") +
    scale_size_continuous(range=c(1.5,9),limits=c(0,100)) +
    coord_cartesian(xlim=c(0,10.4),ylim=c(0.55,1.95),clip="off") +
    theme_void() + theme(legend.position="none",plot.margin=ggplot2::margin(t=0,r=20,b=10,l=20))
}

df_metab_main <- make_bubble_data(pool_metab_main,"Metabolomics")
df_micro_main <- make_bubble_data(pool_micro_main,"Microbiome")
p_fig5c <- patchwork::wrap_plots(
  patchwork::wrap_plots(
    make_bubble_panel(df_metab_main,sprintf("Metabolomics  (Top%d input features, freq \u226550%% shown)",N_SHOW_METAB),"#2980B9"),
    make_bubble_panel(df_micro_main,sprintf("Microbiome  (Top%d input features, freq \u226550%% shown)",N_SHOW_MICRO),"#E67E22"),
    ncol=1, heights=c(n_distinct(df_metab_main$Display_Name),n_distinct(df_micro_main$Display_Name))) +
    patchwork::plot_annotation(
      title=sprintf("Cross-tissue consistency of top-ranked features (Top%d metabolites + Top%d species)",N_SHOW_METAB,N_SHOW_MICRO),
      subtitle="Bubble size=selection frequency across 50 nested CV repeats  |  N/5=number of tissues meeting \u226550% threshold  |  Ordered by cross-tissue consistency",
      theme=theme(plot.title=element_text(size=14,face="bold",hjust=0.5,family="Arial"),
                  plot.subtitle=element_text(size=10,color="grey35",lineheight=1.3,family="Arial"))),
  make_legend_panel(), ncol=1, heights=c(20,2.4))
fig_h_main3 <- min((n_distinct(df_metab_main$Display_Name)+n_distinct(df_micro_main$Display_Name))*0.30+4.8,15)
ggsave(file.path(output_dir,"Fig5c_cross_tissue_stable_features.tiff"),p_fig5c,width=9.5,height=fig_h_main3,device="tiff",dpi=300,compression="lzw",limitsize=FALSE,bg="white")
cat("  [OK] Fig5c saved\n")

# Tissue-specific features supplementary
df_metab_single <- make_bubble_data(pool_metab_single,"Metabolomics")
df_micro_single <- make_bubble_data(pool_micro_single,"Microbiome")
ggsave(file.path(output_dir,"SuppFig_tissue_specific_features.tiff"),
       (make_bubble_panel(df_metab_single,sprintf("Metabolomics - tissue-specific (1/5, Top%d)",N_SHOW_METAB),"#2980B9") +
          make_bubble_panel(df_micro_single,sprintf("Microbiome - tissue-specific (1/5, Top%d)",N_SHOW_MICRO),"#E67E22") +
          patchwork::plot_layout(widths=c(1.8,1))) / make_legend_panel() + patchwork::plot_layout(heights=c(20,2.4)),
       width=14, height=min(max(n_distinct(df_metab_single$Display_Name),n_distinct(df_micro_single$Display_Name))*0.34+4.2,11),
       device="tiff",dpi=300,compression="lzw",limitsize=FALSE,bg="white")

# Per-tissue bar + trend plots
raw_metab_long <- feat_MJ %>% mutate(Organ=meta_J[[COL_TISSUE]],PMI_day=meta_J[[COL_PMI]]) %>%
  pivot_longer(cols=-c(Organ,PMI_day),names_to="Feature",values_to="Value") %>% mutate(Value=log2(Value+1))
raw_micro_long <- as.data.frame(clr_mat(feat_CJ)) %>% mutate(Organ=meta_J[[COL_TISSUE]],PMI_day=meta_J[[COL_PMI]]) %>%
  pivot_longer(cols=-c(Organ,PMI_day),names_to="Feature",values_to="Value")

for (tis in TISSUES) {
  feat_bar_m <- topk_feat_summary %>% filter(Tissue==tis,Omics=="Metabolomics",select_freq>=FREQ_CUT_SHOW) %>%
    slice_min(mean_rank,n=N_METAB_BAR,with_ties=FALSE) %>% arrange(desc(select_freq),mean_rank) %>%
    mutate(Plot_Label=factor(Display_Name,levels=rev(Display_Name)))
  feat_bar_c <- topk_feat_summary %>% filter(Tissue==tis,Omics=="Microbiome",select_freq>=FREQ_CUT_SHOW) %>%
    slice_min(mean_rank,n=N_MICRO_BAR,with_ties=FALSE) %>% arrange(desc(select_freq),mean_rank) %>%
    mutate(Plot_Label=factor(Display_Name,levels=rev(Display_Name)))
  
  make_bar_panel <- function(df, omics_type, panel_title)
    ggplot(df,aes(x=Plot_Label,y=select_freq*100)) +
    geom_col(fill=OMICS_BAR_COLORS[omics_type],width=0.72,alpha=0.88) +
    geom_text(aes(label=sprintf("%.0f%%",select_freq*100)),hjust=-0.12,size=3.0,color="grey25",family="Arial") +
    coord_flip() +
    scale_y_continuous(limits=c(0,120),expand=expansion(mult=c(0,0)),breaks=c(0,25,50,75,100),labels=function(x) paste0(x,"%")) +
    labs(title=panel_title,x=NULL,y="Selection frequency (%)") +
    theme_bw(base_size=11,base_family="Arial") +
    theme(plot.title=element_text(face="bold",size=12,hjust=0.5,family="Arial"),
          axis.text=element_text(family="Arial"),panel.grid.major.y=element_blank(),panel.grid.minor=element_blank())
  
  ggsave(file.path(output_dir,sprintf("Fig5c_ImportanceBar_%s.tiff",tis)),
         make_bar_panel(feat_bar_m,"Metabolomics","Metabolomics") + make_bar_panel(feat_bar_c,"Microbiome","Microbiome") +
           plot_layout(ncol=2,widths=c(1.8,1)) +
           plot_annotation(title=sprintf("%s  Stable input features (freq \u226550%%)",tis),
                           theme=theme(plot.title=element_text(size=13,face="bold",family="Arial"))),
         width=14, height=max(5.5,max(nrow(feat_bar_m),nrow(feat_bar_c))*0.44+3.5),
         device="tiff",dpi=300,compression="lzw")
  
  feat_trend_m <- topk_feat_summary %>% filter(Tissue==tis,Omics=="Metabolomics",select_freq>=FREQ_CUT_SHOW) %>%
    slice_min(mean_rank,n=N_METAB_TREND,with_ties=FALSE) %>% arrange(mean_rank) %>%
    mutate(panel_label=sprintf("%s\n(freq=%.0f%%)",Display_Name,select_freq*100))
  df_trend_m <- raw_metab_long %>% filter(Organ==tis,Feature%in%feat_trend_m$Feature) %>%
    left_join(feat_trend_m%>%select(Feature,panel_label),by="Feature") %>%
    mutate(panel_label=factor(panel_label,levels=feat_trend_m$panel_label))
  ggsave(file.path(output_dir,sprintf("SuppFig8_Trend_Metab_%s.tiff",tis)),
         ggplot(df_trend_m,aes(x=PMI_day,y=Value)) +
           geom_point(color=METAB_TREND_COLOR,size=1.8,alpha=0.60) +
           geom_smooth(method="loess",span=0.8,se=TRUE,color=METAB_TREND_COLOR,fill=METAB_TREND_COLOR,alpha=0.15,linewidth=1.0) +
           scale_x_continuous(breaks=pmi_breaks) + facet_wrap(~panel_label,scales="free_y",ncol=5) +
           labs(title=sprintf("%s - PMI trend of stable metabolites",tis),x="PMI (days)",y="Log2 abundance") +
           theme_bw(base_size=10,base_family="Arial") +
           theme(strip.text=element_text(size=7.5,lineheight=0.85,family="Arial"),strip.background=element_rect(fill="grey95"),
                 panel.spacing=unit(0.4,"lines"),plot.title=element_text(size=12,face="bold",family="Arial"),
                 axis.text=element_text(size=8,family="Arial"),panel.grid.minor=element_blank()),
         width=14, height=ceiling(N_METAB_TREND/5)*3.8+2.5, device="tiff",dpi=300,compression="lzw")
  
  feat_trend_c <- topk_feat_summary %>% filter(Tissue==tis,Omics=="Microbiome",select_freq>=FREQ_CUT_SHOW) %>%
    slice_min(mean_rank,n=N_MICRO_TREND,with_ties=FALSE) %>% arrange(mean_rank) %>%
    mutate(panel_label=sprintf("%s\n(freq=%.0f%%)",Display_Name,select_freq*100))
  df_trend_c <- raw_micro_long %>% filter(Organ==tis,Feature%in%feat_trend_c$Feature) %>%
    left_join(feat_trend_c%>%select(Feature,panel_label),by="Feature") %>%
    mutate(panel_label=factor(panel_label,levels=feat_trend_c$panel_label))
  ggsave(file.path(output_dir,sprintf("SuppFig8_Trend_Micro_%s.tiff",tis)),
         ggplot(df_trend_c,aes(x=PMI_day,y=Value)) +
           geom_point(color=MICRO_TREND_COLOR,size=1.8,alpha=0.60) +
           geom_smooth(method="loess",span=0.9,se=TRUE,color=MICRO_TREND_COLOR,fill=MICRO_TREND_COLOR,alpha=0.15,linewidth=1.0) +
           scale_x_continuous(breaks=pmi_breaks) + facet_wrap(~panel_label,scales="free_y",ncol=N_MICRO_TREND,nrow=1) +
           labs(title=sprintf("%s - PMI trend of stable species",tis),x="PMI (days)",y="CLR abundance") +
           theme_bw(base_size=10,base_family="Arial") +
           theme(strip.text=element_text(size=8.5,lineheight=0.85,family="Arial"),strip.background=element_rect(fill="grey95"),
                 panel.spacing=unit(0.5,"lines"),plot.title=element_text(size=12,face="bold",family="Arial"),
                 axis.text=element_text(size=8,family="Arial"),panel.grid.minor=element_blank()),
         width=max(10,N_MICRO_TREND*2.8), height=5.0, device="tiff",dpi=300,compression="lzw")
  
  cat(sprintf("  [OK] %s: bar + trend plots saved\n", tis))
}

# ── Final summary ──────────────────────────────────────────────────────────────
final_tbl <- bind_rows(
  res_metab_tk%>%group_by(Label)%>%summarise(RMSE_mean=round(mean(RMSE_test,na.rm=TRUE),4),RMSE_sd=round(sd(RMSE_test,na.rm=TRUE),4),.groups="drop")%>%mutate(Model=sprintf("Metab_Top%d",KM_SELECTED),RMSE_str=sprintf("%.3f \u00b1 %.3f",RMSE_mean,RMSE_sd)),
  res_micro_tk%>%group_by(Label)%>%summarise(RMSE_mean=round(mean(RMSE_test,na.rm=TRUE),4),RMSE_sd=round(sd(RMSE_test,na.rm=TRUE),4),.groups="drop")%>%mutate(Model=sprintf("Micro_Top%d",KC_SELECTED),RMSE_str=sprintf("%.3f \u00b1 %.3f",RMSE_mean,RMSE_sd)),
  summ_EF_topk%>%select(Label,Model,RMSE_mean,RMSE_sd,RMSE_str),
  summ_ST_topk%>%select(Label,Model,RMSE_mean,RMSE_sd,RMSE_str)) %>% arrange(Label,RMSE_mean)
print(final_tbl%>%select(Label,Model,RMSE_str), n=40)
write.csv(final_tbl, file.path(output_dir,"FINAL_summary_table.csv"), row.names=FALSE)

cat(sprintf("\n[DONE] v7 complete! Metab K=%d, Micro K=%d\nOutput: %s\n", KM_SELECTED, KC_SELECTED, output_dir))

