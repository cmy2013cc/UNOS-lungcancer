#!/usr/bin/env Rscript

local_lib <- Sys.getenv("EJCTS_R_LIB", unset = "E:/UNOS/EJCTS_reanalysis")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

## Section 8: fair machine-learning comparison for post-transplant lung-cancer prediction
## Primary comparison: identical 13-domain transplantation-baseline predictor set.
## Secondary comparison: identical expanded 19-domain baseline predictor set.
## No post-transplant treatment or rejection variable is used as a predictor.

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260804)

suppressPackageStartupMessages({
  library(riskRegression)
  library(prodlim)
  library(cmprsk)
  library(survival)
  library(glmnet)
  library(randomForestSRC)
  library(xgboost)
  library(ggplot2)
  library(patchwork)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
source_file <- if (length(args) >= 1) args[[1]] else file.path("..", "2.characteristic", "2.lung.rds")
output_dir <- if (length(args) >= 2) args[[2]] else "."
fast_mode <- identical(Sys.getenv("UNOS8_FAST", unset = "0"), "1")
if (!file.exists(source_file)) stop("Input file not found: ", source_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat_raw <- readRDS(source_file)
stopifnot(nrow(dat_raw) == 21104L)
stopifnot(sum(dat_raw$FGSTATUS == 1L) == 459L)
stopifnot(sum(dat_raw$FGSTATUS == 2L) == 9603L)
stopifnot(sum(dat_raw$FGSTATUS == 0L) == 11042L)
stopifnot(!anyNA(dat_raw$FGDAY), !anyNA(dat_raw$FGSTATUS))
stopifnot(!"LUNG_DONOR" %in% names(dat_raw))

core_raw <- c(
  "AGE", "GENDER", "GROUPING", "CIG_USE", "TX_TYPE", "ISCHTIME",
  "CREAT_TRR", "DIAB", "CMV_STATUS", "AGE_DON", "GENDER_DON",
  "HIST_CIG_DON", "DIABETES_DON", "COD_CAD_DON"
)
additional_raw <- c(
  "BMI_TCR", "TBILI", "EBV_SEROSTATUS", "BMI_DON_CALC",
  "HIST_HYPERTENS_DON", "ALCOHOL_HEAVY_DON"
)
required <- c(core_raw, additional_raw, "FGDAY", "FGSTATUS", "TX_YEAR_GROUP")
stopifnot(all(required %in% names(dat_raw)))

feature_sets <- c("Core 13", "Expanded 19")
horizon_years <- c(3, 5, 10)
horizon_days <- horizon_years * 365.25
n_outer <- 5L
n_boot <- if (fast_mode) 30L else 300L
rsf_final_trees <- if (fast_mode) 30L else 150L
xgb_max_rounds <- if (fast_mode) 30L else 150L

mode_value <- function(x) {
  x <- as.character(x[!is.na(x) & as.character(x) != ""])
  names(sort(table(x), decreasing = TRUE))[1]
}

derive_imputation <- function(d, feature_set) {
  raw_vars <- if (feature_set == "Core 13") core_raw else c(core_raw, additional_raw)
  continuous <- intersect(raw_vars, c("AGE", "ISCHTIME", "CREAT_TRR", "AGE_DON", "BMI_TCR", "TBILI", "BMI_DON_CALC"))
  categorical <- setdiff(raw_vars, continuous)
  list(
    continuous = setNames(lapply(continuous, function(v) median(d[[v]], na.rm = TRUE)), continuous),
    categorical = setNames(lapply(categorical, function(v) mode_value(d[[v]])), categorical)
  )
}

fill_missing <- function(x, value) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- as.character(value)
  x
}

prepare_data <- function(d, imp, feature_set) {
  for (v in names(imp$continuous)) d[[v]][is.na(d[[v]])] <- imp$continuous[[v]]
  for (v in names(imp$categorical)) d[[v]] <- fill_missing(d[[v]], imp$categorical[[v]])

  d$AGE10 <- d$AGE / 10
  d$AGE_DON10 <- d$AGE_DON / 10
  d$GENDER <- factor(d$GENDER, levels = c("F", "M"), labels = c("Female", "Male"))
  d$DIAG_GROUP <- factor(
    d$GROUPING, levels = c("A", "D", "Other"),
    labels = c("Chronic obstructive pulmonary disease", "Interstitial lung disease", "Other diagnosis")
  )
  smoke_label <- ifelse(d$CIG_USE == "Y", "Smoking history", "No smoking history")
  tx_label <- ifelse(d$TX_TYPE == "S", "single lung", "double lung")
  d$SMOKE_TX_GROUP <- factor(
    paste(smoke_label, tx_label, sep = " + "),
    levels = c(
      "No smoking history + double lung", "Smoking history + double lung",
      "No smoking history + single lung", "Smoking history + single lung"
    )
  )
  d$DIAB <- factor(d$DIAB, levels = c("N", "Y"), labels = c("No", "Yes"))
  d$CMV_STATUS <- factor(d$CMV_STATUS, levels = c("N", "P"), labels = c("Negative", "Positive"))
  d$GENDER_DON <- factor(d$GENDER_DON, levels = c("F", "M"), labels = c("Female", "Male"))
  d$HIST_CIG_DON <- factor(d$HIST_CIG_DON, levels = c("N", "Y"), labels = c("No", "Yes"))
  d$DIABETES_DON <- factor(d$DIABETES_DON, levels = c("N", "Y"), labels = c("No", "Yes"))
  d$COD_CAD_DON <- factor(d$COD_CAD_DON, levels = c("1", "2", "3", "Other"))
  d$TX_YEAR_GROUP <- factor(d$TX_YEAR_GROUP, levels = c("2005-2010", "2011-2015", "2016-2020", "2021-2023"))

  if (feature_set == "Expanded 19") {
    d$BMI_TCR5 <- d$BMI_TCR / 5
    d$BMI_DON5 <- d$BMI_DON_CALC / 5
    d$EBV_SEROSTATUS <- factor(d$EBV_SEROSTATUS, levels = c("N", "P"), labels = c("Negative", "Positive"))
    d$HIST_HYPERTENS_DON <- factor(d$HIST_HYPERTENS_DON, levels = c("N", "Y"), labels = c("No", "Yes"))
    d$ALCOHOL_HEAVY_DON <- factor(d$ALCOHOL_HEAVY_DON, levels = c("N", "Y"), labels = c("No", "Yes"))
  }
  d
}

core_terms <- c(
  "AGE10", "GENDER", "DIAG_GROUP", "SMOKE_TX_GROUP", "ISCHTIME", "CREAT_TRR",
  "DIAB", "CMV_STATUS", "AGE_DON10", "GENDER_DON", "HIST_CIG_DON",
  "DIABETES_DON", "COD_CAD_DON"
)
expanded_terms <- c(
  core_terms, "BMI_TCR5", "TBILI", "EBV_SEROSTATUS", "BMI_DON5",
  "HIST_HYPERTENS_DON", "ALCOHOL_HEAVY_DON"
)

terms_for <- function(feature_set) if (feature_set == "Core 13") core_terms else expanded_terms
rhs_formula <- function(feature_set) as.formula(paste("~", paste(terms_for(feature_set), collapse = " + ")))
fg_formula <- function(feature_set) as.formula(paste(
  "prodlim::Hist(FGDAY, FGSTATUS) ~", paste(terms_for(feature_set), collapse = " + ")
))
rsf_formula <- function(feature_set) as.formula(paste(
  "Surv(FGDAY, FGSTATUS) ~", paste(terms_for(feature_set), collapse = " + ")
))

make_matrix <- function(d, feature_set) {
  model.matrix(rhs_formula(feature_set), data = d)[, -1, drop = FALSE]
}

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  out <- integer(length(y))
  for (val in sort(unique(y))) {
    idx <- which(y == val)
    out[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  out
}

outer_fold <- integer(nrow(dat_raw))
set.seed(20260804)
for (status_value in sort(unique(dat_raw$FGSTATUS))) {
  idx <- which(dat_raw$FGSTATUS == status_value)
  outer_fold[idx] <- sample(rep(seq_len(n_outer), length.out = length(idx)))
}
stopifnot(identical(as.integer(table(outer_fold)), c(4222L, 4222L, 4221L, 4220L, 4219L)))

step_surv <- function(sf, t) {
  idx <- findInterval(t, sf$time)
  out <- rep(1, length(t))
  keep <- idx > 0
  out[keep] <- sf$surv[idx[keep]]
  pmax(out, 1e-06)
}

ipcw_components <- function(d, target_time) {
  sf <- survival::survfit(survival::Surv(FGDAY, FGSTATUS == 0L) ~ 1, data = d)
  known <- d$FGDAY >= target_time | (d$FGDAY < target_time & d$FGSTATUS != 0L)
  eval_time <- pmin(d$FGDAY, target_time)
  weights <- ifelse(known, 1 / step_surv(sf, eval_time), 0)
  y <- as.integer(d$FGSTATUS == 1L & d$FGDAY <= target_time)
  list(y = y, w = weights, known = known)
}

weighted_auc <- function(score, y, w) {
  keep <- is.finite(score) & is.finite(w) & w > 0
  score <- score[keep]; y <- y[keep]; w <- w[keep]
  case_w <- ifelse(y == 1L, w, 0)
  control_w <- ifelse(y == 0L, w, 0)
  if (sum(case_w) <= 0 || sum(control_w) <= 0) return(NA_real_)
  values <- sort(unique(score))
  agg_case <- tapply(case_w, score, sum)
  agg_control <- tapply(control_w, score, sum)
  case_by <- as.numeric(agg_case[as.character(values)])
  control_by <- as.numeric(agg_control[as.character(values)])
  case_by[is.na(case_by)] <- 0; control_by[is.na(control_by)] <- 0
  control_below <- c(0, head(cumsum(control_by), -1))
  sum(case_by * (control_below + 0.5 * control_by)) / (sum(case_by) * sum(control_by))
}

metric_at <- function(score, comp, calibration = TRUE) {
  y <- comp$y; w <- comp$w
  auc <- weighted_auc(score, y, w)
  prevalence <- sum(w * y) / sum(w)
  brier <- mean(w * (y - score)^2)
  brier_null <- mean(w * (y - prevalence)^2)
  out <- c(AUC = auc, Brier = brier, Scaled_Brier = 1 - brier / brier_null)
  if (calibration) {
    p <- pmin(pmax(score, 1e-06), 1 - 1e-06)
    cal <- suppressWarnings(glm(y ~ qlogis(p), family = binomial(), weights = w))
    out <- c(out, Calibration_intercept = unname(coef(cal)[1]), Calibration_slope = unname(coef(cal)[2]))
  }
  out
}

extract_rsf_cif <- function(obj, target_times, oob = FALSE) {
  arr <- if (oob) obj$cif.oob else obj$cif
  vapply(target_times, function(t) {
    idx <- findInterval(t, obj$time.interest)
    if (idx < 1) rep(0, dim(arr)[1]) else arr[, idx, 1]
  }, numeric(dim(arr)[1]))
}

tune_elastic <- function(x, d, target_time, seed) {
  comp <- ipcw_components(d, target_time)
  keep <- which(comp$w > 0)
  y <- comp$y[keep]; w <- comp$w[keep]; xx <- x[keep, , drop = FALSE]
  foldid <- make_stratified_folds(y, 3L, seed)
  alpha_grid <- if (fast_mode) 0.5 else c(0, 0.5, 1)
  candidates <- lapply(alpha_grid, function(alpha) {
    cv <- cv.glmnet(
      xx, y, family = "binomial", alpha = alpha, weights = w,
      foldid = foldid, type.measure = "deviance", standardize = TRUE,
      nlambda = if (fast_mode) 30 else 70
    )
    data.frame(alpha = alpha, lambda = cv$lambda.min, cvm = min(cv$cvm, na.rm = TRUE))
  })
  cand <- do.call(rbind, candidates)
  best <- cand[which.min(cand$cvm), , drop = FALSE]
  fit <- glmnet(
    xx, y, family = "binomial", alpha = best$alpha, lambda = best$lambda,
    weights = w, standardize = TRUE
  )
  list(fit = fit, alpha = best$alpha, lambda = best$lambda, cvm = best$cvm)
}

predict_elastic <- function(obj, x) {
  as.numeric(predict(obj$fit, newx = x, type = "response", s = obj$lambda))
}

tune_xgb <- function(x, d, target_time, seed) {
  comp <- ipcw_components(d, target_time)
  keep <- which(comp$w > 0)
  y <- comp$y[keep]; w <- comp$w[keep]; xx <- x[keep, , drop = FALSE]
  inner <- make_stratified_folds(y, 3L, seed)
  folds <- split(seq_along(y), inner)
  grid <- if (fast_mode) {
    data.frame(max_depth = 2L, eta = 0.08)
  } else {
    expand.grid(max_depth = c(1L, 2L, 3L), eta = c(0.03, 0.08))
  }
  dtrain <- xgb.DMatrix(xx, label = y, weight = w)
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    params <- list(
      objective = "binary:logistic", eval_metric = "logloss",
      max_depth = as.integer(grid$max_depth[i]), eta = grid$eta[i],
      min_child_weight = 20, subsample = 0.85, colsample_bytree = 0.85,
      tree_method = "hist", nthread = 2
    )
    set.seed(seed + i)
    cv <- xgb.cv(
      params = params, data = dtrain, nrounds = xgb_max_rounds,
      folds = folds, verbose = 0, early_stopping_rounds = if (fast_mode) 5 else 15,
      maximize = FALSE
    )
    evaluation_log <- cv$evaluation_log
    metric_col <- grep("^test_.*_mean$", names(evaluation_log), value = TRUE)[1]
    if (is.na(metric_col)) stop("XGBoost cross-validation did not return a test metric column.")
    best_iter <- cv$best_iteration
    if (length(best_iter) == 0 || !is.finite(best_iter[[1]]) || best_iter[[1]] < 1) {
      # xgboost >=3 may omit best_iteration from xgb.cv even when an evaluation
      # log is present. Select the minimum held-out log loss explicitly.
      best_iter <- which.min(evaluation_log[[metric_col]])
    }
    best_iter <- as.integer(best_iter[[1]])
    rows[[i]] <- data.frame(
      max_depth = params$max_depth, eta = params$eta, nrounds = best_iter,
      logloss = evaluation_log[[metric_col]][best_iter]
    )
  }
  candidates <- do.call(rbind, rows)
  best <- candidates[which.min(candidates$logloss), , drop = FALSE]
  params <- list(
    objective = "binary:logistic", eval_metric = "logloss",
    max_depth = as.integer(best$max_depth), eta = best$eta,
    min_child_weight = 20, subsample = 0.85, colsample_bytree = 0.85,
    tree_method = "hist", nthread = 2
  )
  set.seed(seed + 999)
  fit <- xgb.train(params = params, data = dtrain, nrounds = as.integer(best$nrounds), verbose = 0)
  list(fit = fit, max_depth = best$max_depth, eta = best$eta,
       nrounds = best$nrounds, logloss = best$logloss)
}

predict_xgb <- function(obj, x) as.numeric(predict(obj$fit, x))

tune_rsf <- function(d, feature_set, seed) {
  p <- length(terms_for(feature_set))
  data.frame(mtry = max(1L, round(sqrt(p))), nodesize = 30L, mean_brier = NA_real_)
}

fit_rsf <- function(d, feature_set, tuning, seed, importance = FALSE) {
  set.seed(seed)
  rfsrc(
    rsf_formula(feature_set), data = d, ntree = rsf_final_trees,
    mtry = tuning$mtry, nodesize = tuning$nodesize, splitrule = "logrankCR",
    nsplit = 5, ntime = 50, importance = "none",
    split.depth = if (importance) "all.trees" else FALSE,
    samptype = "swor", sampsize = min(6000L, nrow(d)),
    forest = TRUE, block.size = 10,
    na.action = "na.omit"
  )
}

compact_rsf <- function(fit) {
  for (component in c("chf", "chf.oob", "cif", "cif.oob", "predicted", "predicted.oob",
                      "event.info", "xvar", "yvar", "err.rate", "err.block.rate", "split.depth")) {
    fit[[component]] <- NULL
  }
  fit
}

algorithm_names <- c("Fine-Gray", "Elastic net", "Competing-risk RSF", "XGBoost")
model_grid <- expand.grid(Feature_set = feature_sets, Algorithm = algorithm_names, stringsAsFactors = FALSE)
model_grid$Model_id <- paste(model_grid$Feature_set, model_grid$Algorithm, sep = " | ")
model_grid$Model_label <- paste0(model_grid$Algorithm, " (", model_grid$Feature_set, ")")
model_ids <- model_grid$Model_id

oof <- array(
  NA_real_, dim = c(nrow(dat_raw), length(horizon_days), length(model_ids)),
  dimnames = list(NULL, paste0(horizon_years, "y"), model_ids)
)
hyper_rows <- list()
fold_rows <- list()

cat("Starting Section 8 nested five-fold model comparison", if (fast_mode) "[FAST MODE]" else "", "\n")
for (feature_set in feature_sets) {
  for (fold in seq_len(n_outer)) {
    valid_idx <- which(outer_fold == fold)
    train_idx <- which(outer_fold != fold)
    imp <- derive_imputation(dat_raw[train_idx, , drop = FALSE], feature_set)
    train <- prepare_data(dat_raw[train_idx, , drop = FALSE], imp, feature_set)
    valid <- prepare_data(dat_raw[valid_idx, , drop = FALSE], imp, feature_set)
    x_train <- make_matrix(train, feature_set)
    x_valid <- make_matrix(valid, feature_set)
    stopifnot(identical(colnames(x_train), colnames(x_valid)))

    fg <- riskRegression::FGR(fg_formula(feature_set), cause = 1, data = train)
    stopifnot(isTRUE(fg$crrFit$converged))
    fg_id <- paste(feature_set, "Fine-Gray", sep = " | ")
    oof[valid_idx, , fg_id] <- riskRegression::predictRisk(fg, newdata = valid, times = horizon_days)

    rsf_tuning <- tune_rsf(train, feature_set, 820000 + fold + match(feature_set, feature_sets) * 100)
    rsf <- fit_rsf(train, feature_set, rsf_tuning, 830000 + fold, importance = FALSE)
    rsf_pred <- predict(rsf, newdata = valid)
    rsf_id <- paste(feature_set, "Competing-risk RSF", sep = " | ")
    oof[valid_idx, , rsf_id] <- extract_rsf_cif(rsf_pred, horizon_days)
    hyper_rows[[length(hyper_rows) + 1L]] <- data.frame(
      Feature_set = feature_set, Algorithm = "Competing-risk RSF", Fold = fold,
      Horizon_years = NA, Parameter = "mtry / nodesize / trees",
      Value = paste(rsf_tuning$mtry, rsf_tuning$nodesize, rsf_final_trees, sep = " / ")
    )

    for (j in seq_along(horizon_days)) {
      seed_base <- 840000 + fold * 1000 + j * 10 + match(feature_set, feature_sets)
      en <- tune_elastic(x_train, train, horizon_days[j], seed_base)
      en_id <- paste(feature_set, "Elastic net", sep = " | ")
      oof[valid_idx, j, en_id] <- predict_elastic(en, x_valid)
      hyper_rows[[length(hyper_rows) + 1L]] <- data.frame(
        Feature_set = feature_set, Algorithm = "Elastic net", Fold = fold,
        Horizon_years = horizon_years[j], Parameter = "alpha / lambda",
        Value = sprintf("%.2f / %.6g", en$alpha, en$lambda)
      )

      xgb <- tune_xgb(x_train, train, horizon_days[j], seed_base + 500)
      xgb_id <- paste(feature_set, "XGBoost", sep = " | ")
      oof[valid_idx, j, xgb_id] <- predict_xgb(xgb, x_valid)
      hyper_rows[[length(hyper_rows) + 1L]] <- data.frame(
        Feature_set = feature_set, Algorithm = "XGBoost", Fold = fold,
        Horizon_years = horizon_years[j], Parameter = "depth / eta / rounds",
        Value = sprintf("%d / %.2f / %d", xgb$max_depth, xgb$eta, xgb$nrounds)
      )
    }

    fold_rows[[length(fold_rows) + 1L]] <- data.frame(
      Feature_set = feature_set, Fold = fold, Training_N = length(train_idx), Validation_N = length(valid_idx),
      Validation_lung_cancers = sum(dat_raw$FGSTATUS[valid_idx] == 1L),
      Validation_competing_deaths = sum(dat_raw$FGSTATUS[valid_idx] == 2L),
      Validation_censored = sum(dat_raw$FGSTATUS[valid_idx] == 0L),
      All_models_completed = TRUE
    )
    cat("Completed", feature_set, "outer fold", fold, "of", n_outer, "\n")
    flush.console()
  }
}
stopifnot(!anyNA(oof), all(is.finite(oof)), all(oof >= 0 & oof <= 1))
hyperparameters <- do.call(rbind, hyper_rows)
fold_audit <- do.call(rbind, fold_rows)

## Performance, paired bootstrap differences, ROC, calibration, decision curves.
performance_rows <- list(); difference_rows <- list(); roc_rows <- list()
calibration_rows <- list(); dca_rows <- list(); era_rows <- list(); risk_group_rows <- list()
set.seed(20260804)
boot_indices <- lapply(seq_len(n_boot), function(i) sample.int(nrow(dat_raw), replace = TRUE))

cif_at <- function(d, target_time) {
  fit <- cmprsk::cuminc(ftime = d$FGDAY, fstatus = d$FGSTATUS, cencode = 0)
  event_name <- grep(" 1$", names(fit), value = TRUE)[1]
  if (is.na(event_name)) return(c(estimate = 0, lower = 0, upper = 0))
  obj <- fit[[event_name]]
  idx <- max(which(obj$time <= target_time))
  if (!is.finite(idx)) return(c(estimate = 0, lower = 0, upper = 0))
  est <- obj$est[idx]; se <- sqrt(obj$var[idx])
  c(estimate = est, lower = max(0, est - 1.96 * se), upper = min(1, est + 1.96 * se))
}

for (j in seq_along(horizon_days)) {
  comp <- ipcw_components(dat_raw, horizon_days[j])
  point <- matrix(NA_real_, nrow = length(model_ids), ncol = 5,
                  dimnames = list(model_ids, c("AUC", "Brier", "Scaled_Brier", "Calibration_intercept", "Calibration_slope")))
  for (m in model_ids) point[m, ] <- metric_at(oof[, j, m], comp, TRUE)
  boot <- array(NA_real_, dim = c(n_boot, length(model_ids), 3),
                dimnames = list(NULL, model_ids, c("AUC", "Brier", "Scaled_Brier")))
  for (b in seq_len(n_boot)) {
    idx <- boot_indices[[b]]
    comp_b <- list(y = comp$y[idx], w = comp$w[idx])
    for (m in model_ids) boot[b, m, ] <- metric_at(oof[idx, j, m], comp_b, FALSE)
  }
  obs <- cif_at(dat_raw, horizon_days[j])
  for (m in model_ids) {
    info <- model_grid[match(m, model_grid$Model_id), ]
    ci <- apply(boot[, m, , drop = FALSE], 3, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
    performance_rows[[length(performance_rows) + 1L]] <- data.frame(
      Feature_set = info$Feature_set, Algorithm = info$Algorithm, Model_id = m,
      Horizon_years = horizon_years[j], Lung_cancers_by_horizon = sum(dat_raw$FGSTATUS == 1L & dat_raw$FGDAY <= horizon_days[j]),
      At_risk_at_horizon = sum(dat_raw$FGDAY >= horizon_days[j]), Observed_CIF = obs["estimate"],
      AUC = point[m, "AUC"], AUC_lower95 = ci[1, "AUC"], AUC_upper95 = ci[2, "AUC"],
      Brier = point[m, "Brier"], Brier_lower95 = ci[1, "Brier"], Brier_upper95 = ci[2, "Brier"],
      Scaled_Brier = point[m, "Scaled_Brier"], Calibration_intercept = point[m, "Calibration_intercept"],
      Calibration_slope = point[m, "Calibration_slope"]
    )
    thresholds <- unique(c(0, quantile(oof[, j, m], probs = seq(0, 1, length.out = 101)), 1))
    case_w <- ifelse(comp$y == 1L, comp$w, 0); control_w <- ifelse(comp$y == 0L, comp$w, 0)
    roc_rows[[length(roc_rows) + 1L]] <- do.call(rbind, lapply(thresholds, function(th) data.frame(
      Feature_set = info$Feature_set, Algorithm = info$Algorithm, Horizon_years = horizon_years[j],
      Threshold = th, Sensitivity = sum(case_w[oof[, j, m] >= th]) / sum(case_w),
      Specificity = sum(control_w[oof[, j, m] < th]) / sum(control_w)
    )))

    groups <- cut(oof[, j, m], breaks = unique(quantile(oof[, j, m], probs = seq(0, 1, 0.2))),
                  include.lowest = TRUE, labels = FALSE)
    for (g in sort(unique(groups))) {
      idx <- which(groups == g); obs_g <- cif_at(dat_raw[idx, , drop = FALSE], horizon_days[j])
      calibration_rows[[length(calibration_rows) + 1L]] <- data.frame(
        Feature_set = info$Feature_set, Algorithm = info$Algorithm, Horizon_years = horizon_years[j],
        Risk_quintile = g, N = length(idx), Predicted_risk = mean(oof[idx, j, m]),
        Observed_CIF = obs_g["estimate"], Observed_lower95 = obs_g["lower"], Observed_upper95 = obs_g["upper"]
      )
    }
  }

  for (feature_set in feature_sets) {
    ref <- paste(feature_set, "Fine-Gray", sep = " | ")
    for (algorithm in setdiff(algorithm_names, "Fine-Gray")) {
      m <- paste(feature_set, algorithm, sep = " | ")
      auc_diff_boot <- boot[, m, "AUC"] - boot[, ref, "AUC"]
      brier_diff_boot <- boot[, m, "Brier"] - boot[, ref, "Brier"]
      difference_rows[[length(difference_rows) + 1L]] <- data.frame(
        Feature_set = feature_set, Algorithm = algorithm, Reference = "Fine-Gray",
        Horizon_years = horizon_years[j], Delta_AUC = point[m, "AUC"] - point[ref, "AUC"],
        Delta_AUC_lower95 = quantile(auc_diff_boot, 0.025, na.rm = TRUE),
        Delta_AUC_upper95 = quantile(auc_diff_boot, 0.975, na.rm = TRUE),
        Delta_Brier = point[m, "Brier"] - point[ref, "Brier"],
        Delta_Brier_lower95 = quantile(brier_diff_boot, 0.025, na.rm = TRUE),
        Delta_Brier_upper95 = quantile(brier_diff_boot, 0.975, na.rm = TRUE)
      )
    }
  }
}
performance <- do.call(rbind, performance_rows)
paired_differences <- do.call(rbind, difference_rows)
roc_data <- do.call(rbind, roc_rows)
calibration <- do.call(rbind, calibration_rows)

for (j in 1:2) {
  comp <- ipcw_components(dat_raw, horizon_days[j])
  prevalence <- sum(comp$w * comp$y) / sum(comp$w)
  thresholds <- seq(0.0025, ifelse(j == 1, 0.05, 0.08), length.out = 80)
  for (m in model_ids) {
    info <- model_grid[match(m, model_grid$Model_id), ]
    for (th in thresholds) {
      positive <- oof[, j, m] >= th
      nb_model <- mean(comp$w * (positive * comp$y - positive * (1 - comp$y) * th / (1 - th)))
      dca_rows[[length(dca_rows) + 1L]] <- data.frame(
        Feature_set = info$Feature_set, Algorithm = info$Algorithm, Horizon_years = horizon_years[j],
        Threshold_probability = th, Model = nb_model,
        Alert_all = prevalence - (1 - prevalence) * th / (1 - th), Alert_none = 0
      )
    }
  }
}
dca <- do.call(rbind, dca_rows)

for (m in model_ids) {
  info <- model_grid[match(m, model_grid$Model_id), ]
  for (era in levels(dat_raw$TX_YEAR_GROUP)) {
    idx <- which(dat_raw$TX_YEAR_GROUP == era)
    comp <- ipcw_components(dat_raw[idx, , drop = FALSE], horizon_days[1])
    events <- sum(dat_raw$FGSTATUS[idx] == 1L & dat_raw$FGDAY[idx] <= horizon_days[1])
    obs <- cif_at(dat_raw[idx, , drop = FALSE], horizon_days[1])
    era_rows[[length(era_rows) + 1L]] <- data.frame(
      Feature_set = info$Feature_set, Algorithm = info$Algorithm, Transplant_era = era,
      N = length(idx), Lung_cancers_by_3y = events, Mean_predicted_3y_risk = mean(oof[idx, 1, m]),
      Observed_3y_CIF = obs["estimate"], AUC_3y = if (events >= 20) weighted_auc(oof[idx, 1, m], comp$y, comp$w) else NA_real_
    )
  }
  cutoffs <- quantile(oof[, 1, m], probs = c(0.6, 0.9))
  groups <- cut(oof[, 1, m], breaks = c(-Inf, cutoffs, Inf), labels = c("Low risk", "Intermediate risk", "High risk"))
  for (g in levels(groups)) {
    idx <- which(groups == g); obs3 <- cif_at(dat_raw[idx, , drop = FALSE], horizon_days[1]); obs5 <- cif_at(dat_raw[idx, , drop = FALSE], horizon_days[2])
    risk_group_rows[[length(risk_group_rows) + 1L]] <- data.frame(
      Feature_set = info$Feature_set, Algorithm = info$Algorithm, Risk_group = g, N = length(idx),
      Lung_cancers_total_followup = sum(dat_raw$FGSTATUS[idx] == 1L),
      Cancer_capture_percent = 100 * sum(dat_raw$FGSTATUS[idx] == 1L) / sum(dat_raw$FGSTATUS == 1L),
      Mean_predicted_3y_risk = mean(oof[idx, 1, m]), Observed_3y_CIF = obs3["estimate"],
      Mean_predicted_5y_risk = mean(oof[idx, 2, m]), Observed_5y_CIF = obs5["estimate"]
    )
  }
}
era_audit <- do.call(rbind, era_rows)
risk_groups <- do.call(rbind, risk_group_rows)

## Final full-cohort models and machine-learning importance.
final_models <- list(); importance_rows <- list()
for (feature_set in feature_sets) {
  imp <- derive_imputation(dat_raw, feature_set)
  full <- prepare_data(dat_raw, imp, feature_set)
  x_full <- make_matrix(full, feature_set)

  fg <- riskRegression::FGR(fg_formula(feature_set), cause = 1, data = full)
  stopifnot(isTRUE(fg$crrFit$converged))
  rsf_tuning <- tune_rsf(full, feature_set, 910000 + match(feature_set, feature_sets))
  rsf <- fit_rsf(full, feature_set, rsf_tuning, 920000 + match(feature_set, feature_sets), importance = TRUE)

  en_models <- list(); xgb_models <- list()
  for (j in seq_along(horizon_days)) {
    en_models[[paste0(horizon_years[j], "y")]] <- tune_elastic(x_full, full, horizon_days[j], 930000 + j)
    xgb_models[[paste0(horizon_years[j], "y")]] <- tune_xgb(x_full, full, horizon_days[j], 940000 + j)
    en_obj <- en_models[[paste0(horizon_years[j], "y")]]
    xgb_obj <- xgb_models[[paste0(horizon_years[j], "y")]]
    hyperparameters <- rbind(hyperparameters,
      data.frame(Feature_set = feature_set, Algorithm = "Elastic net", Fold = "Full", Horizon_years = horizon_years[j],
                 Parameter = "alpha / lambda", Value = sprintf("%.2f / %.6g", en_obj$alpha, en_obj$lambda)),
      data.frame(Feature_set = feature_set, Algorithm = "XGBoost", Fold = "Full", Horizon_years = horizon_years[j],
                 Parameter = "depth / eta / rounds", Value = sprintf("%d / %.2f / %d", xgb_obj$max_depth, xgb_obj$eta, xgb_obj$nrounds))
    )
  }
  hyperparameters <- rbind(hyperparameters,
    data.frame(Feature_set = feature_set, Algorithm = "Competing-risk RSF", Fold = "Full", Horizon_years = NA,
               Parameter = "mtry / nodesize / trees", Value = paste(rsf_tuning$mtry, rsf_tuning$nodesize, rsf_final_trees, sep = " / "))
  )

  rsf_depth <- colMeans(rsf$split.depth, na.rm = TRUE)
  names(rsf_depth) <- rsf$xvar.names
  rsf_imp <- 1 / pmax(rsf_depth, 1e-06)
  importance_rows[[length(importance_rows) + 1L]] <- data.frame(
    Feature_set = feature_set, Algorithm = "Competing-risk RSF", Predictor = names(rsf_imp),
    Importance = rsf_imp
  )
  xgb_imp <- xgb.importance(model = xgb_models[["3y"]]$fit)
  if (nrow(xgb_imp)) {
    importance_rows[[length(importance_rows) + 1L]] <- data.frame(
      Feature_set = feature_set, Algorithm = "XGBoost", Predictor = xgb_imp$Feature, Importance = xgb_imp$Gain
    )
  }

  final_models[[feature_set]] <- list(
    imputation = imp, fine_gray = fg, random_survival_forest = compact_rsf(rsf),
    elastic_net = en_models, xgboost = xgb_models, predictor_formula = rhs_formula(feature_set)
  )
}
variable_importance <- do.call(rbind, importance_rows)
variable_importance <- variable_importance[is.finite(variable_importance$Importance), ]
variable_importance <- do.call(rbind, lapply(split(variable_importance, interaction(variable_importance$Feature_set, variable_importance$Algorithm, drop = TRUE)), function(x) {
  x$Normalized_importance <- if (max(x$Importance) > 0) 100 * x$Importance / max(x$Importance) else 0
  x[order(x$Normalized_importance, decreasing = TRUE), ]
}))
rownames(variable_importance) <- NULL

model_bundle <- list(
  models = final_models, feature_sets = list(core_13 = core_terms, expanded_19 = expanded_terms),
  horizons_years = horizon_years, horizons_days = horizon_days,
  source_file = source_file, source_md5 = unname(tools::md5sum(source_file)), generated_on = as.character(Sys.Date()),
  outcome = "First post-transplant lung cancer; death before lung cancer is a competing event",
  validation = "Nested five-fold cross-validation; all preprocessing and tuning inside training folds",
  boundary = "Transplantation-baseline prediction only; no post-transplant treatment or rejection predictors"
)
saveRDS(model_bundle, file.path(output_dir, "machine_learning_models.rds"), compress = "xz")

## Predictor and missingness tables.
predictor_table <- data.frame(
  Predictor_domain = c(
    "Recipient age", "Recipient sex", "Primary diagnosis", "Smoking history x transplant type",
    "Ischemic time", "Serum creatinine", "Recipient diabetes", "CMV serostatus",
    "Donor age", "Donor sex", "Donor smoking history", "Donor diabetes", "Donor cause of death",
    "Recipient BMI", "Total bilirubin", "EBV serostatus", "Donor BMI", "Donor hypertension", "Donor heavy alcohol use"
  ),
  Core_13 = c(rep("Yes", 13), rep("No", 6)), Expanded_19 = "Yes",
  Timing = "Available at transplantation", stringsAsFactors = FALSE
)
missingness <- do.call(rbind, lapply(c(core_raw, additional_raw), function(v) data.frame(
  Variable = v, Missing_N = sum(is.na(dat_raw[[v]]) | as.character(dat_raw[[v]]) == ""),
  Missing_percent = 100 * mean(is.na(dat_raw[[v]]) | as.character(dat_raw[[v]]) == ""),
  Included_core_13 = ifelse(v %in% core_raw, "Yes", "No"), Included_expanded_19 = "Yes"
)))

## Figures.
algorithm_colors <- c("Fine-Gray" = "#17365D", "Elastic net" = "#2C7BB6", "Competing-risk RSF" = "#F28E2B", "XGBoost" = "#C83E1D")
feature_shapes <- c("Core 13" = 16, "Expanded 19" = 17)
theme_pub <- theme_minimal(base_size = 11) + theme(
  plot.title = element_text(face = "bold", color = "#17365D", size = 12),
  plot.subtitle = element_text(color = "#5B677A", size = 9), panel.grid.minor = element_blank(),
  legend.position = "bottom", legend.title = element_blank(), strip.text = element_text(face = "bold", color = "#17365D")
)

perf_plot <- performance[performance$Horizon_years %in% c(3, 5), ]
perf_plot$Horizon_label <- factor(paste0(perf_plot$Horizon_years, " years"), levels = c("3 years", "5 years"))
p_auc <- ggplot(perf_plot, aes(x = AUC, y = Algorithm, color = Algorithm, shape = Feature_set)) +
  geom_errorbar(aes(xmin = AUC_lower95, xmax = AUC_upper95), orientation = "y", width = 0.18) + geom_point(size = 2.8) +
  facet_wrap(~Horizon_label) + scale_color_manual(values = algorithm_colors) + scale_shape_manual(values = feature_shapes) +
  coord_cartesian(xlim = c(min(perf_plot$AUC_lower95) - 0.01, max(perf_plot$AUC_upper95) + 0.01)) +
  labs(title = "A. Cross-validated discrimination", x = "IPCW cumulative/dynamic AUC", y = NULL) + theme_pub
p_brier <- ggplot(perf_plot, aes(x = 100 * Brier, y = Algorithm, color = Algorithm, shape = Feature_set)) +
  geom_errorbar(aes(xmin = 100 * Brier_lower95, xmax = 100 * Brier_upper95), orientation = "y", width = 0.18) + geom_point(size = 2.8) +
  facet_wrap(~Horizon_label, scales = "free_x") + scale_color_manual(values = algorithm_colors) + scale_shape_manual(values = feature_shapes) +
  labs(title = "B. Overall prediction error", x = "IPCW Brier score, %", y = NULL) + theme_pub
p_slope <- ggplot(perf_plot, aes(x = Calibration_slope, y = Algorithm, color = Algorithm, shape = Feature_set)) +
  geom_vline(xintercept = 1, linetype = 2, color = "#888888") + geom_point(size = 2.8) +
  facet_wrap(~Horizon_label) + scale_color_manual(values = algorithm_colors) + scale_shape_manual(values = feature_shapes) +
  labs(title = "C. Calibration slope", x = "Ideal value = 1", y = NULL) + theme_pub
diff_plot <- paired_differences[paired_differences$Horizon_years %in% c(3, 5), ]
diff_plot$Horizon_label <- factor(paste0(diff_plot$Horizon_years, " years"), levels = c("3 years", "5 years"))
p_diff <- ggplot(diff_plot, aes(x = Delta_AUC, y = Algorithm, color = Algorithm, shape = Feature_set)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#888888") +
  geom_errorbar(aes(xmin = Delta_AUC_lower95, xmax = Delta_AUC_upper95), orientation = "y", width = 0.18) + geom_point(size = 2.8) +
  facet_wrap(~Horizon_label) + scale_color_manual(values = algorithm_colors) + scale_shape_manual(values = feature_shapes) +
  labs(title = "D. Paired AUC difference vs Fine-Gray", x = "Delta AUC", y = NULL) + theme_pub
main_figure <- ((p_auc | p_brier) / (p_slope | p_diff)) + plot_annotation(
  title = "Fair comparison of transplantation-baseline prediction algorithms",
  subtitle = "Same cohort, outcomes, outer folds, and preprocessing; circles = core 13, triangles = expanded 19",
  theme = theme(plot.title = element_text(face = "bold", size = 16, color = "#17365D"), plot.subtitle = element_text(size = 10, color = "#5B677A"))
) & theme(legend.position = "none")
ggsave(file.path(output_dir, "Figure8_machine_learning_performance.png"), main_figure, width = 15, height = 10, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "Figure8_machine_learning_performance.pdf"), main_figure, width = 15, height = 10, device = "pdf", bg = "white")

roc_data$Horizon_label <- factor(paste0(roc_data$Horizon_years, " years"), levels = c("3 years", "5 years", "10 years"))
roc_figure <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity, color = Algorithm)) +
  geom_line(linewidth = 0.8) + geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#888888") +
  facet_grid(Feature_set ~ Horizon_label) + scale_color_manual(values = algorithm_colors) +
  coord_equal() + labs(title = "Cross-validated ROC curves", x = "1 - specificity", y = "Sensitivity") + theme_pub
ggsave(file.path(output_dir, "FigureS8_ROC_curves.png"), roc_figure, width = 13, height = 8, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "FigureS8_ROC_curves.pdf"), roc_figure, width = 13, height = 8, device = "pdf", bg = "white")

cal_plot_data <- calibration[calibration$Horizon_years %in% c(3, 5), ]
cal_plot_data$Horizon_label <- factor(paste0(cal_plot_data$Horizon_years, " years"), levels = c("3 years", "5 years"))
cal_figure <- ggplot(cal_plot_data, aes(x = Predicted_risk, y = Observed_CIF, color = Algorithm)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#888888") +
  geom_errorbar(aes(ymin = Observed_lower95, ymax = Observed_upper95), width = 0, alpha = 0.5) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) + facet_grid(Feature_set ~ Horizon_label, scales = "free") +
  scale_color_manual(values = algorithm_colors) + scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "Cross-validated calibration by predicted-risk quintile", x = "Mean predicted risk", y = "Observed cumulative incidence") + theme_pub
ggsave(file.path(output_dir, "FigureS8_calibration.png"), cal_figure, width = 13, height = 8, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "FigureS8_calibration.pdf"), cal_figure, width = 13, height = 8, device = "pdf", bg = "white")

dca_long <- rbind(
  data.frame(dca[, c("Feature_set", "Algorithm", "Horizon_years", "Threshold_probability")], Strategy = "Model", Net_benefit = dca$Model),
  data.frame(dca[, c("Feature_set", "Algorithm", "Horizon_years", "Threshold_probability")], Strategy = "Alert all", Net_benefit = dca$Alert_all),
  data.frame(dca[, c("Feature_set", "Algorithm", "Horizon_years", "Threshold_probability")], Strategy = "Alert none", Net_benefit = dca$Alert_none)
)
dca_models <- dca_long[dca_long$Strategy == "Model", ]
dca_refs <- unique(dca_long[dca_long$Strategy != "Model", c("Feature_set", "Horizon_years", "Threshold_probability", "Strategy", "Net_benefit")])
dca_models$Horizon_label <- factor(paste0(dca_models$Horizon_years, " years"), levels = c("3 years", "5 years"))
dca_refs$Horizon_label <- factor(paste0(dca_refs$Horizon_years, " years"), levels = c("3 years", "5 years"))
dca_figure <- ggplot() +
  geom_line(data = dca_models, aes(Threshold_probability, Net_benefit, color = Algorithm), linewidth = 0.8) +
  geom_line(data = dca_refs, aes(Threshold_probability, Net_benefit, linetype = Strategy), color = "#555555", linewidth = 0.7) +
  facet_grid(Feature_set ~ Horizon_label, scales = "free_x") + scale_color_manual(values = algorithm_colors) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.5)) +
  labs(title = "Exploratory decision-curve analysis", x = "Risk threshold", y = "Net benefit", linetype = NULL) + theme_pub
ggsave(file.path(output_dir, "FigureS8_decision_curves.png"), dca_figure, width = 13, height = 8, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "FigureS8_decision_curves.pdf"), dca_figure, width = 13, height = 8, device = "pdf", bg = "white")

importance_plot_data <- do.call(rbind, lapply(split(variable_importance, interaction(variable_importance$Feature_set, variable_importance$Algorithm, drop = TRUE)), function(x) head(x, 10)))
pretty_predictor <- function(x) {
  exact <- c(
    AGE10 = "Recipient age (per 10 years)", SMOKE_TX_GROUP = "Smoking x transplant group",
    CREAT_TRR = "Recipient creatinine", AGE_DON10 = "Donor age (per 10 years)", BMI_TCR5 = "Recipient BMI (per 5 kg/m2)",
    DIAG_GROUP = "Primary diagnosis", TBILI = "Total bilirubin", ISCHTIME = "Ischemic time", COD_CAD_DON = "Donor cause of death",
    BMI_DON5 = "Donor BMI (per 5 kg/m2)", HIST_CIG_DON = "Donor smoking history", HIST_HYPERTENS_DON = "Donor hypertension",
    DIABETES_DON = "Donor diabetes", GENDER = "Recipient sex", EBV_SEROSTATUS = "EBV serostatus", DIAB = "Recipient diabetes",
    ALCOHOL_HEAVY_DON = "Donor heavy alcohol use", CMV_STATUS = "CMV serostatus", GENDER_DON = "Donor sex"
  )
  out <- ifelse(x %in% names(exact), unname(exact[x]), x)
  out <- sub("^SMOKE_TX_GROUP", "", out)
  out <- sub("^DIAG_GROUP", "Diagnosis: ", out)
  out <- sub("^DIABYes$", "Recipient diabetes: yes", out)
  out <- sub("^GENDERMale$", "Recipient sex: male", out)
  out <- sub("^COD_CAD_DON", "Donor cause of death: category ", out)
  out
}
importance_plot_data$Predictor_label <- pretty_predictor(importance_plot_data$Predictor)
importance_figure <- ggplot(importance_plot_data, aes(Normalized_importance, reorder(Predictor_label, Normalized_importance), fill = Algorithm)) +
  geom_col(width = 0.72) + facet_wrap(vars(Feature_set, Algorithm), ncol = 2, scales = "free_y") +
  scale_fill_manual(values = algorithm_colors) + labs(title = "Exploratory machine-learning variable importance", subtitle = "Normalized within each model; importance is predictive, not causal", x = "Relative importance (maximum = 100)", y = NULL) + theme_pub
ggsave(file.path(output_dir, "FigureS8_variable_importance.png"), importance_figure, width = 13, height = 9, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "FigureS8_variable_importance.pdf"), importance_figure, width = 13, height = 9, device = "pdf", bg = "white")

high_risk <- risk_groups[risk_groups$Risk_group == "High risk", ]
risk_figure <- ggplot(high_risk, aes(x = Algorithm, y = Cancer_capture_percent, fill = Algorithm)) +
  geom_col(width = 0.72) + facet_wrap(~Feature_set) + scale_fill_manual(values = algorithm_colors) +
  geom_text(aes(label = sprintf("%.1f%%", Cancer_capture_percent)), vjust = -0.3, size = 3.5) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Lung cancers concentrated in each model's highest-risk 10%", x = NULL, y = "All lung cancers captured") +
  theme_pub + theme(axis.text.x = element_text(angle = 20, hjust = 1))
ggsave(file.path(output_dir, "FigureS8_high_risk_capture.png"), risk_figure, width = 11, height = 5.5, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "FigureS8_high_risk_capture.pdf"), risk_figure, width = 11, height = 5.5, device = "pdf", bg = "white")

overview <- data.frame(
  Metric = c("Cohort N", "Lung-cancer events", "Competing deaths", "Censored at last follow-up", "Primary feature set", "Expanded feature set",
             "Algorithms", "Validation", "Primary horizon", "Secondary horizon", "Exploratory horizon", "Post-transplant predictors", "Primary interpretation rule"),
  Value = c(nrow(dat_raw), sum(dat_raw$FGSTATUS == 1), sum(dat_raw$FGSTATUS == 2), sum(dat_raw$FGSTATUS == 0),
            "13 transplantation-baseline predictor domains", "19 transplantation-baseline predictor domains",
            paste(algorithm_names, collapse = "; "), "Nested stratified 5-fold cross-validation", "3 years", "5 years", "10 years",
            "None", "Prefer Fine-Gray unless machine learning improves discrimination without worsening calibration or Brier score"),
  stringsAsFactors = FALSE
)
methods <- data.frame(
  Field = c("Source file", "Source MD5", "Generated on", "Outcome", "Competing event", "Core comparison", "Expanded comparison",
            "Outer validation", "Inner tuning", "Imputation", "Elastic net", "Random survival forest", "XGBoost", "Performance", "Calibration", "Decision curve", "Limitation"),
  Value = c(source_file, unname(tools::md5sum(source_file)), as.character(Sys.Date()), "First post-transplant lung cancer", "Death before lung cancer",
            "All algorithms use the identical 13-domain predictor set", "All algorithms use the identical 19-domain predictor set",
            "Stratified five-fold cross-validation with one out-of-fold prediction per recipient", "Three-fold tuning restricted to each outer training set",
            "Training-fold medians for continuous variables and modes for categorical variables; outcomes never impute predictors",
            "Horizon-specific IPCW weighted logistic elastic net", "Competing-risk split rule with prespecified mtry and terminal-node size", "Horizon-specific IPCW weighted gradient boosting",
            "IPCW cumulative/dynamic AUC, Brier score, scaled Brier score; 300 paired bootstrap samples in the full run",
            "Out-of-fold quintiles plus calibration intercept and slope", "Exploratory IPCW net benefit at 3 and 5 years",
            "Internal validation only; no external or prospective validation"), stringsAsFactors = FALSE
)

result_payload <- list(
  Overview = overview, Performance = performance, Paired_differences = paired_differences,
  Calibration = calibration, Risk_groups = risk_groups, Decision_curve = dca,
  Era_audit = era_audit, Hyperparameters = hyperparameters, Fold_audit = fold_audit,
  Variable_importance = variable_importance, Predictor_sets = predictor_table,
  Missingness = missingness, Methods = methods
)
write_json(result_payload, file.path(output_dir, "section8_results.json"), pretty = TRUE, auto_unbox = TRUE, digits = 10, na = "null")
saveRDS(list(Fold = outer_fold, Predictions = oof), file.path(output_dir, "section8_oof_predictions.rds"))

cat("Section 8 machine-learning comparison complete\n")
print(performance[performance$Horizon_years %in% c(3, 5), c("Feature_set", "Algorithm", "Horizon_years", "AUC", "Brier", "Calibration_slope")])
print(paired_differences[paired_differences$Horizon_years %in% c(3, 5), ])
