local_lib <- Sys.getenv("EJCTS_R_LIB", unset = "E:/UNOS/EJCTS_reanalysis")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

## Section 7: clinically usable Fine-Gray risk prediction with internal validation
## Source cohort: corrected new/2.characteristic/2.lung.rds
## Primary horizon: 3 years; key secondary horizon: 5 years; 10 years exploratory

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260804)

suppressPackageStartupMessages({
  library(riskRegression)
  library(prodlim)
  library(cmprsk)
  library(survival)
  library(ggplot2)
  library(patchwork)
  library(openxlsx)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
source_file <- if (length(args) >= 1) args[[1]] else file.path("..", "2.characteristic", "2.lung.rds")
output_dir <- if (length(args) >= 2) args[[2]] else "."
if (!file.exists(source_file)) stop("Input file not found: ", source_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat_raw <- readRDS(source_file)
stopifnot(nrow(dat_raw) == 21104L)
stopifnot(sum(dat_raw$FGSTATUS == 1L) == 459L)
stopifnot(sum(dat_raw$FGSTATUS == 2L) == 9603L)
stopifnot(sum(dat_raw$FGSTATUS == 0L) == 11042L)
stopifnot(!anyNA(dat_raw$FGDAY), !anyNA(dat_raw$FGSTATUS))
stopifnot(!"LUNG_DONOR" %in% names(dat_raw))

predictor_source <- c(
  "AGE", "GENDER", "GROUPING", "CIG_USE", "TX_TYPE", "ISCHTIME",
  "CREAT_TRR", "DIAB", "CMV_STATUS", "AGE_DON", "GENDER_DON",
  "HIST_CIG_DON", "DIABETES_DON", "COD_CAD_DON"
)
required <- c(predictor_source, "FGDAY", "FGSTATUS", "TX_YEAR_GROUP")
stopifnot(all(required %in% names(dat_raw)))

missingness <- data.frame(
  Variable = predictor_source,
  Missing_N = vapply(dat_raw[predictor_source], function(x) sum(is.na(x)), integer(1)),
  Missing_percent = 100 * vapply(dat_raw[predictor_source], function(x) mean(is.na(x)), numeric(1)),
  stringsAsFactors = FALSE
)

mode_value <- function(x) {
  x <- as.character(x[!is.na(x) & as.character(x) != ""])
  names(sort(table(x), decreasing = TRUE))[1]
}

derive_imputation <- function(d) {
  list(
    ISCHTIME = median(d$ISCHTIME, na.rm = TRUE),
    CREAT_TRR = median(d$CREAT_TRR, na.rm = TRUE),
    CIG_USE = mode_value(d$CIG_USE),
    DIAB = mode_value(d$DIAB),
    CMV_STATUS = mode_value(d$CMV_STATUS),
    HIST_CIG_DON = mode_value(d$HIST_CIG_DON),
    DIABETES_DON = mode_value(d$DIABETES_DON),
    COD_CAD_DON = mode_value(d$COD_CAD_DON)
  )
}

fill_missing <- function(x, value) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- as.character(value)
  x
}

prepare_data <- function(d, imp) {
  d$ISCHTIME[is.na(d$ISCHTIME)] <- imp$ISCHTIME
  d$CREAT_TRR[is.na(d$CREAT_TRR)] <- imp$CREAT_TRR
  d$CIG_USE <- fill_missing(d$CIG_USE, imp$CIG_USE)
  d$DIAB <- fill_missing(d$DIAB, imp$DIAB)
  d$CMV_STATUS <- fill_missing(d$CMV_STATUS, imp$CMV_STATUS)
  d$HIST_CIG_DON <- fill_missing(d$HIST_CIG_DON, imp$HIST_CIG_DON)
  d$DIABETES_DON <- fill_missing(d$DIABETES_DON, imp$DIABETES_DON)
  d$COD_CAD_DON <- fill_missing(d$COD_CAD_DON, imp$COD_CAD_DON)

  d$AGE10 <- d$AGE / 10
  d$AGE_DON10 <- d$AGE_DON / 10
  d$GENDER <- factor(d$GENDER, levels = c("F", "M"), labels = c("Female", "Male"))
  d$DIAG_GROUP <- factor(
    d$GROUPING,
    levels = c("A", "D", "Other"),
    labels = c("Chronic obstructive pulmonary disease", "Interstitial lung disease", "Other diagnosis")
  )
  smoke_label <- ifelse(d$CIG_USE == "Y", "Smoking history", "No smoking history")
  tx_label <- ifelse(d$TX_TYPE == "S", "single lung", "double lung")
  d$SMOKE_TX_GROUP <- factor(
    paste(smoke_label, tx_label, sep = " + "),
    levels = c(
      "No smoking history + double lung",
      "Smoking history + double lung",
      "No smoking history + single lung",
      "Smoking history + single lung"
    )
  )
  d$DIAB <- factor(d$DIAB, levels = c("N", "Y"), labels = c("No", "Yes"))
  d$CMV_STATUS <- factor(d$CMV_STATUS, levels = c("N", "P"), labels = c("Negative", "Positive"))
  d$GENDER_DON <- factor(d$GENDER_DON, levels = c("F", "M"), labels = c("Female", "Male"))
  d$HIST_CIG_DON <- factor(d$HIST_CIG_DON, levels = c("N", "Y"), labels = c("No", "Yes"))
  d$DIABETES_DON <- factor(d$DIABETES_DON, levels = c("N", "Y"), labels = c("No", "Yes"))
  d$COD_CAD_DON <- factor(d$COD_CAD_DON, levels = c("1", "2", "3", "Other"))
  d$TX_YEAR_GROUP <- factor(
    d$TX_YEAR_GROUP,
    levels = c("2005-2010", "2011-2015", "2016-2020", "2021-2023")
  )
  d
}

model_formula <- prodlim::Hist(FGDAY, FGSTATUS) ~
  AGE10 + GENDER + DIAG_GROUP + SMOKE_TX_GROUP + ISCHTIME + CREAT_TRR +
  DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON + HIST_CIG_DON +
  DIABETES_DON + COD_CAD_DON

horizon_years <- c(3, 5, 10)
horizon_days <- horizon_years * 365.25
n_fold <- 5L

fold_id <- integer(nrow(dat_raw))
for (status_value in sort(unique(dat_raw$FGSTATUS))) {
  idx <- which(dat_raw$FGSTATUS == status_value)
  fold_id[idx] <- sample(rep(seq_len(n_fold), length.out = length(idx)))
}

oof_risk <- matrix(NA_real_, nrow(dat_raw), length(horizon_days))
colnames(oof_risk) <- paste0("Risk_", horizon_years, "y")
fold_audit <- vector("list", n_fold)

cat("Starting 5-fold cross-validation on", nrow(dat_raw), "recipients\n")
for (fold in seq_len(n_fold)) {
  valid_idx <- which(fold_id == fold)
  train_idx <- which(fold_id != fold)
  imp <- derive_imputation(dat_raw[train_idx, , drop = FALSE])
  train_dat <- prepare_data(dat_raw[train_idx, , drop = FALSE], imp)
  valid_dat <- prepare_data(dat_raw[valid_idx, , drop = FALSE], imp)

  fit <- riskRegression::FGR(model_formula, cause = 1, data = train_dat)
  stopifnot(isTRUE(fit$crrFit$converged))
  oof_risk[valid_idx, ] <- riskRegression::predictRisk(
    fit,
    newdata = valid_dat,
    times = horizon_days
  )

  fold_audit[[fold]] <- data.frame(
    Fold = fold,
    Training_N = length(train_idx),
    Validation_N = length(valid_idx),
    Validation_lung_cancers = sum(dat_raw$FGSTATUS[valid_idx] == 1L),
    Validation_competing_deaths = sum(dat_raw$FGSTATUS[valid_idx] == 2L),
    Validation_censored = sum(dat_raw$FGSTATUS[valid_idx] == 0L),
    Converged = TRUE
  )
  cat("Completed fold", fold, "of", n_fold, "\n")
}
stopifnot(!anyNA(oof_risk))
stopifnot(all(oof_risk >= 0 & oof_risk <= 1))
fold_audit <- do.call(rbind, fold_audit)

imp_final <- derive_imputation(dat_raw)
dat_model <- prepare_data(dat_raw, imp_final)
final_fit <- riskRegression::FGR(model_formula, cause = 1, data = dat_model)
stopifnot(isTRUE(final_fit$crrFit$converged))
full_risk <- riskRegression::predictRisk(final_fit, newdata = dat_model, times = horizon_days)
colnames(full_risk) <- colnames(oof_risk)

cif_at <- function(d, target_time) {
  fit <- cmprsk::cuminc(ftime = d$FGDAY, fstatus = d$FGSTATUS, cencode = 0)
  event_name <- grep(" 1$", names(fit), value = TRUE)[1]
  if (is.na(event_name)) return(c(estimate = 0, lower = 0, upper = 0))
  obj <- fit[[event_name]]
  idx <- max(which(obj$time <= target_time))
  if (!is.finite(idx)) return(c(estimate = 0, lower = 0, upper = 0))
  est <- obj$est[idx]
  se <- sqrt(obj$var[idx])
  c(estimate = est, lower = max(0, est - 1.96 * se), upper = min(1, est + 1.96 * se))
}

step_surv <- function(sf, t) {
  idx <- findInterval(t, sf$time)
  out <- rep(1, length(t))
  keep <- idx > 0
  out[keep] <- sf$surv[idx[keep]]
  pmax(out, 1e-6)
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
  keep <- is.finite(score) & w > 0
  score <- score[keep]
  y <- y[keep]
  w <- w[keep]
  case_w <- ifelse(y == 1L, w, 0)
  control_w <- ifelse(y == 0L, w, 0)
  if (sum(case_w) <= 0 || sum(control_w) <= 0) return(NA_real_)
  agg_case <- tapply(case_w, score, sum)
  agg_control <- tapply(control_w, score, sum)
  values <- sort(unique(score))
  case_by <- as.numeric(agg_case[as.character(values)])
  control_by <- as.numeric(agg_control[as.character(values)])
  case_by[is.na(case_by)] <- 0
  control_by[is.na(control_by)] <- 0
  control_below <- c(0, head(cumsum(control_by), -1))
  numerator <- sum(case_by * (control_below + 0.5 * control_by))
  numerator / (sum(case_by) * sum(control_by))
}

metric_at <- function(score, comp) {
  y <- comp$y
  w <- comp$w
  auc <- weighted_auc(score, y, w)
  prevalence <- sum(w * y) / sum(w)
  brier <- mean(w * (y - score)^2)
  brier_null <- mean(w * (y - prevalence)^2)
  scaled_brier <- 1 - brier / brier_null
  p <- pmin(pmax(score, 1e-6), 1 - 1e-6)
  cal <- suppressWarnings(glm(y ~ qlogis(p), family = binomial(), weights = w))
  c(
    AUC = auc,
    Brier = brier,
    Scaled_Brier = scaled_brier,
    Calibration_intercept = unname(coef(cal)[1]),
    Calibration_slope = unname(coef(cal)[2])
  )
}

set.seed(20260804)
n_boot <- 300L
performance_rows <- vector("list", length(horizon_days))
roc_rows <- vector("list", length(horizon_days))
for (j in seq_along(horizon_days)) {
  comp <- ipcw_components(dat_model, horizon_days[j])
  score <- oof_risk[, j]
  point <- metric_at(score, comp)
  boot <- matrix(NA_real_, n_boot, 3L)
  colnames(boot) <- c("AUC", "Brier", "Scaled_Brier")
  for (b in seq_len(n_boot)) {
    idx <- sample.int(nrow(dat_model), replace = TRUE)
    tmp <- metric_at(score[idx], list(y = comp$y[idx], w = comp$w[idx]))
    boot[b, ] <- tmp[c("AUC", "Brier", "Scaled_Brier")]
  }
  ci <- apply(boot, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  obs <- cif_at(dat_model, horizon_days[j])
  performance_rows[[j]] <- data.frame(
    Horizon_years = horizon_years[j],
    Lung_cancers_by_horizon = sum(dat_model$FGSTATUS == 1L & dat_model$FGDAY <= horizon_days[j]),
    At_risk_at_horizon = sum(dat_model$FGDAY >= horizon_days[j]),
    Observed_CIF = obs["estimate"],
    AUC = point["AUC"],
    AUC_lower95 = ci[1, "AUC"],
    AUC_upper95 = ci[2, "AUC"],
    Brier = point["Brier"],
    Brier_lower95 = ci[1, "Brier"],
    Brier_upper95 = ci[2, "Brier"],
    Scaled_Brier = point["Scaled_Brier"],
    Scaled_Brier_lower95 = ci[1, "Scaled_Brier"],
    Scaled_Brier_upper95 = ci[2, "Scaled_Brier"],
    Calibration_intercept = point["Calibration_intercept"],
    Calibration_slope = point["Calibration_slope"]
  )

  thresholds <- unique(c(0, quantile(score, probs = seq(0, 1, length.out = 151)), 1))
  case_w <- ifelse(comp$y == 1L, comp$w, 0)
  control_w <- ifelse(comp$y == 0L, comp$w, 0)
  roc_rows[[j]] <- do.call(rbind, lapply(thresholds, function(th) {
    data.frame(
      Horizon_years = horizon_years[j],
      Threshold = th,
      Sensitivity = sum(case_w[score >= th]) / sum(case_w),
      Specificity = sum(control_w[score < th]) / sum(control_w)
    )
  }))
}
performance <- do.call(rbind, performance_rows)
roc_data <- do.call(rbind, roc_rows)
rownames(performance) <- NULL
rownames(roc_data) <- NULL

calibration_rows <- list()
for (j in seq_along(horizon_days)) {
  groups <- cut(
    oof_risk[, j],
    breaks = unique(quantile(oof_risk[, j], probs = seq(0, 1, 0.2), na.rm = TRUE)),
    include.lowest = TRUE,
    labels = FALSE
  )
  for (g in sort(unique(groups))) {
    idx <- which(groups == g)
    obs <- cif_at(dat_model[idx, , drop = FALSE], horizon_days[j])
    calibration_rows[[length(calibration_rows) + 1L]] <- data.frame(
      Horizon_years = horizon_years[j],
      Risk_quintile = g,
      N = length(idx),
      Predicted_risk = mean(oof_risk[idx, j]),
      Observed_CIF = obs["estimate"],
      Observed_lower95 = obs["lower"],
      Observed_upper95 = obs["upper"]
    )
  }
}
calibration <- do.call(rbind, calibration_rows)
rownames(calibration) <- NULL

primary_cutoffs <- as.numeric(quantile(oof_risk[, 1], probs = c(0.60, 0.90), na.rm = TRUE))
dat_model$RISK_GROUP_OOF <- cut(
  oof_risk[, 1],
  breaks = c(-Inf, primary_cutoffs, Inf),
  labels = c("Low risk", "Intermediate risk", "High risk")
)
risk_group_rows <- list()
for (g in levels(dat_model$RISK_GROUP_OOF)) {
  idx <- which(dat_model$RISK_GROUP_OOF == g)
  obs3 <- cif_at(dat_model[idx, , drop = FALSE], horizon_days[1])
  obs5 <- cif_at(dat_model[idx, , drop = FALSE], horizon_days[2])
  risk_group_rows[[length(risk_group_rows) + 1L]] <- data.frame(
    Risk_group = g,
    N = length(idx),
    Lung_cancers_total_followup = sum(dat_model$FGSTATUS[idx] == 1L),
    Mean_predicted_3y_risk = mean(oof_risk[idx, 1]),
    Observed_3y_CIF = obs3["estimate"],
    Observed_3y_lower95 = obs3["lower"],
    Observed_3y_upper95 = obs3["upper"],
    Mean_predicted_5y_risk = mean(oof_risk[idx, 2]),
    Observed_5y_CIF = obs5["estimate"],
    Observed_5y_lower95 = obs5["lower"],
    Observed_5y_upper95 = obs5["upper"]
  )
}
risk_groups <- do.call(rbind, risk_group_rows)
rownames(risk_groups) <- NULL

dca_rows <- list()
for (j in 1:2) {
  comp <- ipcw_components(dat_model, horizon_days[j])
  score <- oof_risk[, j]
  prevalence <- sum(comp$w * comp$y) / sum(comp$w)
  thresholds <- seq(0.0025, ifelse(j == 1, 0.05, 0.08), length.out = 100)
  for (th in thresholds) {
    pred_positive <- score >= th
    nb_model <- mean(comp$w * (pred_positive * comp$y - pred_positive * (1 - comp$y) * th / (1 - th)))
    nb_all <- prevalence - (1 - prevalence) * th / (1 - th)
    dca_rows[[length(dca_rows) + 1L]] <- data.frame(
      Horizon_years = horizon_years[j],
      Threshold_probability = th,
      Model = nb_model,
      Treat_all = nb_all,
      Treat_none = 0
    )
  }
}
dca <- do.call(rbind, dca_rows)
rownames(dca) <- NULL

era_rows <- list()
for (era in levels(dat_model$TX_YEAR_GROUP)) {
  idx <- which(dat_model$TX_YEAR_GROUP == era)
  comp <- ipcw_components(dat_model[idx, , drop = FALSE], horizon_days[1])
  events3 <- sum(dat_model$FGSTATUS[idx] == 1L & dat_model$FGDAY[idx] <= horizon_days[1])
  obs <- cif_at(dat_model[idx, , drop = FALSE], horizon_days[1])
  era_rows[[length(era_rows) + 1L]] <- data.frame(
    Transplant_era = era,
    N = length(idx),
    Lung_cancers_by_3y = events3,
    Mean_predicted_3y_risk = mean(oof_risk[idx, 1]),
    Observed_3y_CIF = obs["estimate"],
    AUC_3y = if (events3 >= 20L) weighted_auc(oof_risk[idx, 1], comp$y, comp$w) else NA_real_,
    Interpretation = if (events3 >= 20L) "Estimable" else "Too few 3-year lung-cancer events for stable era-specific AUC"
  )
}
era_audit <- do.call(rbind, era_rows)
rownames(era_audit) <- NULL

coef_est <- final_fit$crrFit$coef
coef_se <- sqrt(diag(final_fit$crrFit$var))
coefficients <- data.frame(
  Model_term = names(coef_est),
  Coefficient = as.numeric(coef_est),
  SE = as.numeric(coef_se),
  sHR = exp(as.numeric(coef_est)),
  Lower95 = exp(as.numeric(coef_est) - 1.96 * as.numeric(coef_se)),
  Upper95 = exp(as.numeric(coef_est) + 1.96 * as.numeric(coef_se)),
  P_value = 2 * pnorm(-abs(as.numeric(coef_est) / as.numeric(coef_se)))
)

predictor_table <- data.frame(
  Predictor_domain = c(
    "Recipient age", "Recipient sex", "Primary diagnosis", "Smoking history × transplant type",
    "Ischemic time", "Serum creatinine", "Recipient diabetes", "CMV serostatus",
    "Donor age", "Donor sex", "Donor smoking history", "Donor diabetes", "Donor cause of death"
  ),
  Coding = c(
    "Continuous, per 10 years", "Female / male", "COPD / interstitial lung disease / other",
    "Four clinically explicit joint categories", "Continuous, hours", "Continuous, mg/dL",
    "No / yes", "Negative / positive", "Continuous, per 10 years", "Female / male",
    "No / yes", "No / yes", "Category 1 / 2 / 3 / other"
  ),
  Timing = "Available at transplantation",
  Included_in_clinical_model = "Yes",
  stringsAsFactors = FALSE
)

excluded_table <- data.frame(
  Variable_or_domain = c("Transplant era", "Serum albumin", "Steroid/rejection variables", "LUNG_DONOR"),
  Decision = c("Validation/audit only", "Excluded", "Excluded", "Excluded upstream"),
  Rationale = c(
    "Not an intrinsic patient characteristic; era-specific calibration is audited separately",
    "62.4% missing in the corrected cohort",
    "Measured after transplantation and would introduce temporal leakage",
    "Donor-origin lung-cancer indicator was removed from the analytic cohort definition"
  ),
  stringsAsFactors = FALSE
)

risk_cutoffs_final <- as.numeric(quantile(full_risk[, 1], probs = c(0.60, 0.90), na.rm = TRUE))
model_bundle <- list(
  model = final_fit,
  imputation = imp_final,
  predictor_formula = model_formula,
  horizons_days = horizon_days,
  horizons_years = horizon_years,
  exploratory_3y_risk_cutoffs = risk_cutoffs_final,
  source_file = source_file,
  source_md5 = unname(tools::md5sum(source_file)),
  generated_on = as.character(Sys.Date()),
  interpretation_boundary = "Internal cross-validation only; external validation required before bedside use"
)
saveRDS(model_bundle, file.path(output_dir, "FineGray_prediction_model.rds"))

cv_predictions <- data.frame(
  Row_ID = seq_len(nrow(dat_model)),
  Fold = fold_id,
  FGDAY = dat_model$FGDAY,
  FGSTATUS = dat_model$FGSTATUS,
  Transplant_era = dat_model$TX_YEAR_GROUP,
  Risk_3y = oof_risk[, 1],
  Risk_5y = oof_risk[, 2],
  Risk_10y = oof_risk[, 3],
  Risk_group_3y = dat_model$RISK_GROUP_OOF
)
saveRDS(cv_predictions, file.path(output_dir, "FineGray_crossvalidated_predictions.rds"))

cif_groups <- cmprsk::cuminc(
  ftime = dat_model$FGDAY,
  fstatus = dat_model$FGSTATUS,
  group = dat_model$RISK_GROUP_OOF,
  cencode = 0
)
cif_curve <- do.call(rbind, lapply(grep(" 1$", names(cif_groups), value = TRUE), function(nm) {
  obj <- cif_groups[[nm]]
  data.frame(
    Risk_group = sub(" 1$", "", nm),
    Days = obj$time,
    CIF = obj$est
  )
}))

palette <- c("3" = "#2C7BB6", "5" = "#F28E2B", "10" = "#B2182B")
risk_palette <- c("Low risk" = "#2C7BB6", "Intermediate risk" = "#F28E2B", "High risk" = "#C83E1D")
theme_pub <- theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", color = "#17365D", size = 12),
    plot.subtitle = element_text(color = "#5B677A", size = 9),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", color = "#17365D")
  )

p_roc <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity, color = factor(Horizon_years))) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#8A8A8A") +
  coord_equal() +
  scale_color_manual(values = palette, labels = paste0(horizon_years, " years")) +
  labs(
    title = "A. Cross-validated discrimination",
    subtitle = paste0("AUC: ", paste(sprintf("%dy %.3f", performance$Horizon_years, performance$AUC), collapse = " | ")),
    x = "1 - specificity", y = "Sensitivity"
  ) + theme_pub

p_cal <- ggplot(calibration[calibration$Horizon_years %in% c(3, 5), ],
                aes(x = Predicted_risk, y = Observed_CIF, color = factor(Horizon_years))) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "#8A8A8A") +
  geom_errorbar(aes(ymin = Observed_lower95, ymax = Observed_upper95), width = 0, alpha = 0.75) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.5) +
  scale_color_manual(values = palette, labels = c("3 years", "5 years")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(title = "B. Calibration by risk quintile", x = "Mean predicted risk", y = "Observed cumulative incidence") +
  theme_pub

dca_long <- rbind(
  data.frame(Horizon_years = dca$Horizon_years, Threshold_probability = dca$Threshold_probability, Strategy = "Clinical Fine-Gray model", Net_benefit = dca$Model),
  data.frame(Horizon_years = dca$Horizon_years, Threshold_probability = dca$Threshold_probability, Strategy = "Alert all", Net_benefit = dca$Treat_all),
  data.frame(Horizon_years = dca$Horizon_years, Threshold_probability = dca$Threshold_probability, Strategy = "Alert none", Net_benefit = dca$Treat_none)
)
p_dca <- ggplot(dca_long, aes(x = Threshold_probability, y = Net_benefit, color = Strategy)) +
  geom_line(linewidth = 0.85) +
  facet_wrap(~ paste0(Horizon_years, " years"), scales = "free_x") +
  scale_color_manual(values = c("Clinical Fine-Gray model" = "#2C7BB6", "Alert all" = "#C83E1D", "Alert none" = "#555555")) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.5)) +
  labs(title = "C. Decision-curve analysis", x = "Risk threshold", y = "Net benefit") + theme_pub

p_cif <- ggplot(cif_curve[cif_curve$Days <= 365.25 * 10, ], aes(x = Days / 365.25, y = CIF, color = Risk_group)) +
  geom_step(linewidth = 0.9) +
  scale_color_manual(values = risk_palette) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(title = "D. Cross-validated 3-year risk strata", x = "Years after transplantation", y = "Observed lung-cancer incidence") +
  theme_pub

figure <- (p_roc | p_cal) / (p_dca | p_cif) +
  plot_annotation(
    title = "Internally validated clinical prediction of post-transplant lung cancer",
    subtitle = "Corrected competing-risk outcome; 5-fold out-of-fold predictions; N = 21,104, lung cancers = 459",
    theme = theme(plot.title = element_text(face = "bold", size = 16, color = "#17365D"), plot.subtitle = element_text(size = 10, color = "#5B677A"))
  )
ggsave(file.path(output_dir, "Figure7_FineGray_prediction_performance.png"), figure, width = 13, height = 10, dpi = 400, bg = "white")
ggsave(file.path(output_dir, "Figure7_FineGray_prediction_performance.pdf"), figure, width = 13, height = 10, device = "pdf", bg = "white")

overview <- data.frame(
  Metric = c(
    "Source cohort N", "Lung-cancer events", "Competing deaths", "Censored at last follow-up",
    "Validation design", "Primary prediction horizon", "Key secondary horizon", "Exploratory horizon",
    "Predictor domains", "Post-transplant predictors", "Transplant era in clinical model", "Maximum predictor missingness",
    "Interpretation boundary"
  ),
  Value = c(
    nrow(dat_model), sum(dat_model$FGSTATUS == 1L), sum(dat_model$FGSTATUS == 2L), sum(dat_model$FGSTATUS == 0L),
    "Stratified 5-fold cross-validation with out-of-fold prediction for every recipient",
    "3 years", "5 years", "10 years", nrow(predictor_table), "None", "No; audited separately",
    sprintf("%.2f%%", max(missingness$Missing_percent)),
    "Internal validation only; external validation and clinical-impact evaluation are required before bedside use"
  )
)

result_payload <- list(
  Overview = overview,
  Predictors = predictor_table,
  Excluded_variables = excluded_table,
  Model_coefficients = coefficients,
  Cross_validated_performance = performance,
  Calibration = calibration,
  Risk_groups = risk_groups,
  Decision_curve = dca,
  Era_audit = era_audit,
  Fold_audit = fold_audit,
  Missingness = missingness,
  Methods = data.frame(
    Field = c(
      "Source file", "Source MD5", "Generated on", "Outcome", "Competing event", "Censor code",
      "Variable timing", "Imputation", "Validation", "Performance", "Calibration", "Decision curve",
      "Risk groups", "Era handling", "Primary limitation"
    ),
    Value = c(
      source_file, unname(tools::md5sum(source_file)), as.character(Sys.Date()),
      "First post-transplant lung cancer (FGSTATUS = 1)", "Death before lung cancer (FGSTATUS = 2)",
      "Alive/event-free at last follow-up (FGSTATUS = 0)",
      "All clinical predictors were available at transplantation",
      "Training-fold medians for continuous variables and modes for categorical variables; outcomes were never used for imputation",
      "Five-fold stratified cross-validation; every reported prediction performance estimate uses out-of-fold predictions",
      "IPCW cumulative/dynamic AUC and Brier score; 300 nonparametric bootstrap samples for metric confidence intervals",
      "Observed cumulative incidence within cross-validated risk quintiles",
      "IPCW net benefit at 3 and 5 years; intended as decision-analytic exploration, not a clinical threshold recommendation",
      "Exploratory 60th and 90th percentiles of cross-validated 3-year predicted risk",
      "Not used as a clinical predictor; era-specific 3-year calibration/discrimination audited for transportability",
      "Single-registry internal validation only; no external or prospective validation"
    ),
    stringsAsFactors = FALSE
  )
)
jsonlite::write_json(result_payload, file.path(output_dir, "section7_results.json"), pretty = TRUE, auto_unbox = TRUE, digits = 10, na = "null")

write.csv(roc_data, file.path(output_dir, "roc_plot_data.csv"), row.names = FALSE)
write.csv(cif_curve, file.path(output_dir, "risk_group_cif_data.csv"), row.names = FALSE)

cat("Section 7 analysis complete\n")
print(overview)
print(performance)
print(risk_groups)
print(era_audit)
