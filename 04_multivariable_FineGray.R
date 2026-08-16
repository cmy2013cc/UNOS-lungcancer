#!/usr/bin/env Rscript

local_lib <- Sys.getenv("EJCTS_R_LIB", unset = "E:/UNOS/EJCTS_reanalysis")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

# Clinically interpretable multivariable Fine-Gray model for post-transplant lung cancer.
# The expanded primary model retains clinically interpretable baseline variables from the
# original analysis, while excluding high-missingness and post-transplant variables.
# The smoking/transplant-type effect and transplant-era adjustment are allowed to change with
# log(1 + time). Transplant era is not presented as a clinical risk factor or in the forest plot.

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else file.path("..", "2.characteristic", "2.lung.rds")
output_dir <- if (length(args) >= 2) args[[2]] else "."

required_packages <- c("cmprsk", "dplyr", "ggplot2", "patchwork", "car", "jsonlite", "splines", "mice")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(cmprsk)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(car)
  library(jsonlite)
  library(mice)
})

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data_raw <- readRDS(input_file)
required_vars <- c(
  "FGDAY", "FGSTATUS", "TX_YEAR", "AGE", "GENDER", "CIG_USE", "GROUPING", "TX_TYPE",
  "ISCHTIME", "CREAT_TRR", "DIAB", "CMV_STATUS", "AGE_DON", "GENDER_DON",
  "HIST_CIG_DON", "DIABETES_DON", "COD_CAD_DON"
)
missing_vars <- setdiff(required_vars, names(data_raw))
if (length(missing_vars) > 0) stop("Required variables missing: ", paste(missing_vars, collapse = ", "))
if (any(is.na(data_raw$FGDAY)) || any(is.na(data_raw$FGSTATUS))) stop("FGDAY and FGSTATUS must be complete.")
if (any(!data_raw$FGSTATUS %in% c(0, 1, 2))) stop("FGSTATUS must be coded 0/1/2.")
if (any(data_raw$FGDAY < 0)) stop("FGDAY contains negative values.")

factor_with_labels <- function(x, levels, labels) {
  factor(as.character(x), levels = levels, labels = labels)
}

data_analysis <- data_raw %>%
  mutate(
    FGDAY_YEAR = FGDAY / 365.25,
    AGE10 = AGE / 10,
    AGE_DON10 = AGE_DON / 10,
    GENDER = factor_with_labels(GENDER, c("F", "M"), c("Female", "Male")),
    CIG_USE = factor_with_labels(CIG_USE, c("N", "Y"), c("No", "Yes")),
    GROUPING = factor_with_labels(
      GROUPING, c("A", "D", "Other"),
      c("COPD", "Interstitial lung disease", "Other diagnosis")
    ),
    TX_TYPE = factor_with_labels(TX_TYPE, c("S", "D"), c("Single lung", "Double lung")),
    DIAB = factor_with_labels(DIAB, c("N", "Y"), c("No", "Yes")),
    DIABETES_DON = factor_with_labels(DIABETES_DON, c("N", "Y"), c("No", "Yes")),
    CMV_STATUS = factor_with_labels(CMV_STATUS, c("N", "P"), c("Negative", "Positive")),
    GENDER_DON = factor_with_labels(GENDER_DON, c("F", "M"), c("Female", "Male")),
    HIST_CIG_DON = factor_with_labels(HIST_CIG_DON, c("N", "Y"), c("No", "Yes")),
    COD_CAD_DON = factor_with_labels(
      COD_CAD_DON, c("1", "2", "3", "Other"),
      c("Category 1", "Category 2", "Category 3", "Other")
    ),
    AGE_GROUP = cut(
      AGE, breaks = c(-Inf, 55, 60, 65, 70, Inf), right = FALSE,
      labels = c("<55", "55-59", "60-64", "65-69", ">=70")
    ),
    TX_YEAR_GROUP = cut(
      TX_YEAR, breaks = c(2004, 2010, 2015, 2020, Inf), right = TRUE,
      labels = c("2005-2010", "2011-2015", "2016-2020", "2021-2023")
    ),
    SMOKE_TX_GROUP = factor(
      paste(CIG_USE, TX_TYPE, sep = " + "),
      levels = c(
        "No + Double lung", "Yes + Double lung",
        "No + Single lung", "Yes + Single lung"
      )
    )
  )

model_vars <- c(
  "AGE10", "GENDER", "CIG_USE", "GROUPING", "TX_TYPE", "TX_YEAR_GROUP",
  "ISCHTIME", "CREAT_TRR", "DIAB", "CMV_STATUS", "AGE_DON10", "GENDER_DON",
  "HIST_CIG_DON", "DIABETES_DON", "COD_CAD_DON",
  "FGDAY_YEAR", "FGSTATUS"
)

# Multiple imputation preserves the full cohort and all lung-cancer events. Outcomes are used
# as predictors but are never imputed. Missingness is <=2.1% for every included variable.
imputation_data <- droplevels(data_analysis[, model_vars])
imputation_methods <- mice::make.method(imputation_data)
imputation_methods[vapply(imputation_data, function(x) !anyNA(x), logical(1))] <- ""
imputation_methods[c("FGDAY_YEAR", "FGSTATUS")] <- ""
predictor_matrix <- mice::make.predictorMatrix(imputation_data)
diag(predictor_matrix) <- 0
predictor_matrix[c("FGDAY_YEAR", "FGSTATUS"), ] <- 0

n_imputations <- 20
imputed <- mice::mice(
  imputation_data,
  m = n_imputations,
  maxit = 10,
  method = imputation_methods,
  predictorMatrix = predictor_matrix,
  seed = 20260803,
  printFlag = FALSE
)

completed_data <- lapply(seq_len(n_imputations), function(i) {
  dat <- mice::complete(imputed, i)
  dat$SMOKE_TX_GROUP <- factor(
    paste(dat$CIG_USE, dat$TX_TYPE, sep = " + "),
    levels = c(
      "No + Double lung", "Yes + Double lung",
      "No + Single lung", "Yes + Single lung"
    )
  )
  droplevels(dat)
})
cat("Multiple imputation complete:", n_imputations, "data sets\n")
model_data <- completed_data[[1]]

clinical_formula <- ~ AGE10 + GENDER + GROUPING + SMOKE_TX_GROUP +
  ISCHTIME + CREAT_TRR + DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON +
  HIST_CIG_DON + DIABETES_DON + COD_CAD_DON + TX_YEAR_GROUP
clinical_matrix <- model.matrix(clinical_formula, data = model_data)[, -1, drop = FALSE]
time_varying_columns <- grep("^(SMOKE_TX_GROUP|TX_YEAR_GROUP)", colnames(clinical_matrix))
time_varying_matrix <- clinical_matrix[, time_varying_columns, drop = FALSE]

fit_one_imputation <- function(dat) {
  x <- model.matrix(clinical_formula, data = dat)[, -1, drop = FALSE]
  tv <- grep("^(SMOKE_TX_GROUP|TX_YEAR_GROUP)", colnames(x))
  crr(
    ftime = dat$FGDAY_YEAR,
    fstatus = dat$FGSTATUS,
    cov1 = x,
    cov2 = x[, tv, drop = FALSE],
    tf = function(u) matrix(log1p(u), nrow = length(u), ncol = length(tv)),
    failcode = 1,
    cencode = 0,
    maxiter = 30
  )
}

fit_list <- lapply(seq_along(completed_data), function(i) {
  cat("Fitting imputed Fine-Gray model", i, "of", n_imputations, "\n")
  flush.console()
  fit_one_imputation(completed_data[[i]])
})
if (!all(vapply(fit_list, function(x) isTRUE(x$converged), logical(1)))) {
  stop("At least one imputed Fine-Gray model did not converge.")
}
term_template <- names(fit_list[[1]]$coef)
if (!all(vapply(fit_list, function(x) identical(names(x$coef), term_template), logical(1)))) {
  stop("Coefficient names differ across imputed Fine-Gray models.")
}

coef_matrix <- do.call(cbind, lapply(fit_list, function(x) x$coef))
pooled_beta <- rowMeans(coef_matrix)
within_variance <- Reduce("+", lapply(fit_list, function(x) x$var)) / n_imputations
between_variance <- stats::cov(t(coef_matrix))
total_variance <- within_variance + (1 + 1 / n_imputations) * between_variance
fit_final <- list(
  coef = setNames(as.numeric(pooled_beta), term_template),
  var = total_variance,
  converged = TRUE,
  pooling = list(
    m = n_imputations,
    within_variance = within_variance,
    between_variance = between_variance,
    total_variance = total_variance
  )
)

fit_proportional_joint <- crr(
  ftime = model_data$FGDAY_YEAR,
  fstatus = model_data$FGSTATUS,
  cov1 = clinical_matrix,
  failcode = 1,
  cencode = 0
)

coef_stats <- data.frame(
  Term = names(fit_final$coef),
  Beta = as.numeric(fit_final$coef),
  SE = sqrt(diag(fit_final$var)),
  stringsAsFactors = FALSE
) %>% mutate(
    Within_variance = diag(within_variance),
    Between_variance = diag(between_variance),
    Rubin_df = ifelse(
      Between_variance <= .Machine$double.eps,
      Inf,
      (n_imputations - 1) *
        (1 + Within_variance / ((1 + 1 / n_imputations) * Between_variance))^2
    ),
    sHR = exp(Beta),
    Critical_value = ifelse(is.finite(Rubin_df), qt(0.975, Rubin_df), qnorm(0.975)),
    Lower_95CI = exp(Beta - Critical_value * SE),
    Upper_95CI = exp(Beta + Critical_value * SE),
    P_value = ifelse(
      is.finite(Rubin_df),
      2 * pt(abs(Beta / SE), Rubin_df, lower.tail = FALSE),
      2 * pnorm(abs(Beta / SE), lower.tail = FALSE)
    )
  )

wald_statistic <- function(fit, indices) {
  beta <- fit$coef[indices]
  variance <- fit$var[indices, indices, drop = FALSE]
  drop(t(beta) %*% solve(variance, beta))
}

wald_p <- function(fit, indices) {
  pchisq(wald_statistic(fit, indices), df = length(indices), lower.tail = FALSE)
}

term_names <- names(fit_final$coef)
block_indices <- list(
  "Recipient age" = grep("^AGE10$", term_names),
  "Recipient sex" = grep("^GENDERMale$", term_names),
  "Primary diagnosis" = grep("^GROUPING", term_names),
  "Smoking history and transplant type" = grep("^SMOKE_TX_GROUP", term_names),
  "Ischemic time" = grep("^ISCHTIME$", term_names),
  "Serum creatinine" = grep("^CREAT_TRR$", term_names),
  "Recipient diabetes" = grep("^DIABYes$", term_names),
  "CMV serostatus" = grep("^CMV_STATUS", term_names),
  "Donor age" = grep("^AGE_DON10$", term_names),
  "Donor sex" = grep("^GENDER_DON", term_names),
  "Donor smoking history" = grep("^HIST_CIG_DON", term_names),
  "Donor diabetes" = grep("^DIABETES_DON", term_names),
  "Donor cause of death" = grep("^COD_CAD_DON", term_names),
  "Transplant era (adjustment only)" = grep("^TX_YEAR_GROUP", term_names)
)
global_tests <- data.frame(
  Variable = names(block_indices),
  Degrees_of_freedom = vapply(block_indices, length, integer(1)),
  Overall_P_value = vapply(block_indices, function(idx) wald_p(fit_final, idx), numeric(1)),
  Role = c(rep("Clinical variable", 13), "Nuisance adjustment; not a clinical risk item"),
  stringsAsFactors = FALSE
)

format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

get_stat <- function(term) {
  row <- coef_stats[coef_stats$Term == term, ]
  if (nrow(row) != 1) stop("Expected exactly one coefficient for term: ", term)
  row
}

time_specific_effect <- function(main_term, time_years) {
  main_index <- match(main_term, colnames(clinical_matrix))
  tv_position <- match(main_index, time_varying_columns)
  if (is.na(main_index) || is.na(tv_position)) stop("Time-varying term not found: ", main_term)
  gamma_index <- ncol(clinical_matrix) + tv_position
  log_time <- log1p(time_years)
  beta_time <- fit_final$coef[[main_index]] + log_time * fit_final$coef[[gamma_index]]
  combination_variance <- function(v) {
    v[main_index, main_index] + log_time^2 * v[gamma_index, gamma_index] +
      2 * log_time * v[main_index, gamma_index]
  }
  within_time <- combination_variance(within_variance)
  between_time <- combination_variance(between_variance)
  variance_time <- combination_variance(total_variance)
  se_time <- sqrt(variance_time)
  rubin_df <- if (between_time <= .Machine$double.eps) Inf else
    (n_imputations - 1) *
      (1 + within_time / ((1 + 1 / n_imputations) * between_time))^2
  critical_value <- if (is.finite(rubin_df)) qt(0.975, rubin_df) else qnorm(0.975)
  p_value <- if (is.finite(rubin_df)) {
    2 * pt(abs(beta_time / se_time), rubin_df, lower.tail = FALSE)
  } else {
    2 * pnorm(abs(beta_time / se_time), lower.tail = FALSE)
  }
  data.frame(
    Term = main_term,
    Time_years = time_years,
    Beta = beta_time,
    SE = se_time,
    sHR = exp(beta_time),
    Lower_95CI = exp(beta_time - critical_value * se_time),
    Upper_95CI = exp(beta_time + critical_value * se_time),
    P_value = p_value,
    stringsAsFactors = FALSE
  )
}

joint_term_map <- c(
  "Smoking history + double lung" = "SMOKE_TX_GROUPYes + Double lung",
  "No smoking history + single lung" = "SMOKE_TX_GROUPNo + Single lung",
  "Smoking history + single lung" = "SMOKE_TX_GROUPYes + Single lung"
)
joint_time_effects <- bind_rows(lapply(c(1, 3, 5, 10), function(time_years) {
  bind_rows(lapply(unname(joint_term_map), time_specific_effect, time_years = time_years)) %>%
    mutate(Category = names(joint_term_map)[match(Term, joint_term_map)])
})) %>%
  select(Category, Term, Time_years, Beta, SE, sHR, Lower_95CI, Upper_95CI, P_value)

row_spec <- data.frame(
  Variable = c(
    "Recipient age", rep("Recipient sex", 2), rep("Primary diagnosis", 3),
    rep("Smoking history and transplant type", 4),
    "Ischemic time", "Serum creatinine", rep("Recipient diabetes", 2),
    rep("CMV serostatus", 2), "Donor age", rep("Donor sex", 2),
    rep("Donor smoking history", 2), rep("Donor diabetes", 2),
    rep("Donor cause of death", 4)
  ),
  Category = c(
    "Per 10-year increase",
    "Female", "Male",
    "COPD", "Interstitial lung disease", "Other diagnosis",
    "No smoking history + double lung", "Smoking history + double lung",
    "No smoking history + single lung", "Smoking history + single lung",
    "Per 1-hour increase", "Per 1-mg/dL increase",
    "No", "Yes", "Negative", "Positive", "Per 10-year increase",
    "Female", "Male", "No", "Yes", "No", "Yes",
    "Category 1", "Category 2", "Category 3", "Other"
  ),
  Term = c(
    "AGE10",
    NA, "GENDERMale",
    NA, "GROUPINGInterstitial lung disease", "GROUPINGOther diagnosis",
    NA, "SMOKE_TX_GROUPYes + Double lung", "SMOKE_TX_GROUPNo + Single lung",
    "SMOKE_TX_GROUPYes + Single lung",
    "ISCHTIME", "CREAT_TRR", NA, "DIABYes", NA, "CMV_STATUSPositive",
    "AGE_DON10", NA, "GENDER_DONMale", NA, "HIST_CIG_DONYes",
    NA, "DIABETES_DONYes", NA, "COD_CAD_DONCategory 2",
    "COD_CAD_DONCategory 3", "COD_CAD_DONOther"
  ),
  Is_reference = c(
    FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE,
    FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE,
    FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)

table4 <- row_spec %>%
  mutate(
    sHR = ifelse(Is_reference, 1, NA_real_),
    Lower_95CI = NA_real_,
    Upper_95CI = NA_real_,
    P_value = NA_real_
  )

for (i in seq_len(nrow(table4))) {
  if (!table4$Is_reference[[i]]) {
    if (grepl("^SMOKE_TX_GROUP", table4$Term[[i]])) {
      stat <- joint_time_effects %>%
        filter(Term == table4$Term[[i]], Time_years == 5)
    } else {
      stat <- get_stat(table4$Term[[i]])
    }
    table4$sHR[[i]] <- stat$sHR
    table4$Lower_95CI[[i]] <- stat$Lower_95CI
    table4$Upper_95CI[[i]] <- stat$Upper_95CI
    table4$P_value[[i]] <- stat$P_value
  }
}

clinical_rows <- global_tests$Variable != "Transplant era (adjustment only)"
clinical_global <- setNames(global_tests$Overall_P_value[clinical_rows], global_tests$Variable[clinical_rows])
table4$Overall_P_value <- NA_real_
for (variable in unique(table4$Variable)) {
  first_row <- which(table4$Variable == variable)[[1]]
  table4$Overall_P_value[[first_row]] <- clinical_global[[variable]]
}

table4 <- table4 %>%
  mutate(
    Reference = ifelse(Is_reference, "Reference", ""),
    `Adjusted sHR (95% CI)` = ifelse(
      Is_reference,
      "Reference",
      sprintf("%.2f (%.2f-%.2f)", sHR, Lower_95CI, Upper_95CI)
    ),
    Effect_time = ifelse(Variable == "Smoking history and transplant type", "5 years", "Constant effect"),
    `P value` = vapply(P_value, format_p, character(1)),
    `Overall P value` = vapply(Overall_P_value, format_p, character(1))
  )

# Pooled global test of the four-level smoking-history/transplant-type construct,
# including its time-varying coefficients in the expanded primary model.
interaction_indices <- block_indices[["Smoking history and transplant type"]]
interaction_df <- length(interaction_indices)
interaction_stat <- wald_statistic(fit_final, interaction_indices)
interaction_p <- pchisq(interaction_stat, df = interaction_df, lower.tail = FALSE)

# Age functional-form assessment within the same adjusted model.
spline_age_matrix <- model.matrix(
  ~ splines::ns(AGE10, df = 3) + GENDER + GROUPING + SMOKE_TX_GROUP +
    ISCHTIME + CREAT_TRR + DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON +
    HIST_CIG_DON + DIABETES_DON + COD_CAD_DON + TX_YEAR_GROUP,
  data = model_data
)[, -1, drop = FALSE]
spline_tv_columns <- grep("^(SMOKE_TX_GROUP|TX_YEAR_GROUP)", colnames(spline_age_matrix))
spline_tv_matrix <- spline_age_matrix[, spline_tv_columns, drop = FALSE]
fit_spline_age <- crr(
  model_data$FGDAY_YEAR, model_data$FGSTATUS,
  cov1 = spline_age_matrix,
  cov2 = spline_tv_matrix,
  tf = function(u) matrix(log1p(u), nrow = length(u), ncol = ncol(spline_tv_matrix)),
  failcode = 1, cencode = 0, maxiter = 20
)
age_spline_lrt <- 2 * (fit_spline_age$loglik - fit_list[[1]]$loglik)
age_spline_p <- pchisq(age_spline_lrt, df = 2, lower.tail = FALSE)

# Pooled Wald test for the additional baseline adjustment block. These variables are retained
# for confounding control regardless of individual or joint statistical significance.
extended_variable_names <- c(
  "Ischemic time", "Serum creatinine", "Recipient diabetes", "CMV serostatus",
  "Donor age", "Donor sex", "Donor smoking history", "Donor diabetes",
  "Donor cause of death"
)
extended_indices <- unlist(block_indices[extended_variable_names], use.names = FALSE)
extended_df <- length(extended_indices)
extended_stat <- wald_statistic(fit_final, extended_indices)
extended_p <- pchisq(extended_stat, df = extended_df, lower.tail = FALSE)

# Collinearity diagnostic using adjusted GVIF.
vif_model <- lm(
  FGDAY_YEAR ~ AGE10 + GENDER + GROUPING + SMOKE_TX_GROUP + ISCHTIME +
    CREAT_TRR + DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON + HIST_CIG_DON +
    DIABETES_DON + COD_CAD_DON + TX_YEAR_GROUP,
  data = model_data
)
vif_raw <- car::vif(vif_model)
if (is.matrix(vif_raw)) {
  vif_table <- data.frame(
    Variable = rownames(vif_raw),
    GVIF = vif_raw[, "GVIF"],
    Degrees_of_freedom = vif_raw[, "Df"],
    Adjusted_GVIF = vif_raw[, "GVIF^(1/(2*Df))"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
} else {
  vif_table <- data.frame(
    Variable = names(vif_raw), GVIF = as.numeric(vif_raw), Degrees_of_freedom = 1,
    Adjusted_GVIF = sqrt(as.numeric(vif_raw)), stringsAsFactors = FALSE
  )
}

# Proportional subdistribution-hazards diagnostic using log(1 + time) interactions.
p_columns <- ncol(clinical_matrix)
fit_time_varying <- crr(
  ftime = model_data$FGDAY_YEAR,
  fstatus = model_data$FGSTATUS,
  cov1 = clinical_matrix,
  cov2 = clinical_matrix,
  tf = function(u) matrix(log1p(u), nrow = length(u), ncol = p_columns),
  failcode = 1,
  cencode = 0,
  maxiter = 20
)
if (!isTRUE(fit_time_varying$converged)) warning("Time-varying diagnostic model did not converge.")

time_indices <- seq_len(p_columns) + p_columns
time_names <- paste0(colnames(clinical_matrix), " x log(1 + time)")
time_beta <- fit_time_varying$coef[time_indices]
time_se <- sqrt(diag(fit_time_varying$var))[time_indices]
ph_term_tests <- data.frame(
  Term = time_names,
  Time_interaction_Beta = as.numeric(time_beta),
  SE = as.numeric(time_se),
  P_value = 2 * pnorm(abs(time_beta / time_se), lower.tail = FALSE),
  stringsAsFactors = FALSE
)

screen_block_main_indices <- list(
  "Recipient age" = grep("^AGE10$", colnames(clinical_matrix)),
  "Recipient sex" = grep("^GENDERMale$", colnames(clinical_matrix)),
  "Primary diagnosis" = grep("^GROUPING", colnames(clinical_matrix)),
  "Smoking history and transplant type" = grep("^SMOKE_TX_GROUP", colnames(clinical_matrix)),
  "Ischemic time" = grep("^ISCHTIME$", colnames(clinical_matrix)),
  "Serum creatinine" = grep("^CREAT_TRR$", colnames(clinical_matrix)),
  "Recipient diabetes" = grep("^DIABYes$", colnames(clinical_matrix)),
  "CMV serostatus" = grep("^CMV_STATUS", colnames(clinical_matrix)),
  "Donor age" = grep("^AGE_DON10$", colnames(clinical_matrix)),
  "Donor sex" = grep("^GENDER_DON", colnames(clinical_matrix)),
  "Donor smoking history" = grep("^HIST_CIG_DON", colnames(clinical_matrix)),
  "Donor diabetes" = grep("^DIABETES_DON", colnames(clinical_matrix)),
  "Donor cause of death" = grep("^COD_CAD_DON", colnames(clinical_matrix)),
  "Transplant era (adjustment only)" = grep("^TX_YEAR_GROUP", colnames(clinical_matrix))
)
ph_block_indices <- lapply(screen_block_main_indices, function(idx) idx + p_columns)
ph_global_tests <- data.frame(
  Variable = names(ph_block_indices),
  Degrees_of_freedom = vapply(ph_block_indices, length, integer(1)),
  Global_time_interaction_P = vapply(
    ph_block_indices,
    function(idx) wald_p(fit_time_varying, idx),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)

# Internal era-influence check: era is retained as adjustment because clinical coefficients shift.
no_era_matrix <- model.matrix(
  ~ AGE10 + GENDER + GROUPING + SMOKE_TX_GROUP + ISCHTIME + CREAT_TRR +
    DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON + HIST_CIG_DON +
    DIABETES_DON + COD_CAD_DON,
  data = model_data
)[, -1, drop = FALSE]
fit_no_era <- crr(model_data$FGDAY_YEAR, model_data$FGSTATUS, cov1 = no_era_matrix,
                  failcode = 1, cencode = 0)
common_terms <- intersect(names(fit_no_era$coef), names(fit_proportional_joint$coef))
era_influence <- data.frame(
  Term = common_terms,
  sHR_without_era_adjustment = exp(fit_no_era$coef[common_terms]),
  sHR_with_era_adjustment = exp(fit_proportional_joint$coef[common_terms]),
  Percent_change = 100 * (
    exp(fit_proportional_joint$coef[common_terms]) - exp(fit_no_era$coef[common_terms])
  ) / exp(fit_no_era$coef[common_terms]),
  stringsAsFactors = FALSE
)
max_era_change <- max(abs(era_influence$Percent_change), na.rm = TRUE)

n_total <- nrow(data_analysis)
n_model <- nrow(model_data)
n_excluded <- n_total - n_model
n_events <- sum(model_data$FGSTATUS == 1)
n_deaths <- sum(model_data$FGSTATUS == 2)
n_censored <- sum(model_data$FGSTATUS == 0)
n_parameters <- length(fit_final$coef)

model_summary <- data.frame(
  Metric = c(
    "Source cohort", "Analysis population after multiple imputation",
    "Excluded because of missing model variables", "Number of imputations",
    "Lung cancer events", "Competing deaths", "Censored at last follow-up", "Estimated parameters",
    "Events per estimated parameter", "Joint smoking/transplant construct, pooled global P value",
    "Age nonlinearity P value (first imputed data set)",
    "Expanded baseline adjustment block, pooled global P value",
    "Maximum adjusted GVIF", "Maximum absolute coefficient change after era adjustment",
    "Model purpose"
  ),
  Value = c(
    n_total, n_model, n_excluded, n_imputations, n_events, n_deaths, n_censored, n_parameters,
    round(n_events / n_parameters, 1), format_p(interaction_p), format_p(age_spline_p),
    format_p(extended_p), round(max(vif_table$Adjusted_GVIF), 3),
    sprintf("%.1f%%", max_era_change),
    "Adjusted association model for clinical risk awareness; not a validated bedside calculator"
  ),
  stringsAsFactors = FALSE
)

variable_decisions <- data.frame(
  Variable_or_group = c(
    "Recipient age", "Recipient sex", "Primary diagnosis",
    "Recipient smoking history x transplant type", "Transplant era",
    "Ischemic time", "Serum creatinine", "Recipient diabetes", "CMV serostatus",
    "Donor age", "Donor smoking history", "Donor diabetes", "Donor sex", "Donor cause of death",
    "Serum albumin", "Post-transplant rejection and immunosuppression variables"
  ),
  Decision = c(
    "Included linearly per 10-year increase", "Included", "Included",
    "Included as a time-varying four-level joint variable", "Included only as time-varying nuisance adjustment",
    rep("Included as baseline adjustment", 9),
    "Excluded from primary model", "Excluded from primary model"
  ),
  Rationale = c(
    "No material nonlinearity was detected in the diagnostic data set; a per-10-year effect preserves information.",
    "Baseline, complete, clinically interpretable, and independently associated.",
    "Baseline clinical context retained regardless of individual coefficient significance.",
    "Clinically plausible interaction and statistically supported; its effect changes over time and is reported at clinically relevant horizons.",
    "Adjusted for secular changes with a time-varying effect because omission materially changed clinical coefficients; not presented as a risk item.",
    "Retained from the original baseline model as a transplant-procedure adjustment variable; low missingness was handled by multiple imputation.",
    "Retained from the original baseline model as a marker of recipient renal status; low missingness was handled by multiple imputation.",
    "Retained from the original baseline model as a clinically relevant recipient comorbidity.",
    "Retained from the original baseline model as transplant infectious-risk context.",
    "Retained from the original baseline model as a donor characteristic, modeled per 10 years.",
    "Retained from the original baseline model as donor exposure history.",
    "Retained from the original baseline model as a donor comorbidity.",
    "Added from the prespecified baseline candidate set and the prior P<0.20 screen.",
    "Retained from the original baseline model for donor-context adjustment; category 1 is the reference.",
    "62.4% missing; including it would require heavy extrapolation and potentially unstable estimates.",
    "Measured after transplantation and therefore unsuitable as baseline predictors in this model."
  ),
  stringsAsFactors = FALSE
)

references <- data.frame(
  Source = c(
    "Fine JP, Gray RJ. JASA. 1999.",
    "Austin PC, Fine JP. Statistics in Medicine. 2017.",
    "Collins GS et al. BMJ. 2015. TRIPOD statement."
  ),
  URL = c(
    "https://doi.org/10.1080/01621459.1999.10474144",
    "https://doi.org/10.1002/sim.7501",
    "https://doi.org/10.1136/bmj.g7594"
  ),
  Purpose = c(
    "Fine-Gray proportional subdistribution hazards model.",
    "Reporting and interpretation of Fine-Gray analyses.",
    "Transparent reporting and validation expectations for prediction models."
  ),
  stringsAsFactors = FALSE
)

# Forest plot: clinical terms only; era is deliberately hidden from the clinical display.
plot_data <- table4 %>%
  mutate(
    Row = rev(seq_len(n())),
    Group_start = !duplicated(Variable),
    Left_label = ifelse(
      Group_start,
      paste0(Variable, "\n  ", Category, ifelse(Is_reference, " (reference)", "")),
      paste0("  ", Category)
    ),
    Estimate_label = `Adjusted sHR (95% CI)`,
    P_label = ifelse(
      Is_reference,
      ifelse(`Overall P value` == "", "", paste0("Overall ", `Overall P value`)),
      `P value`
    ),
    Point_color = ifelse(
      Category == "Smoking history + single lung", "#B2182B",
      ifelse(
        grepl("single lung", Category, fixed = TRUE) | Category == "Smoking history + double lung",
        "#D6604D", "#2166AC"
      )
    )
  )

left_panel <- ggplot(plot_data, aes(x = 0, y = Row, label = Left_label)) +
  geom_text(aes(fontface = ifelse(Group_start, "bold", "plain")), hjust = 0, size = 3.45) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 4, 28, 4))

forest_panel <- ggplot(plot_data, aes(y = Row)) +
  geom_vline(xintercept = 1, color = "#777777", linewidth = 0.6, linetype = "dashed") +
  geom_segment(
    data = plot_data %>% filter(!Is_reference),
    aes(x = Lower_95CI, xend = Upper_95CI, yend = Row),
    linewidth = 0.8, color = "#4D4D4D"
  ) +
  geom_point(aes(x = sHR, color = Point_color, shape = Is_reference), size = 2.8, stroke = 0.8) +
  scale_color_identity() +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 15)) +
  scale_x_log10(
    limits = c(0.25, 40),
    breaks = c(0.25, 0.5, 1, 2, 5, 10, 20),
    labels = c("0.25", "0.5", "1", "2", "5", "10", "20")
  ) +
  labs(x = "Adjusted subdistribution hazard ratio (log scale)", y = NULL) +
  theme_classic(base_size = 10.5) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "none", plot.margin = margin(8, 8, 28, 8)
  )

right_panel <- ggplot(plot_data, aes(y = Row)) +
  geom_text(aes(x = 0, label = Estimate_label), hjust = 0, size = 3.25) +
  geom_text(aes(x = 1, label = P_label), hjust = 1, size = 3.25) +
  annotate("text", x = 0, y = max(plot_data$Row) + 0.75, label = "Adjusted sHR (95% CI)",
           hjust = 0, fontface = "bold", size = 3.25) +
  annotate("text", x = 1, y = max(plot_data$Row) + 0.75, label = "P / overall P",
           hjust = 1, fontface = "bold", size = 3.25) +
  coord_cartesian(xlim = c(0, 1), ylim = c(min(plot_data$Row) - 0.5, max(plot_data$Row) + 1), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(8, 4, 28, 4))

forest_figure <- left_panel + forest_panel + right_panel +
  plot_layout(widths = c(3.15, 3.6, 2.35)) +
  plot_annotation(
    title = "Independent associations with post-transplant lung cancer",
    subtitle = paste0(
      "Fine-Gray model (N = ", format(n_model, big.mark = ","),
      "; lung cancer events = ", n_events,
      "); death treated as a competing event; joint smoking/transplant effect shown at 5 years"
    ),
    caption = "Transplant era was modeled as a time-varying nuisance adjustment and is not displayed as a clinical risk factor.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, color = "#17365D"),
      plot.subtitle = element_text(size = 10.5, color = "#4D4D4D"),
      plot.caption = element_text(size = 9, color = "#595959", hjust = 0)
    )
  )

ggsave(file.path(output_dir, "Figure_multivariable_FineGray_forest.pdf"), forest_figure,
       width = 11.2, height = 12.5, device = "pdf")
ggsave(file.path(output_dir, "Figure_multivariable_FineGray_forest.png"), forest_figure,
       width = 11.2, height = 12.5, dpi = 320, bg = "white")

model_object <- list(
  fit = fit_final,
  clinical_matrix_columns = colnames(clinical_matrix),
  time_varying_matrix_columns = colnames(time_varying_matrix),
  preprocessing = list(
    age_effect = "Linear per 10-year increase",
    smoking_transplant_reference = "No smoking history + double lung",
    transplant_era_adjustment = c("2005-2010", "2011-2015", "2016-2020", "2021-2023"),
    missing_data = paste0(n_imputations, " multiple imputations with Rubin pooling")
  ),
  sample = list(N = n_model, lung_cancer_events = n_events, competing_deaths = n_deaths,
                censored = n_censored),
  source_file = normalizePath(input_file, winslash = "/", mustWork = TRUE),
  generated_on = as.character(Sys.Date())
)
saveRDS(model_object, file.path(output_dir, "main_FineGray_model.rds"))

metadata <- data.frame(
  Field = c(
    "Source file", "Source MD5", "Generated on", "Cohort restriction", "Outcome",
    "Competing event", "Missing-data strategy", "Age functional form",
    "Smoking/transplant-type handling", "Era handling", "Interpretation boundary"
  ),
  Value = c(
    normalizePath(input_file, winslash = "/", mustWork = TRUE),
    unname(tools::md5sum(input_file)), as.character(Sys.Date()),
    "LAS era; TX_DATE >= 2005-05-04 applied upstream",
    "First post-transplant lung cancer (FGSTATUS = 1)",
    "Death before lung cancer (FGSTATUS = 2)",
    paste0("Multiple imputation with ", n_imputations,
           " data sets and Rubin pooling; full cohort and all lung-cancer events retained"),
    "Linear per 10 years; no residual nonlinearity in the first imputed diagnostic data set",
    "Four-level joint variable with log(1 + time) interactions",
    "Time-varying nuisance adjustment only; omitted from the clinical forest plot",
    "Association model for risk awareness; external validation is required before use as a bedside calculator"
  ),
  stringsAsFactors = FALSE
)

json_payload <- list(
  Table4 = table4 %>%
    select(Variable, Category, Reference, Effect_time, `Adjusted sHR (95% CI)`, `P value`, `Overall P value`,
           sHR, Lower_95CI, Upper_95CI, P_value, Overall_P_value),
  Time_specific_joint_effects = joint_time_effects %>%
    mutate(
      `Adjusted sHR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", sHR, Lower_95CI, Upper_95CI),
      `P value` = vapply(P_value, format_p, character(1))
    ),
  Global_tests = global_tests %>% mutate(Overall_P_value_formatted = vapply(Overall_P_value, format_p, character(1))),
  Model_summary = model_summary,
  Variable_decisions = variable_decisions,
  PH_global_tests = ph_global_tests %>%
    mutate(Global_time_interaction_P_formatted = vapply(Global_time_interaction_P, format_p, character(1))),
  PH_term_tests = ph_term_tests %>% mutate(P_value_formatted = vapply(P_value, format_p, character(1))),
  Collinearity = vif_table,
  Era_adjustment_summary = data.frame(
    Metric = c("Maximum absolute sHR change", "Decision"),
    Value = c(sprintf("%.1f%%", max_era_change),
              "Retain era as nuisance adjustment; do not present it as a clinical risk factor"),
    stringsAsFactors = FALSE
  ),
  Interaction_and_functional_form = data.frame(
    Test = c(
      "Joint smoking history and transplant type (pooled global Wald test)",
      "Age natural spline vs linear",
      "Expanded baseline adjustment block (pooled Wald test)"
    ),
    Degrees_of_freedom = c(interaction_df, 2, extended_df),
    Statistic = c(interaction_stat, age_spline_lrt, extended_stat),
    P_value = c(interaction_p, age_spline_p, extended_p),
    P_value_formatted = vapply(c(interaction_p, age_spline_p, extended_p), format_p, character(1)),
    stringsAsFactors = FALSE
  ),
  Joint_group_counts = as.data.frame(with(model_data, table(SMOKE_TX_GROUP, FGSTATUS))) %>%
    rename(Group = SMOKE_TX_GROUP, Status = FGSTATUS, N = Freq),
  Metadata = metadata,
  References = references
)

write_json(
  json_payload,
  path = file.path(output_dir, "multivariable_results.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  dataframe = "rows",
  na = "null",
  digits = NA
)

cat("Clinically interpretable multivariable Fine-Gray analysis complete\n")
cat("Source:", normalizePath(input_file, winslash = "/", mustWork = TRUE), "\n")
cat("N:", n_model, "| Lung cancer:", n_events, "| Competing deaths:", n_deaths,
    "| Censored:", n_censored, "\n")
cat("Interaction P:", format_p(interaction_p), "| Age nonlinearity P:", format_p(age_spline_p),
    "| Expanded baseline block P:", format_p(extended_p), "\n")
cat("Outputs:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
