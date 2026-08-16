#!/usr/bin/env Rscript

local_lib <- Sys.getenv("EJCTS_R_LIB", unset = "E:/UNOS/EJCTS_reanalysis")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

# Sensitivity analyses for post-transplant lung cancer after lung transplantation.
# The primary analysis is the 90-day exclusion cohort created upstream. This script does not
# repeat that exclusion. Instead, it tests stricter 180- and 365-day exclusions, alternative
# competing-risk estimands, complete-case analysis, and alternative transplant-era handling.

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else file.path("..", "2.characteristic", "2.lung.rds")
output_dir <- if (length(args) >= 2) args[[2]] else "."

required_packages <- c("cmprsk", "survival", "dplyr", "ggplot2", "patchwork", "mice", "jsonlite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(cmprsk)
  library(survival)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(mice)
  library(jsonlite)
})

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
input_md5 <- unname(tools::md5sum(input_file))

group_levels <- c(
  "No smoking history + double lung",
  "Smoking history + double lung",
  "No smoking history + single lung",
  "Smoking history + single lung"
)
group_short <- c(
  "No smoking history + double lung" = "No smoking + double lung",
  "Smoking history + double lung" = "Smoking + double lung",
  "No smoking history + single lung" = "No smoking + single lung",
  "Smoking history + single lung" = "Smoking + single lung"
)
group_colors <- c(
  "No smoking history + double lung" = "#2C7BB6",
  "Smoking history + double lung" = "#00A6A6",
  "No smoking history + single lung" = "#F28E2B",
  "Smoking history + single lung" = "#C73E1D"
)

factor_with_labels <- function(x, levels, labels) {
  factor(as.character(x), levels = levels, labels = labels)
}

make_joint_group <- function(cig, tx) {
  factor(
    dplyr::case_when(
      cig == "No" & tx == "Double lung" ~ group_levels[[1]],
      cig == "Yes" & tx == "Double lung" ~ group_levels[[2]],
      cig == "No" & tx == "Single lung" ~ group_levels[[3]],
      cig == "Yes" & tx == "Single lung" ~ group_levels[[4]],
      TRUE ~ NA_character_
    ),
    levels = group_levels
  )
}

data_raw <- readRDS(input_file)
required_vars <- c(
  "FGDAY", "FGSTATUS", "TX_YEAR", "AGE", "GENDER", "CIG_USE", "GROUPING", "TX_TYPE",
  "ISCHTIME", "CREAT_TRR", "DIAB", "CMV_STATUS", "AGE_DON", "GENDER_DON",
  "HIST_CIG_DON", "DIABETES_DON", "COD_CAD_DON"
)
missing_vars <- setdiff(required_vars, names(data_raw))
if (length(missing_vars) > 0) stop("Required variables missing: ", paste(missing_vars, collapse = ", "))
if (anyNA(data_raw$FGDAY) || anyNA(data_raw$FGSTATUS)) stop("FGDAY and FGSTATUS must be complete.")
if (any(!data_raw$FGSTATUS %in% c(0, 1, 2))) stop("FGSTATUS must be coded 0/1/2.")
if (any(data_raw$FGDAY < 0)) stop("FGDAY contains negative values.")
if (any(data_raw$FGSTATUS == 1 & data_raw$FGDAY <= 90)) {
  stop("The input unexpectedly contains lung-cancer events within 90 days; upstream exclusion must be checked.")
}

data_analysis <- data_raw %>%
  mutate(
    ROW_ID = seq_len(n()),
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
    TX_YEAR_GROUP = cut(
      TX_YEAR, breaks = c(2004, 2010, 2015, 2020, Inf), right = TRUE,
      labels = c("2005-2010", "2011-2015", "2016-2020", "2021-2023")
    ),
    SMOKE_TX_GROUP = make_joint_group(CIG_USE, TX_TYPE)
  )

model_vars <- c(
  "AGE10", "GENDER", "CIG_USE", "GROUPING", "TX_TYPE", "TX_YEAR_GROUP",
  "ISCHTIME", "CREAT_TRR", "DIAB", "CMV_STATUS", "AGE_DON10", "GENDER_DON",
  "HIST_CIG_DON", "DIABETES_DON", "COD_CAD_DON", "FGDAY_YEAR", "FGSTATUS"
)
imputation_data <- droplevels(data_analysis[, model_vars])

missingness <- data.frame(
  Variable = names(imputation_data),
  Missing_n = vapply(imputation_data, function(x) sum(is.na(x)), integer(1)),
  Missing_percent = 100 * vapply(imputation_data, function(x) mean(is.na(x)), numeric(1)),
  stringsAsFactors = FALSE
) %>% arrange(desc(Missing_percent), Variable)

n_imputations <- 20L
mi_cache_file <- Sys.getenv("UNOS6_MI_CACHE_FILE", unset = "")
if (nzchar(mi_cache_file) && file.exists(mi_cache_file)) {
  mi_cache <- readRDS(mi_cache_file)
  if (!identical(mi_cache$input_md5, input_md5) || length(mi_cache$completed_data) != n_imputations) {
    stop("The multiple-imputation cache does not match the input data.")
  }
  completed_data <- mi_cache$completed_data
} else {
  imputation_methods <- mice::make.method(imputation_data)
  imputation_methods[vapply(imputation_data, function(x) !anyNA(x), logical(1))] <- ""
  imputation_methods[c("FGDAY_YEAR", "FGSTATUS")] <- ""
  predictor_matrix <- mice::make.predictorMatrix(imputation_data)
  diag(predictor_matrix) <- 0
  predictor_matrix[c("FGDAY_YEAR", "FGSTATUS"), ] <- 0
  imputed <- mice::mice(
    imputation_data, m = n_imputations, maxit = 10,
    method = imputation_methods, predictorMatrix = predictor_matrix,
    seed = 20260803, printFlag = FALSE
  )
  completed_data <- lapply(seq_len(n_imputations), function(i) {
    dat <- mice::complete(imputed, i)
    dat$SMOKE_TX_GROUP <- make_joint_group(dat$CIG_USE, dat$TX_TYPE)
    droplevels(dat)
  })
}
cat("Multiple imputation available:", length(completed_data), "completed data sets\n")

# Harmonize caches created by Section 5.
completed_data <- lapply(completed_data, function(dat) {
  dat$SMOKE_TX_GROUP <- make_joint_group(dat$CIG_USE, dat$TX_TYPE)
  droplevels(dat)
})

base_formula <- ~ AGE10 + GENDER + GROUPING + SMOKE_TX_GROUP +
  ISCHTIME + CREAT_TRR + DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON +
  HIST_CIG_DON + DIABETES_DON + COD_CAD_DON + TX_YEAR_GROUP
no_era_formula <- ~ AGE10 + GENDER + GROUPING + SMOKE_TX_GROUP +
  ISCHTIME + CREAT_TRR + DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON +
  HIST_CIG_DON + DIABETES_DON + COD_CAD_DON

tv_joint_and_era <- function(nms) which(grepl("^(SMOKE_TX_GROUP|TX_YEAR_GROUP)", nms))
tv_joint_only <- function(nms) which(grepl("^SMOKE_TX_GROUP", nms))

fit_one_fg <- function(dat, formula, tv_selector) {
  x <- model.matrix(formula, data = dat)[, -1, drop = FALSE]
  tv <- tv_selector(colnames(x))
  if (length(tv) < 1) stop("At least one time-varying column is required.")
  fit <- cmprsk::crr(
    ftime = dat$FGDAY_YEAR, fstatus = dat$FGSTATUS,
    cov1 = x, cov2 = x[, tv, drop = FALSE],
    tf = function(u) matrix(log1p(u), nrow = length(u), ncol = length(tv)),
    failcode = 1, cencode = 0, maxiter = 30
  )
  if (!isTRUE(fit$converged)) stop("Fine-Gray model did not converge.")
  list(coef = fit$coef, var = fit$var, converged = fit$converged)
}

fit_one_cscox <- function(dat) {
  x <- model.matrix(base_formula, data = dat)[, -1, drop = FALSE]
  original_names <- colnames(x)
  safe_names <- make.names(original_names, unique = TRUE)
  colnames(x) <- safe_names
  tv_original <- original_names[tv_joint_and_era(original_names)]
  tv_safe <- safe_names[tv_joint_and_era(original_names)]
  dd <- data.frame(time = dat$FGDAY_YEAR, event = as.integer(dat$FGSTATUS == 1), x, check.names = FALSE)
  rhs <- c(safe_names, paste0("tt(", tv_safe, ")"))
  f <- as.formula(paste("survival::Surv(time, event) ~", paste(rhs, collapse = " + ")))
  fit <- survival::coxph(
    f, data = dd, ties = "efron",
    tt = function(x, t, ...) x * log1p(t),
    singular.ok = FALSE, model = FALSE, x = FALSE, y = FALSE
  )
  if (any(!is.finite(fit$coef)) || isFALSE(fit$converged)) stop("Cause-specific Cox model did not converge.")
  expected_safe <- c(safe_names, paste0("tt(", tv_safe, ")"))
  index <- match(expected_safe, names(fit$coef))
  if (anyNA(index)) {
    stop("Could not align cause-specific Cox coefficients: ", paste(setdiff(expected_safe, names(fit$coef)), collapse = ", "))
  }
  b <- fit$coef[index]
  v <- vcov(fit)[index, index, drop = FALSE]
  final_names <- c(original_names, paste0(tv_original, "*tf"))
  names(b) <- final_names
  dimnames(v) <- list(final_names, final_names)
  list(coef = b, var = v, converged = TRUE)
}

pool_models <- function(fits) {
  term_names <- names(fits[[1]]$coef)
  if (!all(vapply(fits, function(x) identical(names(x$coef), term_names), logical(1)))) {
    stop("Coefficient names differ across models being pooled.")
  }
  coef_matrix <- do.call(cbind, lapply(fits, `[[`, "coef"))
  within <- Reduce("+", lapply(fits, `[[`, "var")) / length(fits)
  if (length(fits) > 1) {
    between <- stats::cov(t(coef_matrix))
  } else {
    between <- matrix(0, nrow = length(term_names), ncol = length(term_names),
                      dimnames = list(term_names, term_names))
  }
  total <- within + (1 + 1 / length(fits)) * between
  list(
    coef = setNames(rowMeans(coef_matrix), term_names),
    var = total, within = within, between = between, m = length(fits)
  )
}

detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (is.na(detected_cores) || detected_cores < 1) detected_cores <- 2L
requested_cores <- suppressWarnings(as.integer(Sys.getenv("UNOS6_CORES", unset = "6")))
if (is.na(requested_cores) || requested_cores < 1) requested_cores <- 1L
fit_cores <- if (.Platform$OS.type == "windows") 1L else min(requested_cores, detected_cores, n_imputations)
cat("Parallel model workers:", fit_cores, "\n")

fit_cache_dir <- Sys.getenv("UNOS6_FIT_CACHE_DIR", unset = "")
if (nzchar(fit_cache_dir)) dir.create(fit_cache_dir, recursive = TRUE, showWarnings = FALSE)

run_cached_mi <- function(cache_name, scenario_label, runner) {
  cache_file <- if (nzchar(fit_cache_dir)) file.path(fit_cache_dir, paste0(cache_name, ".rds")) else ""
  if (nzchar(cache_file) && file.exists(cache_file)) {
    cache <- readRDS(cache_file)
    if (identical(cache$input_md5, input_md5) && identical(cache$scenario, scenario_label) &&
        length(cache$fits) == n_imputations) {
      cat("Loaded model cache:", scenario_label, "\n")
      return(cache$fits)
    }
  }
  cat("Fitting scenario:", scenario_label, "\n")
  fits <- parallel::mclapply(seq_len(n_imputations), function(i) {
    cat("  ", scenario_label, "- imputation", i, "of", n_imputations, "\n")
    runner(completed_data[[i]], i)
  }, mc.cores = fit_cores, mc.preschedule = TRUE)
  failed <- vapply(fits, inherits, logical(1), what = "try-error")
  if (any(failed)) stop("At least one model failed for scenario: ", scenario_label)
  if (nzchar(cache_file)) saveRDS(list(input_md5 = input_md5, scenario = scenario_label, fits = fits), cache_file)
  fits
}

# Primary 20-imputation Fine-Gray fits can be reused from Section 5.
primary_fit_cache_file <- Sys.getenv("UNOS6_PRIMARY_FIT_CACHE_FILE", unset = "")
if (nzchar(primary_fit_cache_file) && file.exists(primary_fit_cache_file)) {
  primary_cache <- readRDS(primary_fit_cache_file)
  if (!identical(primary_cache$input_md5, input_md5) || length(primary_cache$joint_fits) != n_imputations) {
    stop("The primary Fine-Gray fit cache does not match the input data.")
  }
  primary_fits <- lapply(primary_cache$joint_fits, function(obj) {
    fit <- if (!is.null(obj$fit)) obj$fit else obj
    list(coef = fit$coef, var = fit$var, converged = fit$converged)
  })
  cat("Loaded primary Fine-Gray fits from Section 5 cache\n")
} else {
  primary_fits <- run_cached_mi(
    "primary_fg", "Primary Fine-Gray: upstream 90-day exclusion, time-varying era",
    function(dat, i) fit_one_fg(dat, base_formula, tv_joint_and_era)
  )
}

keep_180 <- !(data_raw$FGSTATUS == 1 & data_raw$FGDAY <= 180)
keep_365 <- !(data_raw$FGSTATUS == 1 & data_raw$FGDAY <= 365)

cscox_fits <- run_cached_mi(
  "cause_specific_cox", "Cause-specific Cox: death censored",
  function(dat, i) fit_one_cscox(dat)
)
lag180_fits <- run_cached_mi(
  "lag180_fg", "Fine-Gray: exclude lung cancers within 180 days",
  function(dat, i) fit_one_fg(droplevels(dat[keep_180, , drop = FALSE]), base_formula, tv_joint_and_era)
)
lag365_fits <- run_cached_mi(
  "lag365_fg", "Fine-Gray: exclude lung cancers within 365 days",
  function(dat, i) fit_one_fg(droplevels(dat[keep_365, , drop = FALSE]), base_formula, tv_joint_and_era)
)
fixed_era_fits <- run_cached_mi(
  "fixed_era_fg", "Fine-Gray: fixed transplant-era adjustment",
  function(dat, i) fit_one_fg(dat, base_formula, tv_joint_only)
)
no_era_fits <- run_cached_mi(
  "no_era_fg", "Fine-Gray: no transplant-era adjustment",
  function(dat, i) fit_one_fg(dat, no_era_formula, tv_joint_only)
)

complete_case_mask <- complete.cases(imputation_data)
complete_case_data <- droplevels(data_analysis[complete_case_mask, , drop = FALSE])
complete_case_data$SMOKE_TX_GROUP <- make_joint_group(complete_case_data$CIG_USE, complete_case_data$TX_TYPE)
complete_case_fit <- fit_one_fg(complete_case_data, base_formula, tv_joint_and_era)

pooled <- list(
  primary = pool_models(primary_fits),
  lag180 = pool_models(lag180_fits),
  lag365 = pool_models(lag365_fits),
  fixed_era = pool_models(fixed_era_fits),
  no_era = pool_models(no_era_fits),
  complete_case = pool_models(list(complete_case_fit)),
  cause_specific = pool_models(cscox_fits)
)

scenario_info <- data.frame(
  Scenario_key = c("primary", "lag180", "lag365", "fixed_era", "no_era", "complete_case", "cause_specific"),
  Scenario = c(
    "Primary: >90-day cohort; time-varying era",
    "Exclude lung cancers within 180 days",
    "Exclude lung cancers within 365 days",
    "Fixed transplant-era adjustment",
    "No transplant-era adjustment",
    "Complete-case Fine-Gray",
    "Cause-specific Cox; death censored"
  ),
  Estimand = c(rep("Subdistribution hazard ratio", 6), "Cause-specific hazard ratio"),
  Missing_data = c(rep("20 multiple imputations", 5), "Complete cases", "20 multiple imputations"),
  N = c(nrow(data_raw), sum(keep_180), sum(keep_365), nrow(data_raw), nrow(data_raw),
        nrow(complete_case_data), nrow(data_raw)),
  Lung_cancer_events = c(sum(data_raw$FGSTATUS == 1), sum(data_raw$FGSTATUS[keep_180] == 1),
                         sum(data_raw$FGSTATUS[keep_365] == 1), rep(sum(data_raw$FGSTATUS == 1), 2),
                         sum(complete_case_data$FGSTATUS == 1), sum(data_raw$FGSTATUS == 1)),
  Competing_deaths = c(sum(data_raw$FGSTATUS == 2), sum(data_raw$FGSTATUS[keep_180] == 2),
                       sum(data_raw$FGSTATUS[keep_365] == 2), rep(sum(data_raw$FGSTATUS == 2), 2),
                       sum(complete_case_data$FGSTATUS == 2), sum(data_raw$FGSTATUS == 2)),
  Censored = c(sum(data_raw$FGSTATUS == 0), sum(data_raw$FGSTATUS[keep_180] == 0),
               sum(data_raw$FGSTATUS[keep_365] == 0), rep(sum(data_raw$FGSTATUS == 0), 2),
               sum(complete_case_data$FGSTATUS == 0), sum(data_raw$FGSTATUS == 0)),
  Era_handling = c("Time-varying", "Time-varying", "Time-varying", "Fixed", "None", "Time-varying", "Time-varying"),
  stringsAsFactors = FALSE
)

joint_main_terms <- c(
  "SMOKE_TX_GROUPSmoking history + double lung",
  "SMOKE_TX_GROUPNo smoking history + single lung",
  "SMOKE_TX_GROUPSmoking history + single lung"
)
names(joint_main_terms) <- group_levels[-1]
contrast_weights <- c(-1, -1, 1)

find_time_index <- function(fit, main_term) {
  idx <- which(startsWith(names(fit$coef), paste0(main_term, "*tf")))
  if (length(idx) != 1) stop("Could not identify time coefficient for ", main_term)
  idx
}

rubin_scalar <- function(fit, L) {
  estimate <- drop(L %*% fit$coef)
  within <- drop(L %*% fit$within %*% t(L))
  between <- drop(L %*% fit$between %*% t(L))
  total <- drop(L %*% fit$var %*% t(L))
  se <- sqrt(total)
  df <- if (fit$m <= 1 || between <= .Machine$double.eps) Inf else
    (fit$m - 1) * (1 + within / ((1 + 1 / fit$m) * between))^2
  critical <- if (is.finite(df)) qt(0.975, df) else qnorm(0.975)
  p <- if (is.finite(df)) 2 * pt(abs(estimate / se), df, lower.tail = FALSE) else
    2 * pnorm(abs(estimate / se), lower.tail = FALSE)
  list(estimate = estimate, se = se, df = df, critical = critical, p = p)
}

effect_rows <- bind_rows(lapply(seq_len(nrow(scenario_info)), function(i) {
  key <- scenario_info$Scenario_key[[i]]
  fit <- pooled[[key]]
  bind_rows(lapply(c(1, 3, 5, 10), function(time_years) {
    bind_rows(lapply(names(joint_main_terms), function(group) {
      main_term <- joint_main_terms[[group]]
      main_idx <- match(main_term, names(fit$coef))
      time_idx <- find_time_index(fit, main_term)
      if (is.na(main_idx)) stop("Could not identify main coefficient for ", main_term)
      L <- matrix(0, nrow = 1, ncol = length(fit$coef))
      L[1, main_idx] <- 1
      L[1, time_idx] <- log1p(time_years)
      stat <- rubin_scalar(fit, L)
      data.frame(
        Scenario_key = key,
        Scenario = scenario_info$Scenario[[i]],
        Estimand = scenario_info$Estimand[[i]],
        Missing_data = scenario_info$Missing_data[[i]],
        Time_years = time_years,
        Comparison_group = group,
        Reference_group = group_levels[[1]],
        Effect_measure = ifelse(key == "cause_specific", "HR", "sHR"),
        Estimate = exp(stat$estimate),
        Lower_95CI = exp(stat$estimate - stat$critical * stat$se),
        Upper_95CI = exp(stat$estimate + stat$critical * stat$se),
        P_value = stat$p,
        Rubin_df = stat$df,
        stringsAsFactors = FALSE
      )
    }))
  }))
}))

interaction_tests <- bind_rows(lapply(seq_len(nrow(scenario_info)), function(i) {
  key <- scenario_info$Scenario_key[[i]]
  fit <- pooled[[key]]
  main_idx <- match(unname(joint_main_terms), names(fit$coef))
  time_idx <- vapply(unname(joint_main_terms), function(term) find_time_index(fit, term), integer(1))
  if (anyNA(main_idx)) stop("Could not identify all joint-group main coefficients for ", key)
  L <- matrix(0, nrow = 2, ncol = length(fit$coef))
  L[1, main_idx] <- contrast_weights
  L[2, time_idx] <- contrast_weights
  theta <- as.numeric(L %*% fit$coef)
  V <- L %*% fit$var %*% t(L)
  statistic <- drop(t(theta) %*% solve(V, theta))
  p_global <- pchisq(statistic, df = 2, lower.tail = FALSE)
  l5 <- c(1, log1p(5))
  est5 <- sum(l5 * theta)
  se5 <- sqrt(drop(t(l5) %*% V %*% l5))
  data.frame(
    Scenario_key = key,
    Scenario = scenario_info$Scenario[[i]],
    Estimand = scenario_info$Estimand[[i]],
    Degrees_of_freedom = 2,
    Wald_statistic = statistic,
    Global_P_interaction = p_global,
    Interaction_ratio_at_5y = exp(est5),
    Interaction_5y_lower = exp(est5 - qnorm(0.975) * se5),
    Interaction_5y_upper = exp(est5 + qnorm(0.975) * se5),
    Interaction_5y_P = 2 * pnorm(abs(est5 / se5), lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}))

# Descriptive Aalen-Johansen estimates in the first completed data set provide an observed-risk
# check that does not depend on the regression baseline subdistribution hazard.
aj_data <- completed_data[[1]]
aj_fit <- cmprsk::cuminc(
  ftime = aj_data$FGDAY_YEAR, fstatus = aj_data$FGSTATUS,
  group = aj_data$SMOKE_TX_GROUP, cencode = 0
)
observed_cif <- bind_rows(lapply(group_levels, function(group) {
  component_name <- paste0(group, " 1")
  component <- aj_fit[[component_name]]
  if (is.null(component)) stop("Aalen-Johansen component not found: ", component_name)
  bind_rows(lapply(c(1, 3, 5, 10), function(t) {
    idx <- max(which(component$time <= t))
    est <- component$est[[idx]]
    se <- sqrt(component$var[[idx]])
    data.frame(
      Group = group, Time_years = t, Observed_CIF = est,
      Lower_95CI = max(0, est - qnorm(0.975) * se),
      Upper_95CI = min(1, est + qnorm(0.975) * se),
      SE = se, stringsAsFactors = FALSE
    )
  }))
}))

# Full coefficient audit table.
all_coefficients <- bind_rows(lapply(seq_len(nrow(scenario_info)), function(i) {
  key <- scenario_info$Scenario_key[[i]]
  fit <- pooled[[key]]
  se <- sqrt(diag(fit$var))
  data.frame(
    Scenario_key = key,
    Scenario = scenario_info$Scenario[[i]],
    Term = names(fit$coef),
    Beta = as.numeric(fit$coef),
    SE = se,
    Ratio = exp(fit$coef),
    Lower_95CI = exp(fit$coef - qnorm(0.975) * se),
    Upper_95CI = exp(fit$coef + qnorm(0.975) * se),
    stringsAsFactors = FALSE
  )
}))

# Primary adjusted 5-year absolute risks from Section 5, if available, for descriptive comparison.
section5_model_file <- Sys.getenv("UNOS6_SECTION5_MODEL_FILE", unset = "")
adjusted_5y <- NULL
if (nzchar(section5_model_file) && file.exists(section5_model_file)) {
  section5 <- readRDS(section5_model_file)
  adjusted_5y <- section5$horizon_risks %>%
    filter(Time_years == 5) %>%
    transmute(
      Group = as.character(Group), Time_years,
      Adjusted_CIF = Adjusted_CIF,
      Lower_95CI = Lower_95CI,
      Upper_95CI = Upper_95CI
    )
}

primary_5y <- effect_rows %>% filter(Scenario_key == "primary", Time_years == 5)
lag365_5y <- effect_rows %>% filter(Scenario_key == "lag365", Time_years == 5)
cscox_5y <- effect_rows %>% filter(Scenario_key == "cause_specific", Time_years == 5)
direction_stable <- all(primary_5y$Estimate > 1) && all(lag365_5y$Estimate > 1) && all(cscox_5y$Estimate > 1)
rank_stable <- identical(
  primary_5y$Comparison_group[order(primary_5y$Estimate)],
  lag365_5y$Comparison_group[order(lag365_5y$Estimate)]
)

summary_list <- list(
  source_file = normalizePath(input_file),
  source_md5 = input_md5,
  generated_on = as.character(Sys.Date()),
  primary_cohort = list(
    N = nrow(data_raw), lung_cancer = sum(data_raw$FGSTATUS == 1),
    competing_death = sum(data_raw$FGSTATUS == 2), censored = sum(data_raw$FGSTATUS == 0)
  ),
  exclusions = list(
    additional_180_day_cases = sum(data_raw$FGSTATUS == 1 & data_raw$FGDAY <= 180),
    additional_365_day_cases = sum(data_raw$FGSTATUS == 1 & data_raw$FGDAY <= 365)
  ),
  complete_case = list(
    N = nrow(complete_case_data), lung_cancer = sum(complete_case_data$FGSTATUS == 1),
    retained_percent = 100 * nrow(complete_case_data) / nrow(data_raw)
  ),
  key_robustness = list(
    all_three_joint_exposure_effects_above_one_at_5y_in_primary_lag365_and_cscox = direction_stable,
    ordering_of_joint_exposure_effects_preserved_primary_vs_lag365 = rank_stable,
    primary_global_interaction_p = interaction_tests$Global_P_interaction[interaction_tests$Scenario_key == "primary"],
    lag365_global_interaction_p = interaction_tests$Global_P_interaction[interaction_tests$Scenario_key == "lag365"],
    cause_specific_global_interaction_p = interaction_tests$Global_P_interaction[interaction_tests$Scenario_key == "cause_specific"]
  )
)

write_intermediates <- identical(Sys.getenv("UNOS6_WRITE_INTERMEDIATES", unset = "0"), "1")
if (write_intermediates) {
  write.csv(scenario_info, file.path(output_dir, "cohort_summary.csv"), row.names = FALSE, na = "")
  write.csv(effect_rows, file.path(output_dir, "effect_estimates.csv"), row.names = FALSE, na = "")
  write.csv(interaction_tests, file.path(output_dir, "interaction_tests.csv"), row.names = FALSE, na = "")
  write.csv(observed_cif, file.path(output_dir, "observed_cif.csv"), row.names = FALSE, na = "")
  write.csv(missingness, file.path(output_dir, "missingness.csv"), row.names = FALSE, na = "")
  write.csv(all_coefficients, file.path(output_dir, "all_coefficients.csv"), row.names = FALSE, na = "")
  if (!is.null(adjusted_5y)) write.csv(adjusted_5y, file.path(output_dir, "adjusted_5y.csv"), row.names = FALSE, na = "")
  jsonlite::write_json(summary_list, file.path(output_dir, "analysis_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 10)
}

saveRDS(
  list(
    pooled_models = pooled,
    scenario_info = scenario_info,
    effect_estimates = effect_rows,
    interaction_tests = interaction_tests,
    observed_cif = observed_cif,
    adjusted_5y = adjusted_5y,
    missingness = missingness,
    summary = summary_list
  ),
  file.path(output_dir, "sensitivity_models.rds")
)

# Main sensitivity forest plot at 5 years.
plot_data <- effect_rows %>%
  filter(Time_years == 5) %>%
  mutate(
    Scenario = factor(Scenario, levels = rev(scenario_info$Scenario)),
    Comparison_group = factor(Comparison_group, levels = group_levels[-1])
  )

fig6 <- ggplot(plot_data, aes(x = Estimate, y = Scenario, color = Comparison_group)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#7A7A7A", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = Lower_95CI, xmax = Upper_95CI), height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.5) +
  scale_x_log10(breaks = c(0.5, 1, 2, 4, 8, 16, 32)) +
  scale_color_manual(values = group_colors[group_levels[-1]], labels = group_short[group_levels[-1]]) +
  facet_wrap(~ Comparison_group, ncol = 1, scales = "free_y",
             labeller = as_labeller(group_short)) +
  labs(
    title = "The 5-year exposure gradient was robust across sensitivity analyses",
    subtitle = "Time-specific adjusted hazard ratios versus no smoking history + double lung",
    x = "Hazard ratio (log scale)", y = NULL, color = NULL,
    caption = "Fine-Gray models report subdistribution HRs; the cause-specific model reports cause-specific HRs. Bars show 95% CIs."
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", color = "#17365D", size = 15),
    plot.subtitle = element_text(color = "#4D5B66", size = 10.5),
    strip.text = element_text(face = "bold", color = "#17365D", hjust = 0),
    strip.background = element_rect(fill = "#EAF1F8", color = NA),
    panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
    legend.position = "none", plot.caption = element_text(hjust = 0, color = "#666666", size = 8.5),
    plot.margin = margin(12, 18, 10, 12)
  )

ggsave(file.path(output_dir, "Figure6_sensitivity_forest.png"), fig6, width = 10.5, height = 10.5, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "Figure6_sensitivity_forest.pdf"), fig6, width = 10.5, height = 10.5, device = "pdf")

# Supplemental cohort and observed-risk checks.
lag_plot <- data.frame(
  Window = factor(c(">90 days (primary)", ">180 days", ">365 days"), levels = c(">90 days (primary)", ">180 days", ">365 days")),
  Events = c(sum(data_raw$FGSTATUS == 1), sum(data_raw$FGSTATUS[keep_180] == 1), sum(data_raw$FGSTATUS[keep_365] == 1))
)
p_lag <- ggplot(lag_plot, aes(x = Window, y = Events, fill = Window)) +
  geom_col(width = 0.62, show.legend = FALSE) +
  geom_text(aes(label = Events), vjust = -0.45, fontface = "bold", color = "#17365D") +
  scale_fill_manual(values = c("#2C7BB6", "#F28E2B", "#C73E1D")) +
  scale_y_continuous(limits = c(0, max(lag_plot$Events) * 1.14), expand = expansion(mult = c(0, 0))) +
  labs(title = "Lung-cancer events retained", subtitle = "Stricter exclusion windows remove 10 and 53 additional cases", x = NULL, y = "Events") +
  theme_minimal(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold", color = "#17365D"), panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 18, hjust = 1))

if (!is.null(adjusted_5y)) {
  risk_compare <- bind_rows(
    adjusted_5y %>% transmute(Group, Method = "Adjusted Fine-Gray", Estimate = Adjusted_CIF, Lower_95CI, Upper_95CI),
    observed_cif %>% filter(Time_years == 5) %>% transmute(Group, Method = "Observed Aalen-Johansen", Estimate = Observed_CIF, Lower_95CI, Upper_95CI)
  ) %>% mutate(Group = factor(Group, levels = group_levels))
  p_risk <- ggplot(risk_compare, aes(x = Group, y = 100 * Estimate, color = Method, group = Method)) +
    geom_errorbar(aes(ymin = 100 * Lower_95CI, ymax = 100 * Upper_95CI), width = 0.12, position = position_dodge(width = 0.42)) +
    geom_point(size = 2.5, position = position_dodge(width = 0.42)) +
    scale_color_manual(values = c("Adjusted Fine-Gray" = "#17365D", "Observed Aalen-Johansen" = "#C73E1D")) +
    scale_x_discrete(labels = group_short) +
    labs(title = "Adjusted and observed 5-year risks show the same gradient", x = NULL, y = "5-year lung-cancer risk (%)", color = NULL) +
    theme_minimal(base_size = 10.5) +
    theme(plot.title = element_text(face = "bold", color = "#17365D"), axis.text.x = element_text(angle = 22, hjust = 1), legend.position = "top", panel.grid.major.x = element_blank())
} else {
  p_risk <- ggplot(observed_cif %>% filter(Time_years == 5), aes(x = Group, y = 100 * Observed_CIF, color = Group)) +
    geom_point(size = 3) + geom_errorbar(aes(ymin = 100 * Lower_95CI, ymax = 100 * Upper_95CI), width = 0.12) +
    scale_color_manual(values = group_colors) + scale_x_discrete(labels = group_short) +
    labs(title = "Observed 5-year cumulative incidence", x = NULL, y = "Risk (%)") +
    theme_minimal(base_size = 10.5) + theme(legend.position = "none", axis.text.x = element_text(angle = 22, hjust = 1))
}

figs <- p_lag + p_risk + plot_layout(widths = c(0.82, 1.35)) +
  plot_annotation(
    title = "Sensitivity checks on cohort definition and absolute risk",
    theme = theme(plot.title = element_text(face = "bold", color = "#17365D", size = 15))
  )
ggsave(file.path(output_dir, "FigureS6_cohort_absolute_risk.png"), figs, width = 12, height = 5.7, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "FigureS6_cohort_absolute_risk.pdf"), figs, width = 12, height = 5.7, device = "pdf")

if (write_intermediates) writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
cat("Sensitivity analysis complete. Outputs written to:", normalizePath(output_dir), "\n")
