#!/usr/bin/env Rscript

local_lib <- Sys.getenv("EJCTS_R_LIB", unset = "E:/UNOS/EJCTS_reanalysis")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

# Effect modification and adjusted absolute lung-cancer risk after lung transplantation.
# This analysis uses the final LAS-era cohort created in new/2.characteristic and preserves
# all recipients and lung-cancer events in the current upstream cohort through 20-fold multiple imputation.

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else file.path("..", "2.characteristic", "2.lung.rds")
output_dir <- if (length(args) >= 2) args[[2]] else "."

required_packages <- c("cmprsk", "dplyr", "ggplot2", "patchwork", "mice", "jsonlite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(cmprsk)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(mice)
  library(jsonlite)
})

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c(
  "No smoking history + double lung",
  "Smoking history + double lung",
  "No smoking history + single lung",
  "Smoking history + single lung"
)
group_colors <- c(
  "No smoking history + double lung" = "#2C7BB6",
  "Smoking history + double lung" = "#00A6A6",
  "No smoking history + single lung" = "#F28E2B",
  "Smoking history + single lung" = "#C73E1D"
)
group_short <- c(
  "No smoking history + double lung" = "No smoking + double lung",
  "Smoking history + double lung" = "Smoking + double lung",
  "No smoking history + single lung" = "No smoking + single lung",
  "Smoking history + single lung" = "Smoking + single lung"
)

factor_with_labels <- function(x, levels, labels) {
  factor(as.character(x), levels = levels, labels = labels)
}

format_p <- function(p) {
  ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
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
    TX_YEAR_GROUP = cut(
      TX_YEAR, breaks = c(2004, 2010, 2015, 2020, Inf), right = TRUE,
      labels = c("2005-2010", "2011-2015", "2016-2020", "2021-2023")
    ),
    SMOKE_TX_GROUP = factor(
      dplyr::case_when(
        CIG_USE == "No" & TX_TYPE == "Double lung" ~ group_levels[[1]],
        CIG_USE == "Yes" & TX_TYPE == "Double lung" ~ group_levels[[2]],
        CIG_USE == "No" & TX_TYPE == "Single lung" ~ group_levels[[3]],
        CIG_USE == "Yes" & TX_TYPE == "Single lung" ~ group_levels[[4]],
        TRUE ~ NA_character_
      ),
      levels = group_levels
    )
  )

model_vars <- c(
  "AGE10", "GENDER", "CIG_USE", "GROUPING", "TX_TYPE", "TX_YEAR_GROUP",
  "ISCHTIME", "CREAT_TRR", "DIAB", "CMV_STATUS", "AGE_DON10", "GENDER_DON",
  "HIST_CIG_DON", "DIABETES_DON", "COD_CAD_DON", "FGDAY_YEAR", "FGSTATUS"
)

# The outcome is allowed to predict missing covariates but is never itself imputed.
imputation_data <- droplevels(data_analysis[, model_vars])
imputation_methods <- mice::make.method(imputation_data)
imputation_methods[vapply(imputation_data, function(x) !anyNA(x), logical(1))] <- ""
imputation_methods[c("FGDAY_YEAR", "FGSTATUS")] <- ""
predictor_matrix <- mice::make.predictorMatrix(imputation_data)
diag(predictor_matrix) <- 0
predictor_matrix[c("FGDAY_YEAR", "FGSTATUS"), ] <- 0

n_imputations <- 20
cache_file <- Sys.getenv("UNOS5_CACHE_FILE", unset = "")
input_md5 <- unname(tools::md5sum(input_file))
if (nzchar(cache_file) && file.exists(cache_file)) {
  cache <- readRDS(cache_file)
  if (!identical(cache$input_md5, input_md5) || length(cache$completed_data) != n_imputations) {
    stop("The imputation cache does not match the current input data.")
  }
  completed_data <- cache$completed_data
} else {
  imputed <- mice::mice(
    imputation_data, m = n_imputations, maxit = 10,
    method = imputation_methods, predictorMatrix = predictor_matrix,
    seed = 20260803, printFlag = FALSE
  )
  completed_data <- lapply(seq_len(n_imputations), function(i) {
    dat <- mice::complete(imputed, i)
    dat$SMOKE_TX_GROUP <- factor(
      dplyr::case_when(
        dat$CIG_USE == "No" & dat$TX_TYPE == "Double lung" ~ group_levels[[1]],
        dat$CIG_USE == "Yes" & dat$TX_TYPE == "Double lung" ~ group_levels[[2]],
        dat$CIG_USE == "No" & dat$TX_TYPE == "Single lung" ~ group_levels[[3]],
        dat$CIG_USE == "Yes" & dat$TX_TYPE == "Single lung" ~ group_levels[[4]],
        TRUE ~ NA_character_
      ),
      levels = group_levels
    )
    dat$TX_TYPE_INT <- relevel(dat$TX_TYPE, ref = "Double lung")
    droplevels(dat)
  })
  if (nzchar(cache_file)) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(input_md5 = input_md5, completed_data = completed_data), cache_file)
  }
}
cat("Multiple imputation complete:", n_imputations, "data sets\n")

base_formula <- ~ AGE10 + GENDER + GROUPING + SMOKE_TX_GROUP +
  ISCHTIME + CREAT_TRR + DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON +
  HIST_CIG_DON + DIABETES_DON + COD_CAD_DON + TX_YEAR_GROUP

fit_formula <- function(dat, formula, tv_selector) {
  x <- model.matrix(formula, data = dat)[, -1, drop = FALSE]
  tv <- tv_selector(colnames(x))
  fit <- crr(
    ftime = dat$FGDAY_YEAR, fstatus = dat$FGSTATUS,
    cov1 = x, cov2 = x[, tv, drop = FALSE],
    tf = function(u) matrix(log1p(u), nrow = length(u), ncol = length(tv)),
    failcode = 1, cencode = 0, maxiter = 30
  )
  if (!isTRUE(fit$converged)) stop("Fine-Gray model did not converge.")
  list(fit = fit, x_columns = colnames(x), tv = tv, tv_columns = colnames(x)[tv])
}

joint_tv_selector <- function(nms) {
  which((grepl("^SMOKE_TX_GROUP", nms) & !grepl(":", nms)) | grepl("^TX_YEAR_GROUP", nms))
}

detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
if (is.na(detected_cores) || detected_cores < 1) detected_cores <- 2L
fit_cores <- if (.Platform$OS.type == "windows") 1L else max(1L, min(4L, detected_cores))
fit_cache_file <- Sys.getenv("UNOS5_FIT_CACHE_FILE", unset = "")
if (nzchar(fit_cache_file) && file.exists(fit_cache_file)) {
  fit_cache <- readRDS(fit_cache_file)
  if (!identical(fit_cache$input_md5, input_md5) || length(fit_cache$joint_fits) != n_imputations) {
    stop("The model-fit cache does not match the current input data.")
  }
  joint_fits <- fit_cache$joint_fits
} else {
  joint_fits <- parallel::mclapply(seq_along(completed_data), function(i) {
    cat("Fitting primary joint model", i, "of", n_imputations, "\n")
    fit_formula(completed_data[[i]], base_formula, joint_tv_selector)
  }, mc.cores = fit_cores)
  if (nzchar(fit_cache_file)) {
    dir.create(dirname(fit_cache_file), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(input_md5 = input_md5, joint_fits = joint_fits), fit_cache_file)
  }
}

pool_fits <- function(fit_objects) {
  fits <- lapply(fit_objects, `[[`, "fit")
  term_names <- names(fits[[1]]$coef)
  if (!all(vapply(fits, function(x) identical(names(x$coef), term_names), logical(1)))) {
    stop("Coefficient names differ across imputed models.")
  }
  coef_matrix <- do.call(cbind, lapply(fits, function(x) x$coef))
  within <- Reduce("+", lapply(fits, function(x) x$var)) / length(fits)
  between <- stats::cov(t(coef_matrix))
  total <- within + (1 + 1 / length(fits)) * between
  list(
    coef = setNames(rowMeans(coef_matrix), term_names), var = total,
    within = within, between = between, m = length(fits)
  )
}

pooled_joint <- pool_fits(joint_fits)

# The four-level joint factor is exactly equivalent to smoking + transplant type + their
# interaction. Linear contrasts therefore recover the formal interaction without refitting:
# beta(Yes+Single) - beta(Yes+Double) - beta(No+Single), for both the baseline and log-time terms.
joint_main_terms <- c(
  "SMOKE_TX_GROUPSmoking history + double lung",
  "SMOKE_TX_GROUPNo smoking history + single lung",
  "SMOKE_TX_GROUPSmoking history + single lung"
)
joint_main_indices <- match(joint_main_terms, names(pooled_joint$coef))
joint_time_indices <- vapply(
  joint_main_terms,
  function(term) {
    idx <- which(startsWith(names(pooled_joint$coef), paste0(term, "*tf")))
    if (length(idx) != 1) stop("Could not identify time coefficient for ", term)
    idx
  },
  integer(1)
)
if (anyNA(joint_main_indices)) stop("Could not identify joint-group main coefficients.")
contrast_weights <- c(-1, -1, 1)
L_interaction <- matrix(0, nrow = 2, ncol = length(pooled_joint$coef))
L_interaction[1, joint_main_indices] <- contrast_weights
L_interaction[2, joint_time_indices] <- contrast_weights
formal_theta <- as.numeric(L_interaction %*% pooled_joint$coef)
formal_theta_var <- L_interaction %*% pooled_joint$var %*% t(L_interaction)
formal_stat <- drop(t(formal_theta) %*% solve(formal_theta_var, formal_theta))
formal_p <- pchisq(formal_stat, df = 2, lower.tail = FALSE)

time_specific_interaction <- bind_rows(lapply(c(1, 3, 5, 10), function(t) {
  l <- log1p(t)
  a <- c(1, l)
  est <- sum(a * formal_theta)
  se <- sqrt(drop(t(a) %*% formal_theta_var %*% a))
  data.frame(
    Time_years = t,
    Interaction_ratio = exp(est),
    Lower_95CI = exp(est - qnorm(0.975) * se),
    Upper_95CI = exp(est + qnorm(0.975) * se),
    P_value = 2 * pnorm(abs(est / se), lower.tail = FALSE)
  )
}))

# Secondary, prespecified heterogeneity screens. These are deliberately exploratory and use the
# first completed data set in proportional Fine-Gray models; the primary smoking-by-transplant
# interaction and all absolute-risk estimates above remain 20-imputation analyses.
secondary_specs <- list(
  "Recipient age (continuous, per 10 years)" = list(
    formula = ~ AGE10 + GENDER + GROUPING + SMOKE_TX_GROUP +
      ISCHTIME + CREAT_TRR + DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON +
      HIST_CIG_DON + DIABETES_DON + COD_CAD_DON + TX_YEAR_GROUP +
      SMOKE_TX_GROUP:AGE10,
    pattern = "AGE10"
  ),
  "Recipient sex" = list(
    formula = ~ AGE10 + GENDER + GROUPING + SMOKE_TX_GROUP +
      ISCHTIME + CREAT_TRR + DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON +
      HIST_CIG_DON + DIABETES_DON + COD_CAD_DON + TX_YEAR_GROUP +
      SMOKE_TX_GROUP:GENDER,
    pattern = "GENDER"
  ),
  "Primary diagnosis" = list(
    formula = ~ AGE10 + GENDER + GROUPING + SMOKE_TX_GROUP +
      ISCHTIME + CREAT_TRR + DIAB + CMV_STATUS + AGE_DON10 + GENDER_DON +
      HIST_CIG_DON + DIABETES_DON + COD_CAD_DON + TX_YEAR_GROUP +
      SMOKE_TX_GROUP:GROUPING,
    pattern = "GROUPING"
  )
)

secondary_results <- bind_rows(lapply(names(secondary_specs), function(label) {
  spec <- secondary_specs[[label]]
  if (identical(label, "Primary diagnosis")) {
    event_cells <- with(
      completed_data[[1]],
      tapply(FGSTATUS == 1, list(SMOKE_TX_GROUP, GROUPING), sum)
    )
    if (any(event_cells < 5, na.rm = TRUE)) {
      return(data.frame(
        Effect_modifier = label,
        Degrees_of_freedom = 6,
        Wald_statistic = NA_real_,
        P_interaction = NA_real_,
        Method = "Not estimated because multiple joint exposure/diagnosis cells contained fewer than 5 lung-cancer events",
        Reason = paste0("Minimum cell event count = ", min(event_cells, na.rm = TRUE))
      ))
    }
  }
  x <- model.matrix(spec$formula, data = completed_data[[1]])[, -1, drop = FALSE]
  fit <- crr(
    completed_data[[1]]$FGDAY_YEAR, completed_data[[1]]$FGSTATUS,
    cov1 = x, failcode = 1, cencode = 0, maxiter = 30
  )
  if (!isTRUE(fit$converged)) stop("Secondary interaction model did not converge for ", label)
  nms <- names(fit$coef)
  idx <- which(grepl(":", nms) & grepl("SMOKE_TX_GROUP", nms) & grepl(spec$pattern, nms))
  if (length(idx) == 0) stop("No secondary interaction terms found for ", label)
  stat <- drop(t(fit$coef[idx]) %*% solve(fit$var[idx, idx, drop = FALSE], fit$coef[idx]))
  data.frame(
    Effect_modifier = label,
    Degrees_of_freedom = length(idx),
    Wald_statistic = stat,
    P_interaction = pchisq(stat, df = length(idx), lower.tail = FALSE),
    Method = "Exploratory proportional Fine-Gray screen in first completed data set",
    Reason = ""
  )
}))
secondary_results$FDR_q_value <- NA_real_
estimable_secondary <- which(!is.na(secondary_results$P_interaction))
secondary_results$FDR_q_value[estimable_secondary] <-
  p.adjust(secondary_results$P_interaction[estimable_secondary], method = "BH")

# Standardized cumulative incidence under each smoking/transplant intervention. Risk gradients
# permit Rubin pooling of model-based uncertainty; the baseline subdistribution-hazard jumps are
# retained from each imputed Fine-Gray fit.
standardized_one <- function(dat, fit_object, group, times) {
  dat_new <- dat
  dat_new$SMOKE_TX_GROUP <- factor(group, levels = group_levels)
  x <- model.matrix(base_formula, data = dat_new)[, -1, drop = FALSE]
  fit <- fit_object$fit
  p <- ncol(x)
  tv <- fit_object$tv
  q <- length(tv)
  beta <- fit$coef[seq_len(p)]
  gamma <- fit$coef[p + seq_len(q)]
  eta <- as.numeric(x %*% beta)
  z <- x[, tv, drop = FALSE]
  z_key <- apply(z, 1, paste, collapse = "|")
  unique_keys <- unique(z_key)
  event_times <- fit$uftime
  jumps <- fit$bfitj

  rows <- lapply(times, function(horizon) {
    if (horizon <= 0) {
      return(list(time = horizon, risk = 0, gradient = rep(0, p + q), within_var = 0))
    }
    H <- numeric(nrow(x))
    dH_gamma <- matrix(0, nrow = nrow(x), ncol = q)
    keep <- event_times <= horizon
    u <- event_times[keep]
    bj <- jumps[keep]
    logu <- log1p(u)
    for (key in unique_keys) {
      idx <- which(z_key == key)
      zk <- z[idx[[1]], ]
      slope <- sum(zk * gamma)
      weighted <- bj * exp(slope * logu)
      baseline_component <- sum(weighted)
      derivative_scalar <- sum(weighted * logu)
      exp_eta <- exp(eta[idx])
      H[idx] <- exp_eta * baseline_component
      dH_gamma[idx, ] <- exp_eta * matrix(
        rep(derivative_scalar * zk, each = length(idx)),
        nrow = length(idx), byrow = FALSE
      )
    }
    survival_component <- exp(-H)
    risk <- mean(1 - survival_component)
    gradient_main <- colMeans((survival_component * H) * x)
    gradient_gamma <- colMeans(survival_component * dH_gamma)
    gradient <- c(gradient_main, gradient_gamma)
    within_var <- drop(t(gradient) %*% fit$var %*% gradient)
    list(time = horizon, risk = risk, gradient = gradient, within_var = within_var)
  })
  names(rows) <- as.character(times)
  rows
}

pool_scalar <- function(estimates, within_variances, transform = c("identity", "logit", "log")) {
  transform <- match.arg(transform)
  m <- length(estimates)
  qbar <- mean(estimates)
  ubar <- mean(within_variances)
  between <- if (m > 1) var(estimates) else 0
  total <- ubar + (1 + 1 / m) * between
  df <- if (between <= .Machine$double.eps) Inf else
    (m - 1) * (1 + ubar / ((1 + 1 / m) * between))^2
  critical <- if (is.finite(df)) qt(0.975, df) else qnorm(0.975)
  if (transform == "logit") {
    q_safe <- min(max(qbar, 1e-8), 1 - 1e-8)
    se_t <- sqrt(total) / (q_safe * (1 - q_safe))
    ci <- plogis(qlogis(q_safe) + c(-1, 1) * critical * se_t)
  } else if (transform == "log") {
    q_safe <- max(qbar, 1e-8)
    se_t <- sqrt(total) / q_safe
    ci <- exp(log(q_safe) + c(-1, 1) * critical * se_t)
  } else {
    ci <- qbar + c(-1, 1) * critical * sqrt(total)
  }
  c(Estimate = qbar, Lower_95CI = ci[[1]], Upper_95CI = ci[[2]], SE = sqrt(total), DF = df)
}

curve_times <- sort(unique(c(0, seq(0.25, 10, by = 0.25), 1, 3, 5, 10)))
prediction_store <- vector("list", n_imputations)
for (i in seq_len(n_imputations)) {
  cat("Standardizing absolute risks for imputation", i, "of", n_imputations, "\n")
  prediction_store[[i]] <- setNames(lapply(group_levels, function(g) {
    standardized_one(completed_data[[i]], joint_fits[[i]], g, curve_times)
  }), group_levels)
}

standardized_curve <- bind_rows(lapply(group_levels, function(g) {
  bind_rows(lapply(curve_times, function(t) {
    estimates <- vapply(prediction_store, function(x) x[[g]][[as.character(t)]]$risk, numeric(1))
    within <- vapply(prediction_store, function(x) x[[g]][[as.character(t)]]$within_var, numeric(1))
    stat <- if (t == 0) c(Estimate = 0, Lower_95CI = 0, Upper_95CI = 0, SE = 0, DF = Inf) else
      pool_scalar(estimates, within, "logit")
    data.frame(
      Group = g, Time_years = t,
      Adjusted_CIF = stat[["Estimate"]], Lower_95CI = stat[["Lower_95CI"]],
      Upper_95CI = stat[["Upper_95CI"]], SE = stat[["SE"]]
    )
  }))
})) %>% mutate(Group = factor(Group, levels = group_levels))

horizon_risks <- standardized_curve %>%
  filter(Time_years %in% c(1, 3, 5, 10)) %>%
  mutate(
    Risk_percent = 100 * Adjusted_CIF,
    Lower_percent = 100 * Lower_95CI,
    Upper_percent = 100 * Upper_95CI
  )

reference_group <- group_levels[[1]]
risk_contrasts <- bind_rows(lapply(group_levels[-1], function(g) {
  bind_rows(lapply(c(1, 3, 5, 10), function(t) {
    est_diff <- numeric(n_imputations)
    var_diff <- numeric(n_imputations)
    est_logrr <- numeric(n_imputations)
    var_logrr <- numeric(n_imputations)
    for (i in seq_len(n_imputations)) {
      a <- prediction_store[[i]][[g]][[as.character(t)]]
      b <- prediction_store[[i]][[reference_group]][[as.character(t)]]
      est_diff[[i]] <- a$risk - b$risk
      grad_diff <- a$gradient - b$gradient
      var_diff[[i]] <- drop(t(grad_diff) %*% joint_fits[[i]]$fit$var %*% grad_diff)
      est_logrr[[i]] <- log(max(a$risk, 1e-8) / max(b$risk, 1e-8))
      grad_logrr <- a$gradient / max(a$risk, 1e-8) - b$gradient / max(b$risk, 1e-8)
      var_logrr[[i]] <- drop(t(grad_logrr) %*% joint_fits[[i]]$fit$var %*% grad_logrr)
    }
    rd <- pool_scalar(est_diff, var_diff, "identity")
    logrr <- pool_scalar(est_logrr, var_logrr, "identity")
    data.frame(
      Comparison_group = g, Reference_group = reference_group, Time_years = t,
      Risk_difference = rd[["Estimate"]], RD_lower_95CI = rd[["Lower_95CI"]],
      RD_upper_95CI = rd[["Upper_95CI"]],
      Risk_ratio = exp(logrr[["Estimate"]]), RR_lower_95CI = exp(logrr[["Lower_95CI"]]),
      RR_upper_95CI = exp(logrr[["Upper_95CI"]])
    )
  }))
}))

# Observed data counts are shown separately from imputed allocation summaries.
observed_counts <- data_analysis %>%
  mutate(
    Group = as.character(SMOKE_TX_GROUP),
    Group = ifelse(is.na(Group), "Smoking history missing", Group),
    Status = factor(FGSTATUS, levels = c(0, 1, 2), labels = c("Censored at last follow-up", "Lung cancer", "Competing death"))
  ) %>%
  count(Group, Status, name = "N") %>%
  group_by(Group) %>% mutate(Group_total = sum(N)) %>% ungroup()

imputed_count_summary <- bind_rows(lapply(seq_along(completed_data), function(i) {
  completed_data[[i]] %>% count(SMOKE_TX_GROUP, FGSTATUS, name = "N") %>% mutate(Imputation = i)
})) %>%
  group_by(SMOKE_TX_GROUP, FGSTATUS) %>%
  summarise(Mean_N = mean(N), Minimum_N = min(N), Maximum_N = max(N), .groups = "drop") %>%
  rename(Group = SMOKE_TX_GROUP, Status = FGSTATUS)

model_summary <- data.frame(
  Metric = c(
    "Source cohort", "Lung cancer events", "Competing deaths", "Censored at last follow-up",
    "Number of imputations", "Primary exposure", "Reference exposure group",
    "Formal interaction test", "Primary absolute-risk horizons",
    "Missing-data note", "Interpretation boundary"
  ),
  Value = c(
    nrow(data_analysis), sum(data_analysis$FGSTATUS == 1), sum(data_analysis$FGSTATUS == 2),
    sum(data_analysis$FGSTATUS == 0), n_imputations,
    "Four-level joint smoking-history/transplant-type variable",
    reference_group,
    paste0("Smoking history x transplant type, pooled 2-df Wald P = ", format_p(formal_p)),
    "1, 3, 5, and 10 years",
    paste0(sum(is.na(data_analysis$CIG_USE)), " recipients had missing smoking history and were retained by multiple imputation"),
    "Adjusted associations and standardized risks for clinical awareness; not a validated screening rule"
  ),
  stringsAsFactors = FALSE
)

metadata <- data.frame(
  Field = c(
    "Source file", "Source MD5", "Generated on", "Cohort", "Outcome", "Competing event",
    "Early-cancer exclusion", "Missing data", "Primary adjustment", "Absolute-risk method",
    "Confidence intervals", "Formal interaction", "Secondary interaction multiplicity"
  ),
  Value = c(
    normalizePath(input_file, winslash = "/", mustWork = TRUE), unname(tools::md5sum(input_file)),
    as.character(Sys.Date()), sprintf("Final LAS-era new cohort: %s recipients", format(nrow(data_analysis), big.mark = ",")),
    "First post-transplant lung cancer (FGSTATUS = 1)", "Death before lung cancer (FGSTATUS = 2)",
    "The upstream 90-day exclusion was used; no duplicate exclusion was applied in this script",
    "20 multiple imputations; outcomes used only as predictors and never imputed",
    "Same expanded baseline Fine-Gray structure as new/4.multi; joint exposure and era allowed to vary with log(1+time)",
    "Regression standardization over the full cohort under each of four exposure interventions",
    "Rubin-pooled model-based 95% CIs using delta-method gradients; baseline-hazard uncertainty is not separately bootstrapped",
    "Exact two-row linear contrast of the four-level joint model: baseline and log-time interaction jointly tested (2 df)",
    "Age and sex are exploratory proportional Fine-Gray screens with Benjamini-Hochberg FDR q values; diagnosis interaction was not estimated because of sparse event cells"
  ),
  stringsAsFactors = FALSE
)

plot_theme <- theme_minimal(base_size = 11.5) +
  theme(
    plot.title = element_text(face = "bold", size = 15, color = "#17365D"),
    plot.subtitle = element_text(size = 10.5, color = "#4D4D4D"),
    plot.caption = element_text(size = 8.7, color = "#595959", hjust = 0),
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    legend.position = "bottom", legend.title = element_blank(),
    strip.text = element_text(face = "bold", color = "#17365D")
  )

panel_a <- ggplot(standardized_curve, aes(Time_years, Adjusted_CIF, color = Group, fill = Group)) +
  geom_ribbon(aes(ymin = Lower_95CI, ymax = Upper_95CI), alpha = 0.11, linewidth = 0) +
  geom_line(linewidth = 1.15) +
  scale_color_manual(values = group_colors, labels = group_short) +
  scale_fill_manual(values = group_colors, labels = group_short) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(labels = function(x) paste0(sprintf("%.1f", 100 * x), "%"), expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "A. Adjusted cumulative incidence",
    subtitle = "Standardized to the full cohort; death treated as a competing event",
    x = "Years after transplantation", y = "Adjusted lung-cancer cumulative incidence"
  ) + plot_theme

horizon_plot_data <- horizon_risks %>%
  mutate(
    Group_label = factor(unname(group_short[as.character(Group)]), levels = rev(unname(group_short[group_levels]))),
    Time_label = factor(paste0(Time_years, " year", ifelse(Time_years == 1, "", "s")),
                        levels = c("1 year", "3 years", "5 years", "10 years"))
  )
panel_b <- ggplot(horizon_plot_data, aes(Risk_percent, Group_label, color = Group)) +
  geom_errorbar(
    aes(xmin = Lower_percent, xmax = Upper_percent),
    orientation = "y", width = 0.12, linewidth = 0.75
  ) +
  geom_point(size = 2.6) +
  facet_wrap(~ Time_label, ncol = 2, scales = "free_x") +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_x_continuous(labels = function(x) paste0(sprintf("%.1f", x), "%"), expand = expansion(mult = c(0.03, 0.12))) +
  labs(title = "B. Time-specific absolute risk", x = "Adjusted cumulative incidence (95% CI)", y = NULL) +
  plot_theme + theme(axis.text.y = element_text(size = 8.5))

panel_c <- ggplot(time_specific_interaction, aes(Time_years, Interaction_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#737373", linewidth = 0.7) +
  geom_ribbon(aes(ymin = Lower_95CI, ymax = Upper_95CI), fill = "#C73E1D", alpha = 0.13) +
  geom_line(color = "#C73E1D", linewidth = 1.1) +
  geom_point(color = "#C73E1D", size = 2.7) +
  scale_x_continuous(breaks = c(1, 3, 5, 10), limits = c(1, 10)) +
  scale_y_log10() +
  labs(
    title = "C. Smoking-by-transplant interaction over time",
    subtitle = paste0("Overall pooled 2-df interaction P = ", format_p(formal_p)),
    x = "Years after transplantation", y = "Interaction ratio (log scale)"
  ) + plot_theme

main_figure <- panel_a / (panel_b | panel_c) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(
    title = "Smoking history, transplant type, and post-transplant lung-cancer risk",
    subtitle = paste0(
      "Fine-Gray analysis; N = ", format(nrow(data_analysis), big.mark = ","),
      ", lung-cancer events = ", sum(data_analysis$FGSTATUS == 1), ", 20 multiple imputations"
    ),
    caption = "Adjusted for recipient, donor, procedural, diagnostic, and transplant-era covariates. Shaded areas and bars are model-based 95% confidence intervals.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 17, color = "#17365D"),
      plot.subtitle = element_text(size = 10.5, color = "#4D4D4D"),
      plot.caption = element_text(size = 8.8, color = "#595959", hjust = 0)
    )
  )

ggsave(file.path(output_dir, "Figure5_effect_modification_absolute_risk.png"), main_figure,
       width = 12, height = 12.2, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "Figure5_effect_modification_absolute_risk.pdf"), main_figure,
       width = 12, height = 12.2, device = "pdf")

secondary_plot_data <- secondary_results %>%
  mutate(
    Effect_modifier = factor(Effect_modifier, levels = rev(Effect_modifier)),
    Estimable = !is.na(P_interaction),
    Plot_P = ifelse(Estimable, P_interaction, 0),
    Significant = ifelse(Estimable, FDR_q_value < 0.05, FALSE),
    Result_label = ifelse(
      Estimable,
      paste0("P=", format_p(P_interaction), "; q=", format_p(FDR_q_value)),
      "Not estimable: sparse event cells"
    )
  )
secondary_figure <- ggplot(secondary_plot_data, aes(Plot_P, Effect_modifier)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "#9B1C1C", linewidth = 0.7) +
  geom_segment(data = secondary_plot_data %>% filter(Estimable),
               aes(x = 0, xend = Plot_P, yend = Effect_modifier), color = "#B7C9DB", linewidth = 1) +
  geom_point(data = secondary_plot_data %>% filter(Estimable), aes(color = Significant), size = 3.4) +
  geom_point(data = secondary_plot_data %>% filter(!Estimable), shape = 4, size = 3.4, color = "#737373") +
  geom_text(aes(label = Result_label),
            hjust = -0.08, size = 3.5, color = "#303030") +
  scale_color_manual(values = c(`FALSE` = "#2C7BB6", `TRUE` = "#C73E1D"), guide = "none") +
  scale_x_continuous(limits = c(0, max(0.12, max(secondary_results$P_interaction, na.rm = TRUE) * 1.35)),
                     expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Prespecified secondary effect-modification tests",
    subtitle = "P values test the complete joint-group interaction block; q values control FDR across three tests",
    x = "P for interaction", y = NULL,
    caption = "These tests are secondary and are not interpreted from subgroup-specific significance alone."
  ) + plot_theme + theme(legend.position = "none")

ggsave(file.path(output_dir, "FigureS_effect_modification_secondary.png"), secondary_figure,
       width = 9.5, height = 4.8, dpi = 320, bg = "white")
ggsave(file.path(output_dir, "FigureS_effect_modification_secondary.pdf"), secondary_figure,
       width = 9.5, height = 4.8, device = "pdf")

analysis_object <- list(
  pooled_joint_fit = pooled_joint,
  formal_interaction = list(
    contrast_matrix = L_interaction, coefficients = formal_theta,
    variance = formal_theta_var, statistic = formal_stat, df = 2, p = formal_p
  ),
  time_specific_interaction = time_specific_interaction,
  standardized_curve = standardized_curve,
  horizon_risks = horizon_risks,
  risk_contrasts = risk_contrasts,
  secondary_interactions = secondary_results,
  sample = list(
    N = nrow(data_analysis), lung_cancer = sum(data_analysis$FGSTATUS == 1),
    competing_death = sum(data_analysis$FGSTATUS == 2), censored = sum(data_analysis$FGSTATUS == 0)
  ),
  source_file = normalizePath(input_file, winslash = "/", mustWork = TRUE),
  source_md5 = unname(tools::md5sum(input_file)), generated_on = as.character(Sys.Date())
)
saveRDS(analysis_object, file.path(output_dir, "effect_modification_absolute_risk_model.rds"))

json_payload <- list(
  Model_summary = model_summary,
  Formal_interaction = data.frame(
    Test = "Smoking history x transplant type (main + log-time interaction)",
    Degrees_of_freedom = 2, Wald_statistic = formal_stat, P_value = formal_p,
    P_value_formatted = format_p(formal_p)
  ),
  Time_specific_interaction = time_specific_interaction %>%
    mutate(
      `Interaction ratio (95% CI)` = sprintf("%.2f (%.2f-%.2f)", Interaction_ratio, Lower_95CI, Upper_95CI),
      `P value` = format_p(P_value)
    ),
  Adjusted_absolute_risks = horizon_risks %>%
    mutate(`Adjusted risk (95% CI)` = sprintf("%.2f%% (%.2f%%-%.2f%%)", Risk_percent, Lower_percent, Upper_percent)),
  Risk_contrasts = risk_contrasts %>%
    mutate(
      `Risk difference, percentage points (95% CI)` = sprintf(
        "%.2f (%.2f-%.2f)", 100 * Risk_difference, 100 * RD_lower_95CI, 100 * RD_upper_95CI
      ),
      `Risk ratio (95% CI)` = sprintf("%.2f (%.2f-%.2f)", Risk_ratio, RR_lower_95CI, RR_upper_95CI)
    ),
  Secondary_interactions = secondary_results %>%
    mutate(`P for interaction` = format_p(P_interaction), `FDR q value` = format_p(FDR_q_value)),
  Observed_group_counts = observed_counts,
  Imputed_group_count_summary = imputed_count_summary,
  Standardized_curve = standardized_curve,
  Metadata = metadata
)
write_json(
  json_payload, file.path(output_dir, "effect_modification_results.json"),
  pretty = TRUE, auto_unbox = TRUE, dataframe = "rows", na = "null", digits = NA
)

cat("Effect-modification and absolute-risk analysis complete\n")
cat("N:", nrow(data_analysis), "| Lung cancer:", sum(data_analysis$FGSTATUS == 1),
    "| Competing deaths:", sum(data_analysis$FGSTATUS == 2), "\n")
cat("Formal smoking x transplant interaction P:", format_p(formal_p), "\n")
cat("Outputs:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
