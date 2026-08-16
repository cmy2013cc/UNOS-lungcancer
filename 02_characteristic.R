#!/usr/bin/env Rscript

local_lib <- Sys.getenv("EJCTS_R_LIB", unset = "E:/UNOS/EJCTS_reanalysis")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

# Reproducible Table 1 analysis for the revised UNOS lung-cancer cohort.
# The source RDS is never modified. Table 1 uses the observed values as stored;
# potentially implausible values are only flagged in a separate review list.

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else file.path("..", "1.2datanew", "3.lungnew.rds")
output_dir <- if (length(args) >= 2) args[[2]] else "."
json_file <- if (length(args) >= 3) args[[3]] else file.path(output_dir, "Table1_statistics.json")
review_json_file <- if (length(args) >= 4) args[[4]] else file.path(output_dir, "Data_quality_review.json")

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- readRDS(input_file)
required <- c(
  "MALIG_P", "TX_YEAR", "AGE", "GENDER", "BMI_TCR", "DIAB", "CIG_USE",
  "GROUPING", "CREAT_TRR", "TBILI", "TOT_SERUM_ALBUM", "CPRA", "TX_TYPE",
  "AGE_DON", "GENDER_DON", "HIST_CIG_DON", "DIABETES_DON"
)
missing_vars <- setdiff(required, names(data))
if (length(missing_vars) > 0) stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
if (any(!data$MALIG_P %in% c("Y", "N"))) stop("MALIG_P contains values other than Y/N")

data$lung_group <- factor(
  ifelse(data$MALIG_P == "Y", "Lung Cancer", "No Lung Cancer"),
  levels = c("Lung Cancer", "No Lung Cancer")
)
data$TX_YEAR_GROUP <- cut(
  data$TX_YEAR,
  breaks = c(2004, 2010, 2015, 2020, Inf),
  labels = c("2005-2010", "2011-2015", "2016-2020", "2021-2023"),
  right = TRUE
)

data$BMI_TCR_REVIEW_FLAG <- !is.na(data$BMI_TCR) & data$BMI_TCR > 80

group_index <- list(
  overall = rep(TRUE, nrow(data)),
  cancer = data$lung_group == "Lung Cancer",
  control = data$lung_group == "No Lung Cancer"
)

fmt_n <- function(x) formatC(as.integer(x), format = "d", big.mark = ",")
fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  if (p > 0.999) return(">0.999")
  sprintf("%.3f", p)
}
fmt_mean_sd <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(x), stats::sd(x))
}
fmt_median_iqr <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  q <- stats::quantile(x, c(0.25, 0.50, 0.75), names = FALSE, na.rm = TRUE)
  sprintf(
    paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
    q[[2]], q[[1]], q[[3]]
  )
}
fmt_cat <- function(x, value, idx) {
  denom <- sum(!is.na(x[idx]))
  n <- sum(x[idx] == value, na.rm = TRUE)
  if (denom == 0) return(paste0(fmt_n(n), " (NA)"))
  paste0(fmt_n(n), " (", sprintf("%.1f", 100 * n / denom), "%)")
}

t_p <- function(x) {
  ok <- !is.na(x) & !is.na(data$lung_group)
  tryCatch(stats::t.test(x[ok] ~ data$lung_group[ok], var.equal = FALSE)$p.value,
           error = function(e) NA_real_)
}
wilcox_p <- function(x) {
  ok <- !is.na(x) & !is.na(data$lung_group)
  tryCatch(stats::wilcox.test(x[ok] ~ data$lung_group[ok], exact = FALSE)$p.value,
           error = function(e) NA_real_)
}
chisq_p <- function(x) {
  ok <- !is.na(x) & !is.na(data$lung_group)
  tab <- table(x[ok], data$lung_group[ok])
  tryCatch(suppressWarnings(stats::chisq.test(tab, correct = FALSE)$p.value),
           error = function(e) NA_real_)
}

rows <- list()
add_row <- function(variable, label, overall = "", cancer = "", control = "",
                    p_value = "", row_type = "variable", indent = 0,
                    method = "", missing = FALSE) {
  rows[[length(rows) + 1L]] <<- data.frame(
    variable = variable,
    characteristic = label,
    overall = overall,
    lung_cancer = cancer,
    no_lung_cancer = control,
    p_value = p_value,
    row_type = row_type,
    indent = indent,
    method = method,
    missing = missing,
    stringsAsFactors = FALSE
  )
}

add_continuous <- function(variable, label, summary_type = c("mean_sd", "median_iqr"),
                           digits = 1, method_label) {
  summary_type <- match.arg(summary_type)
  x <- data[[variable]]
  formatter <- if (summary_type == "mean_sd") fmt_mean_sd else fmt_median_iqr
  p <- if (summary_type == "mean_sd") t_p(x) else wilcox_p(x)
  add_row(
    variable, label,
    formatter(x[group_index$overall], digits),
    formatter(x[group_index$cancer], digits),
    formatter(x[group_index$control], digits),
    fmt_p(p), "variable", 0, method_label
  )
}

add_categorical <- function(variable, label, levels, labels = levels, method_label = "Pearson chi-square test") {
  x <- data[[variable]]
  add_row(variable, label, p_value = fmt_p(chisq_p(x)), row_type = "variable",
          method = method_label)
  for (i in seq_along(levels)) {
    add_row(
      variable, labels[[i]],
      fmt_cat(x, levels[[i]], group_index$overall),
      fmt_cat(x, levels[[i]], group_index$cancer),
      fmt_cat(x, levels[[i]], group_index$control),
      "", "level", 1, ""
    )
  }
}

add_categorical("TX_YEAR_GROUP", "Transplant year", levels(data$TX_YEAR_GROUP))
add_continuous("AGE", "Recipient age, years", "mean_sd", 1, "Welch two-sample t-test")
add_categorical("GENDER", "Recipient sex", c("F", "M"), c("Female", "Male"))
add_continuous("BMI_TCR", "Body mass index, kg/m^2", "median_iqr", 1, "Wilcoxon rank-sum test")
add_categorical("DIAB", "Recipient diabetes", c("N", "Y"), c("No", "Yes"))
add_categorical("CIG_USE", "Recipient smoking history", c("N", "Y"), c("No", "Yes"))
add_categorical(
  "GROUPING", "Primary diagnosis",
  c("A", "D", "Other"),
  c("Chronic obstructive pulmonary disease (COPD)",
    "Interstitial lung disease (ILD)",
    "Other diagnosis")
)
add_continuous("CREAT_TRR", "Serum creatinine, mg/dL", "median_iqr", 2, "Wilcoxon rank-sum test")
add_continuous("TBILI", "Total bilirubin, mg/dL", "median_iqr", 2, "Wilcoxon rank-sum test")
add_continuous("TOT_SERUM_ALBUM", "Serum albumin, g/dL", "median_iqr", 2, "Wilcoxon rank-sum test")
add_continuous("CPRA", "Calculated panel-reactive antibody, %", "median_iqr", 1, "Wilcoxon rank-sum test")
add_categorical("TX_TYPE", "Transplant type", c("S", "D"), c("Single lung", "Double lung"))
add_continuous("AGE_DON", "Donor age, years", "mean_sd", 1, "Welch two-sample t-test")
add_categorical("GENDER_DON", "Donor sex", c("F", "M"), c("Female", "Male"))
add_categorical("HIST_CIG_DON", "Donor smoking history", c("N", "Y"), c("No", "Yes"))
add_categorical("DIABETES_DON", "Donor diabetes", c("N", "Y"), c("No", "Yes"))

table1 <- do.call(rbind, rows)

review_rules <- data.frame(
  variable = c("BMI_TCR", "CREAT_TRR", "TBILI", "TOT_SERUM_ALBUM"),
  review_rule = c(">80 kg/m^2", ">10 mg/dL", ">30 mg/dL", ">6.0 g/dL"),
  rationale = c(
    "Extremely high BMI; possible decimal or unit issue",
    "Extreme creatinine; verify source value and unit",
    "Extreme bilirubin; may be clinically possible but merits verification",
    "Above usual physiologic range; possible source or unit issue"
  ),
  stringsAsFactors = FALSE
)
review_predicates <- list(
  BMI_TCR = function(x) x > 80,
  CREAT_TRR = function(x) x > 10,
  TBILI = function(x) x > 30,
  TOT_SERUM_ALBUM = function(x) x > 6
)
review_records <- list()
for (i in seq_len(nrow(review_rules))) {
  variable <- review_rules$variable[[i]]
  flagged <- which(!is.na(data[[variable]]) & review_predicates[[variable]](data[[variable]]))
  review_rules$flagged_n[[i]] <- length(flagged)
  if (length(flagged) > 0) {
    review_records[[length(review_records) + 1L]] <- data.frame(
      source_row = flagged,
      variable = variable,
      observed_value = data[[variable]][flagged],
      lung_cancer_status = as.character(data$lung_group[flagged]),
      transplant_year = data$TX_YEAR[flagged],
      recipient_age = data$AGE[flagged],
      recipient_sex = as.character(data$GENDER[flagged]),
      review_rule = review_rules$review_rule[[i]],
      stringsAsFactors = FALSE
    )
  }
}
review_records <- do.call(rbind, review_records)

metadata <- list(
  title = "Table 1. Baseline characteristics of lung transplant recipients by post-transplant lung cancer status",
  source_file = normalizePath(input_file, winslash = "/", mustWork = TRUE),
  source_md5 = unname(tools::md5sum(input_file)),
  generated_on = format(Sys.Date(), "%Y-%m-%d"),
  total_n = nrow(data),
  lung_cancer_n = sum(data$lung_group == "Lung Cancer"),
  no_lung_cancer_n = sum(data$lung_group == "No Lung Cancer"),
  continuous_display = "Mean (SD) for recipient and donor age; median [IQR] for BMI, laboratory values, and CPRA",
  tests = "Welch two-sample t-test for recipient and donor age; Wilcoxon rank-sum test for BMI, laboratory values, and CPRA; Pearson chi-square test without continuity correction for categorical variables",
  missing_rule = "Missing rows are not displayed; summaries and percentages use available observations",
  bmi_rule = "Observed BMI values were used without replacement or exclusion; values >80 kg/m^2 are listed separately for review",
  outcome_rule = "Lung cancer status was taken directly from MALIG_P in the revised new dataset; no additional 90-day exclusion was applied at this step"
)

data_output <- data
data_output$BMI_TCR_REVIEW_FLAG <- NULL
saveRDS(data_output, file.path(output_dir, "2.lung.rds"))
jsonlite::write_json(list(metadata = metadata, rows = table1), json_file,
                     pretty = TRUE, auto_unbox = TRUE, na = "null")
jsonlite::write_json(list(source_file = metadata$source_file,
                          rules = review_rules,
                          records = review_records),
                     review_json_file, pretty = TRUE, auto_unbox = TRUE, na = "null")

cat("Table 1 analysis complete\n")
cat("Source:", metadata$source_file, "\n")
cat("N:", metadata$total_n, "| Lung cancer:", metadata$lung_cancer_n,
    "| No lung cancer:", metadata$no_lung_cancer_n, "\n")
cat("BMI values >80 flagged separately for review:", sum(data$BMI_TCR_REVIEW_FLAG), "\n")
cat("Statistics JSON:", normalizePath(json_file, winslash = "/", mustWork = TRUE), "\n")
