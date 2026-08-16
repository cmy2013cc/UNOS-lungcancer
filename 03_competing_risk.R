#!/usr/bin/env Rscript

local_lib <- Sys.getenv("EJCTS_R_LIB", unset = "E:/UNOS/EJCTS_reanalysis")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

# Revised competing-risk analysis for the final UNOS lung-cancer cohort.
# The source dataset is read-only. All outputs are written to output_dir.

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[[1]] else file.path("..", "2.characteristic", "2.lung.rds")
output_dir <- if (length(args) >= 2) args[[2]] else "."

required_packages <- c("cmprsk", "dplyr", "ggplot2", "openxlsx", "patchwork", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(cmprsk)
  library(dplyr)
  library(ggplot2)
  library(openxlsx)
  library(patchwork)
  library(scales)
})

if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data_raw <- readRDS(input_file)
required_vars <- c("FGDAY", "FGSTATUS", "MALIG_P", "TX_YEAR", "AGE", "GENDER", "CIG_USE",
                   "GROUPING", "TX_TYPE")
missing_vars <- setdiff(required_vars, names(data_raw))
if (length(missing_vars) > 0) stop("Missing required variables: ", paste(missing_vars, collapse = ", "))

if (any(is.na(data_raw$FGDAY)) || any(is.na(data_raw$FGSTATUS))) {
  stop("FGDAY and FGSTATUS must be complete before this analysis.")
}
if (any(!data_raw$FGSTATUS %in% c(0, 1, 2))) stop("FGSTATUS must be coded 0/1/2.")
if (any(data_raw$FGDAY < 0)) stop("FGDAY contains negative values.")

# The source cohort is restricted upstream to the LAS era (TX_DATE >= 2005-05-04).
# These non-overlapping groups match Table 1.
data_raw$TX_YEAR_GROUP <- cut(
  data_raw$TX_YEAR,
  breaks = c(2004, 2010, 2015, 2020, Inf),
  labels = c("2005-2010", "2011-2015", "2016-2020", "2021-2023"),
  right = TRUE
)

factor_with_labels <- function(x, levels, labels) {
  factor(as.character(x), levels = levels, labels = labels)
}

data_analysis <- data_raw %>%
  mutate(
    FGDAY_YEAR = FGDAY / 365.25,
    TX_YEAR_GROUP = factor(TX_YEAR_GROUP,
                           levels = c("2005-2010", "2011-2015", "2016-2020", "2021-2023")),
    AGE_GROUP = cut(AGE, breaks = c(-Inf, 55, 60, 65, 70, Inf), right = FALSE,
                    labels = c("<55", "55-59", "60-64", "65-69", ">=70")),
    GENDER = factor_with_labels(GENDER, c("F", "M"), c("Female", "Male")),
    DIAB = factor_with_labels(DIAB, c("N", "Y"), c("No", "Yes")),
    CIG_USE = factor_with_labels(CIG_USE, c("N", "Y"), c("No", "Yes")),
    GROUPING = factor_with_labels(
      GROUPING, c("A", "D", "Other"),
      c("COPD", "Interstitial lung disease", "Other diagnosis")
    ),
    TX_TYPE = factor_with_labels(TX_TYPE, c("S", "D"), c("Single lung", "Double lung")),
    CMV_STATUS = factor_with_labels(CMV_STATUS, c("N", "P"), c("Negative", "Positive")),
    EBV_SEROSTATUS = factor_with_labels(EBV_SEROSTATUS, c("N", "P"), c("Negative", "Positive")),
    GENDER_DON = factor_with_labels(GENDER_DON, c("F", "M"), c("Female", "Male")),
    HIST_CIG_DON = factor_with_labels(HIST_CIG_DON, c("N", "Y"), c("No", "Yes")),
    DIABETES_DON = factor_with_labels(DIABETES_DON, c("N", "Y"), c("No", "Yes")),
    HIST_HYPERTENS_DON = factor_with_labels(HIST_HYPERTENS_DON, c("N", "Y"), c("No", "Yes")),
    ALCOHOL_HEAVY_DON = factor_with_labels(ALCOHOL_HEAVY_DON, c("N", "Y"), c("No", "Yes")),
    COD_CAD_DON = factor(as.character(COD_CAD_DON), levels = c("1", "2", "3", "Other"),
                         labels = c("Category 1", "Category 2", "Category 3", "Other"))
  )

# -----------------------------------------------------------------------------
# Table 2: incidence and timing
# -----------------------------------------------------------------------------

safe_step_estimate <- function(ci_object, target_time, event_code) {
  curve_names <- names(ci_object)[grepl(paste0("(^| )", event_code, "$"), names(ci_object))]
  if (!length(curve_names)) stop("Event code not found in cuminc object: ", event_code)
  curve <- ci_object[[curve_names[[1]]]]
  idx <- findInterval(target_time, curve$time)
  if (idx == 0) return(0)
  curve$est[[idx]]
}

overall_ci <- cmprsk::cuminc(
  ftime = data_analysis$FGDAY_YEAR,
  fstatus = data_analysis$FGSTATUS,
  cencode = 0
)

evaluation_times <- c(3, 5, 10)
table2_cif <- data.frame(
  Time_years = evaluation_times,
  Lung_cancer_CIF = vapply(evaluation_times, safe_step_estimate, numeric(1),
                           ci_object = overall_ci, event_code = 1),
  Death_CIF = vapply(evaluation_times, safe_step_estimate, numeric(1),
                     ci_object = overall_ci, event_code = 2)
)

followup_q <- quantile(data_analysis$FGDAY_YEAR, c(0.25, 0.50, 0.75), na.rm = TRUE)
cancer_times <- data_analysis$FGDAY_YEAR[data_analysis$FGSTATUS == 1]
cancer_q <- quantile(cancer_times, c(0.25, 0.50, 0.75), na.rm = TRUE)
n_total <- nrow(data_analysis)
n_cancer <- sum(data_analysis$FGSTATUS == 1)
n_death <- sum(data_analysis$FGSTATUS == 2)
n_censored <- sum(data_analysis$FGSTATUS == 0)

table2_final <- data.frame(
  Item = c(
    "Total patients",
    "Post-transplant lung cancer, n (%)",
    "Death before lung cancer, n (%)",
    "Censored at last follow-up, n (%)",
    "Median follow-up, years (IQR)",
    "Median time to lung cancer, years (IQR)",
    "3-year cumulative incidence of lung cancer, %",
    "5-year cumulative incidence of lung cancer, %",
    "10-year cumulative incidence of lung cancer, %",
    "3-year cumulative incidence of death, %",
    "5-year cumulative incidence of death, %",
    "10-year cumulative incidence of death, %"
  ),
  Value = c(
    format(n_total, big.mark = ",", scientific = FALSE),
    sprintf("%s (%.2f%%)", format(n_cancer, big.mark = ","), 100 * n_cancer / n_total),
    sprintf("%s (%.2f%%)", format(n_death, big.mark = ","), 100 * n_death / n_total),
    sprintf("%s (%.2f%%)", format(n_censored, big.mark = ","), 100 * n_censored / n_total),
    sprintf("%.2f (%.2f-%.2f)", followup_q[[2]], followup_q[[1]], followup_q[[3]]),
    sprintf("%.2f (%.2f-%.2f)", cancer_q[[2]], cancer_q[[1]], cancer_q[[3]]),
    sprintf("%.2f", 100 * table2_cif$Lung_cancer_CIF[[1]]),
    sprintf("%.2f", 100 * table2_cif$Lung_cancer_CIF[[2]]),
    sprintf("%.2f", 100 * table2_cif$Lung_cancer_CIF[[3]]),
    sprintf("%.2f", 100 * table2_cif$Death_CIF[[1]]),
    sprintf("%.2f", 100 * table2_cif$Death_CIF[[2]]),
    sprintf("%.2f", 100 * table2_cif$Death_CIF[[3]])
  ),
  stringsAsFactors = FALSE
)

table2_basic <- data.frame(
  Total = n_total,
  Lung_cancer_n = n_cancer,
  Death_before_lung_cancer_n = n_death,
  Censored_n = n_censored,
  Followup_median_years = followup_q[[2]],
  Followup_Q1_years = followup_q[[1]],
  Followup_Q3_years = followup_q[[3]],
  Lung_cancer_time_median_years = cancer_q[[2]],
  Lung_cancer_time_Q1_years = cancer_q[[1]],
  Lung_cancer_time_Q3_years = cancer_q[[3]]
)

# -----------------------------------------------------------------------------
# Table 3: prespecified baseline candidate variables only
# -----------------------------------------------------------------------------

candidate_spec <- data.frame(
  variable = c(
    "TX_YEAR_GROUP", "AGE", "GENDER", "BMI_TCR", "DIAB", "CIG_USE", "GROUPING",
    "CREAT_TRR", "TBILI", "TOT_SERUM_ALBUM", "CPRA", "CMV_STATUS", "EBV_SEROSTATUS",
    "TX_TYPE", "ISCHTIME", "AGE_DON", "GENDER_DON", "BMI_DON_CALC", "HIST_CIG_DON",
    "DIABETES_DON", "HIST_HYPERTENS_DON", "COD_CAD_DON", "ALCOHOL_HEAVY_DON"
  ),
  label = c(
    "Transplant era", "Recipient age", "Recipient sex", "Recipient BMI", "Recipient diabetes",
    "Recipient smoking history", "Primary diagnosis", "Serum creatinine", "Total bilirubin",
    "Serum albumin", "Calculated panel-reactive antibody", "CMV serostatus", "EBV serostatus",
    "Transplant type", "Ischemic time", "Donor age", "Donor sex", "Donor BMI",
    "Donor smoking history", "Donor diabetes", "Donor hypertension", "Donor cause of death",
    "Donor heavy alcohol use"
  ),
  scale = c(NA, 10, NA, 5, NA, NA, NA, 1, 1, 1, 10, NA, NA, NA, 1, 10, NA, 5, NA, NA, NA, NA, NA),
  unit = c(
    "", "per 10 years", "", "per 5 kg/m^2", "", "", "", "per 1 mg/dL", "per 1 mg/dL",
    "per 1 g/dL", "per 10 percentage points", "", "", "", "per 1 hour", "per 10 years",
    "", "per 5 kg/m^2", "", "", "", "", ""
  ),
  stringsAsFactors = FALSE
)

missing_candidate_vars <- setdiff(candidate_spec$variable, names(data_analysis))
if (length(missing_candidate_vars) > 0) {
  stop("Candidate variables missing from dataset: ", paste(missing_candidate_vars, collapse = ", "))
}

format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

global_wald_p <- function(beta, variance) {
  if (length(beta) == 1) {
    z <- beta / sqrt(variance[[1, 1]])
    return(2 * stats::pnorm(abs(z), lower.tail = FALSE))
  }
  statistic <- tryCatch(
    drop(t(beta) %*% solve(variance, beta)),
    error = function(e) NA_real_
  )
  if (is.na(statistic)) return(NA_real_)
  stats::pchisq(statistic, df = length(beta), lower.tail = FALSE)
}

run_univariable_fg <- function(variable, label, scale_value, unit_label) {
  x <- data_analysis[[variable]]
  keep <- !is.na(x)
  dat <- data_analysis[keep, c("FGDAY_YEAR", "FGSTATUS"), drop = FALSE]
  x <- x[keep]
  n <- length(x)
  events <- sum(dat$FGSTATUS == 1)
  missing_pct <- 100 * (1 - n / nrow(data_analysis))
  if (n < 50 || length(unique(x)) < 2 || events < 5) return(NULL)

  if (is.factor(x) || is.character(x)) {
    x <- droplevels(factor(x))
    design <- model.matrix(~ x)[, -1, drop = FALSE]
    comparisons <- paste(levels(x)[-1], "vs", levels(x)[1])
    reference <- levels(x)[1]
  } else {
    if (is.na(scale_value)) scale_value <- 1
    design <- matrix(as.numeric(x) / scale_value, ncol = 1)
    colnames(design) <- variable
    comparisons <- unit_label
    reference <- ""
  }

  fit <- tryCatch(
    cmprsk::crr(
      ftime = dat$FGDAY_YEAR,
      fstatus = dat$FGSTATUS,
      cov1 = design,
      failcode = 1,
      cencode = 0
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  beta <- fit$coef
  variance <- fit$var
  se <- sqrt(diag(variance))
  coefficient_p <- 2 * stats::pnorm(abs(beta / se), lower.tail = FALSE)
  overall_p <- global_wald_p(beta, variance)
  shr <- exp(beta)
  lower <- exp(beta - 1.96 * se)
  upper <- exp(beta + 1.96 * se)

  data.frame(
    variable = variable,
    Variable = label,
    Comparison = comparisons,
    Reference = reference,
    N = n,
    Events = events,
    Missing_pct = missing_pct,
    sHR = shr,
    Lower_95CI = lower,
    Upper_95CI = upper,
    sHR_95CI = sprintf("%.2f (%.2f-%.2f)", shr, lower, upper),
    P_value = coefficient_p,
    P_value_fmt = vapply(coefficient_p, format_p, character(1)),
    Overall_P_value = overall_p,
    Overall_P_value_fmt = c(format_p(overall_p), rep("", max(0, length(beta) - 1))),
    stringsAsFactors = FALSE
  )
}

univariable_results <- bind_rows(lapply(seq_len(nrow(candidate_spec)), function(i) {
  run_univariable_fg(
    candidate_spec$variable[[i]], candidate_spec$label[[i]],
    candidate_spec$scale[[i]], candidate_spec$unit[[i]]
  )
}))

univariable_results$Variable_order <- match(univariable_results$variable, candidate_spec$variable)
univariable_results <- univariable_results %>%
  arrange(Variable_order) %>%
  select(-Variable_order)

candidate_summary <- univariable_results %>%
  group_by(variable, Variable) %>%
  summarise(
    N = first(N),
    Events = first(Events),
    Missing_pct = first(Missing_pct),
    Overall_P_value = first(Overall_P_value),
    Overall_P_value_fmt = format_p(first(Overall_P_value)),
    .groups = "drop"
  ) %>%
  mutate(variable_order = match(variable, candidate_spec$variable)) %>%
  arrange(variable_order) %>%
  select(-variable_order)

candidate_p020 <- candidate_summary %>%
  filter(!is.na(Overall_P_value), Overall_P_value < 0.20) %>%
  arrange(Overall_P_value)

excluded_variables <- data.frame(
  Category = c(
    "Outcome or outcome-derived variables",
    "Post-transplant/follow-up variables",
    "Duplicate time representation",
    "Variables with ambiguous coding or extreme missingness"
  ),
  Variables = c(
    "MALIG_P, lung_group, MALIGDAY, DEATHDAY, FGDAY, FGSTATUS, GSTATUS",
    "ACUTE_REJ_EPI, TRTREJ1Y, current/maintenance immunosuppression variables",
    "TX_YEAR (TX_YEAR_GROUP retained instead)",
    "Variables not prespecified in candidate_spec; review against the UNOS data dictionary before adding"
  ),
  Rationale = c(
    "Prevents direct outcome leakage",
    "Not consistently known at the transplantation baseline; baseline use could create temporal bias",
    "Avoids including the same calendar-time information twice",
    "Avoids unstable or uninterpretable screening"
  ),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# Publication-style CIF figures
# -----------------------------------------------------------------------------

gray_p_label <- function(p) {
  if (is.na(p)) return("Gray's test P not available")
  if (p < 0.001) return("Gray's test P < 0.001")
  paste0("Gray's test P = ", sprintf("%.3f", p))
}

plot_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
plot_linetypes <- c("solid", "dashed", "dotdash", "longdash", "twodash")

make_cif_plot <- function(data, variable, title, output_stem) {
  plot_data <- data %>%
    filter(!is.na(.data[[variable]]))
  plot_data[[variable]] <- droplevels(factor(plot_data[[variable]]))
  groups <- levels(plot_data[[variable]])

  ci <- cmprsk::cuminc(
    ftime = plot_data$FGDAY_YEAR,
    fstatus = plot_data$FGSTATUS,
    group = plot_data[[variable]],
    cencode = 0
  )
  curve_names <- names(ci)[grepl(" 1$", names(ci))]
  curve_data <- bind_rows(lapply(curve_names, function(curve_name) {
    curve <- ci[[curve_name]]
    data.frame(
      time = curve$time,
      estimate = curve$est,
      group = sub(" 1$", "", curve_name),
      stringsAsFactors = FALSE
    )
  }))
  curve_data$group <- factor(curve_data$group, levels = groups)

  p_value <- if (!is.null(ci$Tests)) ci$Tests[1, "pv"] else NA_real_
  max_estimate <- max(curve_data$estimate[curve_data$time <= 10], na.rm = TRUE)
  y_max <- max(0.03, ceiling(max_estimate * 120) / 100)
  colors <- setNames(plot_palette[seq_along(groups)], groups)
  linetypes <- setNames(plot_linetypes[seq_along(groups)], groups)

  main_plot <- ggplot(curve_data, aes(time, estimate, color = group, linetype = group)) +
    geom_step(linewidth = 1.15, direction = "hv") +
    annotate(
      "label", x = 6.4, y = y_max * 0.92,
      label = gray_p_label(p_value), size = 3.8,
      linewidth = 0, fill = scales::alpha("white", 0.88)
    ) +
    scale_color_manual(values = colors, breaks = groups) +
    scale_linetype_manual(values = linetypes, breaks = groups) +
    scale_x_continuous(breaks = c(0, 1, 3, 5, 10), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = pretty(c(0, y_max), n = 5),
      labels = scales::percent_format(accuracy = 0.1),
      expand = c(0, 0)
    ) +
    coord_cartesian(xlim = c(0, 10), ylim = c(0, y_max)) +
    labs(title = title, x = NULL, y = "Cumulative incidence of lung cancer", color = NULL, linetype = NULL) +
    theme_classic(base_size = 12.5) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      axis.title.y = element_text(face = "bold"),
      legend.position = "top",
      legend.justification = "left",
      legend.text = element_text(size = 10.5),
      legend.key.width = grid::unit(1.4, "cm"),
      plot.margin = margin(8, 12, 0, 8)
    )

  risk_times <- c(0, 1, 3, 5, 10)
  risk_table <- expand.grid(time = risk_times, group = groups, stringsAsFactors = FALSE)
  risk_table$n <- mapply(function(time, group) {
    sum(plot_data[[variable]] == group & plot_data$FGDAY_YEAR >= time, na.rm = TRUE)
  }, risk_table$time, risk_table$group)
  risk_table$group <- factor(risk_table$group, levels = rev(groups))

  risk_plot <- ggplot(risk_table, aes(time, group, label = scales::comma(n), color = group)) +
    geom_text(size = 3.4, show.legend = FALSE) +
    scale_color_manual(values = colors) +
    scale_x_continuous(limits = c(-0.45, 10.45), breaks = risk_times, expand = c(0, 0)) +
    labs(x = "Years after lung transplantation", y = NULL, subtitle = "Number at risk") +
    theme_classic(base_size = 10.5) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(color = "#333333", margin = margin(r = 8)),
      plot.subtitle = element_text(face = "bold", size = 10.5),
      plot.margin = margin(0, 12, 8, 20)
    )

  complete_plot <- main_plot / risk_plot + patchwork::plot_layout(heights = c(4.2, 1.35))
  pdf_file <- file.path(output_dir, paste0(output_stem, ".pdf"))
  png_file <- file.path(output_dir, paste0(output_stem, ".png"))
  ggsave(pdf_file, complete_plot, width = 7.6, height = 6.4, device = "pdf")
  ggsave(png_file, complete_plot, width = 7.6, height = 6.4, dpi = 320, bg = "white")

  list(main = main_plot, complete = complete_plot, p_value = p_value)
}

fig_tx <- make_cif_plot(data_analysis, "TX_TYPE", "A. Transplant type", "CIF_TX_TYPE")
fig_smoke <- make_cif_plot(data_analysis, "CIG_USE", "B. Recipient smoking history", "CIF_CIG_USE")
fig_diag <- make_cif_plot(data_analysis, "GROUPING", "C. Primary diagnosis", "CIF_GROUPING")
fig_age <- make_cif_plot(data_analysis, "AGE_GROUP", "D. Recipient age", "CIF_AGE")
fig_era <- make_cif_plot(data_analysis, "TX_YEAR_GROUP", "E. Transplant era", "CIF_ERA")

combined_figure <- (fig_tx$main | fig_smoke$main) /
  (fig_diag$main | fig_age$main) /
  fig_era$main +
  patchwork::plot_layout(heights = c(1, 1, 1)) +
  patchwork::plot_annotation(
  title = "Cumulative incidence of post-transplant lung cancer",
  theme = theme(plot.title = element_text(face = "bold", size = 17, hjust = 0.5))
)
ggsave(file.path(output_dir, "CIF_ALL_PANELS.pdf"), combined_figure,
       width = 13, height = 15, device = "pdf")
ggsave(file.path(output_dir, "CIF_ALL_PANELS.png"), combined_figure,
       width = 13, height = 15, dpi = 320, bg = "white")

# -----------------------------------------------------------------------------
# Excel workbooks
# -----------------------------------------------------------------------------

metadata <- data.frame(
  Field = c(
    "Source file", "Source MD5", "Generated on", "Total patients", "Lung cancer events",
    "Competing deaths", "Censored at last follow-up", "Study era", "Transplant era definition", "Candidate-variable rule",
    "Important exclusion rule"
  ),
  Value = c(
    normalizePath(input_file, winslash = "/", mustWork = TRUE),
    unname(tools::md5sum(input_file)),
    format(Sys.Date(), "%Y-%m-%d"),
    n_total, n_cancer, n_death, n_censored,
    "LAS era; TX_DATE >= 2005-05-04 (applied upstream)",
    "2005-2010; 2011-2015; 2016-2020; 2021-2023",
    "Prespecified transplantation-baseline variables; global Wald P values for categorical predictors",
    "Outcome-derived and post-transplant/follow-up variables were not screened"
  ),
  stringsAsFactors = FALSE
)

header_style <- createStyle(fontName = "Arial", fontColour = "#000000", fgFill = "#FFFFFF",
                            textDecoration = "bold", halign = "center", valign = "center",
                            border = c("Top", "Bottom"), borderColour = "#000000")
title_style <- createStyle(fontName = "Arial", fontColour = "#000000", fgFill = "#FFFFFF",
                           textDecoration = "bold", fontSize = 11, halign = "left", valign = "center")
note_style <- createStyle(fontName = "Arial", fontColour = "#000000", fgFill = "#FFFFFF",
                          wrapText = TRUE, valign = "top")
percent_style <- createStyle(numFmt = "0.00")
p_style <- createStyle(numFmt = "0.000")

write_styled_sheet <- function(wb, sheet_name, title, data, widths = "auto") {
  addWorksheet(wb, sheet_name, gridLines = FALSE)
  mergeCells(wb, sheet_name, cols = 1:ncol(data), rows = 1)
  writeData(wb, sheet_name, title, startRow = 1, startCol = 1)
  addStyle(wb, sheet_name, title_style, rows = 1, cols = 1:ncol(data), gridExpand = TRUE)
  setRowHeights(wb, sheet_name, rows = 1, heights = 26)
  writeData(wb, sheet_name, data, startRow = 3, headerStyle = header_style, withFilter = FALSE)
  freezePane(wb, sheet_name, firstActiveRow = 4)
  setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = widths)
  addStyle(wb, sheet_name, createStyle(fontName = "Arial", fontColour = "#000000", fgFill = "#FFFFFF"),
           rows = 4:(nrow(data) + 3), cols = 1:ncol(data), gridExpand = TRUE, stack = TRUE)
  addStyle(wb, sheet_name, createStyle(border = "Bottom", borderColour = "#000000"),
           rows = nrow(data) + 3, cols = 1:ncol(data), gridExpand = TRUE, stack = TRUE)
}

wb2 <- createWorkbook(creator = "UNOS lung cancer analysis")
write_styled_sheet(wb2, "Table2_final", "Table 2. Incidence and timing of post-transplant lung cancer",
                   table2_final, widths = c(55, 24))
write_styled_sheet(wb2, "CIF_3_5_10_years", "Cumulative incidence estimates",
                   table2_cif, widths = c(16, 24, 18))
addStyle(wb2, "CIF_3_5_10_years", createStyle(numFmt = "0.00%"),
         rows = 4:(nrow(table2_cif) + 3), cols = 2:3, gridExpand = TRUE, stack = TRUE)
write_styled_sheet(wb2, "Basic_summary", "Analysis summary (numeric values)",
                   table2_basic, widths = rep(22, ncol(table2_basic)))
write_styled_sheet(wb2, "Metadata", "Reproducibility metadata", metadata, widths = c(32, 95))
addStyle(wb2, "Metadata", note_style, rows = 4:(nrow(metadata) + 3), cols = 2, gridExpand = TRUE, stack = TRUE)
saveWorkbook(wb2, file.path(output_dir, "Table2_incidence_timing_lung_cancer.xlsx"), overwrite = TRUE)

table3_display <- univariable_results %>%
  select(Variable, Comparison, Reference, N, Events, Missing_pct, sHR_95CI,
         P_value_fmt, Overall_P_value_fmt)
names(table3_display) <- c("Variable", "Comparison", "Reference", "N", "Events", "Missing, %",
                           "sHR (95% CI)", "P value", "Overall P value")

candidate_display <- candidate_p020 %>%
  select(Variable, N, Events, Missing_pct, Overall_P_value_fmt)
names(candidate_display) <- c("Variable", "N", "Events", "Missing, %", "Overall P value")

candidate_dictionary <- candidate_spec %>%
  transmute(
    Variable_code = variable,
    Display_name = label,
    Unit_or_comparison = ifelse(unit == "", "Categorical; first listed level is reference", unit),
    Included = "Yes - prespecified baseline candidate"
  )

wb3 <- createWorkbook(creator = "UNOS lung cancer analysis")
write_styled_sheet(wb3, "Univariable_FineGray", "Table 3. Univariable Fine-Gray regression",
                   table3_display, widths = c(31, 36, 24, 12, 12, 14, 22, 14, 18))
addStyle(wb3, "Univariable_FineGray", percent_style,
         rows = 4:(nrow(table3_display) + 3), cols = 6, gridExpand = TRUE, stack = TRUE)
write_styled_sheet(wb3, "Candidate_p_lt_0.20", "Candidate variables with global P < 0.20",
                   candidate_display, widths = c(35, 14, 14, 16, 20))
write_styled_sheet(wb3, "Numeric_results", "Machine-readable univariable Fine-Gray results",
                   univariable_results, widths = "auto")
addStyle(wb3, "Numeric_results", p_style,
         rows = 4:(nrow(univariable_results) + 3), cols = c(12, 14), gridExpand = TRUE, stack = TRUE)
write_styled_sheet(wb3, "Variable_definitions", "Prespecified baseline candidate variables",
                   candidate_dictionary, widths = c(28, 38, 48, 38))
write_styled_sheet(wb3, "Excluded_variables", "Variables deliberately excluded from screening",
                   excluded_variables, widths = c(40, 85, 75))
write_styled_sheet(wb3, "Metadata", "Reproducibility metadata", metadata, widths = c(32, 95))
addStyle(wb3, "Metadata", note_style, rows = 4:(nrow(metadata) + 3), cols = 2, gridExpand = TRUE, stack = TRUE)
saveWorkbook(wb3, file.path(output_dir, "Table3_univariable_FineGray.xlsx"), overwrite = TRUE)

cat("Revised Gray/Fine-Gray analysis complete\n")
cat("Source:", normalizePath(input_file, winslash = "/", mustWork = TRUE), "\n")
cat("N:", n_total, "| Lung cancer:", n_cancer, "| Competing deaths:", n_death,
    "| Censored:", n_censored, "\n")
cat("Transplant eras: 2005-2010; 2011-2015; 2016-2020; 2021-2023\n")
cat("Candidate variables:", nrow(candidate_summary), "| Global P < 0.20:", nrow(candidate_p020), "\n")
cat("Outputs:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
