#!/usr/bin/env Rscript

local_lib <- Sys.getenv("EJCTS_R_LIB", unset = "E:/UNOS/EJCTS_reanalysis")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
build_root <- if (length(args)) args[[1]] else normalizePath("..", mustWork = TRUE)
figure_dir <- file.path(build_root, "09_figures")
table_dir <- file.path(build_root, "10_tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("cmprsk", "dplyr", "ggplot2", "jsonlite", "patchwork", "scales")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(cmprsk)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(patchwork)
  library(scales)
})

# Restrained, colorblind-safe publication palette (Okabe-Ito/Tol family).
paper_blue <- "#0072B2"
paper_sky <- "#56B4E9"
paper_green <- "#009E73"
paper_orange <- "#E69F00"
paper_vermillion <- "#D55E00"
paper_purple <- "#CC79A7"
paper_gray <- "#6B6B6B"

theme_ejcts <- theme_classic(base_family = "Arial", base_size = 10.5) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
    plot.subtitle = element_text(size = 9, hjust = 0),
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black")
  )

save_main_figure <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, "_600ppi.tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, "_preview.png")), plot, width = width, height = height,
         units = "in", dpi = 180, bg = "white")
}

cohort_file <- file.path(build_root, "02_descriptive", "2.lung.rds")
dat <- readRDS(cohort_file)
dat$time_years <- dat$FGDAY / 365.25

# Figure 2: separate panels preserve the clinically relevant scale for each outcome.
overall <- cmprsk::cuminc(dat$time_years, dat$FGSTATUS, cencode = 0)
curve_df <- function(code, label) {
  nms <- names(overall)
  pick <- nms[nms == as.character(code) | grepl(paste0("(^| )", code, "$"), nms)]
  if (!length(pick)) stop("Cumulative-incidence curve missing for event code ", code)
  obj <- overall[[pick[[1]]]]
  data.frame(Time_years = obj$time, CIF = obj$est, Outcome = label)
}
overall_df <- bind_rows(
  curve_df(1, "Lung cancer"),
  curve_df(2, "Death before lung cancer")
) %>% filter(Time_years <= 10)

step_value <- function(curve_data, t) {
  idx <- findInterval(t, curve_data$Time_years)
  if (idx == 0) 0 else curve_data$CIF[[idx]]
}
lung_curve <- overall_df %>% filter(Outcome == "Lung cancer")
death_curve <- overall_df %>% filter(Outcome == "Death before lung cancer")
followup_q <- quantile(dat$time_years, c(0.25, 0.5, 0.75), na.rm = TRUE)
cancer_q <- quantile(dat$time_years[dat$FGSTATUS == 1], c(0.25, 0.5, 0.75), na.rm = TRUE)
overall_metrics <- data.frame(
  Metric = c(
    "Total recipients", "Lung cancer", "Death before lung cancer", "Censored at last follow-up",
    "Median follow-up, years", "Follow-up Q1, years", "Follow-up Q3, years",
    "Median time to lung cancer, years", "Time to lung cancer Q1, years", "Time to lung cancer Q3, years",
    "Lung-cancer CIF at 3 years", "Lung-cancer CIF at 5 years", "Lung-cancer CIF at 10 years",
    "Competing-death CIF at 3 years", "Competing-death CIF at 5 years", "Competing-death CIF at 10 years"
  ),
  Value = c(
    nrow(dat), sum(dat$FGSTATUS == 1), sum(dat$FGSTATUS == 2), sum(dat$FGSTATUS == 0),
    followup_q[[2]], followup_q[[1]], followup_q[[3]], cancer_q[[2]], cancer_q[[1]], cancer_q[[3]],
    step_value(lung_curve, 3), step_value(lung_curve, 5), step_value(lung_curve, 10),
    step_value(death_curve, 3), step_value(death_curve, 5), step_value(death_curve, 10)
  )
)
write.csv(overall_metrics, file.path(table_dir, "overall_incidence_timing.csv"), row.names = FALSE)

p_lung <- ggplot(lung_curve, aes(Time_years, CIF)) +
  geom_step(linewidth = 0.9, colour = paper_vermillion) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(labels = percent_format(accuracy = 0.5), limits = c(0, 0.04),
                     breaks = seq(0, 0.04, 0.01), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "A  Lung cancer", x = "Years after transplantation", y = "Cumulative incidence") +
  theme_ejcts + theme(legend.position = "none")

p_death <- ggplot(death_curve, aes(Time_years, CIF)) +
  geom_step(linewidth = 0.9, colour = paper_blue, linetype = "dashed") +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.70),
                     breaks = seq(0, 0.70, 0.10), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "B  Death before lung cancer", x = "Years after transplantation", y = "Cumulative incidence") +
  theme_ejcts + theme(legend.position = "none")

fig2 <- p_lung | p_death
save_main_figure(fig2, "Figure_2_Overall_Cumulative_Incidence", 7.4, 4.3)

# Sex-disaggregated outcomes required for the EJCTS supplement.
sex_summary <- bind_rows(lapply(c("F", "M"), function(sex_code) {
  d <- dat[dat$GENDER == sex_code, , drop = FALSE]
  ci <- cmprsk::cuminc(d$time_years, d$FGSTATUS, cencode = 0)
  event_name <- names(ci)[names(ci) == "1" | grepl("(^| )1$", names(ci))][1]
  event <- ci[[event_name]]
  cif_at <- function(t) {
    idx <- findInterval(t, event$time)
    if (idx == 0) 0 else event$est[[idx]]
  }
  data.frame(
    Sex = ifelse(sex_code == "F", "Female", "Male"),
    N = nrow(d),
    Lung_cancer_n = sum(d$FGSTATUS == 1),
    Competing_death_n = sum(d$FGSTATUS == 2),
    Censored_n = sum(d$FGSTATUS == 0),
    Lung_cancer_CIF_3y = cif_at(3),
    Lung_cancer_CIF_5y = cif_at(5),
    Lung_cancer_CIF_10y = cif_at(10)
  )
}))
write.csv(sex_summary, file.path(table_dir, "sex_disaggregated_outcomes.csv"), row.names = FALSE)

# Supplementary descriptive cumulative-incidence panels; colour is redundant with line type.
dat$`Transplant type` <- factor(dat$TX_TYPE, levels = c("D", "S"), labels = c("Double lung", "Single lung"))
dat$`Smoking history` <- factor(dat$CIG_USE, levels = c("N", "Y"), labels = c("No", "Yes"))
dat$`Primary diagnosis` <- factor(dat$GROUPING, levels = c("A", "D", "Other"),
                                  labels = c("COPD", "Interstitial lung disease", "Other"))
dat$`Recipient age` <- cut(dat$AGE, breaks = c(-Inf, 55, 60, 65, 70, Inf), right = FALSE,
                           labels = c("<55", "55-59", "60-64", "65-69", ">=70"))

subgroup_panel <- function(group_name, panel_label) {
  group <- dat[[group_name]]
  keep <- !is.na(group)
  ci <- cmprsk::cuminc(dat$time_years[keep], dat$FGSTATUS[keep], group = droplevels(group[keep]), cencode = 0)
  curve_names <- names(ci)[grepl(" 1$", names(ci))]
  curve_data <- bind_rows(lapply(curve_names, function(nm) {
    obj <- ci[[nm]]
    data.frame(Time_years = obj$time, CIF = obj$est, Group = sub(" 1$", "", nm))
  })) %>% filter(Time_years <= 10)
  group_levels <- levels(droplevels(group[keep]))
  curve_data <- curve_data %>% mutate(Group = factor(Group, levels = group_levels))
  n_groups <- length(group_levels)
  colour_values <- rep(c(paper_blue, paper_vermillion, paper_green, paper_orange, paper_purple),
                       length.out = n_groups)
  line_values <- rep(c("solid", "dashed", "dotdash", "dotted", "longdash"), length.out = n_groups)
  names(colour_values) <- names(line_values) <- group_levels
  ggplot(curve_data, aes(Time_years, CIF, colour = Group, linetype = Group)) +
    geom_step(linewidth = 0.72) +
    scale_colour_manual(values = colour_values, breaks = group_levels) +
    scale_linetype_manual(values = line_values, breaks = group_levels) +
    scale_x_continuous(breaks = c(0, 5, 10), limits = c(0, 10), expand = expansion(mult = c(0, 0.01))) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, NA)) +
    labs(title = paste0(panel_label, "  ", group_name), x = "Years", y = "Lung-cancer cumulative incidence") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_ejcts + theme(legend.position = "bottom", legend.text = element_text(size = 6.7),
                        legend.key.width = grid::unit(0.8, "cm"), legend.spacing.x = grid::unit(0.1, "cm"))
}

fig_s1 <- (subgroup_panel("Transplant type", "A") | subgroup_panel("Smoking history", "B")) /
  (subgroup_panel("Primary diagnosis", "C") | subgroup_panel("Recipient age", "D"))
save_main_figure(fig_s1, "Figure_S1_Subgroup_Cumulative_Incidence", 7.2, 8.8)

# Figure 3: adjusted associations plus standardized 5-year absolute risk.
fg <- fromJSON(file.path(build_root, "04_multivariable", "multivariable_results.json"), simplifyDataFrame = TRUE)
em <- fromJSON(file.path(build_root, "05_effect_modification", "effect_modification_results.json"), simplifyDataFrame = TRUE)

forest <- as.data.frame(fg$Table4) %>%
  filter(!is.na(sHR), Reference != "Reference") %>%
  mutate(
    Label = ifelse(Category == "", Variable, paste0(Variable, ": ", Category)),
    Label = gsub("Smoking history and transplant type: ", "", Label, fixed = TRUE),
    Label = gsub("Per 10-year increase", "per 10 years", Label, fixed = TRUE),
    Label = gsub("Per 1-hour increase", "per hour", Label, fixed = TRUE),
    Label = gsub("Per 1-mg/dL increase", "per 1 mg/dL", Label, fixed = TRUE),
    Label = factor(Label, levels = rev(unique(Label)))
  )

panel_a <- ggplot(forest, aes(sHR, Label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, colour = paper_gray) +
  geom_errorbarh(aes(xmin = Lower_95CI, xmax = Upper_95CI), height = 0, linewidth = 0.55,
                 colour = paper_blue) +
  geom_point(shape = 16, size = 2, colour = paper_blue) +
  scale_x_log10(breaks = c(0.25, 1, 5, 20), limits = c(0.2, 30)) +
  labs(title = "A  Adjusted sHRs", x = "sHR (95% CI)", y = NULL) +
  theme_ejcts + theme(axis.text.y = element_text(size = 8.1), legend.position = "none")

risk5 <- as.data.frame(em$Adjusted_absolute_risks) %>%
  filter(Time_years == 5) %>%
  mutate(
    Group = factor(Group, levels = rev(c(
      "No smoking history + double lung", "Smoking history + double lung",
      "No smoking history + single lung", "Smoking history + single lung"
    ))),
    Risk_percent = 100 * Adjusted_CIF,
    Lower_percent = 100 * Lower_95CI,
    Upper_percent = 100 * Upper_95CI
  )

risk_colors <- c(
  "No smoking history + double lung" = paper_blue,
  "Smoking history + double lung" = paper_sky,
  "No smoking history + single lung" = paper_orange,
  "Smoking history + single lung" = paper_vermillion
)

panel_b <- ggplot(risk5, aes(Risk_percent, Group, colour = Group)) +
  geom_errorbarh(aes(xmin = Lower_percent, xmax = Upper_percent), height = 0, linewidth = 0.7) +
  geom_point(shape = 15, size = 2.7) +
  scale_colour_manual(values = risk_colors) +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, max(risk5$Upper_percent) * 1.08)) +
  labs(title = "B  Adjusted risk", x = "5-year incidence", y = NULL) +
  theme_ejcts + theme(axis.text.y = element_text(size = 8.1), legend.position = "none")

fig3 <- panel_a + panel_b + plot_layout(widths = c(1.55, 1))
save_main_figure(fig3, "Figure_3_Adjusted_Associations_and_Absolute_Risk", 7.8, 8.8)

sensitivity_path <- file.path(build_root, "06_sensitivity", "effect_estimates.csv")
if (file.exists(sensitivity_path)) {
  sensitivity <- read.csv(sensitivity_path, check.names = FALSE) %>%
    filter(Time_years == 5) %>%
    mutate(
      Scenario = factor(Scenario, levels = rev(unique(Scenario))),
      Comparison_group = factor(Comparison_group, levels = c(
        "Smoking history + double lung", "No smoking history + single lung", "Smoking history + single lung"
      ))
    )
  sensitivity_colors <- c(
    "Smoking history + double lung" = paper_sky,
    "No smoking history + single lung" = paper_orange,
    "Smoking history + single lung" = paper_vermillion
  )
  fig_s2 <- ggplot(sensitivity, aes(Estimate, Scenario, colour = Comparison_group)) +
    geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.45, colour = paper_gray) +
    geom_errorbarh(aes(xmin = Lower_95CI, xmax = Upper_95CI), height = 0, linewidth = 0.65) +
    geom_point(shape = 16, size = 2.1) +
    scale_colour_manual(values = sensitivity_colors) +
    scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 40)) +
    facet_wrap(~ Comparison_group, ncol = 1, scales = "free_y") +
    labs(x = "Hazard ratio at 5 years (95% CI)", y = NULL) +
    theme_ejcts + theme(axis.text.y = element_text(size = 7.5), legend.position = "none")
  save_main_figure(fig_s2, "Figure_S2_Sensitivity_Analyses", 7.2, 8.8)
}

# Figure 4 is built once the full machine-learning comparison exists.
ml_path <- file.path(build_root, "08_machine_learning", "section8_results.json")
if (file.exists(ml_path)) {
  ml <- fromJSON(ml_path, simplifyDataFrame = TRUE)
  perf <- as.data.frame(ml$Performance) %>%
    filter(Feature_set == "Core 13", Horizon_years %in% c(3, 5)) %>%
    mutate(
      Algorithm = factor(Algorithm, levels = c("Fine-Gray", "Elastic net", "Competing-risk RSF", "XGBoost")),
      Horizon = factor(paste0(Horizon_years, " years"), levels = c("3 years", "5 years"))
    )
  diffs <- as.data.frame(ml$Paired_differences) %>%
    filter(Feature_set == "Core 13", Horizon_years %in% c(3, 5)) %>%
    mutate(
      Label = paste(Algorithm, paste0("(", Horizon_years, " years)")),
      Label = factor(Label, levels = rev(unique(Label)))
    )

  alg_colors <- c("Fine-Gray" = paper_blue, "Elastic net" = paper_orange,
                  "Competing-risk RSF" = paper_green, "XGBoost" = paper_purple)
  alg_shapes <- c("Fine-Gray" = 16, "Elastic net" = 17, "Competing-risk RSF" = 15, "XGBoost" = 18)
  alg_lines <- c("Fine-Gray" = "solid", "Elastic net" = "dashed", "Competing-risk RSF" = "dotdash", "XGBoost" = "dotted")
  metric_panel <- function(y, title, y_label, percent = FALSE, hline = NULL) {
    p <- ggplot(perf, aes(Horizon, .data[[y]], group = Algorithm, colour = Algorithm,
                          shape = Algorithm, linetype = Algorithm)) +
      geom_line(linewidth = 0.65) + geom_point(size = 2.4) +
      scale_colour_manual(values = alg_colors) + scale_shape_manual(values = alg_shapes) +
      scale_linetype_manual(values = alg_lines) +
      labs(title = title, x = NULL, y = y_label) + theme_ejcts
    if (percent) p <- p + scale_y_continuous(labels = percent_format(accuracy = 0.1))
    if (!is.null(hline)) p <- p + geom_hline(yintercept = hline, linetype = "dashed", linewidth = 0.4, colour = paper_gray)
    p
  }
  p_auc <- metric_panel("AUC", "A  Discrimination", "AUC")
  p_brier <- metric_panel("Brier", "B  Prediction error", "Brier score", TRUE)
  p_slope <- metric_panel("Calibration_slope", "C  Calibration", "Calibration slope", FALSE, 1)
  p_diff <- ggplot(diffs, aes(Delta_AUC, Label, colour = Algorithm)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = paper_gray) +
    geom_errorbarh(aes(xmin = Delta_AUC_lower95, xmax = Delta_AUC_upper95), height = 0, linewidth = 0.65) +
    geom_point(shape = 16, size = 2.1) +
    scale_colour_manual(values = alg_colors) +
    guides(colour = "none") +
    labs(title = expression("D  " * Delta * "AUC vs Fine-Gray"), x = expression(Delta * "AUC (95% CI)"), y = NULL) +
    theme_ejcts + theme(legend.position = "none", axis.text.y = element_text(size = 8))

  fig4 <- (p_auc | p_brier) / (p_slope | p_diff) + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
  save_main_figure(fig4, "Figure_4_Model_Performance", 7.2, 7.8)
}

cat("EJCTS figure generation complete\n")
