#!/usr/bin/env Rscript

local_lib <- Sys.getenv("EJCTS_R_LIB", unset = "E:/UNOS/EJCTS_reanalysis")
legacy_lib <- "E:/UNOS/EJCTS_reanalysis"
available_libs <- c(local_lib, legacy_lib)
available_libs <- available_libs[dir.exists(available_libs)]
if (length(available_libs)) .libPaths(c(available_libs, .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
build_root <- if (length(args)) args[[1]] else normalizePath("..", mustWork = TRUE)
figure_dir <- file.path(build_root, "09_figures_uniform")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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

if (.Platform$OS.type == "windows") {
  grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
}

# One visual system for every statistical figure.
pal <- c(
  navy = "#1F4E79", vermillion = "#C44E52", teal = "#2A9D8F",
  amber = "#E69F00", purple = "#7B6BA8", sky = "#56B4E9", gray = "#6E6E6E"
)
group_cols <- c(
  "No smoking history + double lung" = pal[["navy"]],
  "Smoking history + double lung" = pal[["teal"]],
  "No smoking history + single lung" = pal[["amber"]],
  "Smoking history + single lung" = pal[["vermillion"]]
)
group_lines <- c(
  "No smoking history + double lung" = "solid",
  "Smoking history + double lung" = "dashed",
  "No smoking history + single lung" = "dotdash",
  "Smoking history + single lung" = "longdash"
)
alg_cols <- c(
  "Fine-Gray" = pal[["navy"]], "Elastic net" = pal[["amber"]],
  "Competing-risk RSF" = pal[["teal"]], "XGBoost" = pal[["purple"]]
)
alg_shapes <- c("Fine-Gray" = 16, "Elastic net" = 17, "Competing-risk RSF" = 15, "XGBoost" = 18)
alg_lines <- c("Fine-Gray" = "solid", "Elastic net" = "dashed", "Competing-risk RSF" = "dotdash", "XGBoost" = "dotted")

theme_ejcts <- theme_classic(base_family = "Arial", base_size = 10.5) +
  theme(
    text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 11.2, hjust = 0, margin = margin(b = 6)),
    strip.background = element_blank(), strip.text = element_text(face = "bold", colour = "black", size = 9.3),
    legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 8),
    legend.key.width = grid::unit(0.9, "cm"), panel.spacing = grid::unit(0.65, "lines"),
    plot.margin = margin(6, 7, 6, 7)
  )

save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, "_600ppi.tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(figure_dir, paste0(stem, "_preview.png")), plot, width = width, height = height,
         units = "in", dpi = 180, bg = "white")
}

json_df <- function(path, key) as.data.frame(fromJSON(path, simplifyDataFrame = TRUE)[[key]])
panel_title <- function(letter, title) paste0(letter, "  ", title)
pct1 <- percent_format(accuracy = 0.1)

dat <- readRDS(file.path(build_root, "02_descriptive", "2.lung.rds"))
dat$time_years <- dat$FGDAY / 365.25

# MAIN FIGURE 2: both competing outcomes on one clinically interpretable scale.
overall <- cmprsk::cuminc(dat$time_years, dat$FGSTATUS, cencode = 0)
curve_for <- function(code, label) {
  nms <- names(overall)
  pick <- nms[nms == as.character(code) | grepl(paste0("(^| )", code, "$"), nms)][1]
  x <- overall[[pick]]
  data.frame(Time_years = x$time, CIF = x$est, Outcome = label)
}
lung <- curve_for(1, "Lung cancer") %>% filter(Time_years <= 10)
death <- curve_for(2, "Death before lung cancer") %>% filter(Time_years <= 10)
overall_df <- bind_rows(lung, death) %>%
  mutate(Outcome = factor(Outcome, levels = c("Lung cancer", "Death before lung cancer")))
outcome_cols <- c("Lung cancer" = pal[["vermillion"]], "Death before lung cancer" = pal[["navy"]])
outcome_lines <- c("Lung cancer" = "solid", "Death before lung cancer" = "dashed")
endpoints <- overall_df %>% group_by(Outcome) %>% slice_max(Time_years, n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(Label = paste0(as.character(Outcome), ": ", percent(CIF, accuracy = 0.1)),
         Label_y = CIF + ifelse(Outcome == "Lung cancer", .018, -.018))

p2 <- ggplot(overall_df, aes(Time_years, CIF, colour = Outcome, linetype = Outcome)) +
  geom_step(linewidth = .95) +
  geom_text(data = endpoints, aes(x = 9.75, y = Label_y, label = Label),
            inherit.aes = FALSE, hjust = 1, size = 3.05, colour = "black") +
  scale_colour_manual(values = outcome_cols) + scale_linetype_manual(values = outcome_lines) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10), expand = expansion(mult = c(0, .01))) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, .7, .1),
                     limits = c(0, .7), expand = expansion(mult = c(0, .01))) +
  labs(x = "Years after transplantation", y = "Cumulative incidence") +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE), linetype = guide_legend(nrow = 1, byrow = TRUE)) +
  theme_ejcts
save_figure(p2, "Figure_2_Overall_Cumulative_Incidence", 7.2, 4.8)

# MAIN FIGURE 3: standardized curves, absolute risks, and the time-varying interaction.
em_path <- file.path(build_root, "05_effect_modification", "effect_modification_results.json")
em <- fromJSON(em_path, simplifyDataFrame = TRUE)
group_short <- c(
  "No smoking history + double lung" = "No smoking + double lung",
  "Smoking history + double lung" = "Smoking + double lung",
  "No smoking history + single lung" = "No smoking + single lung",
  "Smoking history + single lung" = "Smoking + single lung"
)
short_levels <- unname(group_short[names(group_cols)])
group_cols_short <- setNames(unname(group_cols), short_levels)
group_lines_short <- setNames(unname(group_lines), short_levels)
curves <- as.data.frame(em$Standardized_curve) %>%
  mutate(Group_short = factor(unname(group_short[Group]), levels = short_levels))
risks <- as.data.frame(em$Adjusted_absolute_risks) %>%
  mutate(Group_short = factor(unname(group_short[Group]), levels = rev(short_levels)),
         Horizon = factor(paste0(Time_years, " years"), levels = paste0(c(1, 3, 5, 10), " years")),
         Risk = 100 * Adjusted_CIF, Lower = 100 * Lower_95CI, Upper = 100 * Upper_95CI)
interaction_time <- as.data.frame(em$Time_specific_interaction)

p3a <- ggplot(curves, aes(Time_years, Adjusted_CIF, colour = Group_short, linetype = Group_short)) +
  geom_ribbon(aes(ymin = Lower_95CI, ymax = Upper_95CI, fill = Group_short), alpha = .10, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = .85) +
  scale_colour_manual(values = group_cols_short) + scale_fill_manual(values = group_cols_short) + scale_linetype_manual(values = group_lines_short) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10), expand = expansion(mult = c(0, .01))) +
  scale_y_continuous(labels = percent_format(accuracy = .5), expand = expansion(mult = c(0, .05))) +
  labs(title = panel_title("A", "Standardized cumulative incidence"), x = "Years after transplantation", y = "Adjusted cumulative incidence") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) + theme_ejcts
p3b <- ggplot(risks, aes(Risk, Group_short, colour = Group_short)) +
  geom_errorbar(aes(xmin = Lower, xmax = Upper), orientation = "y", width = 0, linewidth = .65) +
  geom_point(size = 2.25) + facet_wrap(~Horizon, ncol = 2, scales = "free_x") +
  scale_colour_manual(values = group_cols_short) +
  scale_x_continuous(
    breaks = scales::breaks_pretty(n = 4),
    labels = function(x) paste0(sprintf("%.1f", x), "%"),
    expand = expansion(mult = c(.04, .12))
  ) +
  labs(title = panel_title("B", "Time-specific adjusted absolute risk"), x = "Cumulative incidence (95% CI)", y = NULL) +
  theme_ejcts + theme(
    legend.position = "none",
    axis.text.y = element_text(size = 8.2),
    axis.text.x = element_text(size = 8.2),
    strip.text = element_text(size = 8.8),
    panel.spacing.x = grid::unit(1.6, "lines"),
    panel.spacing.y = grid::unit(0.9, "lines"),
    plot.margin = margin(6, 12, 6, 7)
  )
p3c <- ggplot(interaction_time, aes(Time_years, Interaction_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = .45, colour = pal[["gray"]]) +
  geom_ribbon(aes(ymin = Lower_95CI, ymax = Upper_95CI), fill = pal[["navy"]], alpha = .12) +
  geom_line(linewidth = .8, colour = pal[["navy"]]) + geom_point(size = 2.2, colour = pal[["navy"]]) +
  scale_x_continuous(breaks = c(1, 3, 5, 10), limits = c(1, 10)) +
  scale_y_log10(breaks = c(.5, 1, 2, 5, 10, 20), limits = c(.5, 20)) +
  annotate("text", x = 1.15, y = 16.5, label = paste0("Overall interaction\nP = ", em$Formal_interaction$P_value_formatted),
           hjust = 0, vjust = 1, size = 3.0) +
  labs(title = panel_title("C", "Interaction over time"), x = "Years after transplantation", y = "Interaction ratio (95% CI)") +
  theme_ejcts + theme(legend.position = "none")
# Give the four time-specific risk panels the full page width. Keeping Panel B
# beside Panel C made the free x-axis scales and confidence intervals too dense
# at final publication size.
figure3 <- (p3a / p3b / p3c) + plot_layout(heights = c(1.05, 1.28, 0.82))
save_figure(figure3, "Figure_3_Effect_Modification_Absolute_Risk", 7.6, 10.0)

# MAIN FIGURE 4: Core13 only; calibration and paired differences are retained in Table 3/supplement.
ml_path <- file.path(build_root, "08_machine_learning", "section8_results.json")
ml <- fromJSON(ml_path, simplifyDataFrame = TRUE)
perf <- as.data.frame(ml$Performance) %>%
  filter(Feature_set == "Core 13", Horizon_years %in% c(3, 5)) %>%
  mutate(Algorithm = factor(Algorithm, levels = names(alg_cols)),
         Horizon = factor(paste0(Horizon_years, " years"), levels = c("3 years", "5 years")))
p4a <- ggplot(perf, aes(AUC, Algorithm, colour = Algorithm)) +
  geom_errorbar(aes(xmin = AUC_lower95, xmax = AUC_upper95), orientation = "y", width = 0, linewidth = .65) + geom_point(size = 2.4) +
  facet_wrap(~Horizon, ncol = 1) + scale_colour_manual(values = alg_cols) +
  labs(title = panel_title("A", "Discrimination"), x = "Cross-validated AUC (95% CI)", y = NULL) +
  theme_ejcts + theme(legend.position = "none")
p4b <- ggplot(perf, aes(100 * Brier, Algorithm, colour = Algorithm)) +
  geom_errorbar(aes(xmin = 100 * Brier_lower95, xmax = 100 * Brier_upper95), orientation = "y", width = 0, linewidth = .65) + geom_point(size = 2.4) +
  facet_wrap(~Horizon, ncol = 1, scales = "free_x") + scale_colour_manual(values = alg_cols) +
  labs(title = panel_title("B", "Prediction error"), x = "Cross-validated Brier score, % (95% CI)", y = NULL) +
  theme_ejcts + theme(legend.position = "none")
save_figure(p4a | p4b, "Figure_4_Model_Performance", 7.2, 5.6)

# SUPPLEMENTARY FIGURE S1: five prespecified descriptive CIF panels, including era.
dat$`Transplant type` <- factor(dat$TX_TYPE, levels = c("D", "S"), labels = c("Double lung", "Single lung"))
dat$`Smoking history` <- factor(dat$CIG_USE, levels = c("N", "Y"), labels = c("No", "Yes"))
dat$`Primary diagnosis` <- factor(dat$GROUPING, levels = c("A", "D", "Other"), labels = c("COPD", "Interstitial lung disease", "Other"))
dat$`Recipient age` <- cut(dat$AGE, c(-Inf, 55, 60, 65, 70, Inf), right = FALSE,
                           labels = c("<55", "55–59", "60–64", "65–69", "≥70"))
dat$`Transplant era` <- cut(dat$TX_YEAR, c(2004, 2010, 2015, 2020, Inf), right = TRUE,
                            labels = c("2005–2010", "2011–2015", "2016–2020", "2021–2023"))

subgroup_cif <- function(var, letter) {
  g <- dat[[var]]; keep <- !is.na(g); g <- droplevels(g[keep])
  ci <- cmprsk::cuminc(dat$time_years[keep], dat$FGSTATUS[keep], group = g, cencode = 0)
  nms <- names(ci)[grepl(" 1$", names(ci))]
  dd <- bind_rows(lapply(nms, function(nm) {
    z <- ci[[nm]]; data.frame(Time_years = z$time, CIF = z$est, Group = sub(" 1$", "", nm))
  })) %>% filter(Time_years <= 10) %>% mutate(Group = factor(Group, levels = levels(g)))
  cols <- setNames(rep(c(pal[["navy"]], pal[["vermillion"]], pal[["teal"]], pal[["amber"]], pal[["purple"]]), length.out = nlevels(g)), levels(g))
  ltys <- setNames(rep(c("solid", "dashed", "dotdash", "dotted", "longdash"), length.out = nlevels(g)), levels(g))
  ggplot(dd, aes(Time_years, CIF, colour = Group, linetype = Group)) + geom_step(linewidth = .72) +
    scale_colour_manual(values = cols, breaks = levels(g)) + scale_linetype_manual(values = ltys, breaks = levels(g)) +
    scale_x_continuous(breaks = c(0, 5, 10), limits = c(0, 10), expand = expansion(mult = c(0, .01))) +
    scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, .03))) +
    labs(title = panel_title(letter, var), x = "Years", y = "Lung-cancer cumulative incidence") +
    guides(colour = guide_legend(nrow = 2, byrow = TRUE), linetype = guide_legend(nrow = 2, byrow = TRUE)) +
    theme_ejcts + theme(legend.text = element_text(size = 7), legend.key.width = grid::unit(.75, "cm"))
}
s1 <- (subgroup_cif("Transplant type", "A") | subgroup_cif("Smoking history", "B")) /
      (subgroup_cif("Primary diagnosis", "C") | subgroup_cif("Recipient age", "D")) /
      subgroup_cif("Transplant era", "E")
save_figure(s1, "Figure_S1_Subgroup_Cumulative_Incidence", 7.2, 11.0)

# SUPPLEMENTARY FIGURE S2: complete adjusted Fine-Gray model.
fg <- fromJSON(file.path(build_root, "04_multivariable", "multivariable_results.json"), simplifyDataFrame = TRUE)
forest <- as.data.frame(fg$Table4) %>% filter(!is.na(Lower_95CI), Reference != "Reference") %>%
  mutate(Label = ifelse(Category == "", Variable, paste0(Variable, ": ", Category)),
         Label = gsub("Smoking history and transplant type: ", "", Label, fixed = TRUE),
         Label = gsub("Per 10-year increase", "per 10 years", Label, fixed = TRUE),
         Label = gsub("Per 1-hour increase", "per hour", Label, fixed = TRUE),
         Label = gsub("Per 1-mg/dL increase", "per 1 mg/dL", Label, fixed = TRUE),
         Label = factor(Label, levels = rev(unique(Label))))
s2 <- ggplot(forest, aes(sHR, Label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .45, colour = pal[["gray"]]) +
  geom_errorbar(aes(xmin = Lower_95CI, xmax = Upper_95CI), orientation = "y", width = 0, linewidth = .6, colour = pal[["navy"]]) +
  geom_point(size = 2.1, colour = pal[["navy"]]) + scale_x_log10(breaks = c(.25, .5, 1, 2, 5, 10, 20, 40)) +
  labs(x = "Adjusted subdistribution hazard ratio (95% CI)", y = NULL) +
  theme_ejcts + theme(legend.position = "none", axis.text.y = element_text(size = 7.6))
save_figure(s2, "Figure_S2_Multivariable_FineGray_Forest", 7.2, 8.7)

# SUPPLEMENTARY FIGURE S3: sensitivity estimates plus event retention.
sens <- read.csv(file.path(build_root, "06_sensitivity", "effect_estimates.csv"), check.names = FALSE) %>%
  filter(Time_years == 5) %>%
  mutate(Scenario = factor(Scenario, levels = rev(unique(Scenario))),
         Comparison_group = factor(Comparison_group, levels = names(group_cols)[-1]))
s3a <- ggplot(sens, aes(Estimate, Scenario, colour = Comparison_group)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = .45, colour = pal[["gray"]]) +
  geom_errorbar(aes(xmin = Lower_95CI, xmax = Upper_95CI), orientation = "y", width = 0, linewidth = .6) + geom_point(size = 2) +
  scale_colour_manual(values = group_cols[names(group_cols)[-1]]) + scale_x_log10(breaks = c(.5, 1, 2, 5, 10, 20, 40)) +
  facet_wrap(~Comparison_group, ncol = 1) + labs(title = panel_title("A", "Effect estimates"), x = "Hazard ratio at 5 years (95% CI)", y = NULL) +
  theme_ejcts + theme(legend.position = "none", axis.text.y = element_text(size = 7.1), strip.text = element_text(size = 8.2))
coh <- read.csv(file.path(build_root, "06_sensitivity", "cohort_summary.csv"), check.names = FALSE) %>%
  filter(Scenario_key %in% c("primary", "lag180", "lag365")) %>%
  mutate(Window = factor(Scenario_key, levels = c("lag365", "lag180", "primary"), labels = c(">365 days", ">180 days", ">90 days")))
primary_events <- coh$Lung_cancer_events[coh$Scenario_key == "primary"][1]
coh <- coh %>% mutate(Retention = 100 * Lung_cancer_events / primary_events,
                      Label = sprintf("%d (%.1f%%)", Lung_cancer_events, Retention),
                      Label_x = ifelse(Retention < 93, Retention + .65, Retention - .65),
                      Label_hjust = ifelse(Retention < 93, 0, 1))
s3b <- ggplot(coh, aes(Retention, Window)) +
  geom_vline(xintercept = 100, linetype = "dotted", linewidth = .45, colour = pal[["gray"]]) +
  geom_point(size = 2.8, colour = pal[["navy"]]) +
  geom_text(aes(x = Label_x, label = Label, hjust = Label_hjust), size = 2.9) +
  scale_x_continuous(limits = c(84, 103), breaks = c(85, 90, 95, 100), labels = function(x) paste0(x, "%")) +
  labs(title = panel_title("B", "Events retained"), x = "Primary events retained", y = "Exclusion window") +
  theme_ejcts + theme(legend.position = "none")
s3 <- (s3a | s3b) + plot_layout(widths = c(2.15, 1))
save_figure(s3, "Figure_S3_Sensitivity_and_Event_Retention", 7.6, 8.6)

# SUPPLEMENTARY FIGURE S4: internally validated Fine-Gray prediction performance.
pred <- fromJSON(file.path(build_root, "07_prediction", "section7_results.json"), simplifyDataFrame = TRUE)
roc <- read.csv(file.path(build_root, "07_prediction", "roc_plot_data.csv")) %>%
  filter(Horizon_years %in% c(3, 5)) %>% mutate(Horizon = factor(paste0(Horizon_years, " years"), levels = c("3 years", "5 years")))
cal <- as.data.frame(pred$Calibration) %>% filter(Horizon_years %in% c(3, 5)) %>%
  mutate(Horizon = factor(paste0(Horizon_years, " years"), levels = c("3 years", "5 years")))
dca <- as.data.frame(pred$Decision_curve) %>% filter(Horizon_years %in% c(3, 5))
dca_long <- bind_rows(
  dca %>% transmute(Horizon_years, Threshold_probability, Strategy = "Fine-Gray", Net_benefit = Model),
  dca %>% transmute(Horizon_years, Threshold_probability, Strategy = "Alert all", Net_benefit = Treat_all),
  dca %>% transmute(Horizon_years, Threshold_probability, Strategy = "Alert none", Net_benefit = Treat_none)
) %>% mutate(Horizon = factor(paste0(Horizon_years, " years"), levels = c("3 years", "5 years")))
risk_curve <- read.csv(file.path(build_root, "07_prediction", "risk_group_cif_data.csv")) %>%
  filter(Days <= 365.25 * 10) %>% mutate(Years = Days / 365.25, Risk_group = factor(Risk_group, levels = c("Low risk", "Intermediate risk", "High risk")))
risk_cols <- c("Low risk" = pal[["navy"]], "Intermediate risk" = pal[["amber"]], "High risk" = pal[["vermillion"]])

s4a <- ggplot(roc, aes(1 - Specificity, Sensitivity, colour = Horizon)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = pal[["gray"]]) +
  geom_line(linewidth = .8) + scale_colour_manual(values = c("3 years" = pal[["navy"]], "5 years" = pal[["vermillion"]])) + coord_equal() +
  labs(title = panel_title("A", "Discrimination"), x = "1 − specificity", y = "Sensitivity") + theme_ejcts
s4b <- ggplot(cal, aes(Predicted_risk, Observed_CIF, colour = Horizon)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = pal[["gray"]]) +
  geom_errorbar(aes(ymin = Observed_lower95, ymax = Observed_upper95), width = 0, linewidth = .5) + geom_line(linewidth = .7) + geom_point(size = 1.9) +
  scale_colour_manual(values = c("3 years" = pal[["navy"]], "5 years" = pal[["vermillion"]])) + scale_x_continuous(labels = pct1) + scale_y_continuous(labels = pct1) +
  labs(title = panel_title("B", "Calibration"), x = "Predicted risk", y = "Observed cumulative incidence") + theme_ejcts + theme(legend.position = "none")
s4c <- ggplot(dca_long, aes(Threshold_probability, Net_benefit, colour = Strategy, linetype = Strategy)) + geom_line(linewidth = .7) +
  facet_wrap(~Horizon, scales = "free_x") + scale_colour_manual(values = c("Fine-Gray" = pal[["navy"]], "Alert all" = pal[["gray"]], "Alert none" = "black")) +
  scale_linetype_manual(values = c("Fine-Gray" = "solid", "Alert all" = "dashed", "Alert none" = "dotted")) + scale_x_continuous(labels = pct1) +
  labs(title = panel_title("C", "Decision curves"), x = "Risk threshold", y = "Net benefit") + theme_ejcts
s4d <- ggplot(risk_curve, aes(Years, CIF, colour = Risk_group, linetype = Risk_group)) + geom_step(linewidth = .8) +
  scale_colour_manual(values = risk_cols) + scale_linetype_manual(values = c("Low risk" = "solid", "Intermediate risk" = "dashed", "High risk" = "dotdash")) +
  scale_x_continuous(breaks = seq(0, 10, 2), limits = c(0, 10), expand = expansion(mult = c(0, .01))) + scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = panel_title("D", "Risk strata"), x = "Years", y = "Cumulative incidence") + theme_ejcts
save_figure((s4a | s4b) / (s4c | s4d),
            "Figure_S4_FineGray_Prediction_Performance", 7.5, 7.8)

# SUPPLEMENTARY FIGURE S5: calibration of all algorithms and predictor sets.
mlcal <- as.data.frame(ml$Calibration) %>% filter(Horizon_years %in% c(3, 5)) %>%
  mutate(Algorithm = factor(Algorithm, levels = names(alg_cols)),
         Feature_set = factor(Feature_set, levels = c("Core 13", "Expanded 19")),
         Horizon = factor(paste0(Horizon_years, " years"), levels = c("3 years", "5 years")))
s5 <- ggplot(mlcal, aes(Predicted_risk, Observed_CIF, colour = Algorithm, linetype = Algorithm)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = pal[["gray"]]) +
  geom_errorbar(aes(ymin = Observed_lower95, ymax = Observed_upper95), width = 0, linewidth = .45, alpha = .75) +
  geom_line(linewidth = .68) + geom_point(aes(shape = Algorithm), size = 1.8) +
  facet_grid(Feature_set ~ Horizon, scales = "free") + scale_colour_manual(values = alg_cols) + scale_linetype_manual(values = alg_lines) + scale_shape_manual(values = alg_shapes) +
  scale_x_continuous(labels = pct1) + scale_y_continuous(labels = pct1) +
  labs(x = "Mean predicted risk", y = "Observed cumulative incidence") + theme_ejcts
save_figure(s5, "Figure_S5_Machine_Learning_Calibration", 7.4, 6.8)

# SUPPLEMENTARY FIGURE S6: decision curves and high-risk capture in one non-redundant figure.
mldca <- as.data.frame(ml$Decision_curve) %>% filter(Horizon_years %in% c(3, 5)) %>%
  mutate(Algorithm = factor(Algorithm, levels = names(alg_cols)),
         Feature_set = factor(Feature_set, levels = c("Core 13", "Expanded 19")),
         Horizon = factor(paste0(Horizon_years, " years"), levels = c("3 years", "5 years")))
refs <- mldca %>% select(Feature_set, Horizon, Threshold_probability, Alert_all, Alert_none) %>% distinct()
s6a <- ggplot() +
  geom_line(data = mldca, aes(Threshold_probability, Model, colour = Algorithm, linetype = Algorithm), linewidth = .68) +
  geom_line(data = refs, aes(Threshold_probability, Alert_all), colour = pal[["gray"]], linetype = "longdash", linewidth = .55) +
  geom_line(data = refs, aes(Threshold_probability, Alert_none), colour = "black", linetype = "dotted", linewidth = .55) +
  facet_grid(Feature_set ~ Horizon, scales = "free_x") + scale_colour_manual(values = alg_cols) + scale_linetype_manual(values = alg_lines) +
  scale_x_continuous(labels = pct1) + labs(title = panel_title("A", "Decision curves"), x = "Risk threshold", y = "Net benefit") + theme_ejcts
high <- as.data.frame(ml$Risk_groups) %>% filter(Risk_group == "High risk") %>%
  mutate(Algorithm = factor(Algorithm, levels = names(alg_cols)), Feature_set = factor(Feature_set, levels = c("Core 13", "Expanded 19")))
s6b <- ggplot(high, aes(Algorithm, Cancer_capture_percent, fill = Algorithm)) + geom_col(width = .68) +
  geom_text(aes(label = sprintf("%.1f%%", Cancer_capture_percent)), vjust = -.3, size = 2.8) + facet_wrap(~Feature_set) +
  scale_fill_manual(values = alg_cols) + scale_y_continuous(limits = c(0, max(high$Cancer_capture_percent) * 1.16), labels = function(x) paste0(x, "%")) +
  labs(title = panel_title("B", "Highest-risk 10%"), x = NULL, y = "Lung cancers captured") +
  theme_ejcts + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1, size = 7.5))
save_figure(s6a / s6b + plot_layout(heights = c(1.6, 1)), "Figure_S6_Decision_Curves_and_High_Risk_Capture", 7.5, 9.0)

manifest <- data.frame(
  Figure = c("Figure 2", "Figure 3", "Figure 4", paste0("Figure S", 1:6)),
  File_stem = c(
    "Figure_2_Overall_Cumulative_Incidence", "Figure_3_Effect_Modification_Absolute_Risk",
    "Figure_4_Model_Performance", "Figure_S1_Subgroup_Cumulative_Incidence",
    "Figure_S2_Multivariable_FineGray_Forest", "Figure_S3_Sensitivity_and_Event_Retention",
    "Figure_S4_FineGray_Prediction_Performance", "Figure_S5_Machine_Learning_Calibration",
    "Figure_S6_Decision_Curves_and_High_Risk_Capture"
  ),
  Width_in = c(7.2, 7.5, 7.2, 7.2, 7.2, 7.6, 7.5, 7.4, 7.5),
  Raster_ppi = 600,
  Font = "Arial",
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(figure_dir, "figure_manifest.csv"), row.names = FALSE)
cat("Uniform EJCTS figure set complete: ", figure_dir, "\n", sep = "")
