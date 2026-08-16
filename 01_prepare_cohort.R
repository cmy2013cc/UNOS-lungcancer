local_lib <- Sys.getenv("EJCTS_R_LIB", unset = "E:/UNOS/EJCTS_reanalysis")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: 01_prepare_cohort.R <raw_rds> <variable_list_xlsx> <outcome_exclusion_xlsx> <output_dir>")
}

raw_file <- args[[1]]
variable_file <- args[[2]]
outcome_file <- args[[3]]
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(jsonlite)
})

set.seed(20260805)

if (!file.exists(raw_file)) stop("Raw data not found: ", raw_file)
if (!file.exists(variable_file)) stop("Variable list not found: ", variable_file)
if (!file.exists(outcome_file)) stop("Outcome exclusion list not found: ", outcome_file)

raw <- readRDS(raw_file)
vars_to_keep <- read_excel(variable_file)$Variable
vars_to_keep <- intersect(vars_to_keep[!is.na(vars_to_keep)], names(raw))
if (!length(vars_to_keep)) stop("No requested variables were found in the raw dataset.")
data <- raw[, vars_to_keep, drop = FALSE]

standardize_missing <- function(x) {
  if (is.factor(x) || is.character(x)) {
    y <- as.character(x)
    missing_codes <- c("U", "ND", "UNKNOWN", "", "NA", "NAN", "998")
    y[trimws(toupper(y)) %in% missing_codes] <- NA_character_
    return(if (is.factor(x)) factor(y) else y)
  }
  x
}
data <- data.frame(lapply(data, standardize_missing), check.names = FALSE)

data$DIAB <- as.character(data$DIAB)
data$DIAB[data$DIAB == "1"] <- "N"
data$DIAB[data$DIAB %in% c("2", "3", "4", "5")] <- "Y"
data$DIAB <- factor(data$DIAB)

for (v in intersect(c("ABO", "ABO_DON"), names(data))) {
  y <- as.character(data[[v]])
  y[y %in% c("A1", "A2")] <- "A"
  y[y %in% c("A1B", "A2B")] <- "AB"
  data[[v]] <- factor(y)
}

date_vars <- intersect(c("PX_STAT_DATE.x", "TX_DATE", "DX_DATE25", "COMPOSITE_DEATH_DATE"), names(data))
for (v in date_vars) data[[v]] <- as.Date(data[[v]])

flow <- list()
record_flow <- function(label, before, after) {
  data.frame(step = label, before = before, excluded = before - after, remaining = after)
}

n0 <- nrow(data)
cohort <- data %>% filter(TX_DATE >= as.Date("2005-05-04"), TX_DATE <= as.Date("2023-09-30"))
flow[[1]] <- record_flow("Eligible lung transplant records in the LAS era", n0, nrow(cohort))

sequential_exclude <- function(dat, predicate, label) {
  before <- nrow(dat)
  out <- dat %>% filter(!(!!rlang::enquo(predicate)))
  list(data = out, row = record_flow(label, before, nrow(out)))
}

tmp <- sequential_exclude(cohort, !is.na(MULTIORG) & MULTIORG == "Y", "Multi-organ transplant")
cohort <- tmp$data; flow[[length(flow) + 1L]] <- tmp$row
tmp <- sequential_exclude(cohort, !is.na(PREV_TX) & PREV_TX == "Y", "Previous transplant")
cohort <- tmp$data; flow[[length(flow) + 1L]] <- tmp$row
tmp <- sequential_exclude(cohort,
                          (!is.na(MALIG) & MALIG == "Y") |
                            (!is.na(MALIG_TCR) & MALIG_TCR == "Y") |
                            (!is.na(MALIG_TRR) & MALIG_TRR == "Y"),
                          "Pre-transplant recipient malignancy")
cohort <- tmp$data; flow[[length(flow) + 1L]] <- tmp$row
tmp <- sequential_exclude(cohort, !is.na(HIST_CANCER_DON) & HIST_CANCER_DON == "Y", "Donor history of malignancy")
cohort <- tmp$data; flow[[length(flow) + 1L]] <- tmp$row
tmp <- sequential_exclude(cohort, !is.na(LUNG_DONOR) & as.character(LUNG_DONOR) == "1", "Donor-derived lung malignancy flag")
cohort <- tmp$data; flow[[length(flow) + 1L]] <- tmp$row

cohort$MALIG_P <- as.character(cohort$MALIG_P)
lung_or_none <- cohort %>%
  filter((MALIG_P == "Y" & as.character(LUNG) == "1") | MALIG_P == "N")
flow[[length(flow) + 1L]] <- record_flow("Other post-transplant malignancy without lung cancer", nrow(cohort), nrow(lung_or_none))

early_cancer <- with(lung_or_none,
  MALIG_P == "Y" & !is.na(DX_DATE25) & as.numeric(DX_DATE25 - TX_DATE) < 90)
analysis_cohort <- lung_or_none[!early_cancer, , drop = FALSE]
flow[[length(flow) + 1L]] <- record_flow("Lung cancer diagnosed within 90 days", nrow(lung_or_none), nrow(analysis_cohort))

analysis_cohort$FOLLOWDAY <- as.numeric(analysis_cohort$PX_STAT_DATE.x - analysis_cohort$TX_DATE)
analysis_cohort$MALIGDAY <- as.numeric(analysis_cohort$DX_DATE25 - analysis_cohort$TX_DATE)
analysis_cohort$DEATHDAY <- as.numeric(analysis_cohort$COMPOSITE_DEATH_DATE - analysis_cohort$TX_DATE)

valid_lung_cancer <- with(analysis_cohort,
  MALIG_P == "Y" & !is.na(MALIGDAY) & (is.na(DEATHDAY) | MALIGDAY <= DEATHDAY))
postdeath_lung_cancer <- with(analysis_cohort,
  MALIG_P == "Y" & !is.na(MALIGDAY) & !is.na(DEATHDAY) & MALIGDAY > DEATHDAY)

analysis_cohort$FGSTATUS <- with(analysis_cohort,
  ifelse(valid_lung_cancer, 1L, ifelse(PSTATUS == 1, 2L, 0L)))
analysis_cohort$FGDAY <- with(analysis_cohort,
  ifelse(FGSTATUS == 1L, MALIGDAY, ifelse(FGSTATUS == 2L, DEATHDAY, FOLLOWDAY)))
analysis_cohort$MALIG_P[postdeath_lung_cancer] <- "N"

if (anyNA(analysis_cohort$FGDAY) || anyNA(analysis_cohort$FGSTATUS)) stop("Missing event time or status after outcome derivation.")
if (any(analysis_cohort$FGDAY < 0)) stop("Negative follow-up time detected.")

non_na_counts <- colSums(!is.na(analysis_cohort))
is_date <- vapply(analysis_cohort, inherits, logical(1), what = "Date")
keep_vars <- union(names(analysis_cohort)[is_date | non_na_counts > 1000], "MALIGDAY")
analysis_reduced <- analysis_cohort[, keep_vars, drop = FALSE]
analysis_reduced <- analysis_reduced[, !vapply(analysis_reduced, inherits, logical(1), what = "Date"), drop = FALSE]
analysis_reduced$MALIG_P <- as.character(analysis_reduced$MALIG_P)

factor_vars <- names(analysis_reduced)[vapply(analysis_reduced, is.factor, logical(1))]
multi_vars <- factor_vars[vapply(analysis_reduced[factor_vars], function(x) nlevels(x) > 2, logical(1))]
for (v in multi_vars) {
  tab <- table(analysis_reduced[[v]], analysis_reduced$MALIG_P)
  low_levels <- if ("Y" %in% colnames(tab)) rownames(tab)[tab[, "Y"] < 10] else rownames(tab)
  y <- as.character(analysis_reduced[[v]])
  y[y %in% low_levels] <- "Other"
  analysis_reduced[[v]] <- factor(y)
}

factor_vars <- names(analysis_reduced)[vapply(analysis_reduced, is.factor, logical(1))]
binary_vars <- factor_vars[vapply(analysis_reduced[factor_vars], function(x) nlevels(x) == 2, logical(1))]
remove_binary <- binary_vars[vapply(binary_vars, function(v) {
  tab <- table(analysis_reduced[[v]], analysis_reduced$MALIG_P)
  !"Y" %in% colnames(tab) || any(tab[, "Y"] < 10)
}, logical(1))]
analysis_reduced <- analysis_reduced[, !names(analysis_reduced) %in% remove_binary, drop = FALSE]

outcome_remove <- read_excel(outcome_file, col_names = FALSE)[[1]]
analysis_reduced <- analysis_reduced[, !names(analysis_reduced) %in% outcome_remove, drop = FALSE]

flow_df <- do.call(rbind, flow)
write.csv(flow_df, file.path(output_dir, "cohort_flow.csv"), row.names = FALSE, na = "")
write.csv(data.frame(variable = remove_binary), file.path(output_dir, "rare_binary_variables_removed.csv"), row.names = FALSE)
saveRDS(analysis_cohort, file.path(output_dir, "analysis_cohort_full.rds"))
saveRDS(analysis_reduced, file.path(output_dir, "analysis_cohort_model_ready.rds"))

metadata <- list(
  generated_on = as.character(Sys.Date()),
  seed = 20260805,
  raw_file = normalizePath(raw_file, winslash = "/"),
  raw_md5 = unname(tools::md5sum(raw_file)),
  final_n = nrow(analysis_reduced),
  lung_cancer_events = sum(analysis_reduced$FGSTATUS == 1L),
  competing_deaths = sum(analysis_reduced$FGSTATUS == 2L),
  no_lung_cancer_or_death_observed = sum(analysis_reduced$FGSTATUS == 0L),
  censored_for_time_to_event_analysis = sum(analysis_reduced$FGSTATUS == 0L),
  postdeath_lung_cancer_reclassified = sum(postdeath_lung_cancer),
  columns = ncol(analysis_reduced)
)
write_json(metadata, file.path(output_dir, "cohort_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))

cat("Cohort preparation complete\n")
cat("N =", metadata$final_n, " lung cancers =", metadata$lung_cancer_events,
    " competing deaths =", metadata$competing_deaths,
    " no lung cancer or death observed =", metadata$no_lung_cancer_or_death_observed, "\n")
