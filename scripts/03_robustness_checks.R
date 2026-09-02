##============================================================================##
## 03_robustness_checks.R
##
## Two complementary sensitivity analyses for the meta-analytic models:
##
##  (1) GEARY – sensitivity to small-sample means (Lajeunesse 2015)
##      Re-fits models after restricting to effect sizes for which the
##      small-sample corrected standardised mean of both the disturbed and
##      control group passes Geary's rule (≥ 3). Compares to the full-data
##      model results.
##
##  (2) OUTLIERS – sensitivity to influential observations (Harrer et al. 2021)
##      Flags effect sizes whose 95% CI does not overlap the pooled model CI,
##      re-fits models after removing them, and compares to full-data results.
##
##============================================================================##

 
## !! Requires output from 01_richness.R and 02_abundance.R to run !! ##


source("00_functions.R")


library(metafor)
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)


##============================================================================##
## Settings ----
##============================================================================##

# Input files
FILES <- list(
  richness  = "data/Richness.xlsx",
  abundance = "data/Abundance.xlsx",
  mines     = "data/Mines.xlsx"
)

# Outcome 
OUTCOME <- "abundance"   # Options: "richness" or "abundance"

# Analysis mode
WEIGHTED  <- TRUE       # FALSE => unweighted

METHOD_SD <- "Bracken"  # SD imputation method (only used when WEIGHTED = TRUE)
                        # Options: "Bracken", "Median", "HotDeck", "Poisson"

# Variance-covariance matrix 
USE_VCOV  <- TRUE       # TRUE  => shared-control VCV matrix
                        # FALSE => sampling variance vector only

# Random effects structure 
if (OUTCOME == "richness") {
  RANDOM_SPEC <- list(~ 1 | Study_ID / Mine_ID)
  RE_LABEL    <- if (USE_VCOV) "Study_Mine_VCV" else "Study_Mine_noVCV"
} else {
  RANDOM_SPEC <- list(~ 1 | Study_ID / Mine_ID, ~ 1 | binomial_name)
  RE_LABEL    <- if (USE_VCOV) "Study_Mine_Binomial_VCV" else "Study_Mine_Binomial_noVCV"
}

# Grouping levels
LEVELS <- c("realm", "taxo_group_realm")


RUN_DIR <- run_label(WEIGHTED, METHOD_SD)
INROOT  <- file.path("results", OUTCOME, RE_LABEL, RUN_DIR)
OUTROOT <- file.path("results", "robustness_checks", OUTCOME, RE_LABEL, RUN_DIR)
dir_make(OUTROOT)

# Geary threshold
GEARY_THRESHOLD <- 3    

# Minimum effect sizes to fit a model 
MIN_K <- 1


##============================================================================##
## Prepare data ----
##============================================================================##

t0 <- prep_data(FILES[[OUTCOME]], outcome = OUTCOME, mines_path = FILES$mines)

if (!WEIGHTED) {
  t <- set_sampling_variance(t0, weighted = FALSE)
} else {
  t <- imp_by_realm(t0, method_SD = METHOD_SD)
  t <- set_sampling_variance(t, weighted = TRUE)
}


##============================================================================##
## Run checks per level ----
##============================================================================##

for (level in LEVELS) {
  
  message("\n", OUTCOME, " | ", RE_LABEL, " | ", RUN_DIR, " | ", level)
  
  # Read original summary and fitted models 
  res_dir     <- file.path(INROOT, level)
  orig        <- suppressWarnings(read_csv(
    file.path(res_dir, paste0(OUTCOME, "_", level, "_summary.csv")),
    show_col_types = FALSE))
  model_store <- readRDS(
    file.path(res_dir, paste0(OUTCOME, "_", level, "_models.rds")))
  
  if ("ci_lower" %in% names(orig)) orig <- rename(orig, log_RR_CI_lower = ci_lower)
  if ("ci_upper" %in% names(orig)) orig <- rename(orig, log_RR_CI_upper = ci_upper)
  if (!("stratum" %in% names(orig))) names(orig)[1] <- "stratum"
  
  out_base <- file.path(OUTROOT, level)
  dir_make(out_base)
  
  
  ## Geary --------
  
  out_geary <- file.path(out_base, "geary")
  dir_make(out_geary)
  
  df_g <- add_geary_flags(t, threshold = GEARY_THRESHOLD)
  
  pass_rates <- df_g %>%
    group_by(.data[[level]]) %>%
    summarise(n_total   = n(),
              n_pass    = sum(geary_pass, na.rm = TRUE),
              pass_rate = round(100 * n_pass / n_total, 1),
              .groups   = "drop")
  write_csv(pass_rates, file.path(out_geary, "geary_pass_rates.csv"))
  
  sens_geary <- fit_all_strata_robust(filter(df_g, geary_pass), level,
                                      RANDOM_SPEC, use_vcov = USE_VCOV,
                                      min_k = MIN_K)
  write_csv(sens_geary, file.path(out_geary, "summary_geary.csv"))
  
  comp_geary <- make_comparison_robust(orig, sens_geary, "_geary")
  write_csv(comp_geary, file.path(out_geary, "comparison_full_vs_geary.csv"))
  
  ggsave(file.path(out_geary, "forest_full_vs_geary.jpg"),
         plot_comparison_forest_robust(comp_geary, level,
                                       paste0(OUTCOME, " | ", RUN_DIR, " | ", level,
                                              "\nGeary sensitivity (threshold = ", GEARY_THRESHOLD, ")")),
         width = 10, height = max(4, 0.45 * nrow(orig) * 2 + 2), dpi = 300)
  
  
  
  ## Outliers --------
  
  out_outl <- file.path(out_base, "outliers")
  dir_make(out_outl)
  
  t_o <- t %>%
    mutate(se_study     = sqrt(SamplingVariance),
           CI_lower     = RR_log - 1.96 * se_study,
           CI_upper     = RR_log + 1.96 * se_study,
           outlier_flag = FALSE)
  
  for (s in unique(t_o[[level]])) {
    m <- model_store[[as.character(s)]]
    if (is.null(m)) next
    idx <- which(t_o[[level]] == s)
    t_o$outlier_flag[idx] <- with(t_o[idx, ],
                                  CI_upper < m$ci.lb | CI_lower > m$ci.ub)
  }
  
  outlier_counts <- t_o %>%
    group_by(.data[[level]]) %>%
    summarise(n_total     = n(),
              n_outlier   = sum(outlier_flag),
              pct_outlier = round(100 * n_outlier / n_total, 1),
              .groups     = "drop")
  write_csv(outlier_counts, file.path(out_outl, "outlier_counts.csv"))
  
  sens_outl <- fit_all_strata_robust(filter(t_o, !outlier_flag), level,
                                     RANDOM_SPEC, use_vcov = USE_VCOV,
                                     min_k = MIN_K)
  write_csv(sens_outl, file.path(out_outl, "summary_no_outliers.csv"))
  
  comp_outl <- make_comparison_robust(orig, sens_outl, "_no_outl")
  write_csv(comp_outl, file.path(out_outl, "comparison_full_vs_no_outliers.csv"))
  
  ggsave(file.path(out_outl, "forest_full_vs_no_outliers.jpg"),
         plot_comparison_forest_robust(comp_outl, level,
                                       paste0(OUTCOME, " | ", RUN_DIR, " | ", level,
                                              "\nOutlier sensitivity (Harrer et al.)")),
         width = 10, height = max(4, 0.45 * nrow(orig) * 2 + 2), dpi = 300)
  
  if (level == "realm") {
    ggsave(file.path(out_outl, "forest_outliers_highlighted.jpg"),
           plot_forest_outliers_robust(t_o, orig, level,
                                       paste0(OUTCOME, " | ", RUN_DIR, " - outliers highlighted")),
           width = 12, height = 6, dpi = 300, bg = "white")
  }
}
##============================================================================##