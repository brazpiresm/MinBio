##============================================================================##
## 01_richness.R
## Random-effects meta-analysis for species RICHNESS
##============================================================================##

source("00_functions.R")

##============================================================================##
## Settings ----
##============================================================================##

# Input files
FILES <- list(
  richness = "data/Richness.xlsx",
  mines    = "data/Mines.xlsx"
)

# Analysis mode
WEIGHTED   <- TRUE          # FALSE => unweighted

METHOD_SD  <- "Bracken"     # SD imputation method (only used when WEIGHTED = TRUE)
# Options:
#   "Bracken"  – CV-based (recommended)
#   "Median"   – within-group median
#   "HotDeck"  – random draw from observed values
#   "Poisson"  – sqrt(mean) approximation

# Variance-covariance matrix
USE_VCOV   <- TRUE          # TRUE  => shared-control VCV matrix
                            # FALSE => sampling variance vector only

# Random effects structure
RANDOM_SPEC <- list(~ 1 | Study_ID / Mine_ID)
RE_LABEL    <- if (USE_VCOV) "Study_Mine_VCV" else "Study_Mine_noVCV"

# Grouping levels at which models are fitted
LEVELS <- c("realm", "taxo_group_realm")

# Low-k flag threshold for cautious interpretation
LOW_K       <- 10   # flag if k < LOW_K
LOW_STUDIES <- 2    # flag if n_studies < LOW_STUDIES

# Output root
RUN_DIR <- run_label(WEIGHTED, METHOD_SD)
OUTROOT <- file.path("results", "richness", RE_LABEL, RUN_DIR)
dir_make(OUTROOT)


##============================================================================##
## Prepare data ----
##============================================================================##

t0 <- prep_data(FILES$richness, outcome = "richness", mines_path = FILES$mines)


##============================================================================##
## Impute and set sampling variance ----
##============================================================================##

# Ensure reproducibility of HotDeck imputation
if (WEIGHTED && METHOD_SD == "HotDeck") set.seed(27)

if (!WEIGHTED) {
  t <- set_sampling_variance(t0, weighted = FALSE)
} else {
  t <- imp_by_realm(t0, method_SD = METHOD_SD)
  t <- set_sampling_variance(t, weighted = TRUE)
}

##============================================================================##
## Diagnostics ----
##============================================================================##

if (WEIGHTED) save_imputation_diagnostics(OUTROOT, t)
save_sampling_variance_diagnostic(OUTROOT, t)


##============================================================================##
## Fit RE models ----
##============================================================================##

for (level in LEVELS) {
  
  message("Fitting richness models — level: ", level)
  
  out_dir <- file.path(OUTROOT, level)
  dir_make(out_dir)
  
  strata      <- unique(t[[level]])
  results     <- list()
  model_store <- list()
  
  for (s in strata) {
    dd <- filter(t, .data[[level]] == s)
    
    if (nrow(dd) == 0) next
    
    fit <- fit_rma_mv(dd,
                      formula     = RR_log ~ 1,
                      random_spec = RANDOM_SPEC,
                      use_vcov    = USE_VCOV)
    
    if (is.null(fit)) {
      message("  Model failed for stratum: ", s)
      next
    }
    
    row          <- summarise_fit(fit)
    row$stratum  <- s
    row$low_k    <- row$k < LOW_K | row$n_studies < LOW_STUDIES
    
    results[[s]]     <- row
    model_store[[s]] <- fit
    
    flag <- if (row$low_k) " [LOW K]" else ""
    message("  ", s, ": lnRR = ", round(row$log_RR, 3),
            "  [", round(row$ci_lower, 3), ", ", round(row$ci_upper, 3), "]",
            "  k = ", row$k, "  studies = ", row$n_studies, flag)
  }
  
  if (length(results) == 0) next
  
  results_df <- dplyr::bind_rows(results)
  
  readr::write_csv(results_df,
                   file.path(out_dir, paste0("richness_", level, "_summary.csv")))
  saveRDS(model_store,
          file.path(out_dir, paste0("richness_", level, "_models.rds")))
  
  save_effect_plots(out_dir, results_df, level, "richness")
  
  n_flagged <- sum(results_df$low_k, na.rm = TRUE)
  if (n_flagged > 0)
    message("  Note: ", n_flagged, " stratum/strata flagged as low-k -> interpret with caution.")
  
  message("  Saved: ", out_dir)
}



##============================================================================##
## Fit taxonomic group moderator models ----
##============================================================================##

out_dir_mod <- file.path(OUTROOT, "taxo_group_moderator")
dir_make(out_dir_mod)

realms      <- unique(t$realm)
mod_results <- list()
mod_means   <- list()

for (r in realms) {
  
  message("Fitting taxo moderator — realm: ", r)
  
  res <- fit_taxo_moderator(t, realm = r,
                            random_spec  = RANDOM_SPEC,
                            use_vcov     = USE_VCOV,
                            qm_threshold = 0.1)
  
  if (is.null(res)) {
    message("  Model failed for realm: ", r)
    next
  }
  
  mod_results[[r]] <- res
  message("  QM = ", res$QM, "  df = ", res$QM_df, "  p = ", res$QM_pval)
  
  readr::write_csv(res$means,
                   file.path(out_dir_mod, paste0("richness_", r, "_taxo_moderator_means.csv")))
  
  if (!is.null(res$contrasts)) {
    readr::write_csv(res$contrasts,
                     file.path(out_dir_mod, paste0("richness_", r, "_taxo_moderator_contrasts.csv")))
    message("  Pairwise contrasts saved (QM p < 0.1).")
  } else {
    message("  No contrasts (QM p >= 0.1).")
  }
  
  mod_means[[r]] <- res$means
}

if (length(mod_means) > 0) {
  means_all <- dplyr::bind_rows(mod_means)
  readr::write_csv(means_all,
                   file.path(out_dir_mod, "richness_taxo_moderator_means_all.csv"))
  save_taxo_moderator_plot(means_all, out_dir_mod, "richness")
  message("Taxo moderator plots saved: ", out_dir_mod)
}

saveRDS(mod_results,
        file.path(out_dir_mod, "richness_taxo_moderator_models.rds"))




##============================================================================##
## Realm moderator model ----
##============================================================================##
out_dir_realm <- file.path(OUTROOT, "realm_moderator")
dir_make(out_dir_realm)

message("Fitting realm moderator (richness, both realms jointly)")

realm_res <- fit_realm_moderator(t,
                                 random_spec  = RANDOM_SPEC,
                                 use_vcov     = USE_VCOV,
                                 qm_threshold = 0.1)

if (is.null(realm_res)) {
  message("  Realm moderator model failed.")
} else {
  message("  QM = ", realm_res$QM, "  df = ", realm_res$QM_df,
          "  p = ", realm_res$QM_pval)
  
  readr::write_csv(realm_res$means,
                   file.path(out_dir_realm, "richness_realm_moderator_means.csv"))
  
  if (!is.null(realm_res$contrast)) {
    readr::write_csv(realm_res$contrast,
                     file.path(out_dir_realm, "richness_realm_moderator_contrast.csv"))
    message("  Realm contrast: delta lnRR = ", realm_res$contrast$delta_log_RR,
            "  p = ", realm_res$contrast$pval)
  }
  
  saveRDS(realm_res, file.path(out_dir_realm, "richness_realm_moderator_model.rds"))
}



##============================================================================##
## Material group moderator models ----
##============================================================================##

out_dir_mat <- file.path(OUTROOT, "material_group_moderator")
dir_make(out_dir_mat)

mat_results    <- list()
mat_means_list <- list()

for (level in LEVELS) {
  
  strata <- unique(t[[level]])
  
  for (s in strata) {
    
    message("Fitting material moderator — ", level, ": ", s)
    
    res <- fit_material_moderator(t,
                                  stratum_col  = level,
                                  stratum_val  = s,
                                  random_spec  = RANDOM_SPEC,
                                  use_vcov     = USE_VCOV,
                                  qm_threshold = 0.1)
    
    if (is.null(res)) {
      message("  Skipped (insufficient data or < 2 material groups).")
      next
    }
    
    key <- paste0(level, "__", s)
    mat_results[[key]]    <- res
    mat_means_list[[key]] <- res$means
    
    message("  QM = ", res$QM, "  df = ", res$QM_df, "  p = ", res$QM_pval,
            "  QE = ", res$QE, "  QEp = ", res$QE_pval)
    
    readr::write_csv(res$means,
                     file.path(out_dir_mat, paste0("richness_", level, "_", s, "_material_moderator_means.csv")))
    
    if (!is.null(res$contrasts)) {
      readr::write_csv(res$contrasts,
                       file.path(out_dir_mat, paste0("richness_", level, "_", s, "_material_moderator_contrasts.csv")))
      message("  Pairwise contrasts saved (QM p < 0.1).")
    } else {
      message("  No contrasts (QM p >= 0.1).")
    }
  }
}

if (length(mat_means_list) > 0) {
  mat_means_all <- dplyr::bind_rows(mat_means_list)
  readr::write_csv(mat_means_all,
                   file.path(out_dir_mat, "richness_material_moderator_means_all.csv"))
  save_material_moderator_plot(mat_means_all, out_dir_mat, "richness")
  message("Material moderator plots saved: ", out_dir_mat)
}

saveRDS(mat_results,
        file.path(out_dir_mat, "richness_material_moderator_models.rds"))







##============================================================================##
## Distance moderator ----
##============================================================================##
# Meta-regression of lnRR on log(distance_disturbed_in_meters + 1).
# log(dist + 1) retains zero-distance records. Only NAs excluded.
# Two model forms fitted for BOTH realm and taxo_group_realm levels:
#   "linear"    : RR_log ~ log_dist1
#   "quadratic" : RR_log ~ log_dist1 + log_dist1^2

out_dir_dist <- file.path(OUTROOT, "distance_moderator")
dir.create(out_dir_dist, recursive = TRUE, showWarnings = FALSE)

t_dist <- dplyr::filter(t, is.finite(distance_disturbed_in_meters)) %>%
  dplyr::mutate(log_dist1   = log(distance_disturbed_in_meters + 1),
                log_dist1_2 = log_dist1^2)

message("Distance moderator: ", nrow(t_dist), " of ", nrow(t),
        " records have distance data.")

for (mform in c("linear", "quadratic")) {
  
  message("\n--- Distance moderator: ", mform, " ---")
  
  dist_summary_list <- list()
  realm_model_list  <- list()
  taxo_model_list   <- list()
  
  fml <- if (mform == "linear") RR_log ~ log_dist1 else
    RR_log ~ log_dist1 + log_dist1_2
  
  ### ---- Realm level ----
  for (r in unique(t_dist$realm)) {
    
    dd <- dplyr::filter(t_dist, realm == r)
    if (nrow(dd) < 3) {
      message("  realm | ", r, ": k < 3, skipped.")
      next
    }
    
    fit <- suppressWarnings(
      fit_rma_mv(dd, formula = fml, random_spec = RANDOM_SPEC, use_vcov = USE_VCOV)
    )
    if (is.null(fit)) {
      message("  realm | ", r, ": model failed, skipped.")
      next
    }
    
    b    <- as.numeric(fit$b)
    se   <- sqrt(diag(fit$vb))
    z    <- b / se
    p    <- round(2 * pnorm(-abs(z)), 4)
    trms <- rownames(fit$b)
    trms[trms == "intrcpt"] <- "intercept"
    
    dist_summary_list[[paste("realm", r, mform, sep = "_")]] <- data.frame(
      level      = "realm",
      stratum    = r,
      model_form = mform,
      term       = trms,
      estimate   = round(b, 6),
      ci_lower   = round(b - 1.96 * se, 6),
      ci_upper   = round(b + 1.96 * se, 6),
      pval       = p,
      k          = nrow(dd),
      n_studies  = length(unique(dd$Study_ID)),
      QM         = round(fit$QM,  3),
      QM_df      = fit$m,
      QM_pval    = round(fit$QMp, 4),
      QE         = round(fit$QE,  0),
      QE_pval    = round(fit$QEp, 4),
      stringsAsFactors = FALSE
    )
    
    realm_model_list[[r]] <- fit
    message("  realm | ", r, ": k = ", nrow(dd), ", QM p = ", round(fit$QMp, 4))
  }
  
  if (length(realm_model_list) > 0) {
    saveRDS(realm_model_list,
            file.path(out_dir_dist,
                      paste0("richness_distance_", mform, "_realm_models.rds")))
    plot_distance_realm(
      model_list   = realm_model_list,
      df_sub       = t_dist,
      out_path     = file.path(out_dir_dist,
                               paste0("richness_distance_", mform, "_realm.jpg")),
      outcome_name = "richness",
      model_form   = mform,
      weighted     = WEIGHTED
    )
  }
  
  ## ---- taxo_group_realm level ----
  for (tgr in unique(t_dist$taxo_group_realm)) {
    
    dd <- dplyr::filter(t_dist, taxo_group_realm == tgr)
    if (nrow(dd) < 3) {
      message("  taxo_group_realm | ", tgr, ": k < 3, skipped.")
      next
    }
    
    fit <- suppressWarnings(
      fit_rma_mv(dd, formula = fml, random_spec = RANDOM_SPEC, use_vcov = USE_VCOV)
    )
    if (is.null(fit)) {
      message("  taxo_group_realm | ", tgr, ": model failed, skipped.")
      next
    }
    
    b    <- as.numeric(fit$b)
    se   <- sqrt(diag(fit$vb))
    z    <- b / se
    p    <- round(2 * pnorm(-abs(z)), 4)
    trms <- rownames(fit$b)
    trms[trms == "intrcpt"] <- "intercept"
    
    dist_summary_list[[paste("taxo", tgr, mform, sep = "_")]] <- data.frame(
      level      = "taxo_group_realm",
      stratum    = tgr,
      model_form = mform,
      term       = trms,
      estimate   = round(b, 6),
      ci_lower   = round(b - 1.96 * se, 6),
      ci_upper   = round(b + 1.96 * se, 6),
      pval       = p,
      k          = nrow(dd),
      n_studies  = length(unique(dd$Study_ID)),
      QM         = round(fit$QM,  3),
      QM_df      = fit$m,
      QM_pval    = round(fit$QMp, 4),
      QE         = round(fit$QE,  0),
      QE_pval    = round(fit$QEp, 4),
      stringsAsFactors = FALSE
    )
    
    taxo_model_list[[tgr]] <- fit
    message("  taxo_group_realm | ", tgr, ": k = ", nrow(dd),
            ", QM p = ", round(fit$QMp, 4))
  }
  
  if (length(taxo_model_list) > 0) {
    saveRDS(taxo_model_list,
            file.path(out_dir_dist,
                      paste0("richness_distance_", mform, "_taxo_models.rds")))
    plot_distance_taxo(
      model_list   = taxo_model_list,
      df_sub       = t_dist,
      out_dir      = out_dir_dist,
      outcome_name = "richness",
      model_form   = mform,
      weighted     = WEIGHTED
    )
  }
  
  if (length(dist_summary_list) > 0) {
    dist_all <- dplyr::bind_rows(dist_summary_list)
    readr::write_csv(dist_all,
                     file.path(out_dir_dist,
                               paste0("richness_distance_", mform, "_summary.csv")))
    message("Saved: richness_distance_", mform, "_summary.csv")
  }
}





##============================================================================##
## Publication bias: Egger's test + funnel plots ----
##============================================================================##

for (level in LEVELS) {
    
    out_dir    <- file.path(OUTROOT, level)
    model_path <- file.path(
      out_dir,
      paste0("richness_", level, "_models.rds")
    )
    
    model_store <- readRDS(model_path)
    
    run_publication_bias(
      model_store,
      out_dir,
      "richness",
      level,
      random_spec = RANDOM_SPEC,
      use_vcov = USE_VCOV
    )
  }
