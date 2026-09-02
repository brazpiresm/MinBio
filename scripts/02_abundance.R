##============================================================================##
## 02_abundance.R
## Random-effects meta-analysis for species ABUNDANCE
##============================================================================##

source("00_functions.R")


##============================================================================##
## Settings ----
##============================================================================##

# Input files
FILES <- list(
  abundance = "data/Abundance.xlsx",
  mines     = "data/Mines.xlsx"
)

# Analysis mode
WEIGHTED  <- TRUE           # FALSE => unweighted

METHOD_SD <- "Bracken"      # SD imputation method (only used when WEIGHTED = TRUE)
# Options:
#   "Bracken"  – CV-based (recommended)
#   "Median"   – within-group median
#   "HotDeck"  – random draw from observed values
#   "Poisson"  – sqrt(mean) approximation

# Variance-covariance matrix
USE_VCOV  <- TRUE           # TRUE  => shared-control VCV matrix
                            # FALSE => sampling variance vector only

# Random effects structure
# Uncomment ONE candidate. The output folder is named automatically from
# RE_LABEL so results from different structures are kept separate.
#
## 1: study/mine 
# RANDOM_SPEC <- list(~ 1 | Study_ID / Mine_ID)
# RE_LABEL    <- "Study_Mine"
#
## 2: study/mine + genus
# RANDOM_SPEC <- list(~ 1 | Study_ID / Mine_ID, ~ 1 | genus)
# RE_LABEL    <- "Study_Mine_Genus"
#
# 3: study/mine + binomial name
RANDOM_SPEC <- list(~ 1 | Study_ID / Mine_ID, ~ 1 | binomial_name)
RE_LABEL    <- "Study_Mine_Binomial"



# Append VCV status to folder name automatically
RE_LABEL <- paste0(RE_LABEL, if (USE_VCOV) "_VCV" else "_noVCV")

# Grouping levels at which models are fitted
LEVELS <- c("realm", "taxo_group_realm")

# Low-k flag threshold for cautious interpretation 
LOW_K       <- 10  # flag if k < LOW_K
LOW_STUDIES <- 2   # flag if n_studies < LOW_STUDIES

# Output root
RUN_DIR <- run_label(WEIGHTED, METHOD_SD)
OUTROOT <- file.path("results", "abundance", RE_LABEL, RUN_DIR)
dir_make(OUTROOT)


##============================================================================##
## Prepare data ----
##============================================================================##

t0 <- prep_data(FILES$abundance, outcome = "abundance", mines_path = FILES$mines)


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
  
  message("Fitting abundance models — level: ", level)
  
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
      message("  Model failed for stratum: ", s, " (k = ", nrow(dd), ")")
      
      row <- data.frame(
        k         = nrow(dd),
        n_studies = NA_integer_,
        n_mines   = NA_integer_,
        n_sites   = NA_integer_,
        log_RR    = NA_real_,
        ci_lower  = NA_real_,
        ci_upper  = NA_real_,
        pval      = NA_real_,
        QM        = NA_real_,
        QM_pval   = NA_real_,
        QE        = NA_real_,
        QE_pval   = NA_real_,
        I2        = NA_real_,
        stringsAsFactors = FALSE
      )
      row$stratum <- s
      row$low_k   <- TRUE
      
      results[[s]] <- row
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
                   file.path(out_dir, paste0("abundance_", level, "_summary.csv")))
  saveRDS(model_store,
          file.path(out_dir, paste0("abundance_", level, "_models.rds")))
  
  save_effect_plots(out_dir, results_df, level, "abundance")
  
  n_flagged <- sum(results_df$low_k, na.rm = TRUE)
  if (n_flagged > 0)
    message("  Note: ", n_flagged, " stratum/strata flagged as low-k -> interpret with caution.")
  
  message("  Saved: ", out_dir)
}
##============================================================================##


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
                   file.path(out_dir_mod, paste0("abundance_", r, "_taxo_moderator_means.csv")))
  
  if (!is.null(res$contrasts)) {
    readr::write_csv(res$contrasts,
                     file.path(out_dir_mod, paste0("abundance_", r, "_taxo_moderator_contrasts.csv")))
    message("  Pairwise contrasts saved (QM p < 0.1).")
  } else {
    message("  No contrasts (QM p >= 0.1).")
  }
  
  mod_means[[r]] <- res$means
}

if (length(mod_means) > 0) {
  means_all <- dplyr::bind_rows(mod_means)
  readr::write_csv(means_all,
                   file.path(out_dir_mod, "abundance_taxo_moderator_means_all.csv"))
  save_taxo_moderator_plot(means_all, out_dir_mod, "abundance")
  message("Taxo moderator plots saved: ", out_dir_mod)
}

saveRDS(mod_results,
        file.path(out_dir_mod, "abundance_taxo_moderator_models.rds"))


##============================================================================##
## Realm moderator model ----
##============================================================================##
out_dir_realm <- file.path(OUTROOT, "realm_moderator")
dir_make(out_dir_realm)

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
                   file.path(out_dir_realm, "abundance_realm_moderator_means.csv"))
  
  if (!is.null(realm_res$contrast)) {
    readr::write_csv(realm_res$contrast,
                     file.path(out_dir_realm, "abundance_realm_moderator_contrast.csv"))
    message("  Realm contrast: delta lnRR = ", realm_res$contrast$delta_log_RR,
            "  p = ", realm_res$contrast$pval)
  }
  
  saveRDS(realm_res, file.path(out_dir_realm, "abundance_realm_moderator_model.rds"))
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
                     file.path(out_dir_mat, paste0("abundance_", level, "_", s, "_material_moderator_means.csv")))
    
    if (!is.null(res$contrasts)) {
      readr::write_csv(res$contrasts,
                       file.path(out_dir_mat, paste0("abundance_", level, "_", s, "_material_moderator_contrasts.csv")))
      message("  Pairwise contrasts saved (QM p < 0.1).")
    } else {
      message("  No contrasts (QM p >= 0.1).")
    }
  }
}

if (length(mat_means_list) > 0) {
  mat_means_all <- dplyr::bind_rows(mat_means_list)
  readr::write_csv(mat_means_all,
                   file.path(out_dir_mat, "abundance_material_moderator_means_all.csv"))
  save_material_moderator_plot(mat_means_all, out_dir_mat, "abundance")
  message("Material moderator plots saved: ", out_dir_mat)
}

saveRDS(mat_results,
        file.path(out_dir_mat, "abundance_material_moderator_models.rds"))





##============================================================================##
## Trait analysis for freshwater invertebrates ----
##============================================================================##
library(vcd)  # for Cramer's V 

# Traits to test based on our hypotheses
TRAITS      <- c("Disp", "Exit", "Drft")
TRAIT_OUT   <- file.path(OUTROOT, "trait_moderator")
dir_make(TRAIT_OUT)

### Load trait database (Poff et al.) ----
traits_raw <- readxl::read_xls(
  "data/poffetal_traitmatrix.xls",
  sheet = 2, skip = 1)

traits_db <- traits_raw %>%
  dplyr::select(Order, Family, Genus, dplyr::all_of(TRAITS)) %>%
  dplyr::filter(!is.na(Genus), Genus != "") %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(TRAITS), as.integer))

### Load taxonomic classification file to get family per genus ----
taxon_class <- readxl::read_excel(
  "data/abundance_taxonomic_classification.xlsx") %>%
  dplyr::select(genus, family) %>%
  dplyr::distinct()

### Subset abundance data to freshwater invertebrates ----
t_fw_inv <- t %>%
  dplyr::filter(realm == "freshwater", taxo_group == "invertebrates") %>%
  dplyr::mutate(genus = sapply(strsplit(binomial_name, " "),
                               function(x) x[1])) %>%
  dplyr::left_join(taxon_class, by = "genus")

### Match traits ----
# genus level first, family level fallback
traits_genus <- traits_db %>%
  dplyr::select(Genus, dplyr::all_of(TRAITS)) %>%
  dplyr::rename(genus = Genus)

traits_family <- traits_db %>%
  dplyr::select(Family, dplyr::all_of(TRAITS)) %>%
  dplyr::rename(family = Family) %>%
  dplyr::group_by(family) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(TRAITS),
                                 ~ as.integer(round(mean(.x, na.rm = TRUE)))),
                   .groups = "drop")

t_traits <- t_fw_inv %>%
  dplyr::left_join(traits_genus, by = "genus") %>%
  dplyr::mutate(match_level = dplyr::if_else(!is.na(Disp), "genus", NA_character_))

unmatched_idx <- is.na(t_traits$match_level)

t_traits_fam <- t_traits[unmatched_idx, ] %>%
  dplyr::select(-dplyr::all_of(TRAITS)) %>%
  dplyr::left_join(traits_family, by = "family") %>%
  dplyr::mutate(match_level = dplyr::if_else(!is.na(Disp), "family", "unmatched"))

t_traits[unmatched_idx, names(t_traits_fam)] <- t_traits_fam

### Matching summary table ----
match_summary <- t_traits %>%
  dplyr::mutate(match_level = dplyr::coalesce(match_level, "unmatched")) %>%
  dplyr::count(match_level) %>%
  dplyr::mutate(pct = round(100 * n / sum(n), 1)) %>%
  dplyr::rename(n_records = n, percent = pct)

readr::write_csv(match_summary,
                 file.path(TRAIT_OUT, "trait_matching_summary.csv"))
print(match_summary)

t_traits_matched <- t_traits %>%
  dplyr::filter(match_level %in% c("genus", "family")) %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(TRAITS), as.factor))

message("Records with trait information: ", nrow(t_traits_matched),
        " (", round(100 * nrow(t_traits_matched) / nrow(t_fw_inv), 1),
        "% of freshwater invertebrate records)")

### Taxa x traits table ----
## One row per unique genus with trait values and number of records (k).
taxa_trait_table <- t_traits_matched %>%
  dplyr::group_by(genus, family, match_level,
                  dplyr::across(dplyr::all_of(TRAITS))) %>%
  dplyr::summarise(k = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(family, genus)

readr::write_csv(taxa_trait_table,
                 file.path(TRAIT_OUT, "taxa_trait_table.csv"))

### Trait-level sample sizes ----
## For each trait, number of records per level.
trait_level_k <- dplyr::bind_rows(lapply(TRAITS, function(tr) {
  t_traits_matched %>%
    dplyr::count(.data[[tr]], name = "k") %>%
    dplyr::rename(level = dplyr::all_of(tr)) %>%
    dplyr::mutate(trait = tr, level = as.character(level)) %>%
    dplyr::select(trait, level, k)
}))

readr::write_csv(trait_level_k,
                 file.path(TRAIT_OUT, "trait_level_sample_sizes.csv"))

### Cramér's V — pairwise trait associations ----
cramer_results <- dplyr::bind_rows(
  utils::combn(TRAITS, 2, simplify = FALSE) %>%
    lapply(function(pair) {
      tbl <- table(t_traits_matched[[pair[1]]],
                   t_traits_matched[[pair[2]]])
      v   <- vcd::assocstats(tbl)$cramer
      p   <- chisq.test(tbl)$p.value
      data.frame(
        trait1   = pair[1],
        trait2   = pair[2],
        cramer_v = round(v, 3),
        chisq_p  = round(p, 4),
        stringsAsFactors = FALSE
      )
    })
)

readr::write_csv(cramer_results,
                 file.path(TRAIT_OUT, "trait_pairwise_cramers_v.csv"))
message("Cramér's V table:")
print(cramer_results)



##============================================================================##
## Trait moderator models ----
##============================================================================##
uni_results    <- list()
uni_means_list <- list()

for (tr in TRAITS) {
  
  message("Fitting univariate trait moderator: ", tr)
  
  res <- fit_trait_univariate(t_traits_matched,
                              trait        = tr,
                              random_spec  = RANDOM_SPEC,
                              use_vcov     = USE_VCOV,
                              qm_threshold = 0.1)
  
  if (is.null(res)) {
    message("  Skipped (fewer than 2 levels).")
    next
  }
  
  uni_results[[tr]]    <- res
  uni_means_list[[tr]] <- res$means
  
  message("  QM = ", res$QM, "  df = ", res$QM_df, "  p = ", res$QM_pval,
          "  QE = ", res$QE, "  QEp = ", res$QE_pval)
  
  readr::write_csv(res$means,
                   file.path(TRAIT_OUT, paste0("abundance_", tr, "_trait_univariate_means.csv")))
  
  if (!is.null(res$contrasts)) {
    readr::write_csv(res$contrasts,
                     file.path(TRAIT_OUT, paste0("abundance_", tr, "_trait_univariate_contrasts.csv")))
    message("  Pairwise contrasts saved (QM p < 0.1).")
  } else {
    message("  No contrasts (QM p >= 0.1).")
  }
}

if (length(uni_means_list) > 0) {
  uni_means_all <- dplyr::bind_rows(uni_means_list)
  readr::write_csv(uni_means_all,
                   file.path(TRAIT_OUT, "abundance_trait_univariate_means_all.csv"))
}

saveRDS(uni_results,
        file.path(TRAIT_OUT, "abundance_trait_univariate_models.rds"))


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
  
  ### -- Realm level --
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
                      paste0("abundance_distance_", mform, "_realm_models.rds")))
    plot_distance_realm(
      model_list   = realm_model_list,
      df_sub       = t_dist,
      out_path     = file.path(out_dir_dist,
                               paste0("abundance_distance_", mform, "_realm.jpg")),
      outcome_name = "abundance",
      model_form   = mform,
      weighted     = WEIGHTED
    )
  }
  
  ## -- taxo_group_realm level --
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
                      paste0("abundance_distance_", mform, "_taxo_models.rds")))
    plot_distance_taxo(
      model_list   = taxo_model_list,
      df_sub       = t_dist,
      out_dir      = out_dir_dist,
      outcome_name = "abundance",
      model_form   = mform,
      weighted     = WEIGHTED
    )
  }
  
  if (length(dist_summary_list) > 0) {
    dist_all <- dplyr::bind_rows(dist_summary_list)
    readr::write_csv(dist_all,
                     file.path(out_dir_dist,
                               paste0("abundance_distance_", mform, "_summary.csv")))
    message("Saved: abundance_distance_", mform, "_summary.csv")
  }
}




##============================================================================##
## Publication bias: Egger's test + funnel plots ----
##============================================================================##

for (level in LEVELS) {
  
  out_dir    <- file.path(OUTROOT, level)
  model_path <- file.path(
    out_dir,
    paste0("abundance_", level, "_models.rds")
  )
  
  model_store <- readRDS(model_path)
  
  run_publication_bias(
    model_store,
    out_dir,
    "abundance",
    level,
    random_spec = RANDOM_SPEC,
    use_vcov = USE_VCOV
  )
}
