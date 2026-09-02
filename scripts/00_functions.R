##============================================================================##
## 00_functions.R
## (Required to run the other scripts)
##============================================================================##

library(dplyr)
library(tidyr)
library(readxl)
library(metafor)
library(ggplot2)
library(metagear)
library(stringr)
library(RColorBrewer)


##============================================================================##
## Data preparation ----
##============================================================================##

prep_data <- function(path, outcome = c("richness", "abundance"), mines_path) {
  outcome <- match.arg(outcome)

  t <- readxl::read_excel(path) %>%
    left_join(readxl::read_excel(mines_path)) %>%
    as_tibble() %>%
    mutate(across(where(is.character), ~ na_if(.x, "NA")))

  if (outcome == "abundance" && "binomial_name" %in% names(t)) {
    t$genus <- sapply(strsplit(t$binomial_name, " "), function(x) x[1])
    t$genus <- as.factor(t$genus)
  }

  t$taxo_group_realm <- paste(t$realm, t$taxo_group, sep = "_")

  # Remove non-numeric characters (e.g. '<', '>') from 'distance' column
  t$distance_disturbed_in_meters <- as.numeric(
    gsub("[^0-9.]", "", t$distance_disturbed_in_meters)
  )

  numeric_columns <- c("R_disturbed", "R_control",
                       "N_disturbed", "N_control",
                       "SD_disturbed", "SD_control")
  t[numeric_columns] <- lapply(
    t[numeric_columns],
    function(x) as.numeric(gsub("[^0-9.eE+-]", "", as.character(x)))
  )

  t$Study_Site_ID <- paste(t$Study_ID, t$Mine_ID, sep = "_")
  t$Effect_ID     <- seq_len(nrow(t))

  t$RR     <- t$R_disturbed / t$R_control
  t$RR_log <- log(t$RR)
}


##============================================================================##
## Imputation ----
##============================================================================##

# Impute missing SD or N values within a grouping variable using a chosen method
impute_within_group <- function(tab, var, type, grouping_var,
                                method = "Bracken") {
  tab       <- as.data.frame(tab)
  var_x     <- paste0(var, "_", type)
  imp_col   <- paste0("imputed_", var_x)
  
  if (!imp_col %in% names(tab)) tab[[imp_col]] <- "no"
  
  pieces <- lapply(split(tab, tab[[grouping_var]]), function(x) {
    n_complete <- sum(!is.na(x[[var_x]]))
    n_missing  <- sum(is.na(x[[var_x]]))
    
    # Skip if nothing to impute or no complete records to impute from
    if (n_missing == 0 || n_complete == 0) return(x)
    
    x[[imp_col]][is.na(x[[var_x]])] <- grouping_var
    
    if (method == "Median") {
      x[[var_x]][is.na(x[[var_x]])] <- median(x[[var_x]], na.rm = TRUE)
      
    } else if (method == "Bracken") {
      x <- metagear::impute_SD(
        x,
        columnSDnames = var_x,
        columnXnames  = c(paste0("R_", type), paste0("N_", type)),
        method        = "Bracken1992"
      )
      
    } else if (method == "HotDeck") {
      x <- metagear::impute_SD(
        x,
        columnSDnames = var_x,
        columnXnames  = c(paste0("R_", type), paste0("N_", type)),
        method        = "HotDeck"
      )
      
    } else if (method == "Poisson") {
      missing_idx <- is.na(x[[var_x]])
      x[[var_x]][missing_idx] <- sqrt(
        as.numeric(x[[paste0("R_", type)]])[missing_idx]
      )
    }
    x
  })
  
  result <- do.call(rbind, pieces)
  n_imp  <- sum(result[[imp_col]] == grouping_var, na.rm = TRUE)
  cat("  [", grouping_var, "] imputed", n_imp, var_x, "records\n")
  
  tibble::as_tibble(result)
}


## Impute missing SDs and sample sizes separately within each realm.
# Available SD imputation methods (set method_SD in the analysis script):
#   "Bracken"  – CV-based from complete cases (recommended, minimises bias;
#                Bracken 1992; Lajeunesse 2016)
#   "Median"   – simple within-realm median
#   "HotDeck"  – random draw from observed values in the same realm
#   "Poisson"  – sqrt(mean) approximation
imp_by_realm <- function(t, method_SD = "Bracken") {
  t$imputed_N_disturbed  <- "no"
  t$imputed_N_control    <- "no"
  t$imputed_SD_disturbed <- "no"
  t$imputed_SD_control   <- "no"
  
  t <- impute_within_group(t, "N",  "disturbed", "realm", method = "Median")
  t <- impute_within_group(t, "N",  "control",   "realm", method = "Median")
  t <- impute_within_group(t, "SD", "disturbed", "realm", method = method_SD)
  t <- impute_within_group(t, "SD", "control",   "realm", method = method_SD)
  
  # Stop if any values are still missing after imputation
  cols_to_check <- c("SD_disturbed", "SD_control", "N_disturbed", "N_control")
  for (col in cols_to_check) {
    n_missing <- sum(is.na(t[[col]]))
    if (n_missing > 0) {
      stop(
        n_missing, " record(s) still have missing ", col,
        " after realm-level imputation. ",
        "Check for realms where all records have missing ", col, ".",
        call. = FALSE
      )
    }
  }
  t
}

##============================================================================##
## Sampling variance ----
##============================================================================##

# Set sampling variance 
# When weighted = FALSE every record receives variance = 1 (unweighted analysis).
set_sampling_variance <- function(t, weighted = TRUE) {
  if (weighted) {
    t$SamplingVariance <-
      t$SD_disturbed^2 / (t$N_disturbed * t$R_disturbed^2) +
      t$SD_control^2   / (t$N_control   * t$R_control^2)
  } else {
    t$SamplingVariance <- 1
  }
  t
}


##============================================================================##
## Variance-covariance matrix ----
##============================================================================##
# Ensure a matrix is positive-definite. If not, compute the nearest PD matrix.
# Adapted from https://github.com/anabenlop/Island_Rule/blob/master/Scripts/000_Functions.R
PDfunc <- function(m) {
  if (corpcor::is.positive.definite(m)) {
    list(mat = as.matrix(m))$mat
  } else {
    as.matrix(Matrix::nearPD(m)$mat)
  }
}


# Build the block-diagonal sampling-variance–covariance matrix for a modelling subset.
#
# Diagonal elements  : per-effect sampling variances.
# Off-diagonal elements: shared-control covariance for effects that share the
# same Control_ID within the subset (= SD_control^2 / (N_control * R_control^2))
#
# When use_vcov = FALSE an identity matrix is returned, effectively treating
# all effects as independent
#
# Adapted from https://github.com/anabenlop/Island_Rule/blob/master/Scripts/000_Functions.R
build_vcov <- function(df,
                       var_col        = "SamplingVariance",
                       control_id_col = "Control_ID",
                       mean_c         = "R_control",
                       sd_c           = "SD_control",
                       n_c            = "N_control",
                       use_vcov       = TRUE) {

  df <- df[is.finite(df[[var_col]]), , drop = FALSE]
  if (nrow(df) == 0L) return(list(df = df, V = matrix(0, 0, 0)))

  if (!use_vcov) {
    return(list(df = df, V = diag(nrow(df))))
  }

  ids <- unique(df[[control_id_col]])
  # Reorder rows so they align with block-diagonal construction
  idx <- unlist(lapply(ids, function(id) which(df[[control_id_col]] == id)),
                use.names = FALSE)
  df  <- df[idx, , drop = FALSE]

  blocks <- lapply(ids, function(id) {
    blk <- df[df[[control_id_col]] == id, , drop = FALSE]
    m   <- blk[[mean_c]][1]
    s   <- blk[[sd_c]][1]
    n   <- blk[[n_c]][1]
    V   <- matrix(s^2 / (n * m^2), nrow = nrow(blk), ncol = nrow(blk))
    diag(V) <- blk[[var_col]]
    V
  })

  V <- PDfunc(as.matrix(Matrix::bdiag(blocks)))
  list(df = df, V = V)
}


##============================================================================##
## Model fitting ----
##============================================================================##

# Fit a single rma.mv model with the shared-control VCV matrix
# or using the sampling variance vector if use_vcov = F
fit_rma_mv <- function(df, formula, random_spec, use_vcov = TRUE) {
  
  # Unweighted: all variances = 1
  if (!use_vcov || all(df$SamplingVariance == 1)) {
    df2 <- df[is.finite(df$SamplingVariance), , drop = FALSE]
    V   <- df2$SamplingVariance
  } else {
    prep <- build_vcov(df)
    df2  <- prep$df
    V    <- prep$V
  }
  
  if (nrow(df2) == 0L) return(NULL)
  
  fit <- try(
    metafor::rma.mv(formula, V = V, data = df2,
                    random  = random_spec, method = "REML",
                    control = list(optimizer = "optim",
                                   optmethod = "Nelder-Mead",
                                   maxit     = 2000),
                    verbose = FALSE),
    silent = TRUE
  )
  
  if (inherits(fit, "try-error")) return(NULL)
  attr(fit, "._df_used") <- df2
  fit
}


# Extract a tidy summary row from a fitted rma.mv object
summarise_fit <- function(fit) {
  if (is.null(fit)) return(NULL)
  df_used <- attr(fit, "._df_used")
  data.frame(
    k         = fit$k,
    n_studies = dplyr::n_distinct(df_used$Study_ID),
    n_mines   = dplyr::n_distinct(paste(df_used$Study_ID, df_used$Mine_ID)),
    n_sites   = dplyr::n_distinct(
      paste(df_used$Study_ID, df_used$Mine_ID, df_used$Disturbed_site_ID)
    ),
    log_RR    = as.numeric(fit$b[1]),
    ci_lower  = fit$ci.lb,
    ci_upper  = fit$ci.ub,
    pval      = round(fit$pval, 3),
    QM        = round(fit$QM, 3),
    QM_pval   = fit$QMp,
    QE        = round(fit$QE, 0),
    QE_pval   = fit$QEp,
    I2        = 100 * fit$sigma2[1] /
      (fit$sigma2[1] + mean(df_used$SamplingVariance, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

##============================================================================##
## Output helpers ----
##============================================================================##

dir_make <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  normalizePath(path, mustWork = FALSE)
}

# Derive the results sub-folder name from the analysis settings
run_label <- function(weighted, method_SD) {
  if (!weighted) "unweighted" else method_SD
}


##============================================================================##
## Diagnostic plots ----
##============================================================================##

save_imputation_diagnostics <- function(out_dir, t) {
  vars <- list(
    list(col = "SD_disturbed", imp = "imputed_SD_disturbed"),
    list(col = "SD_control",   imp = "imputed_SD_control"),
    list(col = "N_disturbed",  imp = "imputed_N_disturbed"),
    list(col = "N_control",    imp = "imputed_N_control")
  )
  for (v in vars) {
    p <- ggplot(t %>% filter(.data[[v$imp]] == "no")) +
      geom_histogram(aes(x = log10(.data[[v$col]]), fill = "original"),
                     alpha = 1, binwidth = 0.5) +
      geom_histogram(aes(x = log10(.data[[v$col]]), fill = "incl_imputed"),
                     data = t, alpha = 0.3, binwidth = 0.5) +
      scale_fill_manual(values = c(original = "grey40", incl_imputed = "grey70")) +
      theme_minimal() +
      labs(x = paste0("log10(", v$col, ")"), y = "Count", fill = NULL)
    ggsave(file.path(out_dir, paste0("diag_", v$col, ".pdf")), p, width = 7, height = 5)
  }
}


save_sampling_variance_diagnostic <- function(out_dir, t) {
  p <- ggplot(t) +
    geom_histogram(aes(x = log10(SamplingVariance)), binwidth = 0.5) +
    theme_minimal() +
    labs(x = "log10(SamplingVariance)", y = "Count")
  ggsave(file.path(out_dir, "diag_sampling_variance.pdf"), p, width = 7, height = 5)
  readr::write_csv(t, file.path(out_dir, "data_with_variance.csv"))
}


##============================================================================##
## Effect-size plots (for quick checks) ----
##============================================================================##

# Back-transform lnRR to percentage change
pct_change <- function(log_rr) (exp(log_rr) - 1) * 100


save_effect_plots <- function(out_dir, results_df, level, outcome_name) {
  pal       <- RColorBrewer::brewer.pal(9, "Set1")[c(2, 4, 3)]
  cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))

  df <- results_df %>%
    filter(!is.na(log_RR)) %>%
    mutate(
      label_p  = ifelse(pval < 0.01, "p < 0.01",
                        paste0("p = ", formatC(pval, format = "f", digits = 2))),
      pct      = pct_change(log_RR),
      pct_lo   = pct_change(ci_lower),
      pct_hi   = pct_change(ci_upper),
      y_label  = paste0(cap_first(gsub("_", " ", stratum)), " (k=", k, ")")
    )

  # --- lnRR ---
  p1 <- ggplot(df, aes(y = y_label, x = log_RR,
                       xmin = ci_lower, xmax = ci_upper)) +
    theme_bw(base_size = 22) +
    geom_vline(xintercept = 0, lty = 2, colour = "grey60", linewidth = 1.2) +
    geom_errorbarh(height = 0, linewidth = 1.2) +
    geom_point(size = 4) +
    geom_text(aes(label = label_p), vjust = -0.8, size = 4.5, colour = "black") +
    labs(x = "lnRR", y = NULL) +
    theme(axis.title.y = element_blank())

  ggsave(file.path(out_dir, paste0(outcome_name, "_", level, "_lnRR.jpg")),
         p1, width = 180, height = 120, units = "mm", dpi = 300)

  # --- % change ---
  p2 <- ggplot(df, aes(y = y_label, x = pct,
                       xmin = pct_lo, xmax = pct_hi)) +
    theme_bw(base_size = 22) +
    geom_vline(xintercept = 0, lty = 2, colour = "grey60", linewidth = 1.2) +
    geom_errorbarh(height = 0, linewidth = 1.2) +
    geom_point(size = 4) +
    geom_text(aes(label = label_p), vjust = -0.8, size = 4.5, colour = "black") +
    labs(x = "\u0394 (%)", y = NULL) +
    theme(axis.title.y = element_blank())

  ggsave(file.path(out_dir, paste0(outcome_name, "_", level, "_pct_change.jpg")),
         p2, width = 180, height = 120, units = "mm", dpi = 300)
}
##============================================================================##


##============================================================================##
## Robustness checks functions ----
##============================================================================##

# --- Geary's small-sample corrected standardised mean (Lajeunesse 2015) ---
geary_z <- function(mean, sd, n) {
  n  <- as.numeric(n)
  sd <- as.numeric(sd)
  mu <- as.numeric(mean)
  ifelse(
    is.na(mu) | is.na(sd) | is.na(n) | sd <= 0 | n <= 0,
    NA_real_,
    (mu / sd) * (4 * n^(3/2)) / (1 + 4 * n)
  )
}


# --- Flag effect sizes that fail Geary's rule ---
add_geary_flags <- function(df, threshold = 3) {
  df %>%
    mutate(
      geary_disturbed = geary_z(R_disturbed, SD_disturbed, N_disturbed),
      geary_control   = geary_z(R_control,   SD_control,   N_control),
      geary_min       = pmin(geary_disturbed, geary_control, na.rm = FALSE),
      geary_pass      = !is.na(geary_min) & geary_min >= threshold
    )
}


# --- Fit sensitivity models ---
fit_all_strata_robust <- function(df, level, random_spec,
                                  use_vcov = TRUE, min_k = 1) {
  strata <- unique(df[[level]])
  rows <- lapply(strata, function(s) {
    dd <- df[df[[level]] == s, , drop = FALSE]
    if (nrow(dd) < min_k) {
      row <- data.frame(
        k = nrow(dd), n_studies = NA_integer_, n_mines = NA_integer_,
        n_sites = NA_integer_, log_RR = NA_real_, ci_lower = NA_real_,
        ci_upper = NA_real_, pval = NA_real_, QM = NA_real_,
        QM_pval = NA_real_, QE = NA_real_, QE_pval = NA_real_, I2 = NA_real_
      )
    } else {
      fit <- fit_rma_mv(dd, formula = RR_log ~ 1,
                        random_spec = random_spec, use_vcov = use_vcov)
      row <- summarise_fit(fit)
      if (is.null(row)) {
        row <- data.frame(
          k = nrow(dd), n_studies = NA_integer_, n_mines = NA_integer_,
          n_sites = NA_integer_, log_RR = NA_real_, ci_lower = NA_real_,
          ci_upper = NA_real_, pval = NA_real_, QM = NA_real_,
          QM_pval = NA_real_, QE = NA_real_, QE_pval = NA_real_, I2 = NA_real_
        )
      }
    }
    row$stratum <- as.character(s)
    row
  })
  out <- dplyr::bind_rows(rows)
  out[, c("stratum", setdiff(names(out), "stratum"))]
}


# --- Build a comparison table: original results vs. sensitivity subset ---
make_comparison_robust <- function(orig, sens, suffix) {
  # normalise CI column names from summarise_fit() to match the orig convention
  if ("ci_lower" %in% names(sens)) sens <- dplyr::rename(sens, log_RR_CI_lower = ci_lower)
  if ("ci_upper" %in% names(sens)) sens <- dplyr::rename(sens, log_RR_CI_upper = ci_upper)
  
  rename_cols <- function(df, sfx) {
    cols <- setdiff(names(df), "stratum")
    names(df)[names(df) %in% cols] <- paste0(cols, sfx)
    df
  }
  dplyr::left_join(
    rename_cols(orig, "_orig"),
    rename_cols(sens, suffix),
    by = "stratum"
  ) %>%
    dplyr::mutate(
      d_log_RR      = .data[[paste0("log_RR", suffix)]] - log_RR_orig,
      pct_change_RR = 100 * (.data[[paste0("log_RR", suffix)]] - log_RR_orig) /
        abs(log_RR_orig)
    )
}


# --- Forest plot: full dataset vs. sensitivity subset ---
plot_comparison_forest_robust <- function(comp, level, title_str) {
  long <- comp %>%
    dplyr::select(stratum,
                  log_RR_orig, log_RR_CI_lower_orig, log_RR_CI_upper_orig,
                  dplyr::matches("^log_RR(?!_CI).*(?<!orig)$", perl = TRUE),
                  dplyr::matches("^log_RR_CI_lower.*(?<!orig)$", perl = TRUE),
                  dplyr::matches("^log_RR_CI_upper.*(?<!orig)$", perl = TRUE)) %>%
    tidyr::pivot_longer(
      -stratum,
      names_to      = c(".value", "set"),
      names_pattern = "^(log_RR(?:_CI_lower|_CI_upper)?)_(orig|.+)$"
    ) %>%
    dplyr::filter(!is.na(log_RR)) %>%
    dplyr::mutate(set = ifelse(set == "orig", "Full dataset", "Sensitivity"))
  
  off  <- 0.22
  long <- long %>%
    dplyr::group_by(stratum) %>%
    dplyr::mutate(idx = dplyr::cur_group_id()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(ypos = ifelse(set == "Full dataset", idx + off, idx - off))
  
  ord <- unique(long$stratum[order(long$idx)])
  
  ggplot2::ggplot(long, ggplot2::aes(y = ypos, x = log_RR,
                                     xmin = log_RR_CI_lower,
                                     xmax = log_RR_CI_upper,
                                     shape = set)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
    ggplot2::geom_errorbarh(height = 0, linewidth = 0.8) +
    ggplot2::geom_point(size = 2.8, fill = "white", stroke = 0.7) +
    ggplot2::scale_shape_manual(
      values = c("Full dataset" = 16, "Sensitivity" = 17), name = NULL) +
    ggplot2::scale_y_continuous(
      breaks = c(seq_along(ord) + off, seq_along(ord) - off),
      labels = c(paste0(ord, " - full"), paste0(ord, " - sensitivity"))
    ) +
    ggplot2::labs(x = "lnRR (95% CI)", y = NULL, title = title_str) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = element_blank(),
                   panel.grid.minor   = element_blank(),
                   legend.position    = "right")
}


# --- Forest plot with outlier points highlighted in red (realm level) ---
plot_forest_outliers_robust <- function(df_flagged, orig_summary, level,
                                        title_str) {
  df_plot <- df_flagged %>%
    dplyr::arrange(.data[[level]], RR_log) %>%
    dplyr::group_by(.data[[level]]) %>%
    dplyr::mutate(
      plottingID = dplyr::row_number(),
      se         = sqrt(SamplingVariance),
      lowerCI    = RR_log - 1.96 * se,
      upperCI    = RR_log + 1.96 * se,
      status     = ifelse(outlier_flag, "Outlier", "Included")
    ) %>%
    dplyr::ungroup()
  
  model_pts <- orig_summary %>%
    dplyr::mutate(
      realm      = stratum,
      plottingID = -1,
      lowerCI    = log_RR_CI_lower,
      upperCI    = log_RR_CI_upper,
      logRR      = log_RR
    )
  
  pal_realm <- c(terrestrial = "#1b9e77", freshwater = "#7570b3")
  
  plots <- lapply(unique(df_plot[[level]]), function(s) {
    dd  <- dplyr::filter(df_plot, .data[[level]] == s)
    mp  <- dplyr::filter(model_pts, realm == s)
    col <- if (s %in% names(pal_realm)) pal_realm[[s]] else "steelblue"
    
    x_vals <- c(dd$RR_log, mp$logRR)
    x_buf  <- diff(range(x_vals, na.rm = TRUE)) * 0.1
    xlim_v <- range(x_vals, na.rm = TRUE) + c(-x_buf, x_buf)
    
    ggplot2::ggplot(dd, ggplot2::aes(y = plottingID, x = RR_log)) +
      ggplot2::geom_linerange(
        ggplot2::aes(xmin = lowerCI, xmax = upperCI, colour = status),
        linewidth = 0.4, alpha = 0.5) +
      ggplot2::geom_point(ggplot2::aes(colour = status), size = 1.6) +
      ggplot2::scale_colour_manual(
        values = c(Included = "black", Outlier = "red"), name = NULL) +
      ggplot2::geom_vline(xintercept = 0, colour = "black") +
      ggplot2::geom_pointrange(
        data = mp,
        ggplot2::aes(y = plottingID, x = logRR, xmin = lowerCI, xmax = upperCI),
        colour = col, shape = 18, size = 1.2, linewidth = 0.6) +
      ggplot2::geom_vline(data = mp, ggplot2::aes(xintercept = logRR),
                          colour = col, linewidth = 0.9) +
      ggplot2::geom_vline(data = mp, ggplot2::aes(xintercept = lowerCI),
                          colour = col, linetype = "dashed") +
      ggplot2::geom_vline(data = mp, ggplot2::aes(xintercept = upperCI),
                          colour = col, linetype = "dashed") +
      ggplot2::scale_x_continuous(limits = xlim_v) +
      ggplot2::labs(
        x        = "lnRR",
        y        = NULL,
        subtitle = paste0(toupper(substring(s, 1, 1)), substring(s, 2))) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(axis.text.y      = element_blank(),
                     axis.ticks.y     = element_blank(),
                     legend.position  = "right")
  })
  
  cowplot::plot_grid(plotlist = plots, nrow = 1,
                     labels    = paste0("(", letters[seq_along(plots)], ")"),
                     label_size = 12, align = "hv")
}
##============================================================================##



##============================================================================##
## Taxonomic group moderator models ----
##============================================================================##

# Fit one model per realm with taxo_group as a fixed-effect moderator.
fit_taxo_moderator <- function(df, realm, random_spec,
                               use_vcov     = TRUE,
                               qm_threshold = 0.1) {
  
  dd <- dplyr::filter(df, realm == !!realm)
  if (nrow(dd) == 0) return(NULL)
  
  dd$taxo_group <- as.factor(dd$taxo_group)
  lvls          <- levels(dd$taxo_group)
  
  # --- fit WITH intercept to get QM omnibus test ---
  fit_int <- fit_rma_mv(dd, formula = RR_log ~ taxo_group,
                        random_spec = random_spec, use_vcov = use_vcov)
  if (is.null(fit_int)) return(NULL)
  
  QM   <- round(fit_int$QM,  3)
  QMdf <- fit_int$m
  QMp  <- fit_int$QMp
  QE   <- round(fit_int$QE,  0)
  QEp  <- fit_int$QEp
  
  # --- fit WITHOUT intercept to get per-level mean estimates ---
  fit_noint <- fit_rma_mv(dd, formula = RR_log ~ 0 + taxo_group,
                          random_spec = random_spec, use_vcov = use_vcov)
  if (is.null(fit_noint)) return(NULL)
  
  coefs    <- as.numeric(fit_noint$b)
  se_coefs <- sqrt(diag(fit_noint$vb))
  zvals    <- coefs / se_coefs
  pvals <- 2 * pnorm(-abs(zvals))
  
  means <- data.frame(
    taxo_group       = lvls,
    log_RR           = coefs,
    ci_lower         = coefs - 1.96 * se_coefs,
    ci_upper         = coefs + 1.96 * se_coefs,
    pval             = pvals,
    k                = as.integer(table(dd$taxo_group)[lvls]),
    realm            = realm,
    QM               = QM,
    QM_df            = QMdf,
    QM_pval          = QMp,
    QE               = QE,
    QE_pval          = QEp,
    pct_change       = round((exp(coefs) - 1) * 100, 1),
    pct_change_lower = round((exp(coefs - 1.96 * se_coefs) - 1) * 100, 1),
    pct_change_upper = round((exp(coefs + 1.96 * se_coefs) - 1) * 100, 1),
    stringsAsFactors = FALSE
  )
  
  # --- pairwise contrasts from no-intercept fit (only if QM p < threshold) ---
  contrasts_df <- NULL
  if (QMp < qm_threshold && length(lvls) > 1) {
    pairs <- utils::combn(lvls, 2, simplify = FALSE)
    contrasts_df <- dplyr::bind_rows(lapply(pairs, function(p) {
      i <- which(lvls == p[1])
      j <- which(lvls == p[2])
      cvec    <- rep(0, length(coefs))
      cvec[i] <-  1
      cvec[j] <- -1
      delta   <- as.numeric(t(cvec) %*% coefs)
      se_d    <- as.numeric(sqrt(t(cvec) %*% fit_noint$vb %*% cvec))
      zval    <- delta / se_d
      data.frame(
        group1       = p[1],
        group2       = p[2],
        delta_log_RR = round(delta, 4),
        se           = round(se_d,  4),
        z            = round(zval,  3),
        pval = 2 * pnorm(-abs(zval)),
        stringsAsFactors = FALSE
      )
    }))
    contrasts_df$pval_holm <- p.adjust(contrasts_df$pval, method = "holm")
  }
  
  list(model_int   = fit_int,
       model_noint = fit_noint,
       means       = means,
       contrasts   = contrasts_df,
       QM          = QM,
       QM_df       = QMdf,
       QM_pval     = QMp,
       QE          = QE,
       QE_pval     = QEp,
       realm       = realm)
}

# Forest plot of marginal means per taxo_group from the moderator model
save_taxo_moderator_plot <- function(means_df, out_dir, outcome_name) {
  cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  
  means_df <- means_df %>%
    dplyr::filter(!is.na(log_RR)) %>%
    dplyr::mutate(
      label_p  = ifelse(pval < 0.01, "p < 0.01",
                        paste0("p = ", formatC(pval, format = "f", digits = 2))),
      label_QM = paste0("QM p = ", formatC(QM_pval, format = "f", digits = 3)),
      y_label  = paste0(cap_first(taxo_group), " (k=", k, ")")
    )
  
  realms <- unique(means_df$realm)
  
  for (r in realms) {
    df_r     <- dplyr::filter(means_df, realm == r)
    qm_label <- unique(df_r$label_QM)
    
    # lnRR plot
    p1 <- ggplot2::ggplot(df_r,
                          ggplot2::aes(y = y_label, x = log_RR,
                                       xmin = ci_lower, xmax = ci_upper)) +
      ggplot2::theme_bw(base_size = 22) +
      ggplot2::geom_vline(xintercept = 0, lty = 2,
                          colour = "grey60", linewidth = 1.2) +
      ggplot2::geom_errorbarh(height = 0, linewidth = 1.2) +
      ggplot2::geom_point(size = 4) +
      ggplot2::geom_text(ggplot2::aes(label = label_p),
                         vjust = -0.8, size = 4.5, colour = "black") +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = qm_label, hjust = 1.1, vjust = 1.5,
                        size = 5, fontface = "italic") +
      ggplot2::labs(x = "lnRR", y = NULL,
                    title = cap_first(r)) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank())
    
    ggplot2::ggsave(
      file.path(out_dir, paste0(outcome_name, "_", r, "_taxo_moderator_lnRR.jpg")),
      p1, width = 180, height = 120, units = "mm", dpi = 300)
    
    # % change plot
    p2 <- ggplot2::ggplot(df_r,
                          ggplot2::aes(y = y_label, x = pct_change,
                                       xmin = pct_change_lower, xmax = pct_change_upper)) +
      ggplot2::theme_bw(base_size = 22) +
      ggplot2::geom_vline(xintercept = 0, lty = 2,
                          colour = "grey60", linewidth = 1.2) +
      ggplot2::geom_errorbarh(height = 0, linewidth = 1.2) +
      ggplot2::geom_point(size = 4) +
      ggplot2::geom_text(ggplot2::aes(label = label_p),
                         vjust = -0.8, size = 4.5, colour = "black") +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = qm_label, hjust = 1.1, vjust = 1.5,
                        size = 5, fontface = "italic") +
      ggplot2::labs(x = "\u0394 (%)", y = NULL,
                    title = cap_first(r)) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank())
    
    ggplot2::ggsave(
      file.path(out_dir, paste0(outcome_name, "_", r, "_taxo_moderator_pct.jpg")),
      p2, width = 180, height = 120, units = "mm", dpi = 300)
  }
}



##============================================================================##
## Realm moderator model ----
##============================================================================##
# Fits ONE joint model across both realms with realm as a fixed-effect
# moderator, to compare the two reals

fit_realm_moderator <- function(df, random_spec,
                                use_vcov     = TRUE,
                                qm_threshold = 0.1) {
  
  dd <- df
  dd$realm <- as.factor(dd$realm)
  lvls     <- levels(dd$realm)
  
  if (length(lvls) < 2) return(NULL)
  
  if (!use_vcov || all(dd$SamplingVariance == 1)) {
    dd2 <- dd[is.finite(dd$SamplingVariance), , drop = FALSE]
    V   <- dd2$SamplingVariance
  } else {
    prep <- build_vcov(dd)
    dd2  <- prep$df
    V    <- prep$V
  }
  if (nrow(dd2) == 0L) return(NULL)
  
  fit_int <- try(
    metafor::rma.mv(RR_log ~ realm, V = V, data = dd2,
                    random = ~ realm | Study_Site_ID,
                    struct = "DIAG", method = "REML",
                    control = list(optimizer = "optim",
                                   optmethod = "Nelder-Mead", maxit = 2000)),
    silent = TRUE
  )
  if (inherits(fit_int, "try-error")) {
    message("  fit_int failed: ", attr(fit_int, "condition")$message)
    return(NULL)
  }
  
  QM   <- round(fit_int$QM,  3)
  QMdf <- fit_int$m
  QMp  <- fit_int$QMp
  QE   <- round(fit_int$QE,  0)
  QEp  <- fit_int$QEp
  
  fit_noint <- try(
    metafor::rma.mv(RR_log ~ 0 + realm, V = V, data = dd2,
                    random = ~ realm | Study_Site_ID,
                    struct = "DIAG", method = "REML",
                    control = list(optimizer = "optim",
                                   optmethod = "Nelder-Mead", maxit = 2000)),
    silent = TRUE
  )
  if (inherits(fit_noint, "try-error")) {
    message("  fit_noint failed: ", attr(fit_noint, "condition")$message)
    return(NULL)
  }
  
  coef_names  <- rownames(fit_noint$b)
  coef_levels <- sub("^realm", "", coef_names)
  match_idx   <- match(lvls, coef_levels)
  
  cat("  [fit_realm_moderator] coefficient <-> level mapping check:\n")
  print(data.frame(level = lvls, coef_name = coef_names[match_idx]))
  
  if (anyNA(match_idx) || !setequal(coef_levels, lvls)) {
    stop("fit_realm_moderator: coefficient names don't match realm levels.\n",
         "  levels found: ", paste(lvls, collapse = ", "), "\n",
         "  coefs found : ", paste(coef_names, collapse = ", "), call. = FALSE)
  }
  
  if (!is.null(fit_noint$tau2)) {
    cat("  [fit_realm_moderator] tau2 by realm (compare vs separate-model sigma2):\n")
    print(data.frame(tau2 = fit_noint$tau2, sqrt_tau2 = sqrt(fit_noint$tau2)))
  }
  
  coefs    <- as.numeric(fit_noint$b)[match_idx]
  se_coefs <- sqrt(diag(fit_noint$vb))[match_idx]
  pvals    <- 2 * pnorm(-abs(coefs / se_coefs))
  
  means <- data.frame(
    realm            = lvls,
    log_RR           = coefs,
    ci_lower         = coefs - 1.96 * se_coefs,
    ci_upper         = coefs + 1.96 * se_coefs,
    pval             = pvals,
    k                = as.integer(table(dd2$realm)[lvls]),
    QM = QM, QM_df = QMdf, QM_pval = QMp, QE = QE, QE_pval = QEp,
    pct_change       = round((exp(coefs) - 1) * 100, 1),
    pct_change_lower = round((exp(coefs - 1.96 * se_coefs) - 1) * 100, 1),
    pct_change_upper = round((exp(coefs + 1.96 * se_coefs) - 1) * 100, 1),
    stringsAsFactors = FALSE
  )
  
  contrast_df <- NULL
  if (length(lvls) == 2) {
    vb_r  <- fit_noint$vb[match_idx, match_idx, drop = FALSE]
    cvec  <- c(1, -1)
    delta <- as.numeric(t(cvec) %*% coefs)
    se_d  <- as.numeric(sqrt(t(cvec) %*% vb_r %*% cvec))
    zval  <- delta / se_d
    contrast_df <- data.frame(
      group1 = lvls[1], group2 = lvls[2],
      delta_log_RR = round(delta, 4), se = round(se_d, 4),
      z = round(zval, 3), pval = 2 * pnorm(-abs(zval)), 3,
      stringsAsFactors = FALSE
    )
  }
  
  list(model_int = fit_int, model_noint = fit_noint,
       means = means, contrast = contrast_df,
       QM = QM, QM_df = QMdf, QM_pval = QMp, QE = QE, QE_pval = QEp)
}




##============================================================================##
## Material group moderator models ----
##============================================================================##
# Fit one model per stratum (realm or taxo_group_realm) with material_group
# as a fixed-effect moderator. The "others" category is dropped before fitting
# because it pools heterogeneous commodities with no expected common effect.
fit_material_moderator <- function(df, stratum_col, stratum_val, random_spec,
                                   use_vcov     = TRUE,
                                   qm_threshold = 0.1) {
  
  dd <- dplyr::filter(df, .data[[stratum_col]] == stratum_val)
  dd <- dplyr::filter(dd, !is.na(material_group), material_group != "others")
  if (nrow(dd) == 0) return(NULL)
  
  dd$material_group <- droplevels(as.factor(dd$material_group))
  lvls              <- levels(dd$material_group)
  if (length(lvls) < 2) return(NULL)   # need at least 2 groups to test moderator
  
  # --- fit WITH intercept to get QM omnibus test ---
  fit_int <- fit_rma_mv(dd, formula = RR_log ~ material_group,
                        random_spec = random_spec, use_vcov = use_vcov)
  if (is.null(fit_int)) return(NULL)
  
  QM   <- round(fit_int$QM,  3)
  QMdf <- fit_int$m
  QMp  <- fit_int$QMp
  QE   <- round(fit_int$QE,  0)
  QEp  <- fit_int$QEp
  
  # --- fit WITHOUT intercept to get per-level mean estimates ---
  fit_noint <- fit_rma_mv(dd, formula = RR_log ~ 0 + material_group,
                          random_spec = random_spec, use_vcov = use_vcov)
  if (is.null(fit_noint)) return(NULL)
  
  coefs    <- as.numeric(fit_noint$b)
  se_coefs <- sqrt(diag(fit_noint$vb))
  zvals    <- coefs / se_coefs
  pvals    <- 2 * pnorm(-abs(zvals))
  
  means <- data.frame(
    material_group   = lvls,
    log_RR           = coefs,
    ci_lower         = coefs - 1.96 * se_coefs,
    ci_upper         = coefs + 1.96 * se_coefs,
    pval             = pvals,
    k                = as.integer(table(dd$material_group)[lvls]),
    stratum          = stratum_val,
    stratum_col      = stratum_col,
    QM               = QM,
    QM_df            = QMdf,
    QM_pval          = QMp,
    QE               = QE,
    QE_pval          = QEp,
    pct_change       = round((exp(coefs) - 1) * 100, 1),
    pct_change_lower = round((exp(coefs - 1.96 * se_coefs) - 1) * 100, 1),
    pct_change_upper = round((exp(coefs + 1.96 * se_coefs) - 1) * 100, 1),
    stringsAsFactors = FALSE
  )
  
  # --- pairwise contrasts from no-intercept fit (only if QM p < threshold) ---
  contrasts_df <- NULL
  if (QMp < qm_threshold && length(lvls) > 1) {
    pairs <- utils::combn(lvls, 2, simplify = FALSE)
    contrasts_df <- dplyr::bind_rows(lapply(pairs, function(p) {
      i <- which(lvls == p[1])
      j <- which(lvls == p[2])
      cvec    <- rep(0, length(coefs))
      cvec[i] <-  1
      cvec[j] <- -1
      delta   <- as.numeric(t(cvec) %*% coefs)
      se_d    <- as.numeric(sqrt(t(cvec) %*% fit_noint$vb %*% cvec))
      zval    <- delta / se_d
      data.frame(
        group1       = p[1],
        group2       = p[2],
        delta_log_RR = round(delta, 4),
        se           = round(se_d,  4),
        z            = round(zval,  3),
        pval         = 2 * pnorm(-abs(zval)),
        stratum      = stratum_val,
        stratum_col  = stratum_col,
        stringsAsFactors = FALSE
      )
    }))
    contrasts_df$pval_holm <- p.adjust(contrasts_df$pval, method = "holm")
  }
  
  list(model_int   = fit_int,
       model_noint = fit_noint,
       means       = means,
       contrasts   = contrasts_df,
       QM          = QM,
       QM_df       = QMdf,
       QM_pval     = QMp,
       QE          = QE,
       QE_pval     = QEp,
       stratum     = stratum_val,
       stratum_col = stratum_col)
}

# Forest plot of marginal means per material_group from the moderator model.
save_material_moderator_plot <- function(means_df, out_dir, outcome_name) {
  cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  
  means_df <- means_df %>%
    dplyr::filter(!is.na(log_RR)) %>%
    dplyr::mutate(
      label_p  = ifelse(pval < 0.01, "p < 0.01",
                        paste0("p = ", formatC(pval, format = "f", digits = 2))),
      label_QM = paste0("QM p = ", formatC(QM_pval, format = "f", digits = 3)),
      y_label  = paste0(cap_first(material_group), " (k=", k, ")")
    )
  
  strata <- unique(means_df$stratum)
  
  for (s in strata) {
    df_s     <- dplyr::filter(means_df, stratum == s)
    qm_label <- unique(df_s$label_QM)
    
    # lnRR plot
    p1 <- ggplot2::ggplot(df_s,
                          ggplot2::aes(y = y_label, x = log_RR,
                                       xmin = ci_lower, xmax = ci_upper)) +
      ggplot2::theme_bw(base_size = 22) +
      ggplot2::geom_vline(xintercept = 0, lty = 2,
                          colour = "grey60", linewidth = 1.2) +
      ggplot2::geom_errorbarh(height = 0, linewidth = 1.2) +
      ggplot2::geom_point(size = 4) +
      ggplot2::geom_text(ggplot2::aes(label = label_p),
                         vjust = -0.8, size = 4.5, colour = "black") +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = qm_label, hjust = 1.1, vjust = 1.5,
                        size = 5, fontface = "italic") +
      ggplot2::labs(x = "lnRR", y = NULL,
                    title = cap_first(s)) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank())
    
    ggplot2::ggsave(
      file.path(out_dir, paste0(outcome_name, "_", s, "_material_moderator_lnRR.jpg")),
      p1, width = 180, height = 120, units = "mm", dpi = 300)
    
    # % change plot
    p2 <- ggplot2::ggplot(df_s,
                          ggplot2::aes(y = y_label, x = pct_change,
                                       xmin = pct_change_lower, xmax = pct_change_upper)) +
      ggplot2::theme_bw(base_size = 22) +
      ggplot2::geom_vline(xintercept = 0, lty = 2,
                          colour = "grey60", linewidth = 1.2) +
      ggplot2::geom_errorbarh(height = 0, linewidth = 1.2) +
      ggplot2::geom_point(size = 4) +
      ggplot2::geom_text(ggplot2::aes(label = label_p),
                         vjust = -0.8, size = 4.5, colour = "black") +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = qm_label, hjust = 1.1, vjust = 1.5,
                        size = 5, fontface = "italic") +
      ggplot2::labs(x = "\u0394 (%)", y = NULL,
                    title = cap_first(s)) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank())
    
    ggplot2::ggsave(
      file.path(out_dir, paste0(outcome_name, "_", s, "_material_moderator_pct.jpg")),
      p2, width = 180, height = 120, units = "mm", dpi = 300)
  }
}




##============================================================================##
## Trait moderator models ----
##============================================================================##
fit_trait_univariate <- function(df, trait, random_spec,
                                 use_vcov     = TRUE,
                                 qm_threshold = 0.1) {
  
  dd <- df[!is.na(df[[trait]]), ]
  dd[[trait]] <- droplevels(as.factor(dd[[trait]]))
  lvls        <- levels(dd[[trait]])
  if (length(lvls) < 2) return(NULL)
  
  fml_int   <- as.formula(paste("RR_log ~", trait))
  fml_noint <- as.formula(paste("RR_log ~ 0 +", trait))
  
  fit_int <- fit_rma_mv(dd, formula = fml_int,
                        random_spec = random_spec, use_vcov = use_vcov)
  if (is.null(fit_int)) return(NULL)
  
  QM   <- round(fit_int$QM,  3)
  QMdf <- fit_int$m
  QMp  <- fit_int$QMp
  QE   <- round(fit_int$QE,  0)
  QEp  <- fit_int$QEp
  
  fit_noint <- fit_rma_mv(dd, formula = fml_noint,
                          random_spec = random_spec, use_vcov = use_vcov)
  if (is.null(fit_noint)) return(NULL)
  
  coefs    <- as.numeric(fit_noint$b)
  se_coefs <- sqrt(diag(fit_noint$vb))
  zvals    <- coefs / se_coefs
  pvals    <- 2 * pnorm(-abs(zvals))
  
  means <- data.frame(
    trait            = trait,
    level            = lvls,
    log_RR           = coefs,
    ci_lower         = coefs - 1.96 * se_coefs,
    ci_upper         = coefs + 1.96 * se_coefs,
    pval             = pvals,
    k                = as.integer(table(dd[[trait]])[lvls]),
    QM               = QM,
    QM_df            = QMdf,
    QM_pval          = QMp,
    QE               = QE,
    QE_pval          = QEp,
    pct_change       = round((exp(coefs) - 1) * 100, 1),
    pct_change_lower = round((exp(coefs - 1.96 * se_coefs) - 1) * 100, 1),
    pct_change_upper = round((exp(coefs + 1.96 * se_coefs) - 1) * 100, 1),
    stringsAsFactors = FALSE
  )
  
  contrasts_df <- NULL
  if (QMp < qm_threshold && length(lvls) > 1) {
    pairs <- utils::combn(lvls, 2, simplify = FALSE)
    contrasts_df <- dplyr::bind_rows(lapply(pairs, function(p) {
      i <- which(lvls == p[1])
      j <- which(lvls == p[2])
      cvec    <- rep(0, length(coefs))
      cvec[i] <-  1
      cvec[j] <- -1
      delta   <- as.numeric(t(cvec) %*% coefs)
      se_d    <- as.numeric(sqrt(t(cvec) %*% fit_noint$vb %*% cvec))
      zval    <- delta / se_d
      data.frame(
        trait        = trait,
        level1       = p[1],
        level2       = p[2],
        delta_log_RR = round(delta, 4),
        se           = round(se_d,  4),
        z            = round(zval,  3),
        pval         = 2 * pnorm(-abs(zval)), 
        stringsAsFactors = FALSE
      )
    }))
    contrasts_df$pval_holm <- p.adjust(contrasts_df$pval, method = "holm")
  }
  
  list(model_int   = fit_int,
       model_noint = fit_noint,
       means       = means,
       contrasts   = contrasts_df,
       QM          = QM,
       QM_df       = QMdf,
       QM_pval     = QMp,
       QE          = QE,
       QE_pval     = QEp,
       trait       = trait)
}


##============================================================================##
## Distance moderator models ----
##============================================================================##
# Meta-regression of lnRR on log(distance_disturbed_in_meters + 1).
# log(dist + 1) retains zero-distance records; only NAs excluded.
# Two model forms:
#   "linear"    : RR_log ~ log_dist1
#   "quadratic" : RR_log ~ log_dist1 + log_dist1^2
fit_distance_moderator <- function(df,
                                   stratum_col,
                                   stratum_val,
                                   random_spec,
                                   use_vcov   = TRUE,
                                   model_form = c("linear", "quadratic")) {
  
  model_form <- match.arg(model_form)
  
  dd <- dplyr::filter(df,
                      .data[[stratum_col]] == stratum_val,
                      is.finite(distance_disturbed_in_meters))
  if (nrow(dd) < 3) return(NULL)
  
  dd$log_dist1   <- log(dd$distance_disturbed_in_meters + 1)
  dd$log_dist1_2 <- dd$log_dist1^2
  
  formula_use <- if (model_form == "linear") RR_log ~ log_dist1 else
    RR_log ~ log_dist1 + log_dist1_2
  
  fit <- suppressWarnings(
    fit_rma_mv(dd, formula = formula_use,
               random_spec = random_spec, use_vcov = use_vcov)
  )
  if (is.null(fit)) return(NULL)
  
  b    <- as.numeric(fit$b)
  se   <- sqrt(diag(fit$vb))
  z    <- b / se
  p    <- 2 * pnorm(-abs(z))
  trms <- rownames(fit$b)
  trms[trms == "intrcpt"] <- "intercept"
  
  data.frame(
    level        = stratum_col,
    stratum      = stratum_val,
    model_form   = model_form,
    term         = trms,
    estimate     = round(b, 6),
    ci_lower     = round(b - 1.96 * se, 6),
    ci_upper     = round(b + 1.96 * se, 6),
    pval         = p,
    k            = nrow(dd),
    n_studies    = length(unique(dd$Study_ID)),
    QM           = round(fit$QM,  3),
    QM_df        = fit$m,
    QM_pval      = fit$QMp, 
    QE           = round(fit$QE,  0),
    QE_pval      = fit$QEp, 
    stringsAsFactors = FALSE
  )
}


# Realm-level distance plot: 
plot_distance_realm <- function(model_list,
                                df_sub,
                                out_path,
                                outcome_name,
                                model_form,
                                weighted = TRUE,
                                n_grid   = 300) {
  
  pal       <- c(terrestrial = "#1b9e77", freshwater = "#7570b3")
  cap_first <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  
  global_max_km <- max(df_sub$distance_disturbed_in_meters, na.rm = TRUE) / 1000
  dist_seq_m    <- seq(0, max(df_sub$distance_disturbed_in_meters, na.rm = TRUE),
                       length.out = n_grid)
  
  curve_list <- list()
  cross_list <- list()
  point_list <- list()
  annot_list <- list()
  
  for (realm_name in names(model_list)) {
    fit <- model_list[[realm_name]]
    if (is.null(fit)) next
    
    dd <- dplyr::filter(df_sub, realm == realm_name)
    if (nrow(dd) == 0) next
    
    ld1_seq <- log(dist_seq_m + 1)
    is_quad <- "log_dist1_2" %in% rownames(fit$b)
    newdat  <- if (is_quad) data.frame(log_dist1   = ld1_seq,
                                       log_dist1_2 = ld1_seq^2) else
                                         data.frame(log_dist1   = ld1_seq)
    
    preds <- tryCatch(predict(fit, newmods = as.matrix(newdat)),
                      error = function(e) NULL)
    if (is.null(preds)) next
    
    curve_list[[realm_name]] <- data.frame(
      dist_km  = dist_seq_m / 1000,
      pred     = as.numeric(preds$pred),
      ci_lower = as.numeric(preds$ci.lb),
      ci_upper = as.numeric(preds$ci.ub),
      realm    = realm_name
    )
    
    # Zero crossings via linear interpolation
    cr    <- curve_list[[realm_name]]
    signs <- sign(cr$pred)
    cidx  <- which(diff(signs) != 0)
    if (length(cidx) > 0) {
      xc <- sapply(cidx, function(i) {
        x1 <- cr$dist_km[i];   y1 <- cr$pred[i]
        x2 <- cr$dist_km[i+1]; y2 <- cr$pred[i+1]
        x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
      })
      cross_list[[realm_name]] <- data.frame(dist_km = xc, realm = realm_name)
    }
    
    # Points sized by weight
    wt <- if (weighted &&
              "SamplingVariance" %in% names(dd) &&
              any(is.finite(dd$SamplingVariance) & dd$SamplingVariance > 0)) {
      1 / dd$SamplingVariance
    } else {
      rep(1, nrow(dd))
    }
    wt_range  <- range(wt, na.rm = TRUE)
    wt_scaled <- if (diff(wt_range) > 0) {
      (wt - wt_range[1]) / (wt_range[2] - wt_range[1])
    } else {
      rep(0.5, length(wt))
    }
    
    point_list[[realm_name]] <- data.frame(
      dist_km = dd$distance_disturbed_in_meters / 1000,
      lnRR    = dd$RR_log,
      wt_size = 1.5 + 4 * wt_scaled,
      realm   = realm_name
    )
    
    qm_label <- if (fit$QMp < 0.001) "QM p < 0.001" else
      paste0("QM p = ", formatC(fit$QMp, format = "f", digits = 3))
    annot_list[[realm_name]] <- data.frame(
      realm = realm_name,
      label = paste0(cap_first(realm_name), ": ", qm_label),
      stringsAsFactors = FALSE
    )
  }
  
  if (length(curve_list) == 0) return(invisible(NULL))
  
  curves <- dplyr::bind_rows(curve_list)
  points <- dplyr::bind_rows(point_list)
  annots <- dplyr::bind_rows(annot_list)
  
  y_range <- range(c(curves$ci_upper, curves$ci_lower, points$lnRR),
                   na.rm = TRUE)
  y_step  <- diff(y_range) * 0.08
  annots$y_pos <- seq(from = y_range[2] - y_step * 0.3,
                      by   = -y_step,
                      length.out = nrow(annots))
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = 2,
                        colour = "grey50", linewidth = 0.8) +
    ggplot2::geom_point(data  = points,
                        ggplot2::aes(x = dist_km, y = lnRR,
                                     colour = realm, size = wt_size),
                        alpha = 0.30, show.legend = FALSE) +
    ggplot2::geom_ribbon(data = curves,
                         ggplot2::aes(x = dist_km, ymin = ci_lower,
                                      ymax = ci_upper, fill = realm),
                         alpha = 0.20, colour = NA) +
    ggplot2::geom_line(data   = curves,
                       ggplot2::aes(x = dist_km, y = pred, colour = realm),
                       linewidth = 1.3)
  
  if (length(cross_list) > 0) {
    crosses <- dplyr::bind_rows(cross_list)
    p <- p + ggplot2::geom_vline(
      data = crosses,
      ggplot2::aes(xintercept = dist_km, colour = realm),
      linetype = "dashed", linewidth = 0.9, show.legend = FALSE)
  }
  
  p <- p +
    ggplot2::scale_x_continuous(
      name   = "Distance from mine (km)",
      limits = c(0, global_max_km),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_size_identity() +
    ggplot2::scale_colour_manual(values = pal, labels = cap_first) +
    ggplot2::scale_fill_manual(  values = pal, labels = cap_first) +
    ggplot2::labs(y      = "lnRR (95% CI)",
                  colour = "Realm", fill = "Realm",
                  title  = paste0(cap_first(outcome_name),
                                  " \u2014 distance moderator (", model_form, ")")) +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.margin     = ggplot2::margin(5, 130, 5, 5, "pt"),
      plot.clip       = "off"
    )
  
  for (i in seq_len(nrow(annots))) {
    p <- p + ggplot2::annotate("text",
                               x        = global_max_km * 1.02,
                               y        = annots$y_pos[i],
                               label    = annots$label[i],
                               colour   = pal[annots$realm[i]],
                               hjust    = 0, size = 4.5, fontface = "italic")
  }
  
  ggplot2::ggsave(out_path, p, width = 220, height = 150,
                  units = "mm", dpi = 300)
  invisible(p)
}


# taxo_group_realm-level distance plot
plot_distance_taxo <- function(model_list,
                               df_sub,
                               out_dir,
                               outcome_name,
                               model_form,
                               weighted = TRUE,
                               n_grid   = 300) {
  
  cap_first   <- function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  strip_realm <- function(x) sub("^(terrestrial|freshwater)_", "", x)
  
  global_max_km <- max(df_sub$distance_disturbed_in_meters, na.rm = TRUE) / 1000
  dist_seq_m    <- seq(0, max(df_sub$distance_disturbed_in_meters, na.rm = TRUE),
                       length.out = n_grid)
  
  for (realm_label in c("terrestrial", "freshwater")) {
    
    groups <- grep(paste0("^", realm_label), names(model_list), value = TRUE)
    if (length(groups) == 0) next
    
    n_grp    <- length(groups)
    pal_vals <- setNames(
      RColorBrewer::brewer.pal(max(3, n_grp), "Set2")[seq_len(n_grp)],
      groups
    )
    pal_lbl <- setNames(pal_vals, sapply(names(pal_vals), strip_realm))
    
    curve_list <- list()
    cross_list <- list()
    point_list <- list()
    annot_list <- list()
    
    for (g in groups) {
      fit <- model_list[[g]]
      if (is.null(fit)) next
      
      dd <- dplyr::filter(df_sub, taxo_group_realm == g)
      if (nrow(dd) == 0) next
      
      lbl     <- strip_realm(g)
      ld1_seq <- log(dist_seq_m + 1)
      is_quad <- "log_dist1_2" %in% rownames(fit$b)
      newdat  <- if (is_quad) data.frame(log_dist1   = ld1_seq,
                                         log_dist1_2 = ld1_seq^2) else
                                           data.frame(log_dist1   = ld1_seq)
      
      preds <- tryCatch(predict(fit, newmods = as.matrix(newdat)),
                        error = function(e) NULL)
      if (is.null(preds)) next
      
      curve_list[[g]] <- data.frame(
        dist_km  = dist_seq_m / 1000,
        pred     = as.numeric(preds$pred),
        ci_lower = as.numeric(preds$ci.lb),
        ci_upper = as.numeric(preds$ci.ub),
        group    = lbl
      )
      
      cr    <- curve_list[[g]]
      signs <- sign(cr$pred)
      cidx  <- which(diff(signs) != 0)
      if (length(cidx) > 0) {
        xc <- sapply(cidx, function(i) {
          x1 <- cr$dist_km[i];   y1 <- cr$pred[i]
          x2 <- cr$dist_km[i+1]; y2 <- cr$pred[i+1]
          x1 + (0 - y1) * (x2 - x1) / (y2 - y1)
        })
        cross_list[[g]] <- data.frame(dist_km = xc, group = lbl)
      }
      
      wt <- if (weighted &&
                "SamplingVariance" %in% names(dd) &&
                any(is.finite(dd$SamplingVariance) & dd$SamplingVariance > 0)) {
        1 / dd$SamplingVariance
      } else {
        rep(1, nrow(dd))
      }
      wt_range  <- range(wt, na.rm = TRUE)
      wt_scaled <- if (diff(wt_range) > 0) {
        (wt - wt_range[1]) / (wt_range[2] - wt_range[1])
      } else {
        rep(0.5, length(wt))
      }
      
      point_list[[g]] <- data.frame(
        dist_km = dd$distance_disturbed_in_meters / 1000,
        lnRR    = dd$RR_log,
        wt_size = 1.5 + 4 * wt_scaled,
        group   = lbl
      )
      
      qm_label <- if (fit$QMp < 0.001) "QM p < 0.001" else
        paste0("QM p = ", formatC(fit$QMp, format = "f", digits = 3))
      annot_list[[g]] <- data.frame(
        group = lbl,
        label = paste0(cap_first(lbl), ": ", qm_label),
        stringsAsFactors = FALSE
      )
    }
    
    if (length(curve_list) == 0) next
    
    curves <- dplyr::bind_rows(curve_list)
    points <- dplyr::bind_rows(point_list)
    annots <- dplyr::bind_rows(annot_list)
    
    y_range <- range(c(curves$ci_upper, curves$ci_lower, points$lnRR),
                     na.rm = TRUE)
    y_step  <- diff(y_range) * 0.08
    annots$y_pos <- seq(from = y_range[2] - y_step * 0.3,
                        by   = -y_step,
                        length.out = nrow(annots))
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, linetype = 2,
                          colour = "grey50", linewidth = 0.8) +
      ggplot2::geom_point(data  = points,
                          ggplot2::aes(x = dist_km, y = lnRR,
                                       colour = group, size = wt_size),
                          alpha = 0.30, show.legend = FALSE) +
      ggplot2::geom_ribbon(data = curves,
                           ggplot2::aes(x = dist_km, ymin = ci_lower,
                                        ymax = ci_upper, fill = group),
                           alpha = 0.15, colour = NA) +
      ggplot2::geom_line(data   = curves,
                         ggplot2::aes(x = dist_km, y = pred, colour = group),
                         linewidth = 1.2)
    
    if (length(cross_list) > 0) {
      crosses <- dplyr::bind_rows(cross_list)
      p <- p + ggplot2::geom_vline(
        data = crosses,
        ggplot2::aes(xintercept = dist_km, colour = group),
        linetype = "dashed", linewidth = 0.85, show.legend = FALSE)
    }
    
    p <- p +
      ggplot2::scale_x_continuous(
        name   = "Distance from mine (km)",
        limits = c(0, global_max_km),
        expand = ggplot2::expansion(mult = c(0, 0))
      ) +
      ggplot2::scale_size_identity() +
      ggplot2::scale_colour_manual(values = pal_lbl, labels = cap_first) +
      ggplot2::scale_fill_manual(  values = pal_lbl, labels = cap_first) +
      ggplot2::labs(y      = "lnRR (95% CI)",
                    colour = NULL, fill = NULL,
                    title  = paste0(cap_first(outcome_name), " \u2014 ",
                                    realm_label, " (", model_form, ")")) +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.margin     = ggplot2::margin(5, 160, 5, 5, "pt")
      ) +
      ggplot2::coord_cartesian(clip = "off")
    
    
    for (i in seq_len(nrow(annots))) {
      p <- p + ggplot2::annotate("text",
                                 x        = global_max_km * 1.02,
                                 y        = annots$y_pos[i],
                                 label    = annots$label[i],
                                 colour   = pal_lbl[annots$group[i]],
                                 hjust    = 0, size = 3.8, fontface = "italic")
    }
    
    ggplot2::ggsave(
      file.path(out_dir,
                paste0(outcome_name, "_distance_", model_form,
                       "_", realm_label, "_taxo.jpg")),
      p, width = 230, height = 160, units = "mm", dpi = 300)
  }
}


##============================================================================##
## Publication bias: Egger's tests ----
##============================================================================##

egger_test_mv <- function(fit, random_spec, use_vcov = TRUE) {
  
  if (is.null(fit)) return(NULL)
  
  df_used <- attr(fit, "._df_used")
  
  # Effective sample-size predictor following Nakagawa et al. (2022)
  df_used$inv_n_tilda <- with(
    df_used,
    (N_control + N_disturbed) /
      (N_control * N_disturbed)
  )
  
  df_used$sqrt_inv_n_tilda <- sqrt(df_used$inv_n_tilda)
  
  df_used <- df_used[
    is.finite(df_used$RR_log) &
      is.finite(df_used$sqrt_inv_n_tilda) &
      df_used$N_control > 0 &
      df_used$N_disturbed > 0,
    , drop = FALSE
  ]
  
  if (nrow(df_used) < 2) return(NULL)
  
  # Unique effect-size ID: within-study/effect-size random effect
  # in the Nakagawa et al. (2022) multilevel publication-bias model
  df_used$RowID <- seq_len(nrow(df_used))
  
  random_spec_pb <- c(
    random_spec,
    list(~ 1 | RowID)
  )
  
  egger_fit <- fit_rma_mv(
    df_used,
    formula     = RR_log ~ sqrt_inv_n_tilda,
    random_spec = random_spec_pb,
    use_vcov    = use_vcov
  )
  
  if (is.null(egger_fit)) return(NULL)
  
  slope_idx <- which(
    rownames(egger_fit$b) == "sqrt_inv_n_tilda"
  )
  
  if (length(slope_idx) != 1) return(NULL)
  
  data.frame(
    term     = "sqrt_inv_n_tilda",
    estimate = round(as.numeric(egger_fit$b[slope_idx]), 4),
    ci_lower = round(as.numeric(egger_fit$ci.lb[slope_idx]), 4),
    ci_upper = round(as.numeric(egger_fit$ci.ub[slope_idx]), 4),
    pval     = as.numeric(egger_fit$pval[slope_idx]), 
    stringsAsFactors = FALSE
  )
}


## Extract funnel-plot data: meta-analytic residual + effective sample size per effect
extract_funnel_data <- function(fit) {
  
  if (is.null(fit)) return(NULL)
  
  df_used <- attr(fit, "._df_used")
  
  df_used$effective_n <- with(
    df_used,
    (4 * N_control * N_disturbed) /
      (N_control + N_disturbed)
  )
  
  keep <- is.finite(df_used$effective_n) &
    df_used$N_control > 0 &
    df_used$N_disturbed > 0
  
  data.frame(
    effective_n = df_used$effective_n[keep],
    resid       = as.numeric(residuals(fit))[keep],
    overall_est = as.numeric(fit$b[1]),
    stringsAsFactors = FALSE
  )
}


run_publication_bias <- function(model_store, out_dir, outcome_name, level,
                                 random_spec, use_vcov = TRUE) {
  
  dir_make(out_dir)
  egger_rows  <- list()
  funnel_rows <- list()
  
  for (s in names(model_store)) {
    
    fit <- model_store[[s]]
    if (is.null(fit)) next
    
    egger_row <- egger_test_mv(
      fit,
      random_spec = random_spec,
      use_vcov = use_vcov
    )
    
    if (is.null(egger_row)) {
      message("  Egger test failed for stratum: ", s)
      next
    }
    
    egger_row$stratum <- s
    egger_rows[[s]] <- egger_row
    
    funnel_row <- extract_funnel_data(fit)
    funnel_row$stratum <- s
    funnel_rows[[s]] <- funnel_row
  }
  
  if (length(egger_rows) == 0) return(NULL)
  
  egger_out <- dplyr::bind_rows(egger_rows)
  egger_out <- egger_out[
    , c("stratum", setdiff(names(egger_out), "stratum"))
  ]
  
  funnel_out <- dplyr::bind_rows(funnel_rows)
  funnel_out <- funnel_out[
    , c("stratum", setdiff(names(funnel_out), "stratum"))
  ]
  
  readr::write_csv(
    egger_out,
    file.path(
      out_dir,
      paste0(outcome_name, "_", level, "_egger_summary.csv")
    )
  )
  
  readr::write_csv(
    funnel_out,
    file.path(
      out_dir,
      paste0(outcome_name, "_", level, "_funnel_data.csv")
    )
  )
  
  message("  Saved Egger summary + funnel data: ", out_dir)
  
  list(
    egger = egger_out,
    funnel = funnel_out
  )
}