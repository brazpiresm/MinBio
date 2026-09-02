##============================================================================##
## 04_abundance_RE_comparison.R
##============================================================================##

source("00_functions.R")


##============================================================================##
## Settings ----
##============================================================================##

FILES <- list(
  abundance = "data/Abundance.xlsx",
  mines     = "data/Mines.xlsx"
)

METHOD_SD <- "Bracken"
USE_VCOV  <- TRUE   

RE_CANDIDATES <- list(
  "Study_Mine" = list(
    ~ 1 | Study_ID / Mine_ID
  ),
  "Study_Mine_Genus" = list(
    ~ 1 | Study_ID / Mine_ID,
    ~ 1 | genus
  ),
  "Study_Mine_Binomial" = list(
    ~ 1 | Study_ID / Mine_ID,
    ~ 1 | binomial_name
  )
)

OUTROOT <- file.path("results", "abundance", "RE_comparison")
dir_make(OUTROOT)


##============================================================================##
## Prepare and impute data ----
##============================================================================##

t0 <- prep_data(FILES$abundance, outcome = "abundance", mines_path = FILES$mines)
t  <- imp_by_realm(t0, method_SD = METHOD_SD)
t  <- set_sampling_variance(t, weighted = TRUE)


##============================================================================##
## BIC comparison, fitted separately per realm ----
##============================================================================##

all_results <- list()

for (r in c("freshwater", "terrestrial")) {
  
  message("\n--- RE comparison for realm: ", r, " ---")
  dd <- dplyr::filter(t, realm == r)
  
  for (re_name in names(RE_CANDIDATES)) {
    
    fit <- fit_rma_mv(dd,
                      formula     = RR_log ~ 1,
                      random_spec = RE_CANDIDATES[[re_name]],
                      use_vcov    = USE_VCOV)
    
    if (is.null(fit)) {
      message("  ", re_name, ": model failed")
      all_results[[paste(r, re_name)]] <- data.frame(
        realm = r, re_structure = re_name,
        AICc = NA_real_, AIC = NA_real_, BIC = NA_real_, logLik = NA_real_,
        log_RR = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    bic_val    <- BIC(fit)
    aic_val    <- AIC(fit)
    aicc_val   <- AIC(fit, correct = TRUE)
    loglik_val <- as.numeric(logLik(fit))
    
    message("  ", re_name, ": AICc = ", round(aicc_val, 1),
            "  AIC = ", round(aic_val, 1),
            "  BIC = ", round(bic_val, 1),
            "  lnRR = ", round(as.numeric(fit$b[1]), 3),
            "  [", round(fit$ci.lb, 3), ", ", round(fit$ci.ub, 3), "]")
    
    all_results[[paste(r, re_name)]] <- data.frame(
      realm        = r,
      re_structure = re_name,
      AICc         = aicc_val,
      AIC          = aic_val,
      BIC          = bic_val,
      logLik       = loglik_val,
      log_RR       = as.numeric(fit$b[1]),
      ci_lower     = fit$ci.lb,
      ci_upper     = fit$ci.ub,
      stringsAsFactors = FALSE
    )
  }
}

comparison_table <- dplyr::bind_rows(all_results) %>%
  dplyr::group_by(realm) %>%
  dplyr::arrange(AICc, .by_group = TRUE) %>%
  dplyr::mutate(
    delta_AICc = AICc - min(AICc, na.rm = TRUE),
    delta_AIC  = AIC  - min(AIC,  na.rm = TRUE),
    delta_BIC  = BIC  - min(BIC,  na.rm = TRUE)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(realm, re_structure, AICc, delta_AICc, AIC, delta_AIC,
                BIC, delta_BIC, logLik, log_RR, ci_lower, ci_upper)

message("\n--- Full comparison table ---")
print(comparison_table)

readr::write_csv(comparison_table,
                 file.path(OUTROOT, "abundance_RE_comparison_BIC.csv"))
message("\nSaved: ", file.path(OUTROOT, "abundance_RE_comparison_BIC.csv"))
##============================================================================##