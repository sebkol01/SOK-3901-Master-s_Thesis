# =============================================================================
# THESIS MODELLING PIPELINE — Norwegian inflation: Walker TVP + XGBoost SHAP
# =============================================================================
# Authors:  Seb & Amund
# Replaces: Walker_estimation.r, Walker_robustness_USD.R,
#           walker_harmonize_from_outputs.R, Walker_harmonize_robustness.R.
#
# Inputs assumed already on disk:
#   master_data.csv                                       (built upstream)
#   xgb_v8_harmonized_shap_group_signed.csv               (xgb_pipeline.py main)
#   xgb_v8_harmonized_shap_group_abs.csv
#   xgb_harmonized_robustness_shap_group_signed.csv       (xgb_pipeline.py USD)
#   xgb_harmonized_robustness_shap_group_abs.csv
#
# Stages (toggle via the RUN_* flags in section 0):
#   1. Walker backtest for any subset of {main, no_intercept, usd_external}
#   2. Harmonise Walker contributions per spec into block-level objects/shares
#   3. Harmonise XGB SHAP outputs (main + USD robustness) into the same schema
#   4. SHAP <-> Walker comparison: joint table, residual decomposition Δ_kt
#   5. Regime-cutoff sensitivity (±3 months)
#   6. AR(1) and RW baselines + accuracy table
#
# Design invariants preserved verbatim from the original code:
#   - real-time (training-window) standardisation at each forecast origin
#   - filtered (last-time-point) coefficients, NOT smoothed
#   - c_jt = β_jt^filtered · x_jt^scaled; intercept contribution = filtered β_0t
#   - lag-1 information set, h=3, TEST_START = 2020-01-01
#   - 4 chains × 2000 iter (1000 warmup), adapt_delta=0.95, max_treedepth=12
# =============================================================================


# -----------------------------------------------------------------------------
# 0. CONFIGURATION
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(walker)
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(rlang)
  library(purrr)
  library(ggplot2)
  library(lubridate)
  library(knitr)
  library(kableExtra)
})

# --- Working directory: try Amund's path, then Seb's; whichever exists wins.
.set_wd <- function() {
  paths <- c(
    "/Users/amundbech1/Library/CloudStorage/OneDrive-UiTOffice365/Data_master/MasterOppgave",
    "~/Library/CloudStorage/OneDrive-UiTOffice365/Amund Skjalgsønn Bech's files - Data_master/MasterOppgave"
  )
  for (p in paths) if (dir.exists(p)) { setwd(p); return(invisible(p)) }
  warning("None of the predefined working directories exist; using cwd: ", getwd())
}
.set_wd()

# --- Run flags (flip to FALSE to skip a stage) -------------------------------
# Walker backtests are expensive: ~1.5 h per spec. Skip stages you've already
# run by setting their flag to FALSE; outputs from earlier runs are read from
# disk by the downstream stages.
RUN_WALKER         <- list(main = FALSE, usd_external = TRUE)
RUN_HARMONIZE_W    <- TRUE      # build block_signed / shares per Walker spec
RUN_HARMONIZE_XGB  <- TRUE      # build block_signed / shares per XGB SHAP file
RUN_COMPARE        <- TRUE      # SHAP vs Walker: joint table, residual, plots
RUN_REGIME_ROB     <- TRUE      # ±3m regime cutoff sensitivity
RUN_BASELINES      <- TRUE      # AR(1) + RW accuracy table

# --- Estimation hyperparameters (UNCHANGED from the original Walker scripts) -
H            <- 3                         # forecast horizon (months)
TEST_START   <- as.Date("2020-01-01")
# Minimum training observations required before a forecast origin enters the
# backtest. Does not bind under the current TEST_START but guards any
# robustness check that moves the start earlier or shrinks the sample.
MIN_TRAIN_N  <- 50
N_CHAINS     <- 4
N_ITER       <- 2000
N_WARMUP     <- 1000
N_CORES      <- min(4, parallel::detectCores())
SEED         <- 42

set.seed(SEED)
options(mc.cores = N_CORES)
rstan_options(auto_write = TRUE)

# --- Feature set + economic block map (single source of truth) ---------------
FEATURE_COLS <- c(
  "kpi_yoy_lag1", "ppi_yoy_lag1", "oil_yoy_lag1",
  "usd_nok_lag1", "eur_nok_lag1",
  "import_yoy_lag1", "eksport_yoy_lag1",
  "unemp_lag1", "rente_lag1"
)

BLOCK_MAP <- tibble::tribble(
  ~feature,            ~block,
  "kpi_yoy_lag1",     "AR",
  "ppi_yoy_lag1",     "PPI",
  "oil_yoy_lag1",     "Oil",
  "usd_nok_lag1",     "FX",
  "eur_nok_lag1",     "FX",
  "import_yoy_lag1",  "Trade",
  "eksport_yoy_lag1", "Trade",
  "unemp_lag1",       "Labour",
  "rente_lag1",       "Monetary"
)
BLOCK_NAMES <- c("AR", "PPI", "Monetary", "FX", "Oil", "Trade", "Labour")

# --- Regime definitions (consistent with thesis Tabell 5) --------------------
# NB: walker_harmonize_from_outputs.R used a different scheme
# (Pre-shock/Pandemic/Inflation shock/Disinflation/Recent). That scheme is NOT
# the thesis design — discarded here.
assign_regime <- function(Date, shift_months = 0) {
  cut1 <- as.Date("2021-06-01") %m+% months(shift_months)
  cut2 <- as.Date("2023-01-01") %m+% months(shift_months)
  cut3 <- as.Date("2024-06-01") %m+% months(shift_months)
  case_when(
    Date < cut1 ~ "COVID (2020-2021)",
    Date < cut2 ~ "Energy Crisis",
    Date < cut3 ~ "Disinflation",
    TRUE        ~ "Normalization"
  )
}

# --- Walker spec catalogue ---------------------------------------------------
WALKER_SPECS <- list(
  main = list(
    has_intercept = TRUE,
    cols          = list(oil = "oil_price_nok", import = "import", eksport = "eksport"),
    out_prefix    = "walker_"
  ),
  usd_external = list(
    has_intercept = TRUE,
    cols          = list(oil = "oljepris_USD", import = "import_USD", eksport = "eksport_USD"),
    out_prefix    = "walker_usd_external_"
  )
)


# -----------------------------------------------------------------------------
# 1. SHARED HELPERS
# -----------------------------------------------------------------------------

#' Standardise `cols` of `train_df` and `new_df` using means/sds from `train_df`.
#' Vectorised over columns (replaces the original per-column for-loop).
scale_with_train_stats <- function(train_df, new_df, cols) {
  means <- vapply(train_df[cols], mean, numeric(1), na.rm = TRUE)
  sds   <- vapply(train_df[cols], sd,   numeric(1), na.rm = TRUE)
  bad   <- is.na(sds) | sds == 0
  if (any(bad)) {
    stop(sprintf("Zero/NA SD in training data for: %s",
                 paste(names(sds)[bad], collapse = ", ")))
  }
  apply_scale <- function(x_df) {
    x_df[cols] <- mapply(function(x, m, s) (x - m) / s,
                         x_df[cols], means, sds, SIMPLIFY = FALSE)
    x_df
  }
  list(train_scaled = apply_scale(train_df),
       new_scaled   = apply_scale(new_df),
       means        = means, sds = sds)
}

#' Diagnostics for one Walker fit: divergent transitions, max R-hat, min ESS.
extract_fit_diagnostics <- function(fit) {
  sp    <- get_sampler_params(fit$stanfit, inc_warmup = FALSE)
  n_div <- sum(vapply(sp, function(x) sum(x[, "divergent__"]), numeric(1)))
  summ  <- summary(fit$stanfit)$summary
  list(n_div    = n_div,
       max_rhat = max(summ[, "Rhat"], na.rm = TRUE),
       min_ess  = min(summ[, "n_eff"], na.rm = TRUE))
}

#' Real-time prediction at the forecast origin from filtered (last-time)
#' posterior-mean coefficients:  ŷ = β_0(t) + Σ_j β_j(t)·x_j(t).
#' has_intercept = FALSE handles the no_intercept robustness spec.
predict_last_state <- function(fit, new_x, feature_cols, has_intercept = TRUE) {
  coefs_last <- coef(fit) %>%
    ungroup() %>%
    filter(time == max(time)) %>%
    select(beta, mean)
  
  intercept_val <- if (has_intercept) {
    if (!"(Intercept)" %in% coefs_last$beta) {
      stop("has_intercept=TRUE but '(Intercept)' missing from coef(fit). ",
           "Check that the formula uses 'rw1(~ 1 + ...)'.")
    }
    coefs_last$mean[coefs_last$beta == "(Intercept)"]
  } else 0
  
  beta_tbl <- coefs_last %>% filter(beta %in% feature_cols)
  missing  <- setdiff(feature_cols, beta_tbl$beta)
  if (length(missing) > 0) {
    stop(sprintf("Missing coefficients for: %s", paste(missing, collapse = ", ")))
  }
  beta_vals <- setNames(beta_tbl$mean, beta_tbl$beta)
  x_vals    <- as.numeric(new_x[1, names(beta_vals), drop = TRUE])
  as.numeric(intercept_val + sum(beta_vals * x_vals))
}

#' Filtered (last-time) posterior summary at a single forecast origin.
extract_filtered_coefs <- function(fit, origin_date) {
  coef(fit) %>%
    ungroup() %>%
    filter(time == max(time)) %>%
    select(beta, mean, sd, `2.5%`, `50%`, `97.5%`) %>%
    mutate(origin = origin_date)
}

#' Build the full feature matrix for one Walker spec (lag1 information set).
#' `cols` lets us swap in USD-denominated oil/trade columns for the robustness
#' spec without touching the rest of the pipeline.
build_features <- function(df, cols, h = H) {
  df %>%
    mutate(
      kpi_yoy_raw     = (kpi / lag(kpi, 12) - 1) * 100,
      ppi_yoy_raw     = (ppi / lag(ppi, 12) - 1) * 100,
      oil_yoy_raw     = (.data[[cols$oil]]     / lag(.data[[cols$oil]],     12) - 1) * 100,
      import_yoy_raw  = (.data[[cols$import]]  / lag(.data[[cols$import]],  12) - 1) * 100,
      eksport_yoy_raw = (.data[[cols$eksport]] / lag(.data[[cols$eksport]], 12) - 1) * 100,
      
      kpi_yoy_lag1     = lag(kpi_yoy_raw, 1),
      ppi_yoy_lag1     = lag(ppi_yoy_raw, 1),
      oil_yoy_lag1     = lag(oil_yoy_raw, 1),
      usd_nok_lag1     = lag(usd_nok, 1),
      eur_nok_lag1     = lag(eur_nok, 1),
      import_yoy_lag1  = lag(import_yoy_raw, 1),
      eksport_yoy_lag1 = lag(eksport_yoy_raw, 1),
      unemp_lag1       = lag(unemployment, 1),
      rente_lag1       = lag(rente, 1),
      
      target = lead(kpi_yoy_raw, h)
    ) %>%
    select(Date, target, kpi_yoy_raw, all_of(FEATURE_COLS)) %>%
    filter(complete.cases(.))
}


# 2. WALKER ESTIMATION ENGINE
# -----------------------------------------------------------------------------

#' rw1 formula: time-varying intercept controlled by the spec.
make_walker_formula <- function(has_intercept) {
  rhs_intercept <- if (has_intercept) "1 + " else "0 + "
  as.formula(
    paste0("target ~ -1 + rw1(~ ", rhs_intercept,
           paste(FEATURE_COLS, collapse = " + "),
           ", beta = c(0, 10), sigma = c(2, 0.01))")
  )
}

#' One forecast origin: scale, fit, predict, store filtered coefs + scaled x.
#' Returns NULL if the origin is unusable (NA target or insufficient training).
fit_one_origin <- function(origin, model_df, formula, has_intercept) {
  train_raw      <- model_df %>% filter(Date <  origin)
  origin_row_raw <- model_df %>% filter(Date == origin)
  if (nrow(origin_row_raw) != 1) return(NULL)
  y_actual <- origin_row_raw$target
  y_rw     <- origin_row_raw$kpi_yoy_lag1
  if (is.na(y_actual) || is.na(y_rw) || nrow(train_raw) < MIN_TRAIN_N) return(NULL)
  
  scaled       <- scale_with_train_stats(train_raw, origin_row_raw, FEATURE_COLS)
  train_scaled <- scaled$train_scaled
  new_x_scaled <- scaled$new_scaled[, FEATURE_COLS, drop = FALSE]
  
  t_start <- Sys.time()
  fit <- walker(
    formula       = formula,
    data          = train_scaled,
    sigma_y_prior = c(2, 0.01),
    chains        = N_CHAINS, iter = N_ITER, warmup = N_WARMUP,
    cores         = N_CORES,  refresh = 0,
    control       = list(adapt_delta = 0.95, max_treedepth = 12)
  )
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
  
  diag <- extract_fit_diagnostics(fit)
  yhat <- predict_last_state(fit, new_x_scaled, FEATURE_COLS, has_intercept)
  
  filt_coefs <- extract_filtered_coefs(fit, origin)
  x_long     <- as.data.frame(new_x_scaled) %>%
    pivot_longer(everything(), names_to = "beta", values_to = "x_scaled") %>%
    mutate(origin = origin)
  
  list(
    forecast_row = data.frame(
      Date = origin, y_actual = y_actual, y_hat = yhat, y_rw = y_rw,
      n_div = diag$n_div, max_rhat = diag$max_rhat, min_ess = diag$min_ess,
      elapsed = elapsed, stringsAsFactors = FALSE
    ),
    filtered_coefs = filt_coefs,
    x_values       = x_long
  )
}

#' Run the expanding-window backtest for one spec; write all four canonical
#' outputs to disk and return an eval-summary tibble.
#' The outer loop is sequential by design (each iteration is an MCMC fit).
run_walker_spec <- function(spec_name, df) {
  spec       <- WALKER_SPECS[[spec_name]]
  model_df   <- build_features(df, spec$cols, h = H)
  origins    <- model_df %>% filter(Date >= TEST_START, !is.na(target)) %>% pull(Date)
  formula    <- make_walker_formula(spec$has_intercept)
  out_prefix <- spec$out_prefix
  log_file   <- paste0(out_prefix, "diagnostics.txt")
  
  cat(sprintf("\n=== WALKER [%s]: %d origins (~%.1f h at 1.5 min/fit) ===\n",
              spec_name, length(origins), length(origins) * 1.5 / 60))
  cat(sprintf("Started %s | spec=%s | H=%d\n",
              Sys.time(), spec_name, H), file = log_file)
  
  results_list <- vector("list", length(origins))
  for (i in seq_along(origins)) {
    origin <- origins[i]
    cat(sprintf("[%d/%d] %s ... ", i, length(origins), origin))
    res <- tryCatch(
      fit_one_origin(origin, model_df, formula, spec$has_intercept),
      error = function(e) {
        msg <- sprintf("[%s] FAILED: %s\n", origin, e$message)
        cat(msg); cat(msg, file = log_file, append = TRUE); NULL
      }
    )
    if (!is.null(res)) {
      results_list[[i]] <- res
      cat(sprintf("yhat=%.2f (%.1f min, div=%d, Rhat=%.3f)\n",
                  res$forecast_row$y_hat, res$forecast_row$elapsed,
                  res$forecast_row$n_div, res$forecast_row$max_rhat))
      cat(sprintf("[%s] yhat=%.2f div=%d Rhat=%.3f ESS=%.0f elapsed=%.1fmin\n",
                  origin, res$forecast_row$y_hat, res$forecast_row$n_div,
                  res$forecast_row$max_rhat, res$forecast_row$min_ess,
                  res$forecast_row$elapsed),
          file = log_file, append = TRUE)
      # Persist forecasts every 6 origins so a crash doesn't lose progress.
      if (i %% 6 == 0) {
        partial <- bind_rows(map(results_list, "forecast_row"))
        write_csv(partial, paste0(out_prefix, "forecasts.csv"))
      }
    }
  }
  
  results_list <- compact(results_list)
  if (length(results_list) == 0) {
    warning(sprintf("Spec %s produced no results.", spec_name))
    return(invisible(NULL))
  }
  
  forecasts      <- bind_rows(map(results_list, "forecast_row"))
  filtered_coefs <- bind_rows(map(results_list, "filtered_coefs"))
  x_values       <- bind_rows(map(results_list, "x_values"))
  
  # Real-time contributions: c_jt = β_jt^filtered · x_jt^scaled.
  # Intercept contribution = filtered β_0(t) directly.
  contributions <- filtered_coefs %>%
    left_join(x_values, by = c("origin", "beta")) %>%
    mutate(contribution = if_else(beta == "(Intercept)", mean, mean * x_scaled)) %>%
    select(origin, beta, contribution) %>%
    pivot_wider(names_from = beta, values_from = contribution) %>%
    rename(Date = origin)
  
  if ("(Intercept)" %in% names(contributions)) {
    contributions <- contributions %>% rename(intercept = `(Intercept)`)
  } else {
    contributions$intercept <- 0   # no_intercept spec
  }
  
  write_csv(forecasts,      paste0(out_prefix, "forecasts.csv"))
  write_csv(filtered_coefs, paste0(out_prefix, "filtered_coefs.csv"))
  write_csv(x_values,       paste0(out_prefix, "x_values_at_origin.csv"))
  write_csv(contributions,  paste0(out_prefix, "contributions.csv"))
  
  rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
  summary_row <- tibble(
    spec        = spec_name,
    n           = nrow(forecasts),
    rmse_walker = rmse(forecasts$y_hat, forecasts$y_actual),
    rmse_rw     = rmse(forecasts$y_rw,  forecasts$y_actual),
    bias_walker = mean(forecasts$y_hat - forecasts$y_actual, na.rm = TRUE),
    total_div   = sum(forecasts$n_div, na.rm = TRUE),
    max_rhat    = max(forecasts$max_rhat, na.rm = TRUE),
    min_ess     = min(forecasts$min_ess, na.rm = TRUE)
  )
  write_csv(summary_row, paste0(out_prefix, "eval_summary.csv"))
  cat("Spec eval:\n"); print(summary_row)
  invisible(summary_row)
}


# -----------------------------------------------------------------------------
# 3. HARMONIZATION (Walker contributions OR XGB SHAP outputs)
# -----------------------------------------------------------------------------

#' Compute signed and absolute shares per block from a wide block_signed table.
#' Loop kept (7 blocks): transparent, fast at this scale, and any vectorisation
#' would obscure the formula c_kt / Σ_m |c_mt|.
compute_block_shares <- function(block_signed) {
  abs_mat   <- abs(as.matrix(block_signed[, BLOCK_NAMES]))
  total_abs <- rowSums(abs_mat, na.rm = TRUE)
  shares    <- block_signed
  shares$total_abs <- total_abs
  for (b in BLOCK_NAMES) {
    shares[[paste0("S_", b)]] <- ifelse(total_abs > 0,
                                        block_signed[[b]]      / total_abs, NA_real_)
    shares[[paste0("A_", b)]] <- ifelse(total_abs > 0,
                                        abs(block_signed[[b]]) / total_abs, NA_real_)
  }
  shares
}

#' Aggregate Walker contributions to block-level signed/abs/share tables.
harmonize_walker <- function(contributions_path, out_prefix,
                             regime_fn = assign_regime) {
  walker_contrib <- read_csv(contributions_path, show_col_types = FALSE) %>%
    mutate(Date = as.Date(Date))
  if (!"intercept" %in% names(walker_contrib)) walker_contrib$intercept <- 0
  
  missing_features <- setdiff(BLOCK_MAP$feature, names(walker_contrib))
  if (length(missing_features) > 0) {
    stop("Missing feature columns in ", contributions_path, ": ",
         paste(missing_features, collapse = ", "))
  }
  
  feature_long <- walker_contrib %>%
    select(Date, all_of(BLOCK_MAP$feature), intercept) %>%
    pivot_longer(cols = all_of(BLOCK_MAP$feature),
                 names_to = "feature", values_to = "contribution") %>%
    left_join(BLOCK_MAP, by = "feature")
  
  block_signed <- feature_long %>%
    group_by(Date, block) %>%
    summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = block, values_from = contribution) %>%
    left_join(walker_contrib %>% select(Date, intercept), by = "Date") %>%
    rename(baseline_WalkerAdj = intercept) %>%
    mutate(regime = regime_fn(Date)) %>%
    arrange(Date)
  
  block_abs <- block_signed %>%
    mutate(across(all_of(BLOCK_NAMES), abs)) %>%
    select(Date, all_of(BLOCK_NAMES), regime)
  
  shares <- compute_block_shares(block_signed)
  
  # Reconstruction check: baseline + Σ block_signed should ≈ y_hat.
  reconstruction <- block_signed %>%
    mutate(sum_blocks        = rowSums(across(all_of(BLOCK_NAMES)), na.rm = TRUE),
           fitted_from_parts = baseline_WalkerAdj + sum_blocks)
  forecast_path <- sub("contributions\\.csv$", "forecasts.csv", contributions_path)
  if (file.exists(forecast_path)) {
    fc <- read_csv(forecast_path, show_col_types = FALSE) %>%
      mutate(Date = as.Date(Date))
    reconstruction <- reconstruction %>%
      left_join(fc %>% select(Date, y_hat), by = "Date") %>%
      mutate(reconstruction_gap_vs_backtest = y_hat - fitted_from_parts)
  }
  
  share_cols <- grep("^(S_|A_)", names(shares), value = TRUE)
  regime_summary <- shares %>%
    group_by(regime) %>%
    summarise(across(all_of(share_cols), \(x) mean(x, na.rm = TRUE)),
              n_obs = n(), .groups = "drop")
  global_summary <- shares %>%
    summarise(across(all_of(share_cols), \(x) mean(x, na.rm = TRUE)))
  
  write_csv(feature_long,    paste0(out_prefix, "harmonized_feature_long.csv"))
  write_csv(block_signed,    paste0(out_prefix, "harmonized_block_signed.csv"))
  write_csv(block_abs,       paste0(out_prefix, "harmonized_block_abs.csv"))
  write_csv(shares,          paste0(out_prefix, "harmonized_shares.csv"))
  write_csv(regime_summary,  paste0(out_prefix, "harmonized_regime_summary.csv"))
  write_csv(global_summary,  paste0(out_prefix, "harmonized_global_summary.csv"))
  write_csv(reconstruction,  paste0(out_prefix, "harmonized_reconstruction_check.csv"))
  write_csv(BLOCK_MAP,       paste0(out_prefix, "feature_block_map.csv"))
  
  invisible(list(block_signed = block_signed, block_abs = block_abs,
                 shares = shares, regime_summary = regime_summary))
}

#' Aggregate XGB SHAP block-level outputs to the same schema as Walker.
#' XGB's blocks are already aggregated in the Python pipeline; we just rename
#' to the short Walker block names and add regime + shares.
harmonize_xgb <- function(signed_path, out_prefix, regime_fn = assign_regime) {
  shap_signed <- read_csv(signed_path, show_col_types = FALSE)
  date_col    <- intersect(c("Date", "date"), names(shap_signed))[[1]]
  shap_signed <- shap_signed %>% mutate(Date = as.Date(.data[[date_col]]))
  
  rename_map <- c(
    "AR (inflation)"  = "AR",
    "PPI / cost-push" = "PPI",
    "Monetary policy" = "Monetary",
    "FX"              = "FX",
    "Oil"             = "Oil",
    "Trade"           = "Trade",
    "Labour market"   = "Labour"
  )
  
  block_signed <- shap_signed %>%
    select(Date, all_of(names(rename_map)), shap_base) %>%
    rename(!!!setNames(names(rename_map), rename_map)) %>%
    rename(baseline_SHAP = shap_base) %>%
    mutate(regime = regime_fn(Date)) %>%
    arrange(Date)
  
  block_abs <- block_signed %>%
    mutate(across(all_of(BLOCK_NAMES), abs)) %>%
    select(Date, all_of(BLOCK_NAMES), regime)
  
  shares <- compute_block_shares(block_signed)
  
  share_cols <- grep("^(S_|A_)", names(shares), value = TRUE)
  regime_summary <- shares %>%
    group_by(regime) %>%
    summarise(across(all_of(share_cols), \(x) mean(x, na.rm = TRUE)),
              n_obs = n(), .groups = "drop")
  global_summary <- shares %>%
    summarise(across(all_of(share_cols), \(x) mean(x, na.rm = TRUE)))
  
  write_csv(block_signed,   paste0(out_prefix, "harmonized_block_signed.csv"))
  write_csv(block_abs,      paste0(out_prefix, "harmonized_block_abs.csv"))
  write_csv(shares,         paste0(out_prefix, "harmonized_shares.csv"))
  write_csv(regime_summary, paste0(out_prefix, "harmonized_regime_summary.csv"))
  write_csv(global_summary, paste0(out_prefix, "harmonized_global_summary.csv"))
  
  invisible(list(block_signed = block_signed, block_abs = block_abs,
                 shares = shares, regime_summary = regime_summary))
}


# -----------------------------------------------------------------------------
# 4. SHAP <-> WALKER COMPARISON  +  RESIDUAL DECOMPOSITION
# -----------------------------------------------------------------------------

#' Inner join of XGB-SHAP and Walker block-signed tables, plus regime label.
build_joint_attributions <- function(walker_signed_path, xgb_signed_path) {
  ws <- read_csv(walker_signed_path, show_col_types = FALSE) %>%
    mutate(Date = as.Date(Date))
  xs_raw <- read_csv(xgb_signed_path, show_col_types = FALSE)
  date_col <- intersect(c("Date", "date"), names(xs_raw))[[1]]
  xs <- xs_raw %>% mutate(Date = as.Date(.data[[date_col]]))
  
  xgb_long  <- c("AR (inflation)", "PPI / cost-push", "Monetary policy",
                 "FX", "Oil", "Trade", "Labour market")
  xgb_short <- c("AR", "PPI", "Monetary", "FX", "Oil", "Trade", "Labour")
  
  xgb_std <- xs %>%
    transmute(
      Date,
      actual         = .data$actual,
      predicted_raw  = .data$predicted_raw,
      baseline_SHAP  = .data$shap_base,
      !!!setNames(syms(xgb_long), paste0(xgb_short, "_SHAP"))
    )
  walker_std <- ws %>%
    transmute(
      Date,
      baseline_WalkerAdj = baseline_WalkerAdj,
      !!!setNames(syms(BLOCK_NAMES), paste0(BLOCK_NAMES, "_WalkerAdj"))
    )
  
  xgb_std %>%
    left_join(walker_std, by = "Date") %>%
    mutate(regime = assign_regime(Date)) %>%
    arrange(Date)
}

#' Block-level residuals Δ_kt = C_kt^XGB - C_kt^TVP, plus regime-conditional
#' variance and mean-absolute-magnitude summaries.
residual_decomposition <- function(joint) {
  res_long <- joint %>%
    transmute(Date, regime,
              AR       = AR_SHAP       - AR_WalkerAdj,
              PPI      = PPI_SHAP      - PPI_WalkerAdj,
              Monetary = Monetary_SHAP - Monetary_WalkerAdj,
              FX       = FX_SHAP       - FX_WalkerAdj,
              Oil      = Oil_SHAP      - Oil_WalkerAdj,
              Trade    = Trade_SHAP    - Trade_WalkerAdj,
              Labour   = Labour_SHAP   - Labour_WalkerAdj) %>%
    pivot_longer(-c(Date, regime), names_to = "block", values_to = "Delta") %>%
    mutate(block = factor(block, levels = BLOCK_NAMES))
  
  res_metrics <- res_long %>%
    group_by(block, regime) %>%
    summarise(var_delta = var(Delta, na.rm = TRUE),
              mean_abs  = mean(abs(Delta), na.rm = TRUE),
              n_obs     = sum(!is.na(Delta)),
              .groups   = "drop")
  list(long = res_long, metrics = res_metrics)
}

#' Sign-agreement and Pearson correlations (signed and absolute) per block.
#' Vectorised over BLOCK_NAMES (replaces 14 pasted cor() calls in the original).
sign_and_corr <- function(joint) {
  bind_rows(lapply(BLOCK_NAMES, function(b) {
    s <- joint[[paste0(b, "_SHAP")]]
    w <- joint[[paste0(b, "_WalkerAdj")]]
    tibble(
      block       = b,
      pct_agree   = mean(sign(s) == sign(w), na.rm = TRUE),
      corr_signed = cor(s,        w,        use = "complete.obs"),
      corr_abs    = cor(abs(s),   abs(w),   use = "complete.obs")
    )
  }))
}


# -----------------------------------------------------------------------------
# 5. REGIME-CUTOFF SENSITIVITY (±3 MONTHS)
# -----------------------------------------------------------------------------

#' Recompute mean absolute shares per regime under a partition shift.
#' Vectorised: replaces the rowwise()+14-column-paste pattern in the original.
regime_cutoff_robustness <- function(joint, shifts = c(-3, 0, 3)) {
  shap_cols   <- paste0(BLOCK_NAMES, "_SHAP")
  walker_cols <- paste0(BLOCK_NAMES, "_WalkerAdj")
  
  one_shift <- function(s) {
    df       <- joint %>% mutate(regime_alt = assign_regime(Date, s))
    abs_shap <- abs(as.matrix(df[, shap_cols]))
    abs_walk <- abs(as.matrix(df[, walker_cols]))
    tot_shap <- rowSums(abs_shap, na.rm = TRUE)
    tot_walk <- rowSums(abs_walk, na.rm = TRUE)
    
    shares <- df %>% select(Date, regime_alt)
    for (b in BLOCK_NAMES) {
      shares[[paste0("A_", b, "_SHAP")]]   <-
        ifelse(tot_shap > 0, abs(df[[paste0(b, "_SHAP")]])      / tot_shap, NA_real_)
      shares[[paste0("A_", b, "_Walker")]] <-
        ifelse(tot_walk > 0, abs(df[[paste0(b, "_WalkerAdj")]]) / tot_walk, NA_real_)
    }
    
    shares %>%
      group_by(regime = regime_alt) %>%
      summarise(across(starts_with("A_"), \(x) mean(x, na.rm = TRUE) * 100),
                n_obs = n(), .groups = "drop") %>%
      mutate(shift_months = s,
             partition    = if (s == 0) "Original" else sprintf("Shift %+dm", s))
  }
  
  bind_rows(lapply(shifts, one_shift)) %>%
    pivot_longer(starts_with("A_"),
                 names_to      = c("block", "model"),
                 names_pattern = "^A_([^_]+)_(SHAP|Walker)$",
                 values_to     = "share") %>%
    mutate(block     = factor(block, levels = BLOCK_NAMES),
           partition = factor(partition,
                              levels = c("Shift -3m", "Original", "Shift +3m")))
}


# -----------------------------------------------------------------------------
# 6. AR(1) AND RW BASELINES + ACCURACY METRICS
# -----------------------------------------------------------------------------
safe_mape <- function(actual, predicted) {
  ok <- is.finite(actual) & is.finite(predicted) & actual != 0
  if (!any(ok)) return(NA_real_)
  mean(abs((actual[ok] - predicted[ok]) / actual[ok]), na.rm = TRUE) * 100
}

forecast_metrics <- function(actual, predicted, model_name) {
  err <- predicted - actual
  tibble(
    Model = model_name,
    RMSE  = sqrt(mean(err^2, na.rm = TRUE)),
    MAE   = mean(abs(err), na.rm = TRUE),
    MAPE  = safe_mape(actual, predicted),
    Bias  = mean(err, na.rm = TRUE),
    N     = sum(is.finite(actual) & is.finite(predicted))
  )
}

#' AR(1) expanding-window backtest. Fixes a bug in the original
#' (referenced a `kpi_yoy` column that doesn't exist in master_data.csv).
#' Loop kept: each origin needs its own lm refit; cheap, sequential, readable.
run_ar1_backtest <- function(df) {
  model_df <- df %>%
    mutate(
      kpi_yoy_raw  = (kpi / lag(kpi, 12) - 1) * 100,
      kpi_yoy_lag1 = lag(kpi_yoy_raw, 1),
      target       = lead(kpi_yoy_raw, H)
    ) %>%
    filter(!is.na(target), !is.na(kpi_yoy_lag1))
  
  origins <- model_df %>% filter(Date >= TEST_START) %>% pull(Date)
  out <- vector("list", length(origins))
  for (i in seq_along(origins)) {
    origin     <- origins[i]
    train_raw  <- model_df %>% filter(Date <  origin)
    origin_row <- model_df %>% filter(Date == origin)
    if (nrow(train_raw) < MIN_TRAIN_N || nrow(origin_row) != 1) next
    fit <- lm(target ~ kpi_yoy_lag1, data = train_raw)
    out[[i]] <- tibble(Date     = origin,
                       y_actual = origin_row$target,
                       y_ar1    = as.numeric(predict(fit, newdata = origin_row)))
  }
  bind_rows(out)
}


# -----------------------------------------------------------------------------
# 7. PLOTS  (only the figures actually referenced in the thesis are kept)
# -----------------------------------------------------------------------------
theme_thesis <- theme_minimal(base_size = 11, base_family = "serif") +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40", margin = margin(b = 10)),
    strip.text       = element_text(face = "italic", size = 10),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 8, color = "grey30"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 9),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey85"),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(0.8, "lines"),
    plot.margin      = margin(10, 12, 10, 10)
  )

plot_walker_vs_actual <- function(forecasts) {
  forecasts %>%
    pivot_longer(c(y_actual, y_hat, y_rw), names_to = "series", values_to = "value") %>%
    mutate(series = factor(series,
                           levels = c("y_actual", "y_hat", "y_rw"),
                           labels = c("Actual", "Walker", "Random walk"))) %>%
    ggplot(aes(x = Date, y = value, color = series, linetype = series)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values    = c("Actual" = "black",
                                     "Walker" = "steelblue",
                                     "Random walk" = "grey50")) +
    scale_linetype_manual(values = c("Actual" = "solid",
                                     "Walker" = "solid",
                                     "Random walk" = "dashed")) +
    labs(title = sprintf("Walker vs random walk (h=%d)", H),
         y = "YoY KPI inflation (%)", x = NULL, color = NULL, linetype = NULL) +
    theme_thesis
}

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})


walker <- read_csv("walker_forecasts.csv", show_col_types = FALSE) %>%
  mutate(Date = as.Date(as.character(Date))) %>%
  select(Date, Actual = y_actual, Walker = y_hat)

xgb <- read_csv("xgb_v8_harmonized_predictions.csv", show_col_types = FALSE)
date_col <- intersect(c("Date", "date"), names(xgb))[[1]]
xgb <- xgb %>%
  mutate(Date = as.Date(as.character(.data[[date_col]]))) %>%
  select(Date, XGBoost = predicted_raw)

panel <- walker %>% inner_join(xgb, by = "Date")

# RMSE i legendforklaringen
rmse_walker <- sqrt(mean((panel$Actual - panel$Walker)^2))
rmse_xgb    <- sqrt(mean((panel$Actual - panel$XGBoost)^2))


plot_data <- panel %>%
  pivot_longer(c(Actual, Walker, XGBoost),
               names_to = "series", values_to = "value") %>%
  mutate(series = factor(series,
                         levels = c("Actual", "Walker", "XGBoost"),
                         labels = c("Actual",
                                    sprintf("Walker (RMSE = %.3f)", rmse_walker),
                                    sprintf("XGBoost (RMSE = %.3f)", rmse_xgb))))

theme_thesis <- theme_minimal(base_size = 11, base_family = "serif") +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40", margin = margin(b = 10)),
    axis.title       = element_text(size = 10),
    axis.text        = element_text(size = 8, color = "grey30"),
    legend.position  = "bottom",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 9),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey85"),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(10, 12, 10, 10)
  )

p <- ggplot(plot_data, aes(x = Date, y = value, color = series, linetype = series)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = setNames(
    c("black", "steelblue", "#D55E00"),
    levels(plot_data$series)
  )) +
  scale_linetype_manual(values = setNames(
    c("solid", "solid", "dashed"),
    levels(plot_data$series)
  )) +
  labs(
    title    = "Walker vs XGBoost — h = 3 inflation forecast",
    subtitle = sprintf("Out-of-sample period: %s – %s",
                       format(min(panel$Date), "%b %Y"),
                       format(max(panel$Date), "%b %Y")),
    x = NULL, y = "YoY KPI inflation (%)",
    color = NULL, linetype = NULL
  ) +
  theme_thesis

p

ggsave("walker_vs_xgboost.png", p, width = 10, height = 5, dpi = 200)

plot_block_contributions <- function(joint) {
  joint %>%
    select(Date, ends_with("_SHAP"), ends_with("_WalkerAdj")) %>%
    pivot_longer(-Date, names_to = "series", values_to = "value") %>%
    mutate(model     = if_else(grepl("_SHAP$", series), "SHAP", "Walker adjusted"),
           component = sub("_(SHAP|WalkerAdj)$", "", series)) %>%
    filter(component %in% BLOCK_NAMES) %>%
    ggplot(aes(x = Date, y = value, linetype = model)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.45) +
    facet_wrap(~ component, scales = "free_y", ncol = 2) +
    scale_linetype_manual(values = c("SHAP" = "solid", "Walker adjusted" = "dashed")) +
    labs(title    = "SHAP versus adjusted Walker decomposition",
         subtitle = "Block-level contributions, 2020–2025",
         x = NULL, y = "Contribution (pp)", linetype = NULL) +
    theme_thesis
}

plot_residual <- function(res_long) {
  regime_colours <- c("COVID (2020-2021)" = "#5E81AC",
                      "Energy Crisis"     = "#BF616A",
                      "Disinflation"      = "#A3BE8C",
                      "Normalization"     = "#B48EAD")
  ggplot(res_long, aes(x = Date, y = Delta, fill = regime)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_col(alpha = 0.85, width = 25) +
    facet_wrap(~ block, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = regime_colours) +
    labs(title    = "Harmonised attribution residual by block",
         subtitle = expression(Delta[kt] == C[kt]^{XGB} - C[kt]^{TVP}),
         x = NULL, y = "Residual (pp)", fill = "Regime") +
    theme_thesis
}

first_non_na_date <- function(data, var) {
  min(data$Date[!is.na(data[[var]])], na.rm = TRUE)
}


df_oil <- df_master %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= min(Date[!is.na(oil_price_nok)], na.rm = TRUE))

oljeplot <- ggplot(df_oil, aes(x = Date, y = oil_price_nok)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(suffix = " kr", decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Oil price (NOK)",
    x = "Date",
    y = "Price per barrel"
  ) +
  theme_thesis

oljeplot

ggsave("oljeplot.png", oljeplot, width = 10, height = 5, dpi = 200)


arbledigplot <- df_master %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= min(Date[!is.na(unemployment)], na.rm = TRUE)) %>%
  ggplot(aes(x = Date, y = unemployment)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(suffix = " %", decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Unemployment as a share of the labor force",
    x = "Date",
    y = "Unemployment (%)"
  ) +
  theme_thesis

arbledigplot

ggsave("unemployment.png", arbledigplot, width = 10, height = 5, dpi = 200)


kpiplot <- df_master %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= min(Date[!is.na(kpi)], na.rm = TRUE)) %>%
  ggplot(aes(x = Date, y = kpi)) +
  geom_line() +
  labs(
    title = "Consumer Price Index over time (Index 2015 = 100)",
    x = "Date",
    y = "CPI"
  ) +
  theme_thesis
kpiplot

ggsave("cpi.png", kpiplot, width = 10, height = 5, dpi = 200)


ppiplot <- df_master %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= min(Date[!is.na(ppi)], na.rm = TRUE)) %>%
  ggplot(aes(x = Date, y = ppi)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Producer price index over time (2021=100)",
    x = "Date",
    y = "PPI"
  ) +
  theme_thesis

ppiplot

ggsave("PPI.png", ppiplot, width = 10, height = 5, dpi = 200)


usdplot <- df_master %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= first_non_na_date(., "usd_nok")) %>%
  ggplot(aes(x = Date, y = usd_nok)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "USD/NOK exchange rate",
    x = "Date",
    y = "NOK per USD"
  ) +
  theme_thesis

usdplot

ggsave("usd_nok.png", usdplot, width = 10, height = 5, dpi = 200)


eurplot <- df_master %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= first_non_na_date(., "eur_nok")) %>%
  ggplot(aes(x = Date, y = eur_nok)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "EUR/NOK exchange rate",
    x = "Date",
    y = "NOK per EUR"
  ) +
  theme_thesis

eurplot

ggsave("eur_nok.png", eurplot, width = 10, height = 5, dpi = 200)


renteplot <- df_master %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= first_non_na_date(., "rente")) %>%
  ggplot(aes(x = Date, y = rente)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(suffix = " %", decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Policy rate over time",
    x = "Date",
    y = "Policy rate (%)"
  ) +
  theme_thesis

renteplot

ggsave("policy_rate.png", renteplot, width = 10, height = 5, dpi = 200)


tradeplot <- df_master %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= min(Date[!is.na(import) | !is.na(eksport)], na.rm = TRUE)) %>%
  select(Date, import, eksport) %>%
  pivot_longer(
    cols = c(import, eksport),
    names_to = "series",
    values_to = "value"
  ) %>%
  mutate(
    series = recode(
      series,
      import = "Import",
      eksport = "Export"
    )
  ) %>%
  ggplot(aes(x = Date, y = value, linetype = series)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(
      scale = 1e-9,
      suffix = " bn NOK",
      decimal.mark = ","
    ),
    breaks = pretty_breaks(n = 6)
  )+
  labs(
    title = "Norwegian foreign trade",
    x = "Date",
    y = "Billion NOK",
    linetype = NULL
  ) +
  theme_thesis

tradeplot

ggsave("trade.png", tradeplot, width = 10, height = 5, dpi = 200)


oljeusdplot <- df_master %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= first_non_na_date(., "oljepris_USD")) %>%
  ggplot(aes(x = Date, y = oljepris_USD)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(suffix = " USD", decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Oil price in USD",
    x = "Date",
    y = "USD per barrel"
  ) +
  theme_thesis

oljeusdplot

ggsave("olje_usd.png", oljeusdplot, width = 10, height = 5, dpi = 200)


tradeusdplot <- df_master %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Date >= min(Date[!is.na(import_USD) | !is.na(eksport_USD)], na.rm = TRUE)) %>%
  select(Date, import_USD, eksport_USD) %>%
  pivot_longer(
    cols = c(import_USD, eksport_USD),
    names_to = "series",
    values_to = "value"
  ) %>%
  mutate(
    series = recode(
      series,
      import_USD = "Import",
      eksport_USD = "Export"
    )
  ) %>%
  ggplot(aes(x = Date, y = value, linetype = series)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(
      scale = 1e-9,
      suffix = " bn NOK",
      decimal.mark = ","
    ),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Norwegian foreign trade in USD",
    x = "Date",
    y = "Billion USD",
    linetype = NULL
  ) +
  theme_thesis

tradeusdplot

ggsave("trade_usd.png", tradeusdplot, width = 10, height = 5, dpi = 200)


# -----------------------------------------------------------------------------
# 8. EXECUTION (driven by the RUN_* flags in section 0)
# -----------------------------------------------------------------------------
df_master <- read_csv("master_data.csv", show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date)) %>% arrange(Date)
cat(sprintf("Loaded master_data.csv: %s to %s (%d rows)\n",
            min(df_master$Date), max(df_master$Date), nrow(df_master)))

# --- Stage 1: Walker backtests ----------------------------------------------
for (spec_name in names(RUN_WALKER)) {
  if (isTRUE(RUN_WALKER[[spec_name]])) run_walker_spec(spec_name, df_master)
}

# --- Stage 2: Walker harmonisation per spec ---------------------------------
walker_handles <- list()
if (isTRUE(RUN_HARMONIZE_W)) {
  for (spec_name in names(RUN_WALKER)) {
    contrib_path <- paste0(WALKER_SPECS[[spec_name]]$out_prefix, "contributions.csv")
    if (file.exists(contrib_path)) {
      walker_handles[[spec_name]] <- harmonize_walker(
        contributions_path = contrib_path,
        out_prefix         = WALKER_SPECS[[spec_name]]$out_prefix
      )
    } else {
      message("Skipping harmonization for ", spec_name,
              ": ", contrib_path, " not found.")
    }
  }
}

# --- Stage 3: XGB SHAP harmonisation ----------------------------------------
xgb_handles <- list()
if (isTRUE(RUN_HARMONIZE_XGB)) {
  xgb_inputs <- list(
    main    = "xgb_v8_harmonized_shap_group_signed.csv",
    usd_rob = "xgb_harmonized_robustness_shap_group_signed.csv"
  )
  xgb_prefixes <- list(main = "xgb_main_", usd_rob = "xgb_robustness_")
  for (n in names(xgb_inputs)) {
    if (file.exists(xgb_inputs[[n]])) {
      xgb_handles[[n]] <- harmonize_xgb(
        signed_path = xgb_inputs[[n]], out_prefix = xgb_prefixes[[n]]
      )
    } else {
      message("Skipping XGB harmonization for ", n,
              ": ", xgb_inputs[[n]], " not found.")
    }
  }
}

# --- Stage 4: SHAP <-> Walker comparison + residual -------------------------
joint <- NULL
if (isTRUE(RUN_COMPARE)) {
  walker_signed_path <- paste0(WALKER_SPECS$main$out_prefix,
                               "harmonized_block_signed.csv")
  xgb_signed_path    <- "xgb_v8_harmonized_shap_group_signed.csv"
  if (file.exists(walker_signed_path) && file.exists(xgb_signed_path)) {
    joint <- build_joint_attributions(walker_signed_path, xgb_signed_path)
    write_csv(joint, "comparison_harmonised_attributions.csv")
    
    diag <- sign_and_corr(joint)
    write_csv(diag, "comparison_routing_compare.csv")
    cat("\nSign agreement and correlations:\n"); print(diag)
    
    res <- residual_decomposition(joint)
    write_csv(res$long,    "residual_long.csv")
    write_csv(res$metrics, "residual_metrics_long.csv")
    
    p_contrib  <- plot_block_contributions(joint)
    p_residual <- plot_residual(res$long)
    ggsave("comparison_plot_contributions.png", p_contrib,
           width = 13, height = 8, dpi = 150)
    ggsave("residual_plot.png", p_residual, width = 12, height = 9, dpi = 150)
  } else {
    message("RUN_COMPARE=TRUE but inputs missing: ",
            walker_signed_path, " or ", xgb_signed_path)
  }
}

# --- Stage 5: Regime-cutoff sensitivity -------------------------------------
if (isTRUE(RUN_REGIME_ROB) && !is.null(joint)) {
  rob <- regime_cutoff_robustness(joint, shifts = c(-3, 0, 3))
  rob_wide <- rob %>%
    pivot_wider(names_from = partition, values_from = share) %>%
    mutate(`Diff (max-min)` = pmax(`Shift -3m`, Original, `Shift +3m`,
                                   na.rm = TRUE) -
             pmin(`Shift -3m`, Original, `Shift +3m`,
                  na.rm = TRUE))
  write_csv(rob,      "robustness_regime_long.csv")
  write_csv(rob_wide, "robustness_regime_wide.csv")
  cat("\nRegime sensitivity (±3m):\n"); print(rob_wide)
}

# --- Stage 6: Baselines + accuracy table ------------------------------------
# Read csv file


fc <- read_csv("walker_usd_external_forecasts.csv", show_col_types = FALSE)
bind_rows(
  forecast_metrics(fc$y_actual, fc$y_rw,  "Random walk (USD spec)"),
  forecast_metrics(fc$y_actual, fc$y_hat, "Walker (USD spec)")
)


if (isTRUE(RUN_BASELINES)) {
  walker_path <- paste0(WALKER_SPECS$main$out_prefix, "forecasts.csv")
  if (file.exists(walker_path)) {
    fc  <- read_csv(walker_path, show_col_types = FALSE)
    ar1 <- run_ar1_backtest(df_master)
    write_csv(ar1, "ar1_forecasts.csv")
    accuracy <- bind_rows(
      forecast_metrics(fc$y_actual,  fc$y_rw,  "Random walk"),
      forecast_metrics(fc$y_actual,  fc$y_hat, "Walker"),
      forecast_metrics(ar1$y_actual, ar1$y_ar1, "AR(1)"),
    )
    write_csv(accuracy, "accuracy_table.csv")
    cat("\nForecast accuracy:\n"); print(accuracy)
  } else {
    message("RUN_BASELINES=TRUE but ", walker_path, " not found.")
  }
}


fc_main <- read_csv("walker_forecasts.csv",              show_col_types = FALSE)
fc_usd  <- read_csv("walker_usd_external_forecasts.csv", show_col_types = FALSE)
ar1     <- read_csv("ar1_forecasts.csv",                 show_col_types = FALSE)

bind_rows(
  forecast_metrics(fc_main$y_actual, fc_main$y_rw,  "Random walk"),
  forecast_metrics(fc_main$y_actual, fc_main$y_hat, "Walker (main)"),
  forecast_metrics(fc_usd$y_actual,  fc_usd$y_hat,  "Walker (USD external)"),
  forecast_metrics(ar1$y_actual,     ar1$y_ar1,     "AR(1)")
)



# =============================================================================
# Generer og lagre alle figurer fra cachede CSV-er
# Kjør etter source("thesis_pipeline.R") når alle outputs ligger på disk.
# =============================================================================

# --- Last inn det vi trenger -----------------------------------------------
forecasts <- read_csv("walker_forecasts.csv", show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date))
joint <- read_csv("comparison_harmonised_attributions.csv", show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date))
res_long <- read_csv("residual_long.csv", show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date))

# --- Bygg plottene ---------------------------------------------------------
p_walker   <- plot_walker_vs_actual(forecasts)
p_contrib  <- plot_block_contributions(joint)
p_residual <- plot_residual(res_long)

# --- Lagre alle ------------------------------------------------------------
ggsave("walker_forecast_plot.png",          p_walker,   width = 10, height = 5, dpi = 200)
ggsave("comparison_plot_contributions.png", p_contrib,  width = 13, height = 8, dpi = 200)
ggsave("residual_plot.png",                 p_residual, width = 12, height = 9, dpi = 200)

# --- Vis i RStudio --------------------------------------------------------
print(p_walker)
print(p_contrib)
print(p_residual)



BLOCK_NAMES <- c("AR", "PPI", "Monetary", "FX", "Oil", "Trade", "Labour")

# =============================================================================
# MAIN SPEC
# =============================================================================

# --- Last inn Walker block_signed (main) ------------------------------------
walker_main <- read_csv("walker_harmonized_block_signed.csv",
                        show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date)) %>%
  rename(AR_WalkerAdj       = AR,
         PPI_WalkerAdj      = PPI,
         Monetary_WalkerAdj = Monetary,
         FX_WalkerAdj       = FX,
         Oil_WalkerAdj      = Oil,
         Trade_WalkerAdj    = Trade,
         Labour_WalkerAdj   = Labour)

# --- Last inn XGB SHAP signed (main) ----------------------------------------
xgb_main_raw <- read_csv("xgb_v8_harmonized_shap_group_signed.csv",
                         show_col_types = FALSE)
xgb_main <- xgb_main_raw %>%
  mutate(Date = as.Date(date)) %>%
  rename(AR_SHAP       = `AR (inflation)`,
         PPI_SHAP      = `PPI / cost-push`,
         Monetary_SHAP = `Monetary policy`,
         FX_SHAP       = `FX`,
         Oil_SHAP      = `Oil`,
         Trade_SHAP    = `Trade`,
         Labour_SHAP   = `Labour market`) %>%
  select(Date, ends_with("_SHAP"), actual, predicted_raw, shap_base)

# --- Joint tabell main ------------------------------------------------------
joint_main <- xgb_main %>%
  inner_join(walker_main, by = "Date") %>%
  arrange(Date)


# =============================================================================
# USD-EXTERNAL SPEC
# =============================================================================

# --- Last inn Walker block_signed (USD external) ----------------------------
walker_usd <- read_csv("walker_usd_external_harmonized_block_signed.csv",
                       show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date)) %>%
  rename(AR_WalkerAdj       = AR,
         PPI_WalkerAdj      = PPI,
         Monetary_WalkerAdj = Monetary,
         FX_WalkerAdj       = FX,
         Oil_WalkerAdj      = Oil,
         Trade_WalkerAdj    = Trade,
         Labour_WalkerAdj   = Labour)

# --- Last inn XGB SHAP signed (USD robustness) ------------------------------
xgb_usd_raw <- read_csv("xgb_harmonized_robustness_shap_group_signed.csv",
                        show_col_types = FALSE)
xgb_usd <- xgb_usd_raw %>%
  rename(Date = date,
         AR_SHAP       = `AR (inflation)`,
         PPI_SHAP      = `PPI / cost-push`,
         Monetary_SHAP = `Monetary policy`,
         FX_SHAP       = `FX`,
         Oil_SHAP      = `Oil`,
         Trade_SHAP    = `Trade`,
         Labour_SHAP   = `Labour market`) %>%
  select(Date, ends_with("_SHAP"), actual, predicted_raw, shap_base)

# --- Joint tabell USD external ----------------------------------------------
joint_usd <- xgb_usd %>%
  inner_join(walker_usd, by = "Date") %>%
  arrange(Date)


# =============================================================================
# TABELL 1: MEAN SIGNED CONTRIBUTION
# =============================================================================
mean_signed <- tibble(
  block               = BLOCK_NAMES,
  SHAP_main           = sapply(BLOCK_NAMES, \(b) mean(joint_main[[paste0(b, "_SHAP")]],      na.rm = TRUE)),
  Walker_main         = sapply(BLOCK_NAMES, \(b) mean(joint_main[[paste0(b, "_WalkerAdj")]], na.rm = TRUE)),
  SHAP_usd_external   = sapply(BLOCK_NAMES, \(b) mean(joint_usd[[paste0(b, "_SHAP")]],       na.rm = TRUE)),
  Walker_usd_external = sapply(BLOCK_NAMES, \(b) mean(joint_usd[[paste0(b, "_WalkerAdj")]],  na.rm = TRUE))
)

write_csv(mean_signed, "comparison_mean_signed_by_spec.csv")
cat("\n=== MEAN SIGNED CONTRIBUTION (pp) ===\n")
print(mean_signed %>% mutate(across(where(is.numeric), \(x) round(x, 3))))


# =============================================================================
# TABELL 2: MEAN ABSOLUTE CONTRIBUTION
# =============================================================================
mean_abs <- tibble(
  block               = BLOCK_NAMES,
  SHAP_main           = sapply(BLOCK_NAMES, \(b) mean(abs(joint_main[[paste0(b, "_SHAP")]]),      na.rm = TRUE)),
  Walker_main         = sapply(BLOCK_NAMES, \(b) mean(abs(joint_main[[paste0(b, "_WalkerAdj")]]), na.rm = TRUE)),
  SHAP_usd_external   = sapply(BLOCK_NAMES, \(b) mean(abs(joint_usd[[paste0(b, "_SHAP")]]),       na.rm = TRUE)),
  Walker_usd_external = sapply(BLOCK_NAMES, \(b) mean(abs(joint_usd[[paste0(b, "_WalkerAdj")]]),  na.rm = TRUE))
)

write_csv(mean_abs, "comparison_mean_abs_by_spec.csv")
cat("\n=== MEAN ABSOLUTE CONTRIBUTION (pp) ===\n")
print(mean_abs %>% mutate(across(where(is.numeric), \(x) round(x, 3))))


# =============================================================================
# TABELL 3: KORRELASJON SHAP vs WALKER
# =============================================================================
corr_table <- tibble(
  block                  = BLOCK_NAMES,
  pct_agree_main         = sapply(BLOCK_NAMES, \(b) mean(sign(joint_main[[paste0(b, "_SHAP")]]) ==
                                                           sign(joint_main[[paste0(b, "_WalkerAdj")]]),
                                                         na.rm = TRUE)),
  corr_signed_main       = sapply(BLOCK_NAMES, \(b) cor(joint_main[[paste0(b, "_SHAP")]],
                                                        joint_main[[paste0(b, "_WalkerAdj")]],
                                                        use = "complete.obs")),
  corr_abs_main          = sapply(BLOCK_NAMES, \(b) cor(abs(joint_main[[paste0(b, "_SHAP")]]),
                                                        abs(joint_main[[paste0(b, "_WalkerAdj")]]),
                                                        use = "complete.obs")),
  pct_agree_usd          = sapply(BLOCK_NAMES, \(b) mean(sign(joint_usd[[paste0(b, "_SHAP")]]) ==
                                                           sign(joint_usd[[paste0(b, "_WalkerAdj")]]),
                                                         na.rm = TRUE)),
  corr_signed_usd        = sapply(BLOCK_NAMES, \(b) cor(joint_usd[[paste0(b, "_SHAP")]],
                                                        joint_usd[[paste0(b, "_WalkerAdj")]],
                                                        use = "complete.obs")),
  corr_abs_usd           = sapply(BLOCK_NAMES, \(b) cor(abs(joint_usd[[paste0(b, "_SHAP")]]),
                                                        abs(joint_usd[[paste0(b, "_WalkerAdj")]]),
                                                        use = "complete.obs"))
)

write_csv(corr_table, "comparison_correlations_by_spec.csv")
cat("\n=== KORRELASJON SHAP vs Walker ===\n")
print(corr_table %>% mutate(across(where(is.numeric), \(x) round(x, 3))))







suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
})

BLOCK_NAMES  <- c("AR", "PPI", "Monetary", "FX", "Oil", "Trade", "Labour")
BLOCK_LABELS <- c(
  "AR"       = "AR (inflation)",
  "PPI"      = "PPI / cost-push",
  "Monetary" = "Monetary policy",
  "FX"       = "FX",
  "Oil"      = "Oil",
  "Trade"    = "Trade",
  "Labour"   = "Labour market"
)
REGIME_LEVELS <- c("COVID (2020-2021)", "Energy Crisis", "Disinflation", "Normalization")
REGIME_LABELS <- c("COVID", "Energy", "Disinfl.", "Normal")

# =============================================================================
# 1. LAST INN OG BYGG JOINT-TABELLER FOR BEGGE SPECS
# =============================================================================
load_spec <- function(walker_file, xgb_file) {
  walker <- read_csv(walker_file, show_col_types = FALSE)
  if (!inherits(walker$Date, "Date")) walker$Date <- as.Date(as.character(walker$Date))
  walker <- walker %>%
    rename(AR_WalkerAdj = AR, PPI_WalkerAdj = PPI, Monetary_WalkerAdj = Monetary,
           FX_WalkerAdj = FX, Oil_WalkerAdj = Oil, Trade_WalkerAdj = Trade,
           Labour_WalkerAdj = Labour)
  
  xgb <- read_csv(xgb_file, show_col_types = FALSE)
  date_col <- intersect(c("Date", "date"), names(xgb))[[1]]
  xgb$Date <- as.Date(as.character(xgb[[date_col]]))
  xgb <- xgb %>%
    rename(AR_SHAP       = `AR (inflation)`,
           PPI_SHAP      = `PPI / cost-push`,
           Monetary_SHAP = `Monetary policy`,
           FX_SHAP       = `FX`,
           Oil_SHAP      = `Oil`,
           Trade_SHAP    = `Trade`,
           Labour_SHAP   = `Labour market`) %>%
    select(Date, ends_with("_SHAP"), actual, predicted_raw, shap_base)
  
  xgb %>%
    inner_join(walker, by = "Date") %>%
    mutate(regime = factor(case_when(
      Date < as.Date("2021-06-01") ~ "COVID (2020-2021)",
      Date < as.Date("2023-01-01") ~ "Energy Crisis",
      Date < as.Date("2024-06-01") ~ "Disinflation",
      TRUE                          ~ "Normalization"
    ), levels = REGIME_LEVELS)) %>%
    arrange(Date)
}

joint_main <- load_spec(
  "walker_harmonized_block_signed.csv",
  "xgb_v8_harmonized_shap_group_signed.csv"
)
joint_usd <- load_spec(
  "walker_usd_external_harmonized_block_signed.csv",
  "xgb_harmonized_robustness_shap_group_signed.csv"
)

# =============================================================================
# 2. BEREGNINGSFUNKSJONER (kalles separat for hver spec)
# =============================================================================

# --- Mean attribution shares (signed S og absolute A), globalt --------------
compute_mean_shares <- function(joint) {
  abs_shap <- as.matrix(joint[, paste0(BLOCK_NAMES, "_SHAP")])     |> abs()
  abs_walk <- as.matrix(joint[, paste0(BLOCK_NAMES, "_WalkerAdj")]) |> abs()
  tot_shap <- rowSums(abs_shap, na.rm = TRUE)
  tot_walk <- rowSums(abs_walk, na.rm = TRUE)
  
  out <- tibble(block = BLOCK_NAMES)
  for (b in BLOCK_NAMES) {
    s <- joint[[paste0(b, "_SHAP")]]
    w <- joint[[paste0(b, "_WalkerAdj")]]
    out[out$block == b, "S_SHAP"]   <- mean(ifelse(tot_shap > 0, s     / tot_shap, NA), na.rm = TRUE)
    out[out$block == b, "A_SHAP"]   <- mean(ifelse(tot_shap > 0, abs(s)/ tot_shap, NA), na.rm = TRUE)
    out[out$block == b, "S_Walker"] <- mean(ifelse(tot_walk > 0, w     / tot_walk, NA), na.rm = TRUE)
    out[out$block == b, "A_Walker"] <- mean(ifelse(tot_walk > 0, abs(w)/ tot_walk, NA), na.rm = TRUE)
  }
  out
}

# --- Correlation per blokk (sign agreement, signed r, absolute r) ----------
compute_correlations <- function(joint) {
  bind_rows(lapply(BLOCK_NAMES, function(b) {
    s <- joint[[paste0(b, "_SHAP")]]
    w <- joint[[paste0(b, "_WalkerAdj")]]
    tibble(
      block       = b,
      pct_agree   = mean(sign(s) == sign(w), na.rm = TRUE),
      corr_signed = cor(s, w, use = "complete.obs"),
      corr_abs    = cor(abs(s), abs(w), use = "complete.obs")
    )
  }))
}

# --- Block importance by regime (mean A-share x 100, per regime) -----------
compute_block_importance <- function(joint) {
  abs_shap <- as.matrix(joint[, paste0(BLOCK_NAMES, "_SHAP")])     |> abs()
  abs_walk <- as.matrix(joint[, paste0(BLOCK_NAMES, "_WalkerAdj")]) |> abs()
  tot_shap <- rowSums(abs_shap, na.rm = TRUE)
  tot_walk <- rowSums(abs_walk, na.rm = TRUE)
  
  per_row <- tibble(regime = joint$regime)
  for (b in BLOCK_NAMES) {
    per_row[[paste0("A_", b, "_SHAP")]]   <-
      ifelse(tot_shap > 0, abs(joint[[paste0(b, "_SHAP")]])      / tot_shap, NA)
    per_row[[paste0("A_", b, "_Walker")]] <-
      ifelse(tot_walk > 0, abs(joint[[paste0(b, "_WalkerAdj")]]) / tot_walk, NA)
  }
  
  per_row %>%
    group_by(regime) %>%
    summarise(across(starts_with("A_"), \(x) mean(x, na.rm = TRUE) * 100),
              .groups = "drop") %>%
    pivot_longer(starts_with("A_"),
                 names_to      = c("block", "model"),
                 names_pattern = "^A_([^_]+)_(SHAP|Walker)$",
                 values_to     = "value") %>%
    pivot_wider(names_from = regime, values_from = value) %>%
    mutate(block = factor(block, levels = BLOCK_NAMES),
           model = factor(model, levels = c("SHAP", "Walker"),
                          labels = c("XGBoost (SHAP)", "Walker"))) %>%
    arrange(model, block)
}

# --- Regime-conditional residual summaries (Var og E|Delta|) ---------------
compute_residual_summary <- function(joint) {
  bind_rows(lapply(BLOCK_NAMES, function(b) {
    delta <- joint[[paste0(b, "_SHAP")]] - joint[[paste0(b, "_WalkerAdj")]]
    tibble(block = b, regime = joint$regime, delta = delta)
  })) %>%
    group_by(block, regime) %>%
    summarise(var_delta = var(delta,      na.rm = TRUE),
              mean_abs  = mean(abs(delta), na.rm = TRUE),
              .groups   = "drop") %>%
    mutate(block = factor(block, levels = BLOCK_NAMES))
}


# =============================================================================
# 3. SKRIV ALLE ÅTTE TABELLER TIL KONSOLLEN
# =============================================================================
cat("% =====================================================================\n")
cat("% MAIN TEXT\n")
cat("% =====================================================================\n\n")

compute_mean_shares(joint_main)

compute_correlations(joint_main)

compute_block_importance(joint_main)


compute_residual_summary(joint_main)
 

cat("% =====================================================================\n")
cat("% APPENDIX: USD EXTERNAL ROBUSTNESS\n")
cat("% =====================================================================\n\n")

compute_mean_shares(joint_usd)

compute_correlations(joint_usd)

compute_block_importance(joint_usd)

compute_residual_summary(joint_usd)

