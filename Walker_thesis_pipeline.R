# =============================================================================
# WALKER TVP ESTIMATION PIPELINE — Norwegian CPI Inflation Forecasting
# Master's Thesis Replication Code
# Authors: Seb & Amund
#
# This script reproduces all Walker-related estimation, harmonisation,
# comparison, and figure outputs. Run top to bottom in a single session.
# Stages 2 and 3 (MCMC backtests) are slow (~1.5 hours each); all downstream
# stages read the outputs they produce from disk, so they can be re-run
# independently if the CSV outputs already exist.
#
# Required input files (must be present in the working directory):
#   master_data.csv                               — macro data (all variables)
#   xgb_v8_harmonized_shap_group_signed.csv       — XGBoost SHAP main spec
#   xgb_harmonized_robustness_shap_group_signed.csv — XGBoost SHAP USD spec
#   xgb_v8_harmonized_predictions.csv             — XGBoost point forecasts
# =============================================================================


# =============================================================================
# PART 0 — CONFIGURATION
# =============================================================================

suppressPackageStartupMessages({
  library(walker)
  library(rstan)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(lubridate)
  library(scales)
  library(knitr)
  library(kableExtra)
  library(corrplot)
})

# Set working directory — tries Amund's path first, then Seb's.
for (.p in c(
  "/Users/amundbech1/Library/CloudStorage/OneDrive-UiTOffice365/Data_master/MasterOppgave",
  "~/Library/CloudStorage/OneDrive-UiTOffice365/Amund Skjalgsønn Bech's files - Data_master/MasterOppgave"
)) {
  if (dir.exists(.p)) { setwd(.p); break }
}

# Estimation hyperparameters — identical to the original Walker scripts.
H           <- 3                       # forecast horizon (months ahead)
TEST_START  <- as.Date("2020-01-01")   # first expanding-window forecast origin
MIN_TRAIN_N <- 50                      # minimum training observations required
N_CHAINS    <- 4
N_ITER      <- 2000
N_WARMUP    <- 1000
N_CORES     <- min(4, parallel::detectCores())
SEED        <- 42

set.seed(SEED)
options(mc.cores = N_CORES)
rstan_options(auto_write = TRUE)

# Feature columns — the nine lag-1 macro predictors.
FEATURE_COLS <- c(
  "kpi_yoy_lag1",   # lagged YoY CPI (AR term)
  "ppi_yoy_lag1",   # lagged YoY PPI
  "oil_yoy_lag1",   # lagged YoY oil price
  "usd_nok_lag1",   # lagged USD/NOK exchange rate
  "eur_nok_lag1",   # lagged EUR/NOK exchange rate
  "import_yoy_lag1",  # lagged YoY imports
  "eksport_yoy_lag1", # lagged YoY exports
  "unemp_lag1",     # lagged unemployment rate
  "rente_lag1"      # lagged policy rate
)

# Economic block map — used for attribution grouping.
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

# Regime cutoffs consistent with thesis Table 5.
assign_regime <- function(Date) {
  dplyr::case_when(
    Date < as.Date("2021-06-01") ~ "COVID (2020–2021)",
    Date < as.Date("2023-01-01") ~ "Energy Crisis",
    Date < as.Date("2024-06-01") ~ "Disinflation",
    TRUE                          ~ "Normalization"
  )
}


# =============================================================================
# PART 0b — LOOP HELPER FUNCTIONS
# These four functions are called at every forecast origin inside the MCMC
# backtest loops (Parts 2 and 3). Keeping them as functions avoids ~60 lines
# of repeated code per loop and makes the loop body readable.
# =============================================================================

# Standardise cols of new_df using means/SDs computed on train_df only.
# This enforces the real-time information constraint: no look-ahead in scaling.
scale_with_train_stats <- function(train_df, new_df, cols) {
  means <- sapply(train_df[, cols, drop = FALSE], mean, na.rm = TRUE)
  sds   <- sapply(train_df[, cols, drop = FALSE], sd,   na.rm = TRUE)
  bad   <- is.na(sds) | sds == 0
  if (any(bad)) {
    stop(sprintf("Zero/NA SD in training data for: %s",
                 paste(names(sds)[bad], collapse = ", ")))
  }
  train_scaled <- train_df
  new_scaled   <- new_df
  for (col in cols) {
    train_scaled[[col]] <- (train_df[[col]] - means[[col]]) / sds[[col]]
    new_scaled[[col]]   <- (new_df[[col]]   - means[[col]]) / sds[[col]]
  }
  list(train_scaled = train_scaled, new_scaled = new_scaled,
       means = means, sds = sds)
}

# Extract MCMC diagnostics: divergent transitions, worst R-hat, lowest ESS.
extract_fit_diagnostics <- function(fit) {
  sp    <- get_sampler_params(fit$stanfit, inc_warmup = FALSE)
  n_div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  summ  <- summary(fit$stanfit)$summary
  list(n_div    = n_div,
       max_rhat = max(summ[, "Rhat"], na.rm = TRUE),
       min_ess  = min(summ[, "n_eff"], na.rm = TRUE))
}

# One-step-ahead prediction using filtered (last time-point) posterior means:
#   y_hat = beta_0(T) + sum_j beta_j(T) * x_j(T)
predict_last_state <- function(fit, new_x, feature_cols) {
  coefs_last <- coef(fit) %>%
    ungroup() %>%
    filter(time == max(time)) %>%
    select(beta, mean)
  if (!"(Intercept)" %in% coefs_last$beta)
    stop("'(Intercept)' not found in coef(fit).")
  intercept_val <- coefs_last$mean[coefs_last$beta == "(Intercept)"]
  beta_tbl <- coefs_last %>% filter(beta %in% feature_cols)
  missing  <- setdiff(feature_cols, beta_tbl$beta)
  if (length(missing) > 0)
    stop(sprintf("Missing coefficients for: %s", paste(missing, collapse = ", ")))
  beta_vals <- setNames(beta_tbl$mean, beta_tbl$beta)
  x_vals    <- as.numeric(new_x[1, names(beta_vals), drop = TRUE])
  as.numeric(intercept_val + sum(beta_vals * x_vals))
}

# Extract the filtered posterior summary at one forecast origin.
extract_filtered_coefs <- function(fit, origin_date) {
  coef(fit) %>%
    ungroup() %>%
    filter(time == max(time)) %>%
    select(beta, mean, sd, `2.5%`, `50%`, `97.5%`) %>%
    mutate(origin = origin_date)
}

# Accuracy metrics helper — used in Part 10 for multiple model comparisons.
safe_mape <- function(actual, predicted) {
  ok <- is.finite(actual) & is.finite(predicted) & actual != 0
  if (!any(ok)) return(NA_real_)
  mean(abs((actual[ok] - predicted[ok]) / actual[ok]), na.rm = TRUE) * 100
}

forecast_metrics <- function(actual, predicted, model_name) {
  err <- predicted - actual
  data.frame(
    Model = model_name,
    RMSE  = sqrt(mean(err^2, na.rm = TRUE)),
    MAE   = mean(abs(err), na.rm = TRUE),
    MAPE  = safe_mape(actual, predicted),
    Bias  = mean(err, na.rm = TRUE),
    N     = sum(is.finite(actual) & is.finite(predicted))
  )
}


# =============================================================================
# PART 1 — DATA LOADING AND FEATURE ENGINEERING (MAIN SPEC)
# =============================================================================

master_df <- read.csv("master_data.csv", stringsAsFactors = FALSE)
master_df$Date <- as.Date(master_df$Date)
master_df <- master_df %>% arrange(Date)

cat("Raw data:", as.character(min(master_df$Date)), "to", as.character(max(master_df$Date)),
    "(", nrow(master_df), "rows)\n")

# Lag-1 information set: all features are lagged by one month so that the
# information set at each forecast origin matches the XGBoost specification.
# YoY growth rates are computed first, then lagged.
master_df <- master_df %>%
  mutate(
    kpi_yoy_raw     = (kpi           / lag(kpi,           12) - 1) * 100,
    ppi_yoy_raw     = (ppi           / lag(ppi,           12) - 1) * 100,
    oil_yoy_raw     = (oil_price_nok / lag(oil_price_nok, 12) - 1) * 100,
    import_yoy_raw  = (import        / lag(import,        12) - 1) * 100,
    eksport_yoy_raw = (eksport       / lag(eksport,       12) - 1) * 100,

    kpi_yoy_lag1     = lag(kpi_yoy_raw, 1),
    ppi_yoy_lag1     = lag(ppi_yoy_raw, 1),
    oil_yoy_lag1     = lag(oil_yoy_raw, 1),
    usd_nok_lag1     = lag(usd_nok, 1),
    eur_nok_lag1     = lag(eur_nok, 1),
    import_yoy_lag1  = lag(import_yoy_raw, 1),
    eksport_yoy_lag1 = lag(eksport_yoy_raw, 1),
    unemp_lag1       = lag(unemployment, 1),
    rente_lag1       = lag(rente, 1),

    # Target: YoY CPI inflation H = 3 months ahead
    target = lead(kpi_yoy_raw, H)
  )

model_df <- master_df %>%
  select(Date, target, kpi_yoy_raw, all_of(FEATURE_COLS)) %>%
  filter(complete.cases(.))

cat("Model data:", as.character(min(model_df$Date)), "to",
    as.character(max(model_df$Date)), "(", nrow(model_df), "rows)\n")

test_origins <- model_df %>%
  filter(Date >= TEST_START, !is.na(target)) %>%
  pull(Date)

cat("Forecast origins:", length(test_origins), "\n")


# =============================================================================
# PART 2 — WALKER BACKTEST: MAIN SPECIFICATION
# Expanding window, train = {t : t < origin}, real-time coefficient extraction.
# Formula: target ~ -1 + rw1(~ 1 + features, beta=c(0,10), sigma=c(2,0.01))
# The time-varying intercept is included via "1 +" inside rw1().
# All specifications include a time-varying intercept.
# =============================================================================

LOG_FILE <- "walker_diagnostics.txt"

rw_formula_main <- as.formula(
  paste0("target ~ -1 + rw1(~ 1 + ",
         paste(FEATURE_COLS, collapse = " + "),
         ", beta = c(0, 10), sigma = c(2, 0.01))")
)

results_main <- data.frame(
  Date     = as.Date(character()),
  y_actual = numeric(),
  y_hat    = numeric(),
  y_rw     = numeric(),
  n_div    = integer(),
  max_rhat = numeric(),
  min_ess  = numeric(),
  elapsed  = numeric(),
  stringsAsFactors = FALSE
)

filtered_coefs_list <- list()
x_values_list       <- list()

cat("Walker Expanding Window Backtesting (main spec)\n", file = LOG_FILE)
cat(paste0("Started: ", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)
cat(paste0("Config: chains=", N_CHAINS, " iter=", N_ITER,
           " warmup=", N_WARMUP, " adapt_delta=0.95 H=", H, "\n\n"),
    file = LOG_FILE, append = TRUE)

cat("\n=== EXPANDING WINDOW BACKTEST (main spec) ===\n")
cat("Total origins:", length(test_origins), "\n")
cat("Estimated time: ~", round(length(test_origins) * 1.5 / 60, 1),
    "hours\n\n")

for (i in seq_along(test_origins)) {
  origin         <- test_origins[i]
  train_raw      <- model_df %>% filter(Date <  origin)
  origin_row_raw <- model_df %>% filter(Date == origin)

  if (nrow(origin_row_raw) != 1) {
    msg <- sprintf("[%s] FAILED: origin row not unique\n", origin)
    cat(msg); cat(msg, file = LOG_FILE, append = TRUE); next
  }

  y_actual <- origin_row_raw$target
  y_rw     <- origin_row_raw$kpi_yoy_lag1  # random walk = last known YoY CPI
  if (is.na(y_actual) || is.na(y_rw) || nrow(train_raw) < MIN_TRAIN_N) next

  cat(sprintf("[%d/%d] Origin: %s (n_train=%d) ... ",
              i, length(test_origins), origin, nrow(train_raw)))

  tryCatch({
    scaled       <- scale_with_train_stats(train_raw, origin_row_raw, FEATURE_COLS)
    train_scaled <- scaled$train_scaled
    new_x_scaled <- scaled$new_scaled[, FEATURE_COLS, drop = FALSE]

    t_start <- Sys.time()
    fit <- walker(
      formula       = rw_formula_main,
      data          = train_scaled,
      sigma_y_prior = c(2, 0.01),
      chains        = N_CHAINS,
      iter          = N_ITER,
      warmup        = N_WARMUP,
      cores         = N_CORES,
      refresh       = 0,
      control       = list(adapt_delta = 0.95, max_treedepth = 12)
    )
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

    diag  <- extract_fit_diagnostics(fit)
    y_hat <- predict_last_state(fit, new_x_scaled, FEATURE_COLS)

    # Store filtered (real-time) coefficients and scaled x-values.
    filtered_coefs_list[[i]] <- extract_filtered_coefs(fit, origin)
    x_values_list[[i]] <- as.data.frame(new_x_scaled) %>%
      pivot_longer(everything(), names_to = "beta", values_to = "x_scaled") %>%
      mutate(origin = origin)

    results_main <- rbind(results_main, data.frame(
      Date     = origin,
      y_actual = y_actual,
      y_hat    = y_hat,
      y_rw     = y_rw,
      n_div    = diag$n_div,
      max_rhat = diag$max_rhat,
      min_ess  = diag$min_ess,
      elapsed  = elapsed,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("done (%.1f min, div=%d, Rhat=%.3f, ESS=%.0f, yhat=%.2f)\n",
                elapsed, diag$n_div, diag$max_rhat, diag$min_ess, y_hat))
    cat(sprintf("[%s] n=%d, %.1fmin, div=%d, Rhat=%.3f, ESS=%.0f, yhat=%.2f\n",
                origin, nrow(train_raw), elapsed, diag$n_div,
                diag$max_rhat, diag$min_ess, y_hat),
        file = LOG_FILE, append = TRUE)

    write.csv(results_main, "walker_forecasts.csv", row.names = FALSE)

  }, error = function(e) {
    msg <- sprintf("[%s] FAILED: %s\n", origin, e$message)
    cat(msg); cat(msg, file = LOG_FILE, append = TRUE)
  })
}

# Build real-time contributions: c_jt = beta_jt^filtered * x_jt^scaled.
# Intercept contribution equals the filtered intercept beta_0t directly.
if (length(filtered_coefs_list) > 0) {
  filtered_coefs_df <- bind_rows(filtered_coefs_list)
  x_values_df       <- bind_rows(x_values_list)

  write.csv(filtered_coefs_df, "walker_filtered_coefs.csv",     row.names = FALSE)
  write.csv(x_values_df,       "walker_x_values_at_origin.csv", row.names = FALSE)

  walker_contributions_main <- filtered_coefs_df %>%
    left_join(x_values_df, by = c("origin", "beta")) %>%
    mutate(contribution = case_when(
      beta == "(Intercept)" ~ mean,
      TRUE                  ~ mean * x_scaled
    )) %>%
    select(origin, beta, contribution) %>%
    pivot_wider(names_from = beta, values_from = contribution) %>%
    rename(Date = origin, intercept = `(Intercept)`)

  write.csv(walker_contributions_main, "walker_contributions.csv", row.names = FALSE)
  cat("Saved: walker_filtered_coefs.csv, walker_x_values_at_origin.csv,",
      "walker_contributions.csv\n")
}

# Evaluation summary for the main spec.
cat("\n=== EVALUATION (main spec) ===\n")
if (nrow(results_main) > 0) {
  results_main <- results_main %>%
    mutate(err_walker = y_hat - y_actual, err_rw = y_rw - y_actual)

  rmse_walker <- sqrt(mean(results_main$err_walker^2, na.rm = TRUE))
  rmse_rw     <- sqrt(mean(results_main$err_rw^2,     na.rm = TRUE))

  cat(sprintf("Walker:  RMSE=%.4f  Bias=%.4f\n",
              rmse_walker, mean(results_main$err_walker, na.rm = TRUE)))
  cat(sprintf("RW:      RMSE=%.4f  Bias=%.4f\n",
              rmse_rw,     mean(results_main$err_rw,     na.rm = TRUE)))
  cat(sprintf("Ratio (Walker/RW): %.4f\n", rmse_walker / rmse_rw))

  # Diebold-Mariano test with Newey-West HAC correction.
  d <- results_main$err_rw^2 - results_main$err_walker^2
  d <- d[is.finite(d)]
  n <- length(d)
  if (n > 5) {
    nw_lag    <- max(H - 1, floor(n^(1/3)))
    d_mean    <- mean(d)
    gamma_0   <- mean((d - d_mean)^2)
    gamma_sum <- 0
    for (k in 1:nw_lag) {
      if ((n - k) <= 0) break
      gamma_k   <- mean((d[(k + 1):n] - d_mean) * (d[1:(n - k)] - d_mean))
      gamma_sum <- gamma_sum + 2 * (1 - k / (nw_lag + 1)) * gamma_k
    }
    var_d    <- (gamma_0 + gamma_sum) / n
    dm_stat  <- d_mean / sqrt(max(var_d, 1e-10))
    dm_pval  <- 2 * (1 - pnorm(abs(dm_stat)))
    cat(sprintf("\nDM test (walker vs RW): stat=%.3f, p=%.4f\n", dm_stat, dm_pval))
    cat("(Positive stat => walker better)\n")
  }

  # Convergence diagnostics summary table
  conv_diag <- data.frame(
    Diagnostic = c(
      "Number of fits",
      "Divergent transitions per fit",
      "Maximum R-hat",
      "Minimum effective sample size",
      "Mean fit time, minutes"
    ),
    Min    = c(nrow(results_main),
               round(min(results_main$n_div,    na.rm = TRUE), 2),
               round(min(results_main$max_rhat, na.rm = TRUE), 2),
               round(min(results_main$min_ess,  na.rm = TRUE), 0),
               round(min(results_main$elapsed,  na.rm = TRUE), 2)),
    Median = c(NA,
               round(median(results_main$n_div,    na.rm = TRUE), 2),
               round(median(results_main$max_rhat, na.rm = TRUE), 2),
               round(median(results_main$min_ess,  na.rm = TRUE), 0),
               round(median(results_main$elapsed,  na.rm = TRUE), 2)),
    Mean   = c(NA,
               round(mean(results_main$n_div,    na.rm = TRUE), 2),
               round(mean(results_main$max_rhat, na.rm = TRUE), 2),
               round(mean(results_main$min_ess,  na.rm = TRUE), 0),
               round(mean(results_main$elapsed,  na.rm = TRUE), 2)),
    Max    = c(NA,
               round(max(results_main$n_div,    na.rm = TRUE), 2),
               round(max(results_main$max_rhat, na.rm = TRUE), 2),
               round(max(results_main$min_ess,  na.rm = TRUE), 0),
               round(max(results_main$elapsed,  na.rm = TRUE), 2)),
    stringsAsFactors = FALSE
  )
  cat("\nConvergence diagnostics summary:\n")
  print(conv_diag, row.names = FALSE, na.print = "--")
  write.csv(conv_diag, "walker_convergence_diagnostics.csv", row.names = FALSE)

  write.csv(results_main, "walker_forecasts.csv", row.names = FALSE)
  cat("Saved: walker_forecasts.csv\n")
}


# =============================================================================
# PART 3 — WALKER BACKTEST: USD EXTERNAL ROBUSTNESS SPECIFICATION
# Same setup as Part 2 but using USD-denominated oil, imports, and exports.
# Requires columns: oljepris_USD, import_USD, eksport_USD in master_data.csv.
# =============================================================================

df_usd <- read.csv("master_data.csv", stringsAsFactors = FALSE)
df_usd$Date <- as.Date(df_usd$Date)
df_usd <- df_usd %>% arrange(Date)

needed_usd <- c("oljepris_USD", "import_USD", "eksport_USD")
missing_usd <- setdiff(needed_usd, names(df_usd))
if (length(missing_usd) > 0) {
  stop(sprintf("USD spec requires columns: %s. Missing from master_data.csv: %s",
               paste(needed_usd, collapse = ", "),
               paste(missing_usd, collapse = ", ")))
}

# Feature engineering with USD-denominated oil and trade.
df_usd <- df_usd %>%
  mutate(
    kpi_yoy_raw     = (kpi          / lag(kpi,          12) - 1) * 100,
    ppi_yoy_raw     = (ppi          / lag(ppi,          12) - 1) * 100,
    oil_yoy_raw     = (oljepris_USD / lag(oljepris_USD, 12) - 1) * 100,
    import_yoy_raw  = (import_USD   / lag(import_USD,   12) - 1) * 100,
    eksport_yoy_raw = (eksport_USD  / lag(eksport_USD,  12) - 1) * 100,

    kpi_yoy_lag1     = lag(kpi_yoy_raw, 1),
    ppi_yoy_lag1     = lag(ppi_yoy_raw, 1),
    oil_yoy_lag1     = lag(oil_yoy_raw, 1),
    usd_nok_lag1     = lag(usd_nok, 1),
    eur_nok_lag1     = lag(eur_nok, 1),
    import_yoy_lag1  = lag(import_yoy_raw, 1),
    eksport_yoy_lag1 = lag(eksport_yoy_raw, 1),
    unemp_lag1       = lag(unemployment, 1),
    rente_lag1       = lag(rente, 1),

    target = lead(kpi_yoy_raw, H)
  )

model_df_usd <- df_usd %>%
  select(Date, target, kpi_yoy_raw, all_of(FEATURE_COLS)) %>%
  filter(complete.cases(.))

test_origins_usd <- model_df_usd %>%
  filter(Date >= TEST_START, !is.na(target)) %>%
  pull(Date)

LOG_FILE_USD <- "walker_usd_external_diagnostics.txt"

rw_formula_usd <- as.formula(
  paste0("target ~ -1 + rw1(~ 1 + ",
         paste(FEATURE_COLS, collapse = " + "),
         ", beta = c(0, 10), sigma = c(2, 0.01))")
)

results_usd          <- data.frame(
  Date = as.Date(character()), y_actual = numeric(), y_hat = numeric(),
  y_rw = numeric(), n_div = integer(), max_rhat = numeric(),
  min_ess = numeric(), elapsed = numeric(), stringsAsFactors = FALSE
)
filtered_coefs_list_usd <- list()
x_values_list_usd       <- list()

cat(sprintf("Walker Robustness Run (usd_external) started: %s\n", Sys.time()),
    file = LOG_FILE_USD)

cat("\n=== EXPANDING WINDOW BACKTEST (USD external spec) ===\n")
cat("Total origins:", length(test_origins_usd), "\n\n")

for (i in seq_along(test_origins_usd)) {
  origin         <- test_origins_usd[i]
  train_raw      <- model_df_usd %>% filter(Date <  origin)
  origin_row_raw <- model_df_usd %>% filter(Date == origin)

  if (nrow(origin_row_raw) != 1) next
  y_actual <- origin_row_raw$target
  y_rw     <- origin_row_raw$kpi_yoy_lag1
  if (is.na(y_actual) || is.na(y_rw) || nrow(train_raw) < MIN_TRAIN_N) next

  cat(sprintf("[%d/%d] %s ... ", i, length(test_origins_usd), origin))

  tryCatch({
    scaled       <- scale_with_train_stats(train_raw, origin_row_raw, FEATURE_COLS)
    train_scaled <- scaled$train_scaled
    new_x_scaled <- scaled$new_scaled[, FEATURE_COLS, drop = FALSE]

    t_start <- Sys.time()
    fit <- walker(
      formula       = rw_formula_usd,
      data          = train_scaled,
      sigma_y_prior = c(2, 0.01),
      chains        = N_CHAINS,
      iter          = N_ITER,
      warmup        = N_WARMUP,
      cores         = N_CORES,
      refresh       = 0,
      control       = list(adapt_delta = 0.95, max_treedepth = 12)
    )
    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

    diag  <- extract_fit_diagnostics(fit)
    y_hat <- predict_last_state(fit, new_x_scaled, FEATURE_COLS)

    filtered_coefs_list_usd[[i]] <- extract_filtered_coefs(fit, origin)
    x_values_list_usd[[i]] <- as.data.frame(new_x_scaled) %>%
      pivot_longer(everything(), names_to = "beta", values_to = "x_scaled") %>%
      mutate(origin = origin)

    results_usd <- rbind(results_usd, data.frame(
      Date = origin, y_actual = y_actual, y_hat = y_hat, y_rw = y_rw,
      n_div = diag$n_div, max_rhat = diag$max_rhat, min_ess = diag$min_ess,
      elapsed = elapsed, stringsAsFactors = FALSE
    ))

    cat(sprintf("done (%.1f min, yhat=%.2f)\n", elapsed, y_hat))
    write.csv(results_usd, "walker_usd_external_forecasts.csv", row.names = FALSE)

  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message))
  })
}

if (length(filtered_coefs_list_usd) > 0) {
  filtered_coefs_df_usd <- bind_rows(filtered_coefs_list_usd)
  x_values_df_usd       <- bind_rows(x_values_list_usd)

  write.csv(filtered_coefs_df_usd,
            "walker_usd_external_filtered_coefs.csv", row.names = FALSE)

  walker_contributions_usd <- filtered_coefs_df_usd %>%
    left_join(x_values_df_usd, by = c("origin", "beta")) %>%
    mutate(contribution = case_when(
      beta == "(Intercept)" ~ mean,
      TRUE                  ~ mean * x_scaled
    )) %>%
    select(origin, beta, contribution) %>%
    pivot_wider(names_from = beta, values_from = contribution) %>%
    rename(Date = origin, intercept = `(Intercept)`)

  write.csv(walker_contributions_usd,
            "walker_usd_external_contributions.csv", row.names = FALSE)
  cat("Saved: walker_usd_external_contributions.csv\n")
}

if (nrow(results_usd) > 0) {
  results_usd <- results_usd %>%
    mutate(err_walker = y_hat - y_actual, err_rw = y_rw - y_actual)
  cat(sprintf("\n=== USD EXTERNAL EVAL ===\nRMSE=%.4f  RW RMSE=%.4f  Ratio=%.4f\n",
              sqrt(mean(results_usd$err_walker^2, na.rm = TRUE)),
              sqrt(mean(results_usd$err_rw^2,     na.rm = TRUE)),
              sqrt(mean(results_usd$err_walker^2, na.rm = TRUE)) /
                sqrt(mean(results_usd$err_rw^2,   na.rm = TRUE))))
  write.csv(results_usd, "walker_usd_external_forecasts.csv", row.names = FALSE)
}


# =============================================================================
# PART 4 — HARMONISATION: WALKER MAIN SPEC
# Aggregates feature-level contributions to economic blocks and computes
# signed/absolute attribution shares. Reads from walker_contributions.csv.
# =============================================================================

walker_contrib_main <- read_csv("walker_contributions.csv",
                                show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date))


# Feature -> block aggregation.
feature_long_main <- walker_contrib_main %>%
  select(Date, all_of(BLOCK_MAP$feature), intercept) %>%
  pivot_longer(cols = all_of(BLOCK_MAP$feature),
               names_to = "feature", values_to = "contribution") %>%
  left_join(BLOCK_MAP, by = "feature")

# Block-level signed contributions (sum of features within each block).
block_signed_main <- feature_long_main %>%
  group_by(Date, block) %>%
  summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = block, values_from = contribution) %>%
  left_join(walker_contrib_main %>% select(Date, intercept), by = "Date") %>%
  rename(baseline_WalkerAdj = intercept) %>%
  mutate(regime = assign_regime(Date)) %>%
  arrange(Date)

# Block-level absolute contributions.
block_abs_main <- feature_long_main %>%
  group_by(Date, block) %>%
  summarise(abs_contribution = sum(abs(contribution), na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = block, values_from = abs_contribution) %>%
  arrange(Date)

# Attribution shares: S_k = C_k / sum_m |C_m|  (signed)
#                     A_k = |C_k| / sum_m |C_m| (absolute)
block_names_main <- setdiff(names(block_signed_main),
                            c("Date", "baseline_WalkerAdj", "regime"))

shares_main <- block_signed_main %>%
  rowwise() %>%
  mutate(total_abs = sum(abs(c_across(all_of(block_names_main))), na.rm = TRUE)) %>%
  ungroup()

for (blk in block_names_main) {
  shares_main[[paste0("S_", blk)]] <- ifelse(
    shares_main$total_abs > 0,
    shares_main[[blk]] / shares_main$total_abs, NA_real_)
  shares_main[[paste0("A_", blk)]] <- ifelse(
    shares_main$total_abs > 0,
    abs(shares_main[[blk]]) / shares_main$total_abs, NA_real_)
}

share_cols_main   <- grep("^(S_|A_)", names(shares_main), value = TRUE)
regime_summary_main <- shares_main %>%
  group_by(regime) %>%
  summarise(across(all_of(share_cols_main), ~mean(.x, na.rm = TRUE)),
            n_obs = n(), .groups = "drop")
global_summary_main <- shares_main %>%
  summarise(across(all_of(share_cols_main), ~mean(.x, na.rm = TRUE)))

# Reconstruction check: baseline + sum of block contributions should equal y_hat.
reconstruction_main <- block_signed_main %>%
  mutate(sum_blocks        = rowSums(across(all_of(block_names_main)), na.rm = TRUE),
         fitted_from_parts = baseline_WalkerAdj + sum_blocks)
if (file.exists("walker_forecasts.csv")) {
  fc_check <- read_csv("walker_forecasts.csv", show_col_types = FALSE) %>%
    mutate(Date = as.Date(Date))
  reconstruction_main <- reconstruction_main %>%
    left_join(fc_check %>% select(Date, y_hat), by = "Date") %>%
    mutate(reconstruction_gap_vs_backtest = y_hat - fitted_from_parts)
}

write_csv(feature_long_main,     "walker_harmonized_feature_long.csv")
write_csv(block_signed_main,     "walker_harmonized_block_signed.csv")
write_csv(block_abs_main,        "walker_harmonized_block_abs.csv")
write_csv(shares_main,           "walker_harmonized_shares.csv")
write_csv(regime_summary_main,   "walker_harmonized_regime_summary.csv")
write_csv(global_summary_main,   "walker_harmonized_global_summary.csv")
write_csv(reconstruction_main,   "walker_harmonized_reconstruction_check.csv")
write_csv(BLOCK_MAP,             "walker_feature_block_map.csv")
cat("Saved: walker_harmonized_*.csv\n")


# =============================================================================
# PART 5 — HARMONISATION: WALKER USD EXTERNAL SPEC
# Identical procedure to Part 4 applied to the USD robustness contributions.
# =============================================================================

walker_contrib_usd <- read_csv("walker_usd_external_contributions.csv",
                               show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date))

if (!"intercept" %in% names(walker_contrib_usd))
  walker_contrib_usd$intercept <- 0

feature_long_usd <- walker_contrib_usd %>%
  select(Date, all_of(BLOCK_MAP$feature), intercept) %>%
  pivot_longer(cols = all_of(BLOCK_MAP$feature),
               names_to = "feature", values_to = "contribution") %>%
  left_join(BLOCK_MAP, by = "feature")

block_signed_usd <- feature_long_usd %>%
  group_by(Date, block) %>%
  summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = block, values_from = contribution) %>%
  left_join(walker_contrib_usd %>% select(Date, intercept), by = "Date") %>%
  rename(baseline_WalkerAdj = intercept) %>%
  mutate(regime = assign_regime(Date)) %>%
  arrange(Date)

block_abs_usd <- feature_long_usd %>%
  group_by(Date, block) %>%
  summarise(abs_contribution = sum(abs(contribution), na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = block, values_from = abs_contribution) %>%
  arrange(Date)

block_names_usd <- setdiff(names(block_signed_usd),
                           c("Date", "baseline_WalkerAdj", "regime"))

shares_usd <- block_signed_usd %>%
  rowwise() %>%
  mutate(total_abs = sum(abs(c_across(all_of(block_names_usd))), na.rm = TRUE)) %>%
  ungroup()

for (blk in block_names_usd) {
  shares_usd[[paste0("S_", blk)]] <- ifelse(
    shares_usd$total_abs > 0,
    shares_usd[[blk]] / shares_usd$total_abs, NA_real_)
  shares_usd[[paste0("A_", blk)]] <- ifelse(
    shares_usd$total_abs > 0,
    abs(shares_usd[[blk]]) / shares_usd$total_abs, NA_real_)
}

share_cols_usd <- grep("^(S_|A_)", names(shares_usd), value = TRUE)
regime_summary_usd <- shares_usd %>%
  group_by(regime) %>%
  summarise(across(all_of(share_cols_usd), ~mean(.x, na.rm = TRUE)),
            n_obs = n(), .groups = "drop")
global_summary_usd <- shares_usd %>%
  summarise(across(all_of(share_cols_usd), ~mean(.x, na.rm = TRUE)))

write_csv(feature_long_usd,   "walker_usd_external_harmonized_feature_long.csv")
write_csv(block_signed_usd,   "walker_usd_external_harmonized_block_signed.csv")
write_csv(block_abs_usd,      "walker_usd_external_harmonized_block_abs.csv")
write_csv(shares_usd,         "walker_usd_external_harmonized_shares.csv")
write_csv(regime_summary_usd, "walker_usd_external_harmonized_regime_summary.csv")
write_csv(global_summary_usd, "walker_usd_external_harmonized_global_summary.csv")
cat("Saved: walker_usd_external_harmonized_*.csv\n")


# =============================================================================
# PART 6 — HARMONISATION: XGBOOST SHAP OUTPUTS
# Reads Python-produced SHAP files and renames block columns to the short
# Walker names for a consistent schema across both models.
# =============================================================================

# Helper: rename XGBoost long block names to short names.
xgb_rename_map <- c(
  "AR (inflation)"  = "AR",
  "PPI / cost-push" = "PPI",
  "Monetary policy" = "Monetary",
  "FX"              = "FX",
  "Oil"             = "Oil",
  "Trade"           = "Trade",
  "Labour market"   = "Labour"
)

# --- Main spec XGB ---
stopifnot(file.exists("xgb_v8_harmonized_shap_group_signed.csv"))
xgb_signed_main_raw <- read_csv("xgb_v8_harmonized_shap_group_signed.csv",
                                show_col_types = FALSE)
date_col_xgb <- intersect(c("Date", "date"), names(xgb_signed_main_raw))[[1]]

xgb_block_signed_main <- xgb_signed_main_raw %>%
  mutate(Date = as.Date(.data[[date_col_xgb]])) %>%
  select(Date, all_of(names(xgb_rename_map)), shap_base,
         actual, predicted_raw) %>%
  rename(!!!setNames(names(xgb_rename_map), xgb_rename_map)) %>%
  rename(baseline_SHAP = shap_base) %>%
  mutate(regime = assign_regime(Date)) %>%
  arrange(Date)

xgb_block_abs_main <- xgb_block_signed_main %>%
  mutate(across(all_of(BLOCK_NAMES), abs)) %>%
  select(Date, all_of(BLOCK_NAMES), regime)

shares_xgb_main <- xgb_block_signed_main %>%
  select(Date, regime, all_of(BLOCK_NAMES)) %>%
  rowwise() %>%
  mutate(total_abs = sum(abs(c_across(all_of(BLOCK_NAMES))), na.rm = TRUE)) %>%
  ungroup()

for (blk in BLOCK_NAMES) {
  shares_xgb_main[[paste0("S_", blk)]] <- ifelse(
    shares_xgb_main$total_abs > 0,
    shares_xgb_main[[blk]] / shares_xgb_main$total_abs, NA_real_)
  shares_xgb_main[[paste0("A_", blk)]] <- ifelse(
    shares_xgb_main$total_abs > 0,
    abs(shares_xgb_main[[blk]]) / shares_xgb_main$total_abs, NA_real_)
}

share_cols_xgb <- grep("^(S_|A_)", names(shares_xgb_main), value = TRUE)
regime_summary_xgb_main <- shares_xgb_main %>%
  group_by(regime) %>%
  summarise(across(all_of(share_cols_xgb), ~mean(.x, na.rm = TRUE)),
            n_obs = n(), .groups = "drop")

write_csv(xgb_block_signed_main,   "xgb_main_harmonized_block_signed.csv")
write_csv(xgb_block_abs_main,      "xgb_main_harmonized_block_abs.csv")
write_csv(shares_xgb_main,         "xgb_main_harmonized_shares.csv")
write_csv(regime_summary_xgb_main, "xgb_main_harmonized_regime_summary.csv")

# --- USD robustness XGB ---
stopifnot(file.exists("xgb_harmonized_robustness_shap_group_signed.csv"))
xgb_signed_usd_raw <- read_csv("xgb_harmonized_robustness_shap_group_signed.csv",
                               show_col_types = FALSE)
date_col_xgb_usd <- intersect(c("Date", "date"), names(xgb_signed_usd_raw))[[1]]

xgb_block_signed_usd <- xgb_signed_usd_raw %>%
  mutate(Date = as.Date(.data[[date_col_xgb_usd]])) %>%
  select(Date, all_of(names(xgb_rename_map)), shap_base,
         actual, predicted_raw) %>%
  rename(!!!setNames(names(xgb_rename_map), xgb_rename_map)) %>%
  rename(baseline_SHAP = shap_base) %>%
  mutate(regime = assign_regime(Date)) %>%
  arrange(Date)

xgb_block_abs_usd <- xgb_block_signed_usd %>%
  mutate(across(all_of(BLOCK_NAMES), abs)) %>%
  select(Date, all_of(BLOCK_NAMES), regime)

shares_xgb_usd <- xgb_block_signed_usd %>%
  select(Date, regime, all_of(BLOCK_NAMES)) %>%
  rowwise() %>%
  mutate(total_abs = sum(abs(c_across(all_of(BLOCK_NAMES))), na.rm = TRUE)) %>%
  ungroup()

for (blk in BLOCK_NAMES) {
  shares_xgb_usd[[paste0("S_", blk)]] <- ifelse(
    shares_xgb_usd$total_abs > 0,
    shares_xgb_usd[[blk]] / shares_xgb_usd$total_abs, NA_real_)
  shares_xgb_usd[[paste0("A_", blk)]] <- ifelse(
    shares_xgb_usd$total_abs > 0,
    abs(shares_xgb_usd[[blk]]) / shares_xgb_usd$total_abs, NA_real_)
}

regime_summary_xgb_usd <- shares_xgb_usd %>%
  group_by(regime) %>%
  summarise(across(all_of(share_cols_xgb), ~mean(.x, na.rm = TRUE)),
            n_obs = n(), .groups = "drop")

write_csv(xgb_block_signed_usd,   "xgb_robustness_harmonized_block_signed.csv")
write_csv(xgb_block_abs_usd,      "xgb_robustness_harmonized_block_abs.csv")
write_csv(shares_xgb_usd,         "xgb_robustness_harmonized_shares.csv")
write_csv(regime_summary_xgb_usd, "xgb_robustness_harmonized_regime_summary.csv")
cat("Saved: xgb_main_harmonized_*.csv, xgb_robustness_harmonized_*.csv\n")


# =============================================================================
# PART 7 — SHAP vs WALKER COMPARISON: JOINT ATTRIBUTION TABLE
# Merges XGBoost SHAP and Walker contributions into one harmonised table.
# Computes sign agreement, Pearson correlations, and attribution shares.
# =============================================================================

# Build joint table for main spec.
joint_main <- xgb_block_signed_main %>%
  transmute(
    Date,
    actual        = actual,
    predicted_raw = predicted_raw,
    baseline_SHAP = baseline_SHAP,
    AR_SHAP        = AR,
    PPI_SHAP       = PPI,
    Monetary_SHAP  = Monetary,
    FX_SHAP        = FX,
    Oil_SHAP       = Oil,
    Trade_SHAP     = Trade,
    Labour_SHAP    = Labour
  ) %>%
  left_join(
    block_signed_main %>%
      transmute(Date,
                baseline_WalkerAdj = baseline_WalkerAdj,
                AR_WalkerAdj        = AR,
                PPI_WalkerAdj       = PPI,
                Monetary_WalkerAdj  = Monetary,
                FX_WalkerAdj        = FX,
                Oil_WalkerAdj       = Oil,
                Trade_WalkerAdj     = Trade,
                Labour_WalkerAdj    = Labour),
    by = "Date"
  ) %>%
  mutate(regime = assign_regime(Date)) %>%
  arrange(Date)

# Build joint table for USD spec.
joint_usd <- xgb_block_signed_usd %>%
  transmute(
    Date,
    actual        = actual,
    predicted_raw = predicted_raw,
    baseline_SHAP = baseline_SHAP,
    AR_SHAP        = AR,
    PPI_SHAP       = PPI,
    Monetary_SHAP  = Monetary,
    FX_SHAP        = FX,
    Oil_SHAP       = Oil,
    Trade_SHAP     = Trade,
    Labour_SHAP    = Labour
  ) %>%
  left_join(
    block_signed_usd %>%
      transmute(Date,
                baseline_WalkerAdj = baseline_WalkerAdj,
                AR_WalkerAdj        = AR,
                PPI_WalkerAdj       = PPI,
                Monetary_WalkerAdj  = Monetary,
                FX_WalkerAdj        = FX,
                Oil_WalkerAdj       = Oil,
                Trade_WalkerAdj     = Trade,
                Labour_WalkerAdj    = Labour),
    by = "Date"
  ) %>%
  mutate(regime = assign_regime(Date)) %>%
  arrange(Date)

write_csv(joint_main, "comparison_harmonised_attributions.csv")

# Sign agreement and correlation by block — main spec.
shap_cols   <- paste0(BLOCK_NAMES, "_SHAP")
walker_cols <- paste0(BLOCK_NAMES, "_WalkerAdj")

routing_compare_main <- data.frame(
  block       = BLOCK_NAMES,
  pct_agree   = sapply(BLOCK_NAMES, function(b)
    mean(sign(joint_main[[paste0(b, "_SHAP")]]) ==
           sign(joint_main[[paste0(b, "_WalkerAdj")]]), na.rm = TRUE)),
  corr_signed = sapply(BLOCK_NAMES, function(b)
    cor(joint_main[[paste0(b, "_SHAP")]],
        joint_main[[paste0(b, "_WalkerAdj")]], use = "complete.obs")),
  corr_abs    = sapply(BLOCK_NAMES, function(b)
    cor(abs(joint_main[[paste0(b, "_SHAP")]]),
        abs(joint_main[[paste0(b, "_WalkerAdj")]]), use = "complete.obs"))
)

routing_compare_usd <- data.frame(
  block       = BLOCK_NAMES,
  pct_agree   = sapply(BLOCK_NAMES, function(b)
    mean(sign(joint_usd[[paste0(b, "_SHAP")]]) ==
           sign(joint_usd[[paste0(b, "_WalkerAdj")]]), na.rm = TRUE)),
  corr_signed = sapply(BLOCK_NAMES, function(b)
    cor(joint_usd[[paste0(b, "_SHAP")]],
        joint_usd[[paste0(b, "_WalkerAdj")]], use = "complete.obs")),
  corr_abs    = sapply(BLOCK_NAMES, function(b)
    cor(abs(joint_usd[[paste0(b, "_SHAP")]]),
        abs(joint_usd[[paste0(b, "_WalkerAdj")]]), use = "complete.obs"))
)

write_csv(routing_compare_main, "comparison_routing_compare.csv")
write_csv(routing_compare_usd,  "comparison_routing_compare_usd.csv")
cat("\nSign agreement and correlations (main spec):\n")
print(routing_compare_main %>% mutate(across(where(is.numeric), ~round(., 3))))
cat("\nSign agreement and correlations (USD robustness spec):\n")
print(routing_compare_usd  %>% mutate(across(where(is.numeric), ~round(., 3))))

# Mean signed and absolute attribution shares — both specs.
compute_mean_shares <- function(joint) {
  abs_s <- abs(as.matrix(joint[, shap_cols]))
  abs_w <- abs(as.matrix(joint[, walker_cols]))
  tot_s <- rowSums(abs_s, na.rm = TRUE)
  tot_w <- rowSums(abs_w, na.rm = TRUE)
  out <- data.frame(block = BLOCK_NAMES)
  for (b in BLOCK_NAMES) {
    s <- joint[[paste0(b, "_SHAP")]]; w <- joint[[paste0(b, "_WalkerAdj")]]
    out[out$block == b, "S_SHAP"]   <- mean(ifelse(tot_s > 0, s / tot_s, NA),     na.rm = TRUE)
    out[out$block == b, "A_SHAP"]   <- mean(ifelse(tot_s > 0, abs(s) / tot_s, NA), na.rm = TRUE)
    out[out$block == b, "S_Walker"] <- mean(ifelse(tot_w > 0, w / tot_w, NA),     na.rm = TRUE)
    out[out$block == b, "A_Walker"] <- mean(ifelse(tot_w > 0, abs(w) / tot_w, NA), na.rm = TRUE)
  }
  out
}

mean_shares_main <- compute_mean_shares(joint_main)
mean_shares_usd  <- compute_mean_shares(joint_usd)
write_csv(mean_shares_main, "comparison_mean_shares_main.csv")
write_csv(mean_shares_usd,  "comparison_mean_shares_usd.csv")
cat("\nMean attribution shares (main spec):\n")
print(mean_shares_main %>% mutate(across(where(is.numeric), ~round(., 3))))
cat("\nMean attribution shares (USD robustness spec):\n")
print(mean_shares_usd  %>% mutate(across(where(is.numeric), ~round(., 3))))

# Block importance by regime (mean absolute share × 100) — both specs.
compute_block_importance <- function(joint) {
  abs_s <- abs(as.matrix(joint[, shap_cols]))
  abs_w <- abs(as.matrix(joint[, walker_cols]))
  tot_s <- rowSums(abs_s, na.rm = TRUE)
  tot_w <- rowSums(abs_w, na.rm = TRUE)
  per_row <- data.frame(regime = joint$regime)
  for (b in BLOCK_NAMES) {
    per_row[[paste0("A_", b, "_SHAP")]]   <-
      ifelse(tot_s > 0, abs(joint[[paste0(b, "_SHAP")]]) / tot_s, NA)
    per_row[[paste0("A_", b, "_Walker")]] <-
      ifelse(tot_w > 0, abs(joint[[paste0(b, "_WalkerAdj")]]) / tot_w, NA)
  }
  per_row %>%
    group_by(regime) %>%
    summarise(across(starts_with("A_"), ~mean(.x, na.rm = TRUE) * 100),
              .groups = "drop") %>%
    pivot_longer(starts_with("A_"),
                 names_to      = c("block", "model"),
                 names_pattern = "^A_([^_]+)_(SHAP|Walker)$",
                 values_to     = "value") %>%
    mutate(block = factor(block, levels = BLOCK_NAMES),
           model = factor(model, levels = c("SHAP", "Walker"),
                          labels = c("XGBoost (SHAP)", "Walker"))) %>%
    arrange(model, block)
}

block_importance_main <- compute_block_importance(joint_main)
block_importance_usd  <- compute_block_importance(joint_usd)
write_csv(block_importance_main, "comparison_block_importance_main.csv")
write_csv(block_importance_usd,  "comparison_block_importance_usd.csv")

cat("\nBlock importance by regime (main spec):\n")
print(block_importance_main %>%
        pivot_wider(names_from = regime, values_from = value) %>%
        mutate(across(where(is.numeric), ~round(.x, 1))))


# =============================================================================
# PART 8 — ATTRIBUTION RESIDUAL DECOMPOSITION
# Delta_kt = C_kt^XGB - C_kt^TVP
# Computes regime-conditional variance and mean absolute magnitude per block.
# =============================================================================

residual_df_main <- joint_main %>%
  mutate(
    Delta_AR       = AR_SHAP       - AR_WalkerAdj,
    Delta_PPI      = PPI_SHAP      - PPI_WalkerAdj,
    Delta_Monetary = Monetary_SHAP - Monetary_WalkerAdj,
    Delta_FX       = FX_SHAP       - FX_WalkerAdj,
    Delta_Oil      = Oil_SHAP      - Oil_WalkerAdj,
    Delta_Trade    = Trade_SHAP    - Trade_WalkerAdj,
    Delta_Labour   = Labour_SHAP   - Labour_WalkerAdj
  ) %>%
  select(Date, regime, starts_with("Delta_"))

residual_long_main <- residual_df_main %>%
  pivot_longer(
    cols         = starts_with("Delta_"),
    names_to     = "block",
    names_prefix = "Delta_",
    values_to    = "Delta"
  ) %>%
  mutate(block = factor(block, levels = BLOCK_NAMES))

residual_metrics_main <- residual_long_main %>%
  group_by(block, regime) %>%
  summarise(
    var_delta = var(Delta, na.rm = TRUE),
    mean_abs  = mean(abs(Delta), na.rm = TRUE),
    n_obs     = sum(!is.na(Delta)),
    .groups   = "drop"
  ) %>%
  mutate(
    regime = factor(regime,
                    levels = c("COVID (2020–2021)", "Energy Crisis",
                               "Disinflation", "Normalization"),
                    labels = c("COVID", "Energy", "Disinfl.", "Normal"))
  ) %>%
  arrange(block, regime)

write_csv(residual_long_main,    "residual_long.csv")
write_csv(residual_metrics_main, "residual_metrics_long.csv")

cat("\nResidual metrics (Var and E|Delta| by block and regime, main spec):\n")
print(residual_metrics_main %>% mutate(across(where(is.numeric), ~round(., 3))))

# USD robustness residual decomposition — identical procedure using joint_usd.
residual_df_usd <- joint_usd %>%
  mutate(
    Delta_AR       = AR_SHAP       - AR_WalkerAdj,
    Delta_PPI      = PPI_SHAP      - PPI_WalkerAdj,
    Delta_Monetary = Monetary_SHAP - Monetary_WalkerAdj,
    Delta_FX       = FX_SHAP       - FX_WalkerAdj,
    Delta_Oil      = Oil_SHAP      - Oil_WalkerAdj,
    Delta_Trade    = Trade_SHAP    - Trade_WalkerAdj,
    Delta_Labour   = Labour_SHAP   - Labour_WalkerAdj
  ) %>%
  select(Date, regime, starts_with("Delta_"))

residual_long_usd <- residual_df_usd %>%
  pivot_longer(
    cols         = starts_with("Delta_"),
    names_to     = "block",
    names_prefix = "Delta_",
    values_to    = "Delta"
  ) %>%
  mutate(block = factor(block, levels = BLOCK_NAMES))

residual_metrics_usd <- residual_long_usd %>%
  group_by(block, regime) %>%
  summarise(
    var_delta = var(Delta, na.rm = TRUE),
    mean_abs  = mean(abs(Delta), na.rm = TRUE),
    n_obs     = sum(!is.na(Delta)),
    .groups   = "drop"
  ) %>%
  mutate(
    regime = factor(regime,
                    levels = c("COVID (2020–2021)", "Energy Crisis",
                               "Disinflation", "Normalization"),
                    labels = c("COVID", "Energy", "Disinfl.", "Normal"))
  ) %>%
  arrange(block, regime)

write_csv(residual_long_usd,    "residual_long_usd.csv")
write_csv(residual_metrics_usd, "residual_metrics_long_usd.csv")

cat("\nResidual metrics (Var and E|Delta| by block and regime, USD robustness spec):\n")
print(n = 28, residual_metrics_usd %>% mutate(across(where(is.numeric), ~round(., 3))))


# =============================================================================
# PART 9 — AR(1) BASELINE AND ACCURACY TABLE
# Expanding-window AR(1) forecast: target ~ kpi_yoy_lag1.
# Combined accuracy table: RW, AR(1), Walker (main), Walker (USD), XGBoost.
# =============================================================================

# AR(1) expanding-window backtest.
model_df_ar1 <- master_df %>%
  mutate(
    kpi_yoy_raw  = (kpi / lag(kpi, 12) - 1) * 100,
    kpi_yoy_lag1 = lag(kpi_yoy_raw, 1),
    target       = lead(kpi_yoy_raw, H)
  ) %>%
  filter(!is.na(target), !is.na(kpi_yoy_lag1))

origins_ar1 <- model_df_ar1 %>% filter(Date >= TEST_START) %>% pull(Date)
ar1_rows    <- vector("list", length(origins_ar1))

for (i in seq_along(origins_ar1)) {
  origin     <- origins_ar1[i]
  train_ar1  <- model_df_ar1 %>% filter(Date <  origin)
  origin_ar1 <- model_df_ar1 %>% filter(Date == origin)
  if (nrow(train_ar1) < MIN_TRAIN_N || nrow(origin_ar1) != 1) next
  fit_ar1    <- lm(target ~ kpi_yoy_lag1, data = train_ar1)
  ar1_rows[[i]] <- data.frame(
    Date     = origin,
    y_actual = origin_ar1$target,
    y_ar1    = as.numeric(predict(fit_ar1, newdata = origin_ar1))
  )
}
ar1_results <- do.call(rbind, ar1_rows)
ar1_results$Date <- as.Date(ar1_results$Date, origin = "1970-01-01")
write_csv(ar1_results, "ar1_forecasts.csv")

# Accuracy table.
fc_main_acc <- read_csv("walker_forecasts.csv",              show_col_types = FALSE)
fc_usd_acc  <- read_csv("walker_usd_external_forecasts.csv", show_col_types = FALSE)

accuracy_table <- bind_rows(
  forecast_metrics(fc_main_acc$y_actual, fc_main_acc$y_rw,  "Random walk"),
  forecast_metrics(fc_main_acc$y_actual, fc_main_acc$y_hat, "Walker (main)"),
  forecast_metrics(fc_usd_acc$y_actual,  fc_usd_acc$y_hat,  "Walker (USD external)"),
  forecast_metrics(ar1_results$y_actual, ar1_results$y_ar1, "AR(1)")
)

# Add XGBoost if predictions file is available.
if (file.exists("xgb_v8_harmonized_predictions.csv")) {
  xgb_pred <- read_csv("xgb_v8_harmonized_predictions.csv",
                       show_col_types = FALSE)
  dc_xgb <- intersect(c("Date", "date"), names(xgb_pred))[[1]]
  xgb_pred$Date <- as.Date(xgb_pred[[dc_xgb]])
  combined_xgb <- fc_main_acc %>%
    mutate(Date = as.Date(Date)) %>%
    inner_join(xgb_pred %>% select(Date, predicted_raw), by = "Date")
  accuracy_table <- bind_rows(
    accuracy_table,
    forecast_metrics(combined_xgb$y_actual, combined_xgb$predicted_raw, "XGBoost")
  )
}

write_csv(accuracy_table, "accuracy_table.csv")
cat("\nForecast accuracy:\n")
print(accuracy_table)

# =============================================================================
# DIEBOLD-MARIANO TESTS OF EQUAL PREDICTIVE ACCURACY
# H0: equal expected squared-error loss (two-sided).
# Newey-West variance correction at lag h-1 = H-1.
# Requires the forecast package: install.packages("forecast")
# =============================================================================

library(forecast)

# Forecast errors on the Walker evaluation sample (actual - predicted).
e_walker <- fc_main_acc$y_actual - fc_main_acc$y_hat
e_rw     <- fc_main_acc$y_actual - fc_main_acc$y_rw

# AR(1) errors — align dates with Walker sample.
ar1_aligned <- ar1_results %>%
  inner_join(fc_main_acc %>% mutate(Date = as.Date(Date)) %>% select(Date),
             by = "Date")
e_ar1 <- ar1_aligned$y_actual - ar1_aligned$y_ar1

run_dm <- function(e1, e2, name1, name2) {
  test  <- dm.test(e1, e2, alternative = "two.sided", h = H, power = 2)
  stars <- ifelse(test$p.value < 0.01, "***",
           ifelse(test$p.value < 0.05, "**",
           ifelse(test$p.value < 0.10, "*", "")))
  data.frame(
    Model1      = name1,
    Model2      = name2,
    RMSE1       = round(sqrt(mean(e1^2, na.rm = TRUE)), 3),
    RMSE2       = round(sqrt(mean(e2^2, na.rm = TRUE)), 3),
    RMSE_ratio  = round(sqrt(mean(e1^2, na.rm = TRUE)) /
                        sqrt(mean(e2^2, na.rm = TRUE)), 3),
    DM_stat     = round(as.numeric(test$statistic), 2),
    p_value     = round(test$p.value, 3),
    sig         = stars,
    stringsAsFactors = FALSE
  )
}

dm_results <- bind_rows(
  run_dm(e_walker, e_rw,  "Walker (main)", "Random walk"),
  run_dm(e_walker, e_ar1, "Walker (main)", "AR(1)")
)

# Add XGBoost comparisons if predictions file was loaded.
if (exists("combined_xgb")) {
  e_xgb_aligned <- combined_xgb$y_actual - combined_xgb$predicted_raw
  e_rw_xgb      <- combined_xgb$y_rw - combined_xgb$y_actual  # align length
  e_rw_xgb      <- combined_xgb$y_actual - combined_xgb$y_rw
  e_walker_xgb  <- combined_xgb$y_actual - combined_xgb$y_hat
  dm_results <- bind_rows(
    dm_results,
    run_dm(e_xgb_aligned, e_rw_xgb,     "XGBoost (main)", "Random walk"),
    run_dm(e_xgb_aligned, e_walker_xgb, "XGBoost (main)", "Walker (main)")
  )
}

cat("\nDiebold-Mariano tests (two-sided, squared error loss, h =", H, "):\n")
print(dm_results, row.names = FALSE)
write_csv(dm_results, "dm_test_results.csv")


# =============================================================================
# PART 10 — DESCRIPTIVE FIGURES
# All figures read from master_data.csv (no API calls required).
# =============================================================================

df_master <- read_csv("master_data.csv", show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date)) %>%
  arrange(Date)

theme_thesis <- theme_minimal(base_size = 11, base_family = "serif") +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, color = "grey40",
                                   margin = margin(b = 10)),
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

# Helper: find the first date with non-NA values for a variable.
first_obs_date <- function(data, var) {
  min(data$Date[!is.na(data[[var]])], na.rm = TRUE)
}

# CPI
p_cpi <- df_master %>%
  filter(Date >= first_obs_date(df_master, "kpi")) %>%
  ggplot(aes(x = Date, y = kpi)) +
  geom_line() +
  labs(title = "Consumer Price Index (index 2015 = 100)",
       x = "Date", y = "CPI") +
  theme_thesis
ggsave("cpi.png", p_cpi, width = 10, height = 5, dpi = 200)

# PPI
p_ppi <- df_master %>%
  filter(Date >= first_obs_date(df_master, "ppi")) %>%
  ggplot(aes(x = Date, y = ppi)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(decimal.mark = ","),
                     breaks = pretty_breaks(n = 6)) +
  labs(title = "Producer Price Index (2021 = 100)", x = "Date", y = "PPI") +
  theme_thesis
ggsave("PPI.png", p_ppi, width = 10, height = 5, dpi = 200)

# Unemployment
p_unemp <- df_master %>%
  filter(Date >= first_obs_date(df_master, "unemployment")) %>%
  ggplot(aes(x = Date, y = unemployment)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(suffix = " %", decimal.mark = ","),
                     breaks = pretty_breaks(n = 6)) +
  labs(title = "Unemployment as a share of the labour force",
       x = "Date", y = "Unemployment (%)") +
  theme_thesis
ggsave("unemployment.png", p_unemp, width = 10, height = 5, dpi = 200)

# Oil price in NOK
p_oil_nok <- df_master %>%
  filter(Date >= first_obs_date(df_master, "oil_price_nok")) %>%
  ggplot(aes(x = Date, y = oil_price_nok)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(suffix = " kr", decimal.mark = ","),
                     breaks = pretty_breaks(n = 6)) +
  labs(title = "Oil price (NOK per barrel)", x = "Date", y = "Price per barrel") +
  theme_thesis
ggsave("oljeplot.png", p_oil_nok, width = 10, height = 5, dpi = 200)

# Oil price in USD
p_oil_usd <- df_master %>%
  filter(Date >= first_obs_date(df_master, "oljepris_USD")) %>%
  ggplot(aes(x = Date, y = oljepris_USD)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(suffix = " USD", decimal.mark = ","),
                     breaks = pretty_breaks(n = 6)) +
  labs(title = "Oil price (USD per barrel)", x = "Date", y = "USD per barrel") +
  theme_thesis
ggsave("olje_usd.png", p_oil_usd, width = 10, height = 5, dpi = 200)

# USD/NOK
p_usd <- df_master %>%
  filter(Date >= first_obs_date(df_master, "usd_nok")) %>%
  ggplot(aes(x = Date, y = usd_nok)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(decimal.mark = ","),
                     breaks = pretty_breaks(n = 6)) +
  labs(title = "USD/NOK exchange rate", x = "Date", y = "NOK per USD") +
  theme_thesis
ggsave("usd_nok.png", p_usd, width = 10, height = 5, dpi = 200)

# EUR/NOK
p_eur <- df_master %>%
  filter(Date >= first_obs_date(df_master, "eur_nok")) %>%
  ggplot(aes(x = Date, y = eur_nok)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(decimal.mark = ","),
                     breaks = pretty_breaks(n = 6)) +
  labs(title = "EUR/NOK exchange rate", x = "Date", y = "NOK per EUR") +
  theme_thesis
ggsave("eur_nok.png", p_eur, width = 10, height = 5, dpi = 200)

# Policy rate
p_rente <- df_master %>%
  filter(Date >= first_obs_date(df_master, "rente")) %>%
  ggplot(aes(x = Date, y = rente)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(suffix = " %", decimal.mark = ","),
                     breaks = pretty_breaks(n = 6)) +
  labs(title = "Policy rate over time", x = "Date", y = "Policy rate (%)") +
  theme_thesis
ggsave("policy_rate.png", p_rente, width = 10, height = 5, dpi = 200)

# Trade (NOK)
p_trade <- df_master %>%
  filter(Date >= min(Date[!is.na(df_master$import) | !is.na(df_master$eksport)],
                     na.rm = TRUE)) %>%
  select(Date, import, eksport) %>%
  pivot_longer(c(import, eksport), names_to = "series", values_to = "value") %>%
  mutate(series = recode(series, import = "Import", eksport = "Export")) %>%
  ggplot(aes(x = Date, y = value, linetype = series)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(scale = 1e-9, suffix = " bn NOK",
                                           decimal.mark = ","),
                     breaks = pretty_breaks(n = 6)) +
  labs(title = "Norwegian foreign trade", x = "Date",
       y = "Billion NOK", linetype = NULL) +
  theme_thesis
ggsave("trade.png", p_trade, width = 10, height = 5, dpi = 200)

# Trade (USD)
p_trade_usd <- df_master %>%
  filter(Date >= min(Date[!is.na(df_master$import_USD) |
                            !is.na(df_master$eksport_USD)], na.rm = TRUE)) %>%
  select(Date, import_USD, eksport_USD) %>%
  pivot_longer(c(import_USD, eksport_USD),
               names_to = "series", values_to = "value") %>%
  mutate(series = recode(series,
                         import_USD  = "Import",
                         eksport_USD = "Export")) %>%
  ggplot(aes(x = Date, y = value, linetype = series)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(labels = label_number(scale = 1e-9, suffix = " bn USD",
                                           decimal.mark = ","),
                     breaks = pretty_breaks(n = 6)) +
  labs(title = "Norwegian foreign trade in USD", x = "Date",
       y = "Billion USD", linetype = NULL) +
  theme_thesis
ggsave("trade_usd.png", p_trade_usd, width = 10, height = 5, dpi = 200)

cat("Saved: descriptive figures (cpi, PPI, unemployment, oil, fx, trade)\n")


# =============================================================================
# PART 11 — MODEL OUTPUT FIGURES
# =============================================================================

# --- Walker forecast vs actual vs random walk ---
fc_plot <- read_csv("walker_forecasts.csv", show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date))

p_forecast <- fc_plot %>%
  pivot_longer(c(y_actual, y_hat, y_rw),
               names_to = "series", values_to = "value") %>%
  mutate(series = factor(series,
                         levels = c("y_actual", "y_hat", "y_rw"),
                         labels = c("Actual", "Walker", "Random walk"))) %>%
  ggplot(aes(x = Date, y = value, color = series, linetype = series)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("Actual"      = "black",
                                "Walker"      = "steelblue",
                                "Random walk" = "grey50")) +
  scale_linetype_manual(values = c("Actual"      = "solid",
                                   "Walker"      = "solid",
                                   "Random walk" = "dashed")) +
  labs(title    = sprintf("Walker vs random walk (h = %d months)", H),
       subtitle = "Out-of-sample backtesting period",
       y        = "YoY CPI inflation (%)",
       x        = NULL, color = NULL, linetype = NULL) +
  theme_thesis

ggsave("walker_forecast_plot.png", p_forecast, width = 10, height = 5, dpi = 200)

# --- Forecast errors ---
p_error <- fc_plot %>%
  mutate(err_walker = y_hat - y_actual, err_rw = y_rw - y_actual) %>%
  pivot_longer(c(err_walker, err_rw),
               names_to = "series", values_to = "error") %>%
  mutate(series = recode(series,
                         err_walker = "Walker",
                         err_rw     = "Random walk")) %>%
  ggplot(aes(x = Date, y = error, color = series, linetype = series)) +
  geom_hline(yintercept = 0, color = "grey40", linetype = "dashed") +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("Walker" = "steelblue", "Random walk" = "grey50")) +
  scale_linetype_manual(values = c("Walker" = "solid", "Random walk" = "dashed")) +
  labs(title = sprintf("Forecast errors (h = %d months)", H),
       y = "Error (pp)", x = NULL, color = NULL, linetype = NULL) +
  theme_thesis

ggsave("walker_error_plot.png", p_error, width = 10, height = 5, dpi = 200)

# --- Walker vs XGBoost ---
if (file.exists("xgb_v8_harmonized_predictions.csv")) {
  xgb_pred_plot <- read_csv("xgb_v8_harmonized_predictions.csv",
                            show_col_types = FALSE)
  dc <- intersect(c("Date", "date"), names(xgb_pred_plot))[[1]]
  xgb_pred_plot$Date <- as.Date(xgb_pred_plot[[dc]])

  panel_wx <- fc_plot %>%
    select(Date, Actual = y_actual, Walker = y_hat) %>%
    inner_join(xgb_pred_plot %>% select(Date, XGBoost = predicted_raw),
               by = "Date")

  rmse_w <- sqrt(mean((panel_wx$Actual - panel_wx$Walker)^2))
  rmse_x <- sqrt(mean((panel_wx$Actual - panel_wx$XGBoost)^2))

  panel_long <- panel_wx %>%
    pivot_longer(c(Actual, Walker, XGBoost),
                 names_to = "series", values_to = "value") %>%
    mutate(series = factor(series,
                           levels = c("Actual", "Walker", "XGBoost"),
                           labels = c("Actual",
                                      sprintf("Walker (RMSE = %.3f)", rmse_w),
                                      sprintf("XGBoost (RMSE = %.3f)", rmse_x))))

  p_wx <- ggplot(panel_long,
                 aes(x = Date, y = value, color = series, linetype = series)) +
    geom_line(linewidth = 0.6) +
    scale_color_manual(values = setNames(
      c("black", "steelblue", "#D55E00"), levels(panel_long$series))) +
    scale_linetype_manual(values = setNames(
      c("solid", "solid", "dashed"), levels(panel_long$series))) +
    labs(title    = sprintf("Walker vs XGBoost — h = %d inflation forecast", H),
         subtitle = sprintf("Out-of-sample period: %s – %s",
                            format(min(panel_wx$Date), "%b %Y"),
                            format(max(panel_wx$Date), "%b %Y")),
         x = NULL, y = "YoY CPI inflation (%)", color = NULL, linetype = NULL) +
    theme_thesis

  ggsave("walker_vs_xgboost.png", p_wx, width = 10, height = 5, dpi = 200)
}

# --- Block contributions: SHAP vs Walker ---
comparison_long <- joint_main %>%
  select(Date, ends_with("_SHAP"), ends_with("_WalkerAdj")) %>%
  pivot_longer(-Date, names_to = "series", values_to = "value") %>%
  mutate(
    model     = if_else(grepl("_SHAP$", series), "SHAP", "Walker adjusted"),
    component = sub("_(SHAP|WalkerAdj)$", "", series)
  ) %>%
  filter(component %in% BLOCK_NAMES)

p_contrib <- ggplot(comparison_long,
                    aes(x = Date, y = value, linetype = model)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.45) +
  facet_wrap(~ component, scales = "free_y", ncol = 2) +
  scale_linetype_manual(
    values = c("SHAP" = "solid", "Walker adjusted" = "dashed")) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(title    = "SHAP versus adjusted Walker decomposition",
       subtitle = "Block-level contributions to predicted inflation",
       x = NULL, y = "Contribution (pp)", linetype = NULL) +
  theme_thesis

ggsave("comparison_plot_contributions.png", p_contrib,
       width = 13, height = 8, dpi = 150)

# --- Attribution residual by block ---
regime_colours <- c(
  "COVID (2020–2021)" = "#5E81AC",
  "Energy Crisis"          = "#BF616A",
  "Disinflation"           = "#A3BE8C",
  "Normalization"          = "#B48EAD"
)

p_residual <- ggplot(residual_long_main,
                     aes(x = Date, y = Delta, fill = regime)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_col(alpha = 0.85, width = 25) +
  facet_wrap(~ block, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = regime_colours) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title    = "Harmonised attribution residual by block",
    subtitle = expression(Delta[kt] == C[kt]^{XGB} - C[kt]^{TVP}),
    x = NULL, y = "Residual (pp)", fill = "Regime"
  ) +
  theme_thesis

ggsave("residual_plot.png", p_residual, width = 12, height = 9, dpi = 150)

cat("Saved: all model output figures\n")
cat("\n=== PIPELINE COMPLETE ===\n")





# =============================================================================
# APPENDIX — PREDICTOR CORRELATION MATRIX
# Pearson correlations between the nine lag-1 predictors over the full
# estimation sample. Used to support the routing-under-correlation discussion
# in Section 5. Outputs:
#   - LaTeX table (kable -> stdout)
#   - PDF correlogram (appendix_correlation_matrix.pdf)
#   - CSV with raw correlations (appendix_predictor_correlations.csv)
# =============================================================================

# Use the same complete-case sample that the models are estimated on.
X_cor <- model_df %>%
  select(all_of(FEATURE_COLS)) %>%
  as.data.frame()

cor_mat <- cor(X_cor, method = "pearson", use = "complete.obs")

# Pretty display names matching the thesis block taxonomy.
pretty_names <- c(
  kpi_yoy_lag1     = "CPI inflation (AR)",
  ppi_yoy_lag1     = "PPI inflation",
  oil_yoy_lag1     = "Oil price",
  usd_nok_lag1     = "USD/NOK",
  eur_nok_lag1     = "EUR/NOK",
  import_yoy_lag1  = "Imports",
  eksport_yoy_lag1 = "Exports",
  unemp_lag1       = "Unemployment",
  rente_lag1       = "Policy rate"
)
rownames(cor_mat) <- pretty_names[rownames(cor_mat)]
colnames(cor_mat) <- pretty_names[colnames(cor_mat)]

# --- LaTeX table -------------------------------------------------------------
cor_tbl <- cor_mat %>%
  round(2) %>%
  as.data.frame()

cor_tex <- cor_tbl %>%
  kable(format     = "latex",
        booktabs   = TRUE,
        caption    = "Pearson correlations between lagged predictors, full estimation sample",
        label      = "predictor_correlations",
        align      = "r") %>%
  kable_styling(latex_options = c("scale_down", "hold_position"))

cat("\n=== Predictor correlation matrix (LaTeX) ===\n")
print(cor_tex)

# --- Visual correlogram (saved to PDF) --------------------------------------
pdf("appendix_correlation_matrix.pdf", width = 7, height = 6)
corrplot(cor_mat,
         method     = "color",
         type       = "upper",
         order      = "original",
         addCoef.col = "black",
         number.cex  = 0.7,
         tl.col      = "black",
         tl.srt      = 45,
         tl.cex      = 0.85,
         diag        = FALSE,
         mar         = c(0, 0, 1, 0))
dev.off()

# --- CSV for reference ------------------------------------------------------
write.csv(cor_mat, "appendix_predictor_correlations.csv")

cat("\nSaved:\n",
    "  appendix_correlation_matrix.pdf\n",
    "  appendix_predictor_correlations.csv\n", sep = "")