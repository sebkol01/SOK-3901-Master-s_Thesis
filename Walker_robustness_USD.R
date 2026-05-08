# =============================================================================
# Walker Robustness Checks
# =============================================================================
# This script runs three Walker specifications:
#   1. BASELINE:      -1 + rw1(~ 1 + ...)    NOK-denominated oil/trade
#   2. NO_INTERCEPT:  -1 + rw1(~ 0 + ...)    NOK-denominated oil/trade
#   3. USD_EXTERNAL:  -1 + rw1(~ 1 + ...)    USD-denominated oil/trade
#
# Each specification produces its own output files with a SPEC-prefixed name,
# so downstream harmonization and comparison scripts can read them separately.
# =============================================================================

# Set working directory (as in Walker_estimation.R)
# setwd("...")

suppressPackageStartupMessages({
  library(walker)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(rstan)
  library(readr)
  library(rlang)
})

# --- Configuration -----------------------------------------------------------
H <- 3
TEST_START <- as.Date("2020-01-01")
MIN_TRAIN_N <- 50
N_CHAINS <- 4
N_ITER <- 2000
N_WARMUP <- 1000
N_CORES <- min(4, parallel::detectCores())
SEED <- 42

# Which specification to run. Choose one at a time, OR run all three
# sequentially at the bottom of the script.
#   "baseline"     -> matches the main paper spec
#   "no_intercept" -> drops time-varying intercept
#   "usd_external" -> uses USD-denominated oil and trade
SPEC <- "usd_external"

set.seed(SEED)
options(mc.cores = N_CORES)
rstan_options(auto_write = TRUE)

# --- Helpers (unchanged from main script) -----------------------------------
scale_with_train_stats <- function(train_df, new_df, cols) {
  means <- sapply(train_df[, cols, drop = FALSE], mean, na.rm = TRUE)
  sds   <- sapply(train_df[, cols, drop = FALSE], sd,   na.rm = TRUE)

  bad_sd <- is.na(sds) | sds == 0
  if (any(bad_sd)) {
    stop(sprintf("Zero/NA SD in training data for: %s",
                 paste(names(sds)[bad_sd], collapse = ", ")))
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

extract_fit_diagnostics <- function(fit) {
  sp <- get_sampler_params(fit$stanfit, inc_warmup = FALSE)
  n_div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  summ <- summary(fit$stanfit)$summary
  max_rhat <- max(summ[, "Rhat"], na.rm = TRUE)
  min_ess  <- min(summ[, "n_eff"], na.rm = TRUE)
  list(n_div = n_div, max_rhat = max_rhat, min_ess = min_ess)
}

predict_last_state <- function(fit, new_x, feature_cols, has_intercept = TRUE) {
  coefs_last <- coef(fit) %>%
    ungroup() %>%
    filter(time == max(time)) %>%
    select(beta, mean)

  if (has_intercept) {
    if (!"(Intercept)" %in% coefs_last$beta) {
      stop("Expected '(Intercept)' in coefficients but not found.")
    }
    intercept_val <- coefs_last$mean[coefs_last$beta == "(Intercept)"]
  } else {
    intercept_val <- 0
  }

  beta_tbl <- coefs_last %>% filter(beta %in% feature_cols)
  beta_vals <- setNames(beta_tbl$mean, beta_tbl$beta)
  x_vals <- as.numeric(new_x[1, names(beta_vals), drop = TRUE])

  as.numeric(intercept_val + sum(beta_vals * x_vals))
}

extract_filtered_coefs <- function(fit, origin_date) {
  coef(fit) %>%
    ungroup() %>%
    filter(time == max(time)) %>%
    select(beta, mean, sd, `2.5%`, `50%`, `97.5%`) %>%
    mutate(origin = origin_date)
}

# --- Load data ---------------------------------------------------------------
df <- read.csv("master_data.csv", stringsAsFactors = FALSE)
df$Date <- as.Date(df$Date)
df <- df %>% arrange(Date)

# --- Feature engineering -----------------------------------------------------
# For USD_EXTERNAL spec, we need oljepris_USD, import_USD, eksport_USD in the
# raw data. These should be added to master_data.csv before running this spec.
# If they do not exist, the script will fall back to NOK-denominated versions
# and warn.

if (SPEC == "usd_external") {
  needed <- c("oljepris_USD", "import_USD", "eksport_USD")
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0) {
    stop(sprintf("USD_EXTERNAL spec requires columns: %s. Missing: %s.
     Add these to master_data.csv before running this spec.",
                 paste(needed, collapse = ", "),
                 paste(missing, collapse = ", ")))
  }
  df <- df %>%
    mutate(
      oil_price_for_model = oljepris_USD,
      import_for_model    = import_USD,
      eksport_for_model   = eksport_USD
    )
  cat("Using USD-denominated oil and trade (oljepris_USD, import_USD, eksport_USD).\n")
} else {
  df <- df %>%
    mutate(
      oil_price_for_model = oil_price_nok,
      import_for_model    = import,
      eksport_for_model   = eksport
    )
  cat("Using NOK-denominated oil and trade (baseline).\n")
}

df <- df %>%
  mutate(
    kpi_yoy_raw     = (kpi / lag(kpi, 12) - 1) * 100,
    ppi_yoy_raw     = (ppi / lag(ppi, 12) - 1) * 100,
    oil_yoy_raw     = (oil_price_for_model / lag(oil_price_for_model, 12) - 1) * 100,
    import_yoy_raw  = (import_for_model / lag(import_for_model, 12) - 1) * 100,
    eksport_yoy_raw = (eksport_for_model / lag(eksport_for_model, 12) - 1) * 100,

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

feature_cols <- c(
  "kpi_yoy_lag1", "ppi_yoy_lag1", "oil_yoy_lag1",
  "usd_nok_lag1", "eur_nok_lag1",
  "import_yoy_lag1", "eksport_yoy_lag1",
  "unemp_lag1", "rente_lag1"
)

model_df <- df %>%
  select(Date, target, kpi_yoy_raw, all_of(feature_cols)) %>%
  filter(complete.cases(.))

test_origins <- model_df %>%
  filter(Date >= TEST_START, !is.na(target)) %>%
  pull(Date)

# --- Build formula based on SPEC --------------------------------------------
if (SPEC == "no_intercept") {
  rw_formula <- as.formula(
    paste0(
      "target ~ -1 + rw1(~ 0 + ",
      paste(feature_cols, collapse = " + "),
      ", beta = c(0, 10), sigma = c(2, 0.01))"
    )
  )
  has_intercept <- FALSE
  cat("Specification: NO_INTERCEPT (rw1 ~ 0 + features)\n")
} else {
  rw_formula <- as.formula(
    paste0(
      "target ~ -1 + rw1(~ 1 + ",
      paste(feature_cols, collapse = " + "),
      ", beta = c(0, 10), sigma = c(2, 0.01))"
    )
  )
  has_intercept <- TRUE
  cat(sprintf("Specification: %s (rw1 ~ 1 + features, with time-varying intercept)\n",
              toupper(SPEC)))
}

# --- Output file prefix ------------------------------------------------------
OUT_PREFIX <- paste0("walker_", SPEC, "_")
LOG_FILE   <- paste0(OUT_PREFIX, "diagnostics.txt")

cat(sprintf("Output prefix: %s\n", OUT_PREFIX))
cat(sprintf("Total forecast origins: %d\n\n", length(test_origins)))

# --- Expanding window backtest ----------------------------------------------
results <- data.frame(
  Date = as.Date(character()), y_actual = numeric(), y_hat = numeric(),
  y_rw = numeric(), n_div = integer(), max_rhat = numeric(),
  min_ess = numeric(), elapsed = numeric(),
  stringsAsFactors = FALSE
)

filtered_coefs_list <- list()
x_values_list <- list()

cat(sprintf("Walker Robustness Run (%s) started: %s\n", SPEC, Sys.time()),
    file = LOG_FILE)

for (i in seq_along(test_origins)) {
  origin <- test_origins[i]
  train_raw <- model_df %>% filter(Date < origin)
  origin_row_raw <- model_df %>% filter(Date == origin)

  if (nrow(origin_row_raw) != 1) next
  y_actual <- origin_row_raw$target
  y_rw <- origin_row_raw$kpi_yoy_lag1
  if (is.na(y_actual) || is.na(y_rw) || nrow(train_raw) < MIN_TRAIN_N) next

  cat(sprintf("[%d/%d] %s ... ", i, length(test_origins), origin))

  tryCatch({
    scaled <- scale_with_train_stats(train_raw, origin_row_raw, feature_cols)
    train_scaled <- scaled$train_scaled
    new_x_scaled <- scaled$new_scaled[, feature_cols, drop = FALSE]

    t_start <- Sys.time()

    fit <- walker(
      formula = rw_formula, data = train_scaled,
      sigma_y_prior = c(2, 0.01),
      chains = N_CHAINS, iter = N_ITER, warmup = N_WARMUP,
      cores = N_CORES, refresh = 0,
      control = list(adapt_delta = 0.95, max_treedepth = 12)
    )

    elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
    diag <- extract_fit_diagnostics(fit)
    y_hat <- predict_last_state(fit, new_x_scaled, feature_cols, has_intercept)

    filtered_coefs_list[[i]] <- extract_filtered_coefs(fit, origin)
    x_values_list[[i]] <- as.data.frame(new_x_scaled) %>%
      pivot_longer(everything(), names_to = "beta", values_to = "x_scaled") %>%
      mutate(origin = origin)

    results <- rbind(results, data.frame(
      Date = origin, y_actual = y_actual, y_hat = y_hat, y_rw = y_rw,
      n_div = diag$n_div, max_rhat = diag$max_rhat, min_ess = diag$min_ess,
      elapsed = elapsed, stringsAsFactors = FALSE
    ))

    cat(sprintf("done (%.1f min, yhat=%.2f)\n", elapsed, y_hat))
    write.csv(results, paste0(OUT_PREFIX, "forecasts.csv"), row.names = FALSE)

  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", e$message))
  })
}

# --- Build real-time contributions ------------------------------------------
if (length(filtered_coefs_list) > 0) {
  filtered_coefs_df <- bind_rows(filtered_coefs_list)
  x_values_df <- bind_rows(x_values_list)

  write.csv(filtered_coefs_df, paste0(OUT_PREFIX, "filtered_coefs.csv"),
            row.names = FALSE)

  realtime_contributions_long <- filtered_coefs_df %>%
    left_join(x_values_df, by = c("origin", "beta")) %>%
    mutate(
      contribution = case_when(
        beta == "(Intercept)" ~ mean,
        TRUE ~ mean * x_scaled
      )
    ) %>%
    select(origin, beta, contribution)

  walker_contributions <- realtime_contributions_long %>%
    pivot_wider(names_from = beta, values_from = contribution) %>%
    rename(Date = origin)

  if ("(Intercept)" %in% names(walker_contributions)) {
    walker_contributions <- walker_contributions %>%
      rename(intercept = `(Intercept)`)
  } else {
    walker_contributions$intercept <- 0
  }

  write.csv(walker_contributions, paste0(OUT_PREFIX, "contributions.csv"),
            row.names = FALSE)
}

# --- Evaluation summary ------------------------------------------------------
if (nrow(results) > 0) {
  results <- results %>%
    mutate(err_walker = y_hat - y_actual, err_rw = y_rw - y_actual)

  rmse_walker <- sqrt(mean(results$err_walker^2, na.rm = TRUE))
  rmse_rw     <- sqrt(mean(results$err_rw^2, na.rm = TRUE))
  mae_walker  <- mean(abs(results$err_walker), na.rm = TRUE)

  summary_row <- data.frame(
    spec = SPEC,
    n_forecasts = nrow(results),
    rmse = rmse_walker,
    mae = mae_walker,
    rmse_rw = rmse_rw,
    rmse_ratio = rmse_walker / rmse_rw,
    total_div = sum(results$n_div, na.rm = TRUE),
    max_rhat = max(results$max_rhat, na.rm = TRUE),
    min_ess = min(results$min_ess, na.rm = TRUE)
  )

  write.csv(summary_row, paste0(OUT_PREFIX, "eval_summary.csv"),
            row.names = FALSE)

  cat(sprintf("\n=== %s SUMMARY ===\n", toupper(SPEC)))
  print(summary_row)
}

cat("\nDone. Output files have prefix:", OUT_PREFIX, "\n")
