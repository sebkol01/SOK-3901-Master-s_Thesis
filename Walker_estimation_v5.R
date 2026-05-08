# =============================================================================
# Walker Benchmark v5: Lag1 Features, leakage-free backtest
# =============================================================================
#
# Main fixes relative to earlier version:
#   - No global scaling before backtest (train-only scaling at each origin)
#   - Old version removed from execution
#   - Safer coefficient extraction by name
#   - Extra MCMC diagnostics (divergences, Rhat, ESS)
#   - More robust file/path handling
#
# Setup:
#   - Target: YoY KPI inflation, h=3 months ahead
#   - Features: lag1 of 9 macro variables (same information set as XGBoost)
#   - Expanding window backtesting from January 2020
#   - train = Date < origin (strict, no leakage)
#   - Evaluation: RMSE, bias, Diebold-Mariano vs random walk
#
# =============================================================================

# --- Working directory (two-user setup) --------------------------------------
path_amund <- "/Users/amundbech1/Library/CloudStorage/OneDrive-UiTOffice365/Data_master/MasterOppgave"
path_seb   <- "~/Library/CloudStorage/OneDrive-UiTOffice365/Amund SkjalgsC8nn Bech's files - Data_master/MasterOppgave"

if (dir.exists(path_amund)) {
  setwd(path_amund)
} else if (dir.exists(path.expand(path_seb))) {
  setwd(path.expand(path_seb))
} else {
  stop("Neither Amund nor Seb working directory was found. Update the paths at the top of the script.")
}

# --- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(walker)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(rstan)
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
LOG_FILE <- "walker_diagnostics.txt"

set.seed(SEED)
options(mc.cores = N_CORES)
rstan_options(auto_write = TRUE)

# --- Helpers -----------------------------------------------------------------
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
  
  list(
    train_scaled = train_scaled,
    new_scaled   = new_scaled,
    means        = means,
    sds          = sds
  )
}

extract_fit_diagnostics <- function(fit) {
  sp <- get_sampler_params(fit$stanfit, inc_warmup = FALSE)
  n_div <- sum(sapply(sp, function(x) sum(x[, "divergent__"])))
  summ <- summary(fit$stanfit)$summary
  max_rhat <- max(summ[, "Rhat"], na.rm = TRUE)
  min_ess  <- min(summ[, "n_eff"], na.rm = TRUE)
  
  list(n_div = n_div, max_rhat = max_rhat, min_ess = min_ess)
}

predict_last_state <- function(fit, new_x, feature_cols) {
  coefs_last <- coef(fit) %>%
    ungroup() %>%
    filter(time == max(time)) %>%
    select(beta, mean)
  
  intercept_val <- if ("(Intercept)" %in% coefs_last$beta) {
    coefs_last$mean[coefs_last$beta == "(Intercept)"]
  } else {
    0
  }
  
  beta_tbl <- coefs_last %>%
    filter(beta %in% feature_cols)
  
  missing_betas <- setdiff(feature_cols, beta_tbl$beta)
  if (length(missing_betas) > 0) {
    stop(sprintf("Missing coefficients for: %s",
                 paste(missing_betas, collapse = ", ")))
  }
  
  beta_vals <- setNames(beta_tbl$mean, beta_tbl$beta)
  x_vals <- as.numeric(new_x[1, names(beta_vals), drop = TRUE])
  
  as.numeric(intercept_val + sum(beta_vals * x_vals))
}

# --- Load and prepare data ---------------------------------------------------
df <- read.csv("master_data.csv", stringsAsFactors = FALSE)
df$Date <- as.Date(df$Date)
df <- df %>% arrange(Date)

cat("Raw data:", as.character(min(df$Date)), "to", as.character(max(df$Date)),
    "(", nrow(df), "rows)\n")

# --- Feature engineering with LAG1 -------------------------------------------
df <- df %>%
  mutate(
    # First compute raw YoY variables at time t
    kpi_yoy_raw     = (kpi / lag(kpi, 12) - 1) * 100,
    ppi_yoy_raw     = (ppi / lag(ppi, 12) - 1) * 100,
    oil_yoy_raw     = (oil_price_nok / lag(oil_price_nok, 12) - 1) * 100,
    import_yoy_raw  = (import / lag(import, 12) - 1) * 100,
    eksport_yoy_raw = (eksport / lag(eksport, 12) - 1) * 100,
    
    # Then lag everything by 1 month to match XGBoost information set
    kpi_yoy_lag1     = lag(kpi_yoy_raw, 1),
    ppi_yoy_lag1     = lag(ppi_yoy_raw, 1),
    oil_yoy_lag1     = lag(oil_yoy_raw, 1),
    usd_nok_lag1     = lag(usd_nok, 1),
    eur_nok_lag1     = lag(eur_nok, 1),
    import_yoy_lag1  = lag(import_yoy_raw, 1),
    eksport_yoy_lag1 = lag(eksport_yoy_raw, 1),
    unemp_lag1       = lag(unemployment, 1),
    rente_lag1       = lag(rente, 1),
    
    # Target is h-step-ahead inflation
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

cat("Model data:", as.character(min(model_df$Date)), "to",
    as.character(max(model_df$Date)), "(", nrow(model_df), "rows)\n")

# --- Forecast origins --------------------------------------------------------
test_origins <- model_df %>%
  filter(Date >= TEST_START, !is.na(target)) %>%
  pull(Date)

cat("Forecast origins:", length(test_origins), "\n")

# --- Walker formula ----------------------------------------------------------
rw_formula <- as.formula(
  paste0(
    "target ~ -1 + rw1(~ ",
    paste(feature_cols, collapse = " + "),
    ", beta = c(0, 10), sigma = c(2, 0.01))"
  )
)

# --- Expanding window backtesting --------------------------------------------
results <- data.frame(
  Date      = as.Date(character()),
  y_actual  = numeric(),
  y_hat     = numeric(),
  y_rw      = numeric(),
  n_div     = integer(),
  max_rhat  = numeric(),
  min_ess   = numeric(),
  elapsed   = numeric(),
  stringsAsFactors = FALSE
)

cat("Walker Expanding Window Backtesting (v5, lag1 features, leakage-free)\n",
    file = LOG_FILE)
cat(paste0("Started: ", Sys.time(), "\n"), file = LOG_FILE, append = TRUE)
cat(paste0(
  "Config: chains=", N_CHAINS,
  " iter=", N_ITER,
  " warmup=", N_WARMUP,
  " adapt_delta=0.95",
  " H=", H,
  "\n\n"
), file = LOG_FILE, append = TRUE)

cat("\n=== EXPANDING WINDOW BACKTESTING ===\n")
cat("Total origins:", length(test_origins), "\n")
cat("Estimated time: ~", round(length(test_origins) * 1.5 / 60, 1),
    "hours (at ~1.5 min/fit)\n\n")

for (i in seq_along(test_origins)) {
  origin <- test_origins[i]
  
  train_raw <- model_df %>% filter(Date < origin)
  origin_row_raw <- model_df %>% filter(Date == origin)
  
  if (nrow(origin_row_raw) != 1) {
    msg <- sprintf("[%s] FAILED: origin row not unique\n", origin)
    cat(msg)
    cat(msg, file = LOG_FILE, append = TRUE)
    next
  }
  
  y_actual <- origin_row_raw$target
  y_rw <- origin_row_raw$kpi_yoy_lag1
  
  if (is.na(y_actual) || is.na(y_rw) || nrow(train_raw) < MIN_TRAIN_N) next
  
  cat(sprintf("[%d/%d] Origin: %s (n_train=%d) ... ",
              i, length(test_origins), origin, nrow(train_raw)))
  
  tryCatch({
    scaled <- scale_with_train_stats(
      train_df = train_raw,
      new_df   = origin_row_raw,
      cols     = feature_cols
    )
    
    train_scaled <- scaled$train_scaled
    new_x_scaled <- scaled$new_scaled[, feature_cols, drop = FALSE]
    
    t_start <- Sys.time()
    
    fit <- walker(
      formula       = rw_formula,
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
    diag <- extract_fit_diagnostics(fit)
    y_hat <- predict_last_state(fit, new_x_scaled, feature_cols)
    
    results <- rbind(results, data.frame(
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
    
    write.csv(results, "walker_forecasts.csv", row.names = FALSE)
    
  }, error = function(e) {
    msg <- sprintf("[%s] FAILED: %s\n", origin, e$message)
    cat(msg)
    cat(msg, file = LOG_FILE, append = TRUE)
  })
}

# --- Evaluation --------------------------------------------------------------
cat("\n\n=== EVALUATION ===\n")

if (nrow(results) > 0) {
  results <- results %>%
    mutate(
      err_walker = y_hat - y_actual,
      err_rw     = y_rw - y_actual
    )
  
  rmse_walker <- sqrt(mean(results$err_walker^2, na.rm = TRUE))
  rmse_rw     <- sqrt(mean(results$err_rw^2, na.rm = TRUE))
  bias_walker <- mean(results$err_walker, na.rm = TRUE)
  bias_rw     <- mean(results$err_rw, na.rm = TRUE)
  
  cat(sprintf("Walker:  RMSE=%.4f  Bias=%.4f\n", rmse_walker, bias_walker))
  cat(sprintf("RW:      RMSE=%.4f  Bias=%.4f\n", rmse_rw, bias_rw))
  cat(sprintf("Ratio (Walker/RW): %.4f\n", rmse_walker / rmse_rw))
  
  # Diebold-Mariano with Newey-West HAC
  d <- results$err_rw^2 - results$err_walker^2
  d <- d[is.finite(d)]
  n <- length(d)
  
  if (n > 5) {
    nw_lag <- max(H - 1, floor(n^(1/3)))
    d_mean <- mean(d)
    
    gamma_0 <- mean((d - d_mean)^2)
    gamma_sum <- 0
    
    for (k in 1:nw_lag) {
      if ((n - k) <= 0) break
      gamma_k <- mean((d[(k + 1):n] - d_mean) * (d[1:(n - k)] - d_mean))
      gamma_sum <- gamma_sum + 2 * (1 - k / (nw_lag + 1)) * gamma_k
    }
    
    var_d <- (gamma_0 + gamma_sum) / n
    dm_stat <- d_mean / sqrt(max(var_d, 1e-10))
    dm_pval <- 2 * (1 - pnorm(abs(dm_stat)))
    
    cat(sprintf("\nDM test (walker vs RW): stat=%.3f, p=%.4f\n", dm_stat, dm_pval))
    cat("(Positive stat => walker better)\n")
  } else {
    cat("\nDM test skipped: too few forecast errors.\n")
  }
  
  cat(sprintf("\nTotal divergent transitions: %d (across %d fits)\n",
              sum(results$n_div, na.rm = TRUE), nrow(results)))
  cat(sprintf("Worst Rhat: %.3f\n", max(results$max_rhat, na.rm = TRUE)))
  cat(sprintf("Lowest ESS: %.0f\n", min(results$min_ess, na.rm = TRUE)))
  cat(sprintf("Total computation time: %.1f hours\n",
              sum(results$elapsed, na.rm = TRUE) / 60))
  
  write.csv(results, "walker_forecasts.csv", row.names = FALSE)
  cat("Saved: walker_forecasts.csv\n")
} else {
  cat("No successful backtest fits were produced.\n")
}

# =============================================================================
# PART 2: Full-sample beta_t paths (for SHAP comparison)
# =============================================================================
cat("\n\n=== FULL-SAMPLE FIT ===\n")

full_data_raw <- model_df %>% filter(!is.na(target))

tryCatch({
  full_scaling <- scale_with_train_stats(
    train_df = full_data_raw,
    new_df   = full_data_raw,
    cols     = feature_cols
  )
  full_data <- full_scaling$train_scaled
  
  t0 <- Sys.time()
  
  fit_full <- walker(
    formula       = rw_formula,
    data          = full_data,
    sigma_y_prior = c(2, 0.01),
    chains        = N_CHAINS,
    iter          = 3000,
    warmup        = 1500,
    cores         = N_CORES,
    refresh       = 500,
    control       = list(adapt_delta = 0.97, max_treedepth = 14)
  )
  
  t_full <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  diag_full <- extract_fit_diagnostics(fit_full)
  
  cat(sprintf("Completed in %.1f minutes\n", t_full))
  cat(sprintf("Divergent transitions: %d\n", diag_full$n_div))
  cat(sprintf("Worst Rhat: %.3f\n", diag_full$max_rhat))
  cat(sprintf("Lowest ESS: %.0f\n", diag_full$min_ess))
  
  coefs_full <- coef(fit_full) %>% ungroup()
  
  beta_paths <- coefs_full %>%
    select(time, beta, mean) %>%
    pivot_wider(names_from = beta, values_from = mean) %>%
    mutate(Date = full_data$Date[time]) %>%
    select(Date, everything(), -time)
  
  beta_ci <- coefs_full %>%
    select(time, beta, `2.5%`, `97.5%`) %>%
    pivot_wider(names_from = beta, values_from = c(`2.5%`, `97.5%`)) %>%
    mutate(Date = full_data$Date[time]) %>%
    select(Date, everything(), -time)
  
  scaling_info <- data.frame(
    feature = feature_cols,
    mean    = full_scaling$means[feature_cols],
    sd      = full_scaling$sds[feature_cols],
    row.names = NULL
  )
  
  write.csv(beta_paths, "walker_beta_paths.csv", row.names = FALSE)
  write.csv(beta_ci, "walker_beta_ci.csv", row.names = FALSE)
  write.csv(scaling_info, "walker_scaling_info.csv", row.names = FALSE)
  cat("Saved: walker_beta_paths.csv, walker_beta_ci.csv, walker_scaling_info.csv\n")
  
  # --- Walker contributions (beta_t * x_t) for SHAP comparison ---
  walker_contributions <- data.frame(Date = full_data$Date)
  
  for (col in feature_cols) {
    if (!col %in% names(beta_paths)) {
      stop(sprintf("Column '%s' missing in beta_paths.", col))
    }
    walker_contributions[[col]] <- beta_paths[[col]] * full_data[[col]]
  }
  
  if ("(Intercept)" %in% names(beta_paths)) {
    walker_contributions[["intercept"]] <- beta_paths[["(Intercept)"]]
  }
  
  write.csv(walker_contributions, "walker_contributions.csv", row.names = FALSE)
  cat("Saved: walker_contributions.csv\n")
  
  # --- Beta path plot ---
  beta_long <- coefs_full %>%
    filter(beta != "(Intercept)") %>%
    mutate(Date = full_data$Date[time])
  
  regime_breaks <- as.Date(c("2020-03-01", "2021-09-01", "2023-06-01", "2024-06-01"))
  
  p_beta <- ggplot(beta_long, aes(x = Date)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "steelblue", alpha = 0.2) +
    geom_line(aes(y = mean), color = "steelblue", linewidth = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = regime_breaks, linetype = "dotted", color = "red", alpha = 0.5) +
    facet_wrap(~ beta, scales = "free_y", ncol = 3) +
    labs(
      title = "Time-Varying Coefficients (Walker RW1, lag1 features)",
      subtitle = "Posterior mean with 95% CI. Red lines = regime boundaries.",
      y = expression(beta[t] ~ "(scaled features)"),
      x = NULL
    ) +
    theme_minimal() +
    theme(strip.text = element_text(size = 9))
  
  ggsave("walker_beta_paths_plot.png", p_beta, width = 12, height = 8, dpi = 150)
  
  p_builtin <- plot_coefs(fit_full, scales = "free") + theme_minimal()
  ggsave("walker_coefs_builtin.png", p_builtin, width = 12, height = 8, dpi = 150)
  
  cat("Saved plots.\n")
  
}, error = function(e) {
  cat(sprintf("Full-sample fit FAILED: %s\n", e$message))
})

# =============================================================================
# PART 3: Forecast plots
# =============================================================================
if (nrow(results) > 0) {
  p_forecast <- ggplot(results, aes(x = Date)) +
    geom_line(aes(y = y_actual, color = "Actual"), linewidth = 0.8) +
    geom_line(aes(y = y_hat, color = "Walker"), linewidth = 0.8) +
    geom_line(aes(y = y_rw, color = "Random Walk"), linewidth = 0.7, linetype = "dashed") +
    scale_color_manual(values = c(
      "Actual" = "black",
      "Walker" = "steelblue",
      "Random Walk" = "grey50"
    )) +
    labs(
      title = "Walker vs Random Walk: h=3 YoY KPI Forecast (lag1 features)",
      y = "YoY KPI Inflation (%)",
      x = NULL,
      color = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave("walker_forecast_plot.png", p_forecast, width = 10, height = 5, dpi = 150)
  
  p_error <- ggplot(results %>% mutate(
    err_walker = y_hat - y_actual,
    err_rw     = y_rw - y_actual
  ), aes(x = Date)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_line(aes(y = err_walker, color = "Walker"), linewidth = 0.7) +
    geom_line(aes(y = err_rw, color = "Random Walk"), linewidth = 0.7, linetype = "dashed") +
    scale_color_manual(values = c(
      "Walker" = "steelblue",
      "Random Walk" = "grey50"
    )) +
    labs(
      title = "Forecast Errors (h=3, lag1 features)",
      y = "Error (pp)",
      x = NULL,
      color = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  ggsave("walker_error_plot.png", p_error, width = 10, height = 5, dpi = 150)
  
  cat("Saved forecast plots.\n")
  print(p_forecast)
} else {
  cat("Skipping forecast plots because results is empty.\n")
}

# =============================================================================
# PART 4: Post-run diagnostics from saved forecasts
# =============================================================================
if (file.exists("walker_forecasts.csv")) {
  res <- read.csv("walker_forecasts.csv", stringsAsFactors = FALSE)
  res$Date <- as.Date(res$Date)
  
  cat("\n=== BASIC STATS ===\n")
  cat("n:", nrow(res), "\n")
  cat("y_actual range:", paste(range(res$y_actual, na.rm = TRUE), collapse = " to "), "\n")
  cat("y_hat range:", paste(range(res$y_hat, na.rm = TRUE), collapse = " to "), "\n")
  cat("y_rw range:", paste(range(res$y_rw, na.rm = TRUE), collapse = " to "), "\n")
  
  res$err <- res$y_hat - res$y_actual
  res$err_rw <- res$y_rw - res$y_actual
  
  cat("\n=== EVALUATION ===\n")
  cat("Walker RMSE:", round(sqrt(mean(res$err^2, na.rm = TRUE)), 4), "\n")
  cat("RW RMSE:", round(sqrt(mean(res$err_rw^2, na.rm = TRUE)), 4), "\n")
  cat("Walker bias:", round(mean(res$err, na.rm = TRUE), 4), "\n")
  cat("RW bias:", round(mean(res$err_rw, na.rm = TRUE), 4), "\n")
  cat("Ratio:", round(
    sqrt(mean(res$err^2, na.rm = TRUE)) / sqrt(mean(res$err_rw^2, na.rm = TRUE)), 4
  ), "\n")
  
  cat("\n=== WORST PREDICTIONS ===\n")
  res_sorted <- res[order(-abs(res$err)), ]
  print(head(res_sorted[, c("Date", "y_actual", "y_hat", "y_rw", "err")], 10))
  
  cat("\n=== CORRELATION ===\n")
  cat("cor(y_hat, y_actual):", round(cor(res$y_hat, res$y_actual, use = "complete.obs"), 4), "\n")
  cat("cor(y_rw, y_actual):", round(cor(res$y_rw, res$y_actual, use = "complete.obs"), 4), "\n")
  
  cat("\n=== DIVERGENCES ===\n")
  cat("Total:", sum(res$n_div, na.rm = TRUE), "\n")
  if ("max_rhat" %in% names(res)) cat("Worst Rhat:", max(res$max_rhat, na.rm = TRUE), "\n")
  if ("min_ess" %in% names(res)) cat("Lowest ESS:", min(res$min_ess, na.rm = TRUE), "\n")
  
  cat("\n=== SLOW ORIGINS (>10 min) ===\n")
  slow <- res[res$elapsed > 10, c("Date", "elapsed", "y_hat")]
  print(slow)
  
  cat("\n=== FIRST FEW RW CHECKS ===\n")
  print(head(res[, c("Date", "y_actual", "y_rw")], 5))
  
  cat("\n=== MATCH AGAINST MODEL DATA ===\n")
  print(
    model_df %>%
      filter(Date >= TEST_START) %>%
      head(5) %>%
      select(Date, kpi_yoy_raw, kpi_yoy_lag1)
  )
} else {
  cat("\nwalker_forecasts.csv not found, so post-run diagnostics were skipped.\n")
}

cat("\n=== ALL DONE ===\n")
