
# setwd for Amund
setwd("/Users/amundbech1/Library/CloudStorage/OneDrive-UiTOffice365/Data_master/MasterOppgave")

# setwd for Seb
setwd("~/Library/CloudStorage/OneDrive-UiTOffice365/Amund Skjalgsønn Bech's files - Data_master/MasterOppgave")

##############
# --- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(walker)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(rstan)
  library(dplyr)
  library(readr)
  library(rlang)
  library(knitr)
  library(kableExtra)
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
  
  if (!"(Intercept)" %in% coefs_last$beta) {
    stop("Expected '(Intercept)' in coefficients but not found. 
          Check that formula uses 'rw1(~ 1 + ...)'")
  }
  intercept_val <- coefs_last$mean[coefs_last$beta == "(Intercept)"]
  
  beta_tbl <- coefs_last %>% filter(beta %in% feature_cols)
  
  missing_betas <- setdiff(feature_cols, beta_tbl$beta)
  if (length(missing_betas) > 0) {
    stop(sprintf("Missing coefficients for: %s",
                 paste(missing_betas, collapse = ", ")))
  }
  
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
# Walker formula - MED tidsvarierende intercept
rw_formula <- as.formula(
  paste0(
    "target ~ -1 + rw1(~ 1 + ",                    # <-- "1 +" lagt til
    paste(feature_cols, collapse = " + "),
    ", beta = c(0, 10), sigma = c(2, 0.01))"
  )
)

# =============================================================================
# PART 1: Expanding window backtesting (real-time, with filtered coefs stored)
# =============================================================================
results <- data.frame(
  Date      = as.Date(character()),
  y_actual  = numeric(),
  y_hat     = numeric(),
  y_rw      = numeric(),
  y_ar1     = numeric(),
  n_div     = integer(),
  max_rhat  = numeric(),
  min_ess   = numeric(),
  elapsed   = numeric(),
  stringsAsFactors = FALSE
)

filtered_coefs_list <- list()
x_values_list <- list()

cat("Walker Expanding Window Backtesting (v6: intercept + real-time decomposition)\n",
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
    
    # Store filtered coefs (real-time, NOT smoothed)
    filtered_coefs_list[[i]] <- extract_filtered_coefs(fit, origin)
    
    # Store x-values at origin (standardized scale) for decomposition
    x_values_list[[i]] <- as.data.frame(new_x_scaled) %>%
      pivot_longer(everything(), names_to = "beta", values_to = "x_scaled") %>%
      mutate(origin = origin)
    
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

# --- Build real-time decomposition + save -----------------------------------
if (length(filtered_coefs_list) > 0) {
  filtered_coefs_df <- bind_rows(filtered_coefs_list)
  x_values_df <- bind_rows(x_values_list)
  
  write.csv(filtered_coefs_df, "walker_filtered_coefs.csv", row.names = FALSE)
  write.csv(x_values_df, "walker_x_values_at_origin.csv", row.names = FALSE)
  
  # Build real-time contributions: c_jt = beta_jt * x_jt (intercept = mean directly)
  realtime_contributions_long <- filtered_coefs_df %>%
    left_join(x_values_df, by = c("origin", "beta")) %>%
    mutate(
      contribution = case_when(
        beta == "(Intercept)" ~ mean,
        TRUE ~ mean * x_scaled
      )
    ) %>%
    select(origin, beta, contribution)
  
  # Wide format - matches old walker_contributions.csv structure so downstream code works
  walker_contributions <- realtime_contributions_long %>%
    pivot_wider(names_from = beta, values_from = contribution) %>%
    rename(Date = origin, intercept = `(Intercept)`)
  
  write.csv(walker_contributions, "walker_contributions.csv", row.names = FALSE)
  cat("Saved: walker_filtered_coefs.csv, walker_x_values_at_origin.csv, walker_contributions.csv (REAL-TIME)\n")
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
# PART 2: Full-sample fit (for coefficient path visualization only)
# NOTE: These SMOOTHED contributions are NOT used for SHAP comparison.
# The real-time contributions above (walker_contributions.csv) are what 
# downstream harmonization code should use.
# =============================================================================
# cat("\n\n=== FULL-SAMPLE FIT (reference only) ===\n")
# 
# full_data_raw <- model_df %>% filter(!is.na(target))
# 
# tryCatch({
#   full_scaling <- scale_with_train_stats(
#     train_df = full_data_raw,
#     new_df   = full_data_raw,
#     cols     = feature_cols
#   )
#   full_data <- full_scaling$train_scaled
# 
#   t0 <- Sys.time()
# 
#   fit_full <- walker(
#     formula       = rw_formula,
#     data          = full_data,
#     sigma_y_prior = c(2, 0.01),
#     chains        = N_CHAINS,
#     iter          = 3000,
#     warmup        = 1500,
#     cores         = N_CORES,
#     refresh       = 500,
#     control       = list(adapt_delta = 0.97, max_treedepth = 14)
#   )
# 
#   t_full <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
#   diag_full <- extract_fit_diagnostics(fit_full)
# 
#   cat(sprintf("Completed in %.1f minutes\n", t_full))
#   cat(sprintf("Divergent transitions: %d\n", diag_full$n_div))
#   cat(sprintf("Worst Rhat: %.3f\n", diag_full$max_rhat))
#   cat(sprintf("Lowest ESS: %.0f\n", diag_full$min_ess))
# 
#   coefs_full <- coef(fit_full) %>% ungroup()
# 
#   beta_paths <- coefs_full %>%
#     select(time, beta, mean) %>%
#     pivot_wider(names_from = beta, values_from = mean) %>%
#     mutate(Date = full_data$Date[time]) %>%
#     select(Date, everything(), -time)
# 
#   beta_ci <- coefs_full %>%
#     select(time, beta, `2.5%`, `97.5%`) %>%
#     pivot_wider(names_from = beta, values_from = c(`2.5%`, `97.5%`)) %>%
#     mutate(Date = full_data$Date[time]) %>%
#     select(Date, everything(), -time)
# 
#   scaling_info <- data.frame(
#     feature = feature_cols,
#     mean    = full_scaling$means[feature_cols],
#     sd      = full_scaling$sds[feature_cols],
#     row.names = NULL
#   )
# 
#   write.csv(beta_paths, "walker_beta_paths.csv", row.names = FALSE)
#   write.csv(beta_ci, "walker_beta_ci.csv", row.names = FALSE)
#   write.csv(scaling_info, "walker_scaling_info.csv", row.names = FALSE)
#   cat("Saved: walker_beta_paths.csv, walker_beta_ci.csv, walker_scaling_info.csv\n")
# 
#   # Full-sample (SMOOTHED) contributions - saved separately for reference
#   walker_contributions_fullsample <- data.frame(Date = full_data$Date)
# 
#   for (col in feature_cols) {
#     if (!col %in% names(beta_paths)) {
#       stop(sprintf("Column '%s' missing in beta_paths.", col))
#     }
#     walker_contributions_fullsample[[col]] <- beta_paths[[col]] * full_data[[col]]
#   }
# 
#   if ("(Intercept)" %in% names(beta_paths)) {
#     walker_contributions_fullsample[["intercept"]] <- beta_paths[["(Intercept)"]]
#   }
# 
#   write.csv(walker_contributions_fullsample,
#             "walker_contributions_fullsample.csv", row.names = FALSE)
#   cat("Saved: walker_contributions_fullsample.csv (smoothed - reference only)\n")
# 
#   # Beta path plot
#   beta_long <- coefs_full %>%
#     mutate(Date = full_data$Date[time])
# 
#   regime_breaks <- as.Date(c("2020-03-01", "2021-09-01", "2023-06-01", "2024-06-01"))
# 
#   p_beta <- ggplot(beta_long, aes(x = Date)) +
#     geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "steelblue", alpha = 0.2) +
#     geom_line(aes(y = mean), color = "steelblue", linewidth = 0.7) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
#     geom_vline(xintercept = regime_breaks, linetype = "dotted", color = "red", alpha = 0.5) +
#     facet_wrap(~ beta, scales = "free_y", ncol = 3) +
#     labs(
#       title = "Time-Varying Coefficients (Walker RW1, lag1 + intercept)",
#       subtitle = "Posterior mean with 95% CI. Red lines = regime boundaries.",
#       y = expression(beta[t] ~ "(scaled features)"),
#       x = NULL
#     ) +
#     theme_minimal() +
#     theme(strip.text = element_text(size = 9))
# 
#   ggsave("walker_beta_paths_plot.png", p_beta, width = 12, height = 8, dpi = 150)
# 
#   p_builtin <- plot_coefs(fit_full, scales = "free") + theme_minimal()
#   ggsave("walker_coefs_builtin.png", p_builtin, width = 12, height = 8, dpi = 150)
# 
#   cat("Saved plots.\n")
# 
# }, error = function(e) {
#   cat(sprintf("Full-sample fit FAILED: %s\n", e$message))
# })
# p_beta
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
p_forecast
p_error
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

# -----------------------------
# Input files produced by Walker run
# -----------------------------
walker_contrib_file <- "walker_contributions.csv"
walker_beta_file    <- "walker_beta_paths.csv"
walker_forecast_file <- "walker_forecasts.csv"

stopifnot(file.exists(walker_contrib_file))
stopifnot(file.exists(walker_beta_file))

walker_contrib <- read_csv(walker_contrib_file, show_col_types = FALSE)
walker_beta    <- read_csv(walker_beta_file, show_col_types = FALSE)

if (!"Date" %in% names(walker_contrib)) stop("walker_contributions.csv must contain Date")
if (!"Date" %in% names(walker_beta)) stop("walker_beta_paths.csv must contain Date")

walker_contrib <- walker_contrib %>% mutate(Date = as.Date(Date))
walker_beta    <- walker_beta %>% mutate(Date = as.Date(Date))

# -----------------------------
# Economic block map
# -----------------------------
block_map <- tibble::tribble(
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

missing_features <- setdiff(block_map$feature, names(walker_contrib))
if (length(missing_features) > 0) {
  stop(sprintf("Missing expected feature columns in walker_contributions.csv: %s",
               paste(missing_features, collapse = ", ")))
}

# -----------------------------
# The real-time Walker contributions are built from the expanding-window
# backtest: at each forecast origin, features are standardized using
# training-window means/sds, so x has mean ~0 in the training sample.
# Therefore beta*x is approximately the centred contribution:
#   c_tilde_jt = (x_jt - xbar_j) * beta_jt, since xbar_j ~= 0 after scaling.
# baseline_WalkerAdj equals the time-varying intercept from the real-time fit.
# This decomposition is on the same informational footing as the SHAP values
# (both use only information available at forecast origin).
# -----------------------------

if (!"intercept" %in% names(walker_contrib)) {
  if ("(Intercept)" %in% names(walker_beta)) {
    walker_contrib <- walker_contrib %>%
      left_join(walker_beta %>% select(Date, `(Intercept)`) %>% rename(intercept = `(Intercept)`), by = "Date")
  } else {
    stop("Could not find intercept in walker_contributions.csv or walker_beta_paths.csv")
  }
}


# -----------------------------
# Feature-level long object
# -----------------------------
walker_feature_long <- walker_contrib %>%
  select(Date, all_of(block_map$feature), intercept) %>%
  pivot_longer(cols = all_of(block_map$feature), names_to = "feature", values_to = "contribution") %>%
  left_join(block_map, by = "feature")

# -----------------------------
# Block-level signed contributions
# -----------------------------
walker_block_signed <- walker_feature_long %>%
  group_by(Date, block) %>%
  summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = block, values_from = contribution) %>%
  left_join(walker_contrib %>% select(Date, intercept), by = "Date") %>%
  rename(baseline_WalkerAdj = intercept) %>%
  arrange(Date)

# -----------------------------
# Block-level absolute contributions
# -----------------------------
walker_block_abs <- walker_feature_long %>%
  group_by(Date, block) %>%
  summarise(abs_contribution = sum(abs(contribution), na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = block, values_from = abs_contribution) %>%
  arrange(Date)

# -----------------------------
# Signed share: C_kt / sum_m |C_mt|
# Absolute share: |C_kt| / sum_m |C_mt|
# -----------------------------
block_names <- setdiff(names(walker_block_signed), c("Date", "baseline_WalkerAdj"))

walker_shares <- walker_block_signed %>%
  rowwise() %>%
  mutate(total_abs = sum(abs(c_across(all_of(block_names))), na.rm = TRUE)) %>%
  ungroup()

for (blk in block_names) {
  walker_shares[[paste0("S_", blk)]] <- ifelse(
    walker_shares$total_abs > 0,
    walker_shares[[blk]] / walker_shares$total_abs,
    NA_real_
  )
  walker_shares[[paste0("A_", blk)]] <- ifelse(
    walker_shares$total_abs > 0,
    abs(walker_shares[[blk]]) / walker_shares$total_abs,
    NA_real_
  )
}

# -----------------------------
# Regimes
# -----------------------------
walker_shares <- walker_shares %>%
  mutate(
    regime = case_when(
      Date < as.Date("2021-06-01") ~ "COVID (2020–2021)",
      Date < as.Date("2023-01-01") ~ "Energy Crisis",
      Date < as.Date("2024-06-01") ~ "Disinflation",
      TRUE ~ "Normalization"
    )
  )
share_cols <- grep("^(S_|A_)", names(walker_shares), value = TRUE)
walker_regime_summary <- walker_shares %>%
  group_by(regime) %>%
  summarise(across(all_of(share_cols), ~mean(.x, na.rm = TRUE)), .groups = "drop")

walker_global_summary <- walker_shares %>%
  summarise(across(all_of(share_cols), ~mean(.x, na.rm = TRUE)))

# -----------------------------
# Long share object for plotting
# -----------------------------
walker_share_long <- walker_shares %>%
  select(Date, regime, all_of(share_cols)) %>%
  pivot_longer(cols = all_of(share_cols), names_to = "series", values_to = "share") %>%
  mutate(
    measure = ifelse(grepl("^S_", series), "Signed", "Absolute"),
    block = sub("^(S_|A_)", "", series),
    model = "Walker adjusted"
  ) %>%
  select(Date, regime, measure, block, model, share)

# -----------------------------
# Diagnostics: do block sums + baseline reconstruct fitted value?
# -----------------------------
walker_reconstruction <- walker_block_signed %>%
  mutate(sum_blocks = rowSums(across(all_of(block_names)), na.rm = TRUE),
         fitted_from_parts = baseline_WalkerAdj + sum_blocks)

if (file.exists(walker_forecast_file)) {
  walker_forecasts <- read_csv(walker_forecast_file, show_col_types = FALSE) %>%
    mutate(Date = as.Date(Date))
  walker_reconstruction <- walker_reconstruction %>%
    left_join(walker_forecasts %>% select(Date, y_hat), by = "Date") %>%
    mutate(reconstruction_gap_vs_backtest = y_hat - fitted_from_parts)
}

# -----------------------------
# Save outputs
# -----------------------------
write_csv(walker_feature_long, "walker_harmonized_feature_long.csv")
write_csv(walker_block_signed, "walker_harmonized_block_signed.csv")
write_csv(walker_block_abs, "walker_harmonized_block_abs.csv")
write_csv(walker_shares, "walker_harmonized_shares.csv")
write_csv(walker_share_long, "walker_harmonized_share_long.csv")
write_csv(walker_regime_summary, "walker_harmonized_regime_summary.csv")
write_csv(walker_global_summary, "walker_harmonized_global_summary.csv")
write_csv(walker_reconstruction, "walker_harmonized_reconstruction_check.csv")
write_csv(block_map, "walker_feature_block_map.csv")

cat("Saved:\n")
cat("- walker_harmonized_feature_long.csv\n")
cat("- walker_harmonized_block_signed.csv\n")
cat("- walker_harmonized_block_abs.csv\n")
cat("- walker_harmonized_shares.csv\n")
cat("- walker_harmonized_share_long.csv\n")
cat("- walker_harmonized_regime_summary.csv\n")
cat("- walker_harmonized_global_summary.csv\n")
cat("- walker_harmonized_reconstruction_check.csv\n")
cat("- walker_feature_block_map.csv\n")

if ("reconstruction_gap_vs_backtest" %in% names(walker_reconstruction)) {
  cat("\nBacktest-vs-full-sample reconstruction gap summary:\n")
  print(summary(walker_reconstruction$reconstruction_gap_vs_backtest))
}

# ------------------------------------------------------------
# Comparing XGBoost SHAP and Walker
# ------------------------------------------------------------



# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
xgb_signed_file    <- "xgb_v8_harmonized_shap_group_signed.csv"
xgb_abs_file       <- "xgb_v8_harmonized_shap_group_abs.csv"
walker_signed_file <- "walker_harmonized_block_signed.csv"
walker_abs_file    <- "walker_harmonized_block_abs.csv"

stopifnot(file.exists(xgb_signed_file))
stopifnot(file.exists(xgb_abs_file))
stopifnot(file.exists(walker_signed_file))
stopifnot(file.exists(walker_abs_file))

# ------------------------------------------------------------
# Load
# ------------------------------------------------------------
xgb_signed    <- read_csv(xgb_signed_file, show_col_types = FALSE)
xgb_abs       <- read_csv(xgb_abs_file, show_col_types = FALSE)
walker_signed <- read_csv(walker_signed_file, show_col_types = FALSE)
walker_abs    <- read_csv(walker_abs_file, show_col_types = FALSE)

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

require_col <- function(df, candidates, label) {
  nm <- first_existing(df, candidates)
  if (is.na(nm)) stop(sprintf("Could not find %s. Tried: %s", label, paste(candidates, collapse = ", ")))
  nm
}

assign_regime <- function(Date) {
  case_when(
    Date < as.Date("2021-06-01") ~ "COVID (2020–2021)",
    Date < as.Date("2023-01-01") ~ "Energy Crisis",
    Date < as.Date("2024-06-01") ~ "Disinflation",
    TRUE ~ "Normalization"
  )
}

# ------------------------------------------------------------
# Date columns
# ------------------------------------------------------------
xgb_date_col <- require_col(xgb_signed, c("Date", "date"), "XGBoost date column")
walker_date_col <- require_col(walker_signed, c("Date", "date"), "Walker date column")

xgb_signed[[xgb_date_col]] <- as.Date(xgb_signed[[xgb_date_col]])
xgb_abs[[require_col(xgb_abs, c("Date", "date"), "XGBoost abs date column")]] <- as.Date(xgb_abs[[require_col(xgb_abs, c("Date", "date"), "XGBoost abs date column")]])
walker_signed[[walker_date_col]] <- as.Date(walker_signed[[walker_date_col]])
walker_abs[[require_col(walker_abs, c("Date", "date"), "Walker abs date column")]] <- as.Date(walker_abs[[require_col(walker_abs, c("Date", "date"), "Walker abs date column")]])

# ------------------------------------------------------------
# Column names from harmonized XGBoost notebook
# ------------------------------------------------------------
xgb_actual_col    <- require_col(xgb_signed, c("actual"), "XGBoost actual")
xgb_pred_col      <- require_col(xgb_signed, c("predicted_raw"), "XGBoost raw prediction")
xgb_base_col      <- require_col(xgb_signed, c("shap_base", "baseline_SHAP"), "XGBoost baseline")
xgb_regime_col    <- first_existing(xgb_signed, c("regime"))

xgb_block_cols <- c(
  AR = require_col(xgb_signed, c("AR (inflation)", "AR"), "XGBoost AR block"),
  PPI = require_col(xgb_signed, c("PPI / cost-push", "PPI"), "XGBoost PPI block"),
  Monetary = require_col(xgb_signed, c("Monetary policy", "Monetary"), "XGBoost Monetary block"),
  FX = require_col(xgb_signed, c("FX"), "XGBoost FX block"),
  Oil = require_col(xgb_signed, c("Oil"), "XGBoost Oil block"),
  Trade = require_col(xgb_signed, c("Trade"), "XGBoost Trade block"),
  Labour = require_col(xgb_signed, c("Labour market", "Labour"), "XGBoost Labour block")
)

# ------------------------------------------------------------
# Column names from Walker harmonization script
# ------------------------------------------------------------
walker_base_col <- require_col(walker_signed, c("baseline_WalkerAdj", "baseline", "intercept"), "Walker baseline")
walker_block_cols <- c(
  AR = require_col(walker_signed, c("AR"), "Walker AR block"),
  PPI = require_col(walker_signed, c("PPI"), "Walker PPI block"),
  Monetary = require_col(walker_signed, c("Monetary"), "Walker Monetary block"),
  FX = require_col(walker_signed, c("FX"), "Walker FX block"),
  Oil = require_col(walker_signed, c("Oil"), "Walker Oil block"),
  Trade = require_col(walker_signed, c("Trade"), "Walker Trade block"),
  Labour = require_col(walker_signed, c("Labour"), "Walker Labour block")
)

# ------------------------------------------------------------
# Standardized combined objects
# ------------------------------------------------------------
xgb_std <- tibble(
  Date = xgb_signed[[xgb_date_col]],
  actual = xgb_signed[[xgb_actual_col]],
  predicted_raw = xgb_signed[[xgb_pred_col]],
  baseline_SHAP = xgb_signed[[xgb_base_col]],
  regime_xgb = if (!is.na(xgb_regime_col)) xgb_signed[[xgb_regime_col]] else NA_character_,
  AR_SHAP = xgb_signed[[xgb_block_cols[["AR"]]]],
  PPI_SHAP = xgb_signed[[xgb_block_cols[["PPI"]]]],
  Monetary_SHAP = xgb_signed[[xgb_block_cols[["Monetary"]]]],
  FX_SHAP = xgb_signed[[xgb_block_cols[["FX"]]]],
  Oil_SHAP = xgb_signed[[xgb_block_cols[["Oil"]]]],
  Trade_SHAP = xgb_signed[[xgb_block_cols[["Trade"]]]],
  Labour_SHAP = xgb_signed[[xgb_block_cols[["Labour"]]]]
)

walker_std <- tibble(
  Date = walker_signed[[walker_date_col]],
  baseline_WalkerAdj = walker_signed[[walker_base_col]],
  AR_WalkerAdj = walker_signed[[walker_block_cols[["AR"]]]],
  PPI_WalkerAdj = walker_signed[[walker_block_cols[["PPI"]]]],
  Monetary_WalkerAdj = walker_signed[[walker_block_cols[["Monetary"]]]],
  FX_WalkerAdj = walker_signed[[walker_block_cols[["FX"]]]],
  Oil_WalkerAdj = walker_signed[[walker_block_cols[["Oil"]]]],
  Trade_WalkerAdj = walker_signed[[walker_block_cols[["Trade"]]]],
  Labour_WalkerAdj = walker_signed[[walker_block_cols[["Labour"]]]]
)

harmonised_attributions <- xgb_std %>%
  left_join(walker_std, by = "Date") %>%
  mutate(regime = assign_regime(Date)) %>%
  arrange(Date)

# ------------------------------------------------------------
# Additivity checks
# ------------------------------------------------------------
shap_block_names <- c("AR_SHAP", "PPI_SHAP", "Monetary_SHAP", "FX_SHAP", "Oil_SHAP", "Trade_SHAP", "Labour_SHAP")
walker_block_names <- c("AR_WalkerAdj", "PPI_WalkerAdj", "Monetary_WalkerAdj", "FX_WalkerAdj", "Oil_WalkerAdj", "Trade_WalkerAdj", "Labour_WalkerAdj")

xgb_check <- harmonised_attributions %>%
  mutate(
    pred_from_shap_blocks = baseline_SHAP + rowSums(across(all_of(shap_block_names)), na.rm = TRUE),
    diff = predicted_raw - pred_from_shap_blocks
  )

walker_check <- harmonised_attributions %>%
  mutate(
    fitted_from_parts = baseline_WalkerAdj + rowSums(across(all_of(walker_block_names)), na.rm = TRUE)
  )
sum(contrib_shares$total_abs_Walker)
sum(contrib_shares$total_abs_SHAP)

# ------------------------------------------------------------
# Signed and absolute shares
# ------------------------------------------------------------
contrib_shares <- harmonised_attributions %>%
  rowwise() %>%
  mutate(
    total_abs_SHAP = sum(abs(c_across(all_of(shap_block_names))), na.rm = TRUE),
    total_abs_Walker = sum(abs(c_across(all_of(walker_block_names))), na.rm = TRUE),

    S_AR_SHAP = if_else(total_abs_SHAP > 0, AR_SHAP / total_abs_SHAP, NA_real_),
    S_PPI_SHAP = if_else(total_abs_SHAP > 0, PPI_SHAP / total_abs_SHAP, NA_real_),
    S_Monetary_SHAP = if_else(total_abs_SHAP > 0, Monetary_SHAP / total_abs_SHAP, NA_real_),
    S_FX_SHAP = if_else(total_abs_SHAP > 0, FX_SHAP / total_abs_SHAP, NA_real_),
    S_Oil_SHAP = if_else(total_abs_SHAP > 0, Oil_SHAP / total_abs_SHAP, NA_real_),
    S_Trade_SHAP = if_else(total_abs_SHAP > 0, Trade_SHAP / total_abs_SHAP, NA_real_),
    S_Labour_SHAP = if_else(total_abs_SHAP > 0, Labour_SHAP / total_abs_SHAP, NA_real_),

    S_AR_Walker = if_else(total_abs_Walker > 0, AR_WalkerAdj / total_abs_Walker, NA_real_),
    S_PPI_Walker = if_else(total_abs_Walker > 0, PPI_WalkerAdj / total_abs_Walker, NA_real_),
    S_Monetary_Walker = if_else(total_abs_Walker > 0, Monetary_WalkerAdj / total_abs_Walker, NA_real_),
    S_FX_Walker = if_else(total_abs_Walker > 0, FX_WalkerAdj / total_abs_Walker, NA_real_),
    S_Oil_Walker = if_else(total_abs_Walker > 0, Oil_WalkerAdj / total_abs_Walker, NA_real_),
    S_Trade_Walker = if_else(total_abs_Walker > 0, Trade_WalkerAdj / total_abs_Walker, NA_real_),
    S_Labour_Walker = if_else(total_abs_Walker > 0, Labour_WalkerAdj / total_abs_Walker, NA_real_),

    A_AR_SHAP = if_else(total_abs_SHAP > 0, abs(AR_SHAP) / total_abs_SHAP, NA_real_),
    A_PPI_SHAP = if_else(total_abs_SHAP > 0, abs(PPI_SHAP) / total_abs_SHAP, NA_real_),
    A_Monetary_SHAP = if_else(total_abs_SHAP > 0, abs(Monetary_SHAP) / total_abs_SHAP, NA_real_),
    A_FX_SHAP = if_else(total_abs_SHAP > 0, abs(FX_SHAP) / total_abs_SHAP, NA_real_),
    A_Oil_SHAP = if_else(total_abs_SHAP > 0, abs(Oil_SHAP) / total_abs_SHAP, NA_real_),
    A_Trade_SHAP = if_else(total_abs_SHAP > 0, abs(Trade_SHAP) / total_abs_SHAP, NA_real_),
    A_Labour_SHAP = if_else(total_abs_SHAP > 0, abs(Labour_SHAP) / total_abs_SHAP, NA_real_),

    A_AR_Walker = if_else(total_abs_Walker > 0, abs(AR_WalkerAdj) / total_abs_Walker, NA_real_),
    A_PPI_Walker = if_else(total_abs_Walker > 0, abs(PPI_WalkerAdj) / total_abs_Walker, NA_real_),
    A_Monetary_Walker = if_else(total_abs_Walker > 0, abs(Monetary_WalkerAdj) / total_abs_Walker, NA_real_),
    A_FX_Walker = if_else(total_abs_Walker > 0, abs(FX_WalkerAdj) / total_abs_Walker, NA_real_),
    A_Oil_Walker = if_else(total_abs_Walker > 0, abs(Oil_WalkerAdj) / total_abs_Walker, NA_real_),
    A_Trade_Walker = if_else(total_abs_Walker > 0, abs(Trade_WalkerAdj) / total_abs_Walker, NA_real_),
    A_Labour_Walker = if_else(total_abs_Walker > 0, abs(Labour_WalkerAdj) / total_abs_Walker, NA_real_)
  ) %>%
  ungroup()

# ------------------------------------------------------------
# Regime summary
# ------------------------------------------------------------
share_cols <- grep("^(S|A)_", names(contrib_shares), value = TRUE)
regime_summary <- contrib_shares %>%
  group_by(regime) %>%
  summarise(across(all_of(share_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# ------------------------------------------------------------
# Sign agreement
# ------------------------------------------------------------
sign_agreement <- harmonised_attributions %>%
  mutate(
    agree_AR = sign(AR_SHAP) == sign(AR_WalkerAdj),
    agree_PPI = sign(PPI_SHAP) == sign(PPI_WalkerAdj),
    agree_Monetary = sign(Monetary_SHAP) == sign(Monetary_WalkerAdj),
    agree_FX = sign(FX_SHAP) == sign(FX_WalkerAdj),
    agree_Oil = sign(Oil_SHAP) == sign(Oil_WalkerAdj),
    agree_Trade = sign(Trade_SHAP) == sign(Trade_WalkerAdj),
    agree_Labour = sign(Labour_SHAP) == sign(Labour_WalkerAdj)
  )

sign_summary <- sign_agreement %>%
  summarise(
    pct_agree_AR = mean(agree_AR, na.rm = TRUE),
    pct_agree_PPI = mean(agree_PPI, na.rm = TRUE),
    pct_agree_Monetary = mean(agree_Monetary, na.rm = TRUE),
    pct_agree_FX = mean(agree_FX, na.rm = TRUE),
    pct_agree_Oil = mean(agree_Oil, na.rm = TRUE),
    pct_agree_Trade = mean(agree_Trade, na.rm = TRUE),
    pct_agree_Labour = mean(agree_Labour, na.rm = TRUE)
  )

# ------------------------------------------------------------
# Routing compare: correlation of signed and absolute contributions
# ------------------------------------------------------------
routing_compare <- tibble(
  block = c("AR", "PPI", "Monetary", "FX", "Oil", "Trade", "Labour"),
  corr_signed = c(
    cor(harmonised_attributions$AR_SHAP, harmonised_attributions$AR_WalkerAdj, use = "complete.obs"),
    cor(harmonised_attributions$PPI_SHAP, harmonised_attributions$PPI_WalkerAdj, use = "complete.obs"),
    cor(harmonised_attributions$Monetary_SHAP, harmonised_attributions$Monetary_WalkerAdj, use = "complete.obs"),
    cor(harmonised_attributions$FX_SHAP, harmonised_attributions$FX_WalkerAdj, use = "complete.obs"),
    cor(harmonised_attributions$Oil_SHAP, harmonised_attributions$Oil_WalkerAdj, use = "complete.obs"),
    cor(harmonised_attributions$Trade_SHAP, harmonised_attributions$Trade_WalkerAdj, use = "complete.obs"),
    cor(harmonised_attributions$Labour_SHAP, harmonised_attributions$Labour_WalkerAdj, use = "complete.obs")
  ),
  corr_abs = c(
    cor(abs(harmonised_attributions$AR_SHAP), abs(harmonised_attributions$AR_WalkerAdj), use = "complete.obs"),
    cor(abs(harmonised_attributions$PPI_SHAP), abs(harmonised_attributions$PPI_WalkerAdj), use = "complete.obs"),
    cor(abs(harmonised_attributions$Monetary_SHAP), abs(harmonised_attributions$Monetary_WalkerAdj), use = "complete.obs"),
    cor(abs(harmonised_attributions$FX_SHAP), abs(harmonised_attributions$FX_WalkerAdj), use = "complete.obs"),
    cor(abs(harmonised_attributions$Oil_SHAP), abs(harmonised_attributions$Oil_WalkerAdj), use = "complete.obs"),
    cor(abs(harmonised_attributions$Trade_SHAP), abs(harmonised_attributions$Trade_WalkerAdj), use = "complete.obs"),
    cor(abs(harmonised_attributions$Labour_SHAP), abs(harmonised_attributions$Labour_WalkerAdj), use = "complete.obs")
  )
)

# ------------------------------------------------------------
# Long objects for plots
# ------------------------------------------------------------
comparison_df <- harmonised_attributions %>%
  select(Date,
         AR_SHAP, PPI_SHAP, Monetary_SHAP, FX_SHAP, Oil_SHAP, Trade_SHAP, Labour_SHAP,
         AR_WalkerAdj, PPI_WalkerAdj, Monetary_WalkerAdj, FX_WalkerAdj, Oil_WalkerAdj, Trade_WalkerAdj, Labour_WalkerAdj) %>%
  pivot_longer(cols = -Date, names_to = "series", values_to = "value") %>%
  mutate(
    model = case_when(
      grepl("_SHAP$", series) ~ "SHAP",
      grepl("_WalkerAdj$", series) ~ "Walker adjusted",
      TRUE ~ NA_character_
    ),
    component = gsub("_(SHAP|WalkerAdj)$", "", series)
  )

share_df <- contrib_shares %>%
  select(Date, regime, starts_with("S_"), starts_with("A_")) %>%
  pivot_longer(
    cols = -c(Date, regime),
    names_to = c("measure", "block", "model"),
    names_pattern = "(S|A)_([^_]+)_(SHAP|Walker)",
    values_to = "share"
  ) %>%
  mutate(
    measure = if_else(measure == "S", "Signed", "Absolute"),
    model = if_else(model == "SHAP", "SHAP", "Walker adjusted")
  )


# --- Thesis theme -----------------------------------------------------------
theme_thesis <- theme_minimal(base_size = 11, base_family = "serif") +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0),
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

# --- p1: Contribution levels ------------------------------------------------
p1 <- ggplot(comparison_df, aes(x = Date, y = value, linetype = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_line(linewidth = 0.45) +
  facet_wrap(~ component, scales = "free_y", ncol = 2) +
  scale_linetype_manual(values = c("SHAP" = "solid", "Walker adjusted" = "dashed")) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title    = "SHAP versus adjusted Walker decomposition",
    subtitle = "Block-level contributions to predicted inflation, 2020–2025",
    x = NULL, y = "Contribution (pp)", linetype = NULL
  ) +
  theme_thesis

# --- p2: Attribution shares --------------------------------------------------
p2 <- ggplot(share_df, aes(x = Date, y = share, linetype = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_line(linewidth = 0.45) +
  facet_grid(measure ~ block, scales = "free_y") +
  scale_linetype_manual(values = c("SHAP" = "solid", "Walker adjusted" = "dashed")) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  labs(
    title    = "Attribution shares: signed and absolute",
    subtitle = "Share of total absolute contribution per block",
    x = NULL, y = "Share", linetype = NULL
  ) +
  theme_thesis

p1
p2
?walker
# ------------------------------------------------------------
# Save
# ------------------------------------------------------------
write_csv(harmonised_attributions, "comparison_harmonised_attributions.csv")
write_csv(contrib_shares, "comparison_contrib_shares.csv")
write_csv(regime_summary, "comparison_regime_summary.csv")
write_csv(sign_summary, "comparison_sign_summary.csv")
write_csv(routing_compare, "comparison_routing_compare.csv")
write_csv(xgb_check, "comparison_xgb_additivity_check.csv")
write_csv(walker_check, "comparison_walker_additivity_check.csv")
write_csv(comparison_df, "comparison_long_contributions.csv")
write_csv(share_df, "comparison_long_shares.csv")

ggsave("comparison_plot_contributions.png", p1, width = 13, height = 8, dpi = 150)
ggsave("comparison_plot_shares.png", p2, width = 14, height = 8, dpi = 150)

cat("Saved comparison outputs.\n")
print(sign_summary)
print(routing_compare)
print(summary(xgb_check$diff))
print(p1)
print(p2)

# --- Table 1: Agreement overview ---
tab1_df <- routing_compare %>%
  arrange(desc(corr_signed)) %>%
  mutate(
    corr_signed = sprintf("%.3f", corr_signed),
    corr_abs = sprintf("%.3f", corr_abs)
  )

colnames(tab1_df) <- c("Block", "Signed r", "Absolute r")

tab1_df %>%
  kbl(
    col.names = c("Block", "Signed correlation", "Absolute correlation"),
    digits = 3, booktabs = TRUE,
    caption = "Correlation between SHAP and Walker adjusted contributions by block"
  ) %>%
  kable_styling(latex_options = c("hold_position"))

# --- Table 2: Mean attribution shares ---
tab2 <- share_df %>%
  group_by(block, measure, model) %>%
  summarise(mean_share = mean(share, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = c(measure, model), values_from = mean_share) %>%
  select(block, starts_with("Absolute"), starts_with("Signed"))

tab2 %>%
  kbl(
    col.names = c("Block", "SHAP", "Walker", "SHAP", "Walker"),
    digits = 3, booktabs = TRUE,
    caption = "Mean attribution shares by block and method"
  ) %>%
  add_header_above(c(" " = 1, "Absolute share" = 2, "Signed share" = 2)) %>%
  kable_styling(latex_options = c("hold_position"))

# Be om Table 5-tilsvarende med ny data
ws <- read.csv("walker_harmonized_regime_summary.csv")
print(ws)

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------
pred_df <- harmonised_attributions %>%
  select(Date, actual, xgb = predicted_raw) %>%
  left_join(
     <-  %>% select(Date, walker = y_hat),
    by = "Date"
  ) %>%
  filter(!is.na(actual), !is.na(xgb), !is.na(walker))

rmse_xgb    <- sqrt(mean((pred_df$actual - pred_df$xgb)^2))
rmse_walker <- sqrt(mean((pred_df$actual - pred_df$walker)^2))

pred_long <- pred_df %>%
  pivot_longer(cols = c(actual, xgb, walker), names_to = "series", values_to = "value") %>%
  mutate(series = factor(series,
                         levels = c("actual", "xgb", "walker"),
                         labels = c("Actual", "XGBoost", "Walker TVP")
  ))

rmse_label <- sprintf("RMSE:  XGBoost = %.3f,  Walker = %.3f", rmse_xgb, rmse_walker)

ggplot(pred_long, aes(x = Date, y = value, color = series, linetype = series)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = c("Actual" = "black", "XGBoost" = "#08306B", "Walker TVP" = "#D62728")) +
  scale_linetype_manual(values = c("Actual" = "solid", "XGBoost" = "dashed", "Walker TVP" = "dotted")) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  annotate("text", x = min(pred_df$Date) + 30, y = max(pred_df$actual) * 0.95,
           label = rmse_label, hjust = 0, size = 3, family = "serif") +
  labs(
    title = "Forecasts: XGBoost versus Walker TVP",
    x = NULL, y = "YoY inflation (%)", color = NULL, linetype = NULL
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 0.8))) +
  theme_thesis

########

# -----------------------------
# Build table data
# -----------------------------
table_abs_kable_data <- contrib_shares %>%
  group_by(regime) %>%
  summarise(
    A_AR_SHAP       = mean(A_AR_SHAP, na.rm = TRUE),
    A_PPI_SHAP      = mean(A_PPI_SHAP, na.rm = TRUE),
    A_Monetary_SHAP = mean(A_Monetary_SHAP, na.rm = TRUE),
    A_FX_SHAP       = mean(A_FX_SHAP, na.rm = TRUE),
    A_Oil_SHAP      = mean(A_Oil_SHAP, na.rm = TRUE),
    A_Trade_SHAP    = mean(A_Trade_SHAP, na.rm = TRUE),
    A_Labour_SHAP   = mean(A_Labour_SHAP, na.rm = TRUE),

    A_AR_Walker       = mean(A_AR_Walker, na.rm = TRUE),
    A_PPI_Walker      = mean(A_PPI_Walker, na.rm = TRUE),
    A_Monetary_Walker = mean(A_Monetary_Walker, na.rm = TRUE),
    A_FX_Walker       = mean(A_FX_Walker, na.rm = TRUE),
    A_Oil_Walker      = mean(A_Oil_Walker, na.rm = TRUE),
    A_Trade_Walker    = mean(A_Trade_Walker, na.rm = TRUE),
    A_Labour_Walker   = mean(A_Labour_Walker, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -regime,
    names_to = c("share_type", "block", "model"),
    names_pattern = "^(A)_(AR|PPI|Monetary|FX|Oil|Trade|Labour)_(SHAP|Walker)$",
    values_to = "value"
  ) %>%
  mutate(
    value = value * 100,
    model = case_when(
      model == "SHAP"   ~ "XGBoost (SHAP)",
      model == "Walker" ~ "Walker"
    ),
    block = case_when(
      block == "AR"       ~ "AR (inflation)",
      block == "PPI"      ~ "PPI / cost-push",
      block == "Monetary" ~ "Monetary policy",
      block == "FX"       ~ "FX",
      block == "Oil"      ~ "Oil",
      block == "Trade"    ~ "Trade",
      block == "Labour"   ~ "Labour market"
    ),
    regime = case_when(
      regime == "COVID (2020–2021)" ~ "COVID",
      regime == "Energy Crisis"     ~ "Energy",
      regime == "Disinflation"      ~ "Disinfl.",
      regime == "Normalization"     ~ "Normal"
    )
  ) %>%
  select(model, block, regime, value) %>%
  pivot_wider(names_from = regime, values_from = value) %>%
  mutate(
    model = factor(model, levels = c("XGBoost (SHAP)", "Walker")),
    block = factor(
      block,
      levels = c(
        "AR (inflation)",
        "PPI / cost-push",
        "Monetary policy",
        "FX",
        "Oil",
        "Trade",
        "Labour market"
      )
    )
  ) %>%
  arrange(model, block)

# Keep only the printed columns
table_abs_print <- table_abs_kable_data %>%
  select(model, block, COVID, Energy, `Disinfl.`, Normal)

# -----------------------------
# KableExtra table
# -----------------------------
kbl(
  table_abs_print %>% select(-model),
  col.names = c("", "COVID", "Energy", "Disinfl.", "Normal"),
  digits = 1,
  align = c("l", "c", "c", "c", "c"),
  caption = "Block importance by regime (% of total absolute contribution)",
  booktabs = TRUE,
  linesep = ""
) %>%
  kable_styling(
    full_width = FALSE,
    position = "center",
    font_size = 13,
    bootstrap_options = c("striped", "condensed")
  ) %>%
  group_rows("XGBoost (SHAP)", 1, 7, italic = TRUE, bold = FALSE) %>%
  group_rows("Walker", 8, 14, italic = TRUE, bold = FALSE) %>%
  row_spec(0, bold = FALSE) %>%
  footnote(
    general = "Entries are mean absolute contribution shares within each regime. Shares sum to 100 within model and regime, up to rounding.",
    general_title = "Note: "
  )


# =============================================================================
# PART 5: Harmonised attribution residual (Øysteins research note)
# Delta_kt = C_kt^XGB - C_kt^TVP
# =============================================================================

# -----------------------------
# Step 1: Compute block-level residuals at each time point
# -----------------------------
residual_df <- harmonised_attributions %>%
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

# -----------------------------
# Step 2: Long format for plotting and summarising
# -----------------------------
residual_long <- residual_df %>%
  pivot_longer(
    cols = starts_with("Delta_"),
    names_to = "block",
    names_prefix = "Delta_",
    values_to = "Delta"
  ) %>%
  mutate(
    block = factor(block, 
                   levels = c("AR", "PPI", "Monetary", "FX", "Oil", "Trade", "Labour"))
  )

# -----------------------------
# Step 3: Regime-conditional summaries
#   var_delta  = Var(Delta_kt | regime)
#   mean_abs   = E(|Delta_kt| | regime)
# -----------------------------
residual_metrics_long <- residual_long %>%
  group_by(block, regime) %>%
  summarise(
    var_delta  = var(Delta, na.rm = TRUE),
    mean_abs   = mean(abs(Delta), na.rm = TRUE),
    n_obs      = sum(!is.na(Delta)),
    .groups    = "drop"
  ) %>%
  mutate(
    regime = factor(regime,
                    levels = c("COVID (2020–2021)", "Energy Crisis",
                               "Disinflation", "Normalization"),
                    labels = c("COVID", "Energy", "Disinfl.", "Normal"))
  ) %>%
  arrange(block, regime)

print(residual_metrics_long)

# -----------------------------
# Step 4: Wide presentation tables (one for variance, one for mean abs)
# -----------------------------
var_wide <- residual_metrics_long %>%
  select(block, regime, var_delta) %>%
  pivot_wider(names_from = regime, values_from = var_delta) %>%
  mutate(statistic = "Var(Delta)") %>%
  select(statistic, block, COVID, Energy, `Disinfl.`, Normal)

mean_abs_wide <- residual_metrics_long %>%
  select(block, regime, mean_abs) %>%
  pivot_wider(names_from = regime, values_from = mean_abs) %>%
  mutate(statistic = "E(|Delta|)") %>%
  select(statistic, block, COVID, Energy, `Disinfl.`, Normal)

residual_table_combined <- bind_rows(var_wide, mean_abs_wide)

print(residual_table_combined)

# -----------------------------
# Step 5: LaTeX table (matches the structure in the thesis)
# -----------------------------
residual_table_combined %>%
  select(-statistic) %>%
  kbl(
    col.names = c("Block", "COVID", "Energy", "Disinfl.", "Normal"),
    digits = 3,
    booktabs = TRUE,
    align = c("l", "c", "c", "c", "c"),
    caption = "Regime-conditional residual summaries, $\\Delta_{kt} = C^{\\text{XGB}}_{kt} - C^{\\text{TVP}}_{kt}$",
    linesep = ""
  ) %>%
  kable_styling(
    full_width = FALSE,
    position = "center",
    font_size = 11,
    latex_options = c("hold_position")
  ) %>%
  pack_rows("$\\mathrm{Var}(\\Delta_{kt} \\mid r)$", 1, 7, 
            italic = TRUE, escape = FALSE) %>%
  pack_rows("$\\mathbb{E}(|\\Delta_{kt}| \\mid r)$", 8, 14, 
            italic = TRUE, escape = FALSE) %>%
  footnote(
    general = "Entries are computed from the harmonised block-level contributions on the out-of-sample backtesting period. Variance captures the volatility of architectural disagreement within each regime; mean absolute magnitude captures its typical size in percentage points.",
    general_title = "Note: ",
    threeparttable = TRUE,
    escape = FALSE
  )

# -----------------------------
# Step 6: Visualisation - residual over time, faceted by block, coloured by regime
# -----------------------------
regime_colours <- c(
  "COVID (2020–2021)" = "#5E81AC",
  "Energy Crisis"     = "#BF616A",
  "Disinflation"      = "#A3BE8C",
  "Normalization"     = "#B48EAD"
)

p_residual <- ggplot(
  residual_long,
  aes(x = Date, y = Delta, fill = regime)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_col(alpha = 0.85, width = 25) +
  facet_wrap(~ block, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = regime_colours) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(
    title = "Harmonised attribution residual by block",
    subtitle = expression(Delta[kt] ~ "= " ~ C[kt]^{XGB} - C[kt]^{TVP}),
    x = NULL,
    y = "Residual (pp)",
    fill = "Regime"
  ) +
  theme_thesis +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "italic", size = 10),
    panel.spacing = unit(0.8, "lines")
  )

print(p_residual)

# -----------------------------
# Step 7: Save outputs
# -----------------------------
write_csv(residual_df, "residual_wide.csv")
write_csv(residual_long, "residual_long.csv")
write_csv(residual_metrics_long, "residual_metrics_long.csv")
write_csv(residual_table_combined, "residual_table_combined.csv")
ggsave("residual_plot.png", p_residual, width = 12, height = 9, dpi = 150)

cat("\nSaved residual outputs:\n")
cat("- residual_wide.csv (Delta per block, wide)\n")
cat("- residual_long.csv (Delta per block, long)\n")
cat("- residual_metrics_long.csv (Var and E|Delta| per block-regime)\n")
cat("- residual_table_combined.csv (presentation table)\n")
cat("- residual_plot.png (visualisation)\n")



# -----------------------------
# Robustness Check
# -----------------------------
# =============================================================================
# Robustness check: Regime partition sensitivity (±3 month shift)
# Tests whether block importance patterns in Tabell 5 are robust to the
# specific choice of regime cutoffs.
# =============================================================================

# -----------------------------
# Step 1: Define three alternative regime partitions
# -----------------------------
# Original cutoffs (from assign_regime function above):
#   COVID:        < 2021-06-01
#   Energy:       2021-06-01 to < 2023-01-01
#   Disinflation: 2023-01-01 to < 2024-06-01
#   Normal:       >= 2024-06-01

assign_regime_shifted <- function(Date, shift_months = 0) {
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

# -----------------------------
# Step 2: Compute block shares per regime under each partition
# -----------------------------
compute_regime_shares <- function(data, partition_label, shift_months) {
  data %>%
    mutate(regime_alt = assign_regime_shifted(Date, shift_months)) %>%
    rowwise() %>%
    mutate(
      total_abs_SHAP = sum(abs(c_across(all_of(shap_block_names))), na.rm = TRUE),
      total_abs_Walker = sum(abs(c_across(all_of(walker_block_names))), na.rm = TRUE),
      
      A_AR_SHAP = abs(AR_SHAP) / total_abs_SHAP,
      A_PPI_SHAP = abs(PPI_SHAP) / total_abs_SHAP,
      A_Monetary_SHAP = abs(Monetary_SHAP) / total_abs_SHAP,
      A_FX_SHAP = abs(FX_SHAP) / total_abs_SHAP,
      A_Oil_SHAP = abs(Oil_SHAP) / total_abs_SHAP,
      A_Trade_SHAP = abs(Trade_SHAP) / total_abs_SHAP,
      A_Labour_SHAP = abs(Labour_SHAP) / total_abs_SHAP,
      
      A_AR_Walker = abs(AR_WalkerAdj) / total_abs_Walker,
      A_PPI_Walker = abs(PPI_WalkerAdj) / total_abs_Walker,
      A_Monetary_Walker = abs(Monetary_WalkerAdj) / total_abs_Walker,
      A_FX_Walker = abs(FX_WalkerAdj) / total_abs_Walker,
      A_Oil_Walker = abs(Oil_WalkerAdj) / total_abs_Walker,
      A_Trade_Walker = abs(Trade_WalkerAdj) / total_abs_Walker,
      A_Labour_Walker = abs(Labour_WalkerAdj) / total_abs_Walker
    ) %>%
    ungroup() %>%
    group_by(regime_alt) %>%
    summarise(
      across(starts_with("A_"), ~ mean(.x, na.rm = TRUE) * 100),
      n_obs = n(),
      .groups = "drop"
    ) %>%
    mutate(partition = partition_label, shift_months = shift_months) %>%
    select(partition, shift_months, regime = regime_alt, n_obs, everything())
}

# -----------------------------
# Step 3: Run for all three partitions
# -----------------------------
robustness_results <- bind_rows(
  compute_regime_shares(harmonised_attributions, "Shift -3m", -3),
  compute_regime_shares(harmonised_attributions, "Original",   0),
  compute_regime_shares(harmonised_attributions, "Shift +3m", +3)
)

# -----------------------------
# Step 4: Pivot to long format for clean comparison table
# -----------------------------
robustness_long <- robustness_results %>%
  pivot_longer(
    cols = starts_with("A_"),
    names_to = c("block", "model"),
    names_pattern = "^A_([^_]+)_(SHAP|Walker)$",
    values_to = "share"
  ) %>%
  mutate(
    block = factor(block, levels = c("AR", "PPI", "Monetary", "FX", "Oil", "Trade", "Labour")),
    regime = factor(regime,
                    levels = c("COVID (2020-2021)", "Energy Crisis", "Disinflation", "Normalization"),
                    labels = c("COVID", "Energy", "Disinfl.", "Normal")),
    partition = factor(partition, levels = c("Shift -3m", "Original", "Shift +3m"))
  )

# -----------------------------
# Step 5: Wide table - one row per (model, block, regime), columns for each partition
# -----------------------------
robustness_table <- robustness_long %>%
  select(model, block, regime, partition, share) %>%
  pivot_wider(names_from = partition, values_from = share) %>%
  mutate(
    `Diff (max-min)` = pmax(`Shift -3m`, `Original`, `Shift +3m`, na.rm = TRUE) -
      pmin(`Shift -3m`, `Original`, `Shift +3m`, na.rm = TRUE)
  ) %>%
  arrange(model, block, regime)

print(robustness_table, n = Inf)

# -----------------------------
# Step 6: Summary - max absolute change across partitions per (model, regime)
# -----------------------------
robustness_summary <- robustness_long %>%
  group_by(model, block, regime) %>%
  summarise(
    max_diff_pp = max(share, na.rm = TRUE) - min(share, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(model, regime) %>%
  summarise(
    max_block_diff_pp = max(max_diff_pp, na.rm = TRUE),
    mean_block_diff_pp = mean(max_diff_pp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(model, regime)

cat("\n=== Sensitivity summary ===\n")
cat("Max difference (in percentage points) in any block share within each regime,\n")
cat("when shifting cutoffs by +/- 3 months:\n\n")
print(robustness_summary)

# -----------------------------
# Step 7: LaTeX-ready table for appendix
# -----------------------------
robustness_table_latex <- robustness_table %>%
  mutate(
    model = case_when(
      model == "SHAP"   ~ "XGBoost (SHAP)",
      model == "Walker" ~ "Walker"
    )
  ) %>%
  select(model, block, regime, `Shift -3m`, Original, `Shift +3m`)

kbl(
  robustness_table_latex %>% select(-model),
  col.names = c("Block", "Regime", "Shift -3m", "Original", "Shift +3m"),
  digits = 1,
  align = c("l", "l", "c", "c", "c"),
  caption = "Robustness check: block importance shares (\\%) under alternative regime partitions (cutoffs shifted by +/- 3 months)",
  booktabs = TRUE,
  linesep = ""
) %>%
  kable_styling(
    full_width = FALSE,
    position = "center",
    font_size = 10,
    latex_options = c("hold_position", "scale_down")
  ) %>%
  pack_rows("XGBoost (SHAP)", 1, 28, italic = TRUE) %>%
  pack_rows("Walker", 29, 56, italic = TRUE) %>%
  footnote(
    general = "Entries are mean absolute contribution shares within each regime, expressed in percentage points. The Original column reproduces the cutoffs used in Table 5; Shift -3m and Shift +3m move all cutoffs three months earlier or later. Hovedfunn: PPI-konvergens i Energy-regimet og Walker monetary-policy-dominans i Normal er kvalitativt uendret.",
    general_title = "Note: ",
    threeparttable = TRUE
  )

# -----------------------------
# Step 8: Save outputs
# -----------------------------
write_csv(robustness_results, "robustness_results_wide.csv")
write_csv(robustness_long,    "robustness_results_long.csv")
write_csv(robustness_table,   "robustness_table.csv")
write_csv(robustness_summary, "robustness_summary.csv")

cat("\nSaved robustness outputs:\n")
cat("- robustness_results_wide.csv\n")
cat("- robustness_results_long.csv\n")
cat("- robustness_table.csv (sammenligning per block/regime)\n")
cat("- robustness_summary.csv (max diff per regime)\n")





######################
# Henter eval metrics
######################
safe_mape <- function(actual, predicted) {
  ok <- is.finite(actual) & is.finite(predicted) & actual != 0
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

walker_results <- read_csv("walker_forecasts.csv", show_col_types = FALSE)

walker_accuracy_table <- bind_rows(
  forecast_metrics(walker_results$y_actual, walker_results$y_rw,  "Random walk"),
  forecast_metrics(walker_results$y_actual, walker_results$y_hat, "Walker")
)

print(walker_accuracy_table)

# AR(1) forecasts
# Prosess som funker hvis man ikke har kjørt hele kodefilen. 

df <- read_csv("master_data.csv", show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date)) %>%
  arrange(Date)

# Same target/lag structure as Walker file
model_df <- df %>%
  arrange(Date) %>%
  mutate(
    kpi_yoy = 100 * (kpi / lag(kpi, 12) - 1),
    kpi_yoy_lag1 = lag(kpi_yoy, 1),
    target = lead(kpi_yoy, H)
  ) %>%
  filter(!is.na(target), !is.na(kpi_yoy_lag1))

test_origins <- model_df %>%
  filter(Date >= TEST_START) %>%
  pull(Date)

ar1_results <- data.frame(
  Date = as.Date(character()),
  y_ar1 = numeric()
)

for (origin in test_origins) {
  train_raw <- model_df %>% filter(Date < origin)
  origin_row <- model_df %>% filter(Date == origin)
  
  if (nrow(train_raw) < MIN_TRAIN_N || nrow(origin_row) != 1) next
  
  ar1_fit <- lm(target ~ kpi_yoy_lag1, data = train_raw)
  y_ar1 <- as.numeric(predict(ar1_fit, newdata = origin_row))
  
  ar1_results <- rbind(
    ar1_results,
    data.frame(Date = as.Date(origin, origin = "1970-01-01"), y_ar1 = y_ar1)  )
}
ar1_results <- ar1_results %>%
  mutate(Date = as.Date(Date, origin = "1970-01-01"))
print(ar1_results)

write_csv(ar1_results, "ar1_forecasts.csv")
