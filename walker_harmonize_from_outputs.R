library(dplyr)
library(readr)
library(tidyr)

# Robust working directory block
try(setwd("/Users/amundbech1/Library/CloudStorage/OneDrive-UiTOffice365/Data_master/MasterOppgave"), silent = TRUE)
try(setwd("~/Library/CloudStorage/OneDrive-UiTOffice365/Amund Skjalgsønn Bech's files - Data_master/MasterOppgave"), silent = TRUE)

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
# Important note:
# The full-sample Walker fit used features scaled with full-sample means/sds.
# That means x has mean ~0 in the fitted sample, so beta*x is already the
# centred contribution object from the research note:
#   c_tilde_jt = (x_jt - xbar_j) * beta_jt
# because xbar_j ~= 0 after scaling.
# Thus baseline_WalkerAdj is approximately the intercept process.
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
# Shares per the research note
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
      Date < as.Date("2020-03-01") ~ "Pre-shock",
      Date < as.Date("2021-09-01") ~ "Pandemic",
      Date < as.Date("2023-06-01") ~ "Inflation shock",
      Date < as.Date("2024-06-01") ~ "Disinflation",
      TRUE ~ "Recent"
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
