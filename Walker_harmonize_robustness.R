# =============================================================================
# Walker Robustness Harmonization
# =============================================================================
# Takes the output from Walker_robustness.R and computes block-level
# attribution objects, shares, and regime summaries for each spec.
#
# Run this once per spec by changing SPEC at the top.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

SPEC <- "usd_external"  # "baseline", "no_intercept", or "usd_external"
OUT_PREFIX <- paste0("walker_", SPEC, "_")

# --- Load inputs -------------------------------------------------------------
walker_contrib <- read_csv(paste0(OUT_PREFIX, "contributions.csv"),
                           show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date))

# --- Block map ---------------------------------------------------------------
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

# Ensure intercept column exists (no_intercept spec sets it to 0)
if (!"intercept" %in% names(walker_contrib)) walker_contrib$intercept <- 0

# --- Feature-level long ------------------------------------------------------
walker_feature_long <- walker_contrib %>%
  select(Date, all_of(block_map$feature), intercept) %>%
  pivot_longer(cols = all_of(block_map$feature),
               names_to = "feature", values_to = "contribution") %>%
  left_join(block_map, by = "feature")

# --- Block-level signed ------------------------------------------------------
walker_block_signed <- walker_feature_long %>%
  group_by(Date, block) %>%
  summarise(contribution = sum(contribution, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = block, values_from = contribution) %>%
  left_join(walker_contrib %>% select(Date, intercept), by = "Date") %>%
  rename(baseline_WalkerAdj = intercept) %>%
  arrange(Date)

# --- Block-level absolute ----------------------------------------------------
walker_block_abs <- walker_feature_long %>%
  group_by(Date, block) %>%
  summarise(abs_contribution = sum(abs(contribution), na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = block, values_from = abs_contribution) %>%
  arrange(Date)

# --- Regime assignment -------------------------------------------------------
# Set of alternative regime partitions. Change here to test sensitivity.
REGIME_VARIANT <- "main"  # "main" or "alt1" or "alt2"

assign_regime <- function(Date, variant = "main") {
  if (variant == "main") {
    case_when(
      Date < as.Date("2021-06-01") ~ "COVID (2020-2021)",
      Date < as.Date("2023-01-01") ~ "Energy Crisis",
      Date < as.Date("2024-06-01") ~ "Disinflation",
      TRUE ~ "Normalization"
    )
  } else if (variant == "alt1") {
    # Earlier cutoff for end of COVID (match Norges Bank rate hike, Sept 2021)
    case_when(
      Date < as.Date("2021-09-01") ~ "COVID (2020-2021)",
      Date < as.Date("2023-06-01") ~ "Energy Crisis",
      Date < as.Date("2024-09-01") ~ "Disinflation",
      TRUE ~ "Normalization"
    )
  } else if (variant == "alt2") {
    # Based on headline CPI threshold (rough)
    case_when(
      Date < as.Date("2021-03-01") ~ "COVID (2020-2021)",
      Date < as.Date("2023-03-01") ~ "Energy Crisis",
      Date < as.Date("2024-03-01") ~ "Disinflation",
      TRUE ~ "Normalization"
    )
  }
}

walker_block_signed <- walker_block_signed %>%
  mutate(regime = assign_regime(Date, REGIME_VARIANT))
walker_block_abs <- walker_block_abs %>%
  mutate(regime = assign_regime(Date, REGIME_VARIANT))

# --- Shares ------------------------------------------------------------------
block_names <- c("AR", "PPI", "Monetary", "FX", "Oil", "Trade", "Labour")

walker_shares <- walker_block_signed %>%
  rowwise() %>%
  mutate(total_abs = sum(abs(c_across(all_of(block_names))), na.rm = TRUE)) %>%
  ungroup()

for (blk in block_names) {
  walker_shares[[paste0("S_", blk)]] <- ifelse(
    walker_shares$total_abs > 0,
    walker_shares[[blk]] / walker_shares$total_abs, NA_real_
  )
  walker_shares[[paste0("A_", blk)]] <- ifelse(
    walker_shares$total_abs > 0,
    abs(walker_shares[[blk]]) / walker_shares$total_abs, NA_real_
  )
}



# --- Regime summary (the key output for comparison) --------------------------
share_cols <- grep("^(S_|A_)", names(walker_shares), value = TRUE)

regime_summary <- walker_shares %>%
  group_by(regime) %>%
  summarise(
    across(all_of(share_cols), ~mean(.x, na.rm = TRUE)),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(spec = SPEC, regime_variant = REGIME_VARIANT)

global_summary <- walker_shares %>%
  summarise(across(all_of(share_cols), ~mean(.x, na.rm = TRUE))) %>%
  mutate(spec = SPEC, regime_variant = REGIME_VARIANT)

# --- Save --------------------------------------------------------------------
write_csv(walker_block_signed,
          paste0(OUT_PREFIX, "block_signed_", REGIME_VARIANT, ".csv"))
write_csv(walker_block_abs,
          paste0(OUT_PREFIX, "block_abs_", REGIME_VARIANT, ".csv"))
write_csv(walker_shares,
          paste0(OUT_PREFIX, "shares_", REGIME_VARIANT, ".csv"))
write_csv(regime_summary,
          paste0(OUT_PREFIX, "regime_summary_", REGIME_VARIANT, ".csv"))
write_csv(global_summary,
          paste0(OUT_PREFIX, "global_summary_", REGIME_VARIANT, ".csv"))

cat(sprintf("Done. SPEC=%s, REGIME=%s.\n", SPEC, REGIME_VARIANT))
cat("Output files:\n")
cat(sprintf("  - %sblock_signed_%s.csv\n", OUT_PREFIX, REGIME_VARIANT))
cat(sprintf("  - %sblock_abs_%s.csv\n", OUT_PREFIX, REGIME_VARIANT))
cat(sprintf("  - %sshares_%s.csv\n", OUT_PREFIX, REGIME_VARIANT))
cat(sprintf("  - %sregime_summary_%s.csv\n", OUT_PREFIX, REGIME_VARIANT))
cat(sprintf("  - %sglobal_summary_%s.csv\n", OUT_PREFIX, REGIME_VARIANT))

sum(walker_shares$total_abs)


# --- Sanity check: absolute shares should sum to 1 across blocks per row ---

abs_walker_cols <- c("A_AR", "A_PPI", "A_Monetary",
                     "A_FX", "A_Oil", "A_Trade", "A_Labour")

share_sums <- walker_shares %>%
  mutate(
    sum_A_Walker = rowSums(across(all_of(abs_walker_cols)), na.rm = TRUE)
  ) %>%
  select(Date, sum_A_Walker)

tol <- 1e-10

cat("=== Absolute shares sum-to-one check ===\n")
cat(sprintf("Walker - max |sum - 1| = %.2e  | rows OK: %d / %d\n",
            max(abs(share_sums$sum_A_Walker - 1), na.rm = TRUE),
            sum(abs(share_sums$sum_A_Walker - 1) < tol, na.rm = TRUE),
            sum(!is.na(share_sums$sum_A_Walker))))

# Vis eventuelle avvikende rader
violations <- share_sums %>%
  filter(abs(sum_A_Walker - 1) > tol)

if (nrow(violations) > 0) {
  cat("\nRader som ikke summerer til 1:\n")
  print(violations)
} else {
  cat("\nAlle rader summerer til 1 innenfor toleranse.\n")
}
mean(walker_shares$A_AR)
mean(walker_shares$A_FX)
mean(walker_shares$A_Oil)


##########---------------------------------------------------------------------
table_abs_kable_data <- walker_shares %>%
  group_by(regime) %>%
  summarise(
    A_AR       = mean(A_AR, na.rm = TRUE),
    A_PPI      = mean(A_PPI, na.rm = TRUE),
    A_Monetary = mean(A_Monetary, na.rm = TRUE),
    A_FX       = mean(A_FX, na.rm = TRUE),
    A_Oil      = mean(A_Oil, na.rm = TRUE),
    A_Trade    = mean(A_Trade, na.rm = TRUE),
    A_Labour   = mean(A_Labour, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = -regime,
    names_to = c("share_type", "block"),
    names_pattern = "^(A)_(AR|PPI|Monetary|FX|Oil|Trade|Labour)$",
    values_to = "value"
  ) %>%
  mutate(
    value = value * 100,
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
      str_detect(regime, "COVID")         ~ "COVID",
      str_detect(regime, "Energy")        ~ "Energy",
      str_detect(regime, "Disinflation")  ~ "Disinfl.",
      str_detect(regime, "Normalization") ~ "Normal",
      TRUE ~ regime
    )
  ) %>%
  select(block, regime, value) %>%
  pivot_wider(names_from = regime, values_from = value) %>%
  mutate(
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
  arrange(block)

table_abs_print <- table_abs_kable_data %>%
  select(block, COVID, Energy, `Disinfl.`, Normal)

kbl(
  table_abs_print,
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
  row_spec(0, bold = FALSE) %>%
  footnote(
    general = "Entries are mean absolute contribution shares within each regime. Shares sum to 100 within each regime, up to rounding.",
    general_title = "Note: "
  )


tab2 <- walker_shares %>%
  pivot_longer(
    cols = matches("^(A|S)_(AR|PPI|Monetary|FX|Oil|Trade|Labour)$"),
    names_to = c("measure", "block"),
    names_pattern = "^(A|S)_(AR|PPI|Monetary|FX|Oil|Trade|Labour)$",
    values_to = "share"
  ) %>%
  mutate(
    measure = recode(measure,
                     A = "Absolute",
                     S = "Signed"
    ),
    block = recode(block,
                   AR = "AR (inflation)",
                   PPI = "PPI / cost-push",
                   Monetary = "Monetary policy",
                   FX = "FX",
                   Oil = "Oil",
                   Trade = "Trade",
                   Labour = "Labour market"
    )
  ) %>%
  group_by(block, measure) %>%
  summarise(mean_share = mean(share, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = measure, values_from = mean_share) %>%
  select(block, Absolute, Signed) %>%
  mutate(
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
  arrange(block)

tab2 %>%
  kbl(
    col.names = c("Block", "Absolute share", "Signed share"),
    digits = 3,
    booktabs = TRUE,
    caption = "Mean attribution shares by block"
  ) %>%
  kable_styling(latex_options = c("hold_position"))
############################
setwd("/Users/amundbech1/Library/CloudStorage/OneDrive-UiTOffice365/Data_master/MasterOppgave")

# --- Last inn XGBoost signed bidrag ------------------------------------------
xgb_signed <- read_csv(paste0("xgb_harmonized_shap_group_signed_robustness.csv"),
                       show_col_types = FALSE) %>%
  mutate(Date = as.Date(date)) %>%
  select(-date)

# --- Map fra XGBoost lange blokk-navn til Walker korte navn ------------------
xgb_block_rename <- c(
  "AR (inflation)"  = "AR",
  "PPI / cost-push" = "PPI",
  "Monetary policy" = "Monetary",
  "FX"              = "FX",
  "Oil"             = "Oil",
  "Trade"           = "Trade",
  "Labour market"   = "Labour"
)

# --- Bygg xgb_block_signed med samme struktur som walker_block_signed --------
# Walker har: Date | AR PPI Monetary FX Oil Trade Labour | baseline_WalkerAdj | regime
xgb_block_signed <- xgb_signed %>%
  select(Date, all_of(names(xgb_block_rename)), shap_base) %>%
  rename(!!!setNames(names(xgb_block_rename), xgb_block_rename)) %>%
  rename(baseline_SHAP = shap_base) %>%
  mutate(regime = assign_regime(Date, REGIME_VARIANT)) %>%
  arrange(Date)

read_csv(walker_shares )

# --- Bygg xgb_block_abs (parallel til walker_block_abs) ----------------------
xgb_block_abs <- xgb_block_signed %>%
  mutate(across(all_of(block_names), abs)) %>%
  select(Date, all_of(block_names), regime)

# --- Shares (samme formel som walker_shares, linje 101-118) ------------------
xgb_shares <- xgb_block_signed %>%
  rowwise() %>%
  mutate(total_abs = sum(abs(c_across(all_of(block_names))), na.rm = TRUE)) %>%
  ungroup()

for (blk in block_names) {
  xgb_shares[[paste0("S_", blk)]] <- ifelse(
    xgb_shares$total_abs > 0,
    xgb_shares[[blk]] / xgb_shares$total_abs, NA_real_
  )
  xgb_shares[[paste0("A_", blk)]] <- ifelse(
    xgb_shares$total_abs > 0,
    abs(xgb_shares[[blk]]) / xgb_shares$total_abs, NA_real_
  )
}

# --- Regime og global summary (samme som Walker, linje 123-134) -------------
xgb_regime_summary <- xgb_shares %>%
  group_by(regime) %>%
  summarise(
    across(all_of(share_cols), ~mean(.x, na.rm = TRUE)),
    n_obs = n(),
    .groups = "drop"
  ) %>%
  mutate(spec = "xgb_robustness", regime_variant = REGIME_VARIANT)

xgb_global_summary <- xgb_shares %>%
  summarise(across(all_of(share_cols), ~mean(.x, na.rm = TRUE))) %>%
  mutate(spec = "xgb_robustness", regime_variant = REGIME_VARIANT)

# --- Lagre med xgb_-prefix og samme suffiks-konvensjon som Walker -----------
XGB_OUT_PREFIX <- "xgb_robustness_"

write_csv(xgb_block_signed,
          paste0(XGB_OUT_PREFIX, "block_signed_", REGIME_VARIANT, ".csv"))
write_csv(xgb_block_abs,
          paste0(XGB_OUT_PREFIX, "block_abs_", REGIME_VARIANT, ".csv"))
write_csv(xgb_shares,
          paste0(XGB_OUT_PREFIX, "shares_", REGIME_VARIANT, ".csv"))
write_csv(xgb_regime_summary,
          paste0(XGB_OUT_PREFIX, "regime_summary_", REGIME_VARIANT, ".csv"))
write_csv(xgb_global_summary,
          paste0(XGB_OUT_PREFIX, "global_summary_", REGIME_VARIANT, ".csv"))

cat(sprintf("\nXGBoost robustness done. REGIME=%s.\n", REGIME_VARIANT))
cat("Output files:\n")
cat(sprintf("  - %sblock_signed_%s.csv\n",    XGB_OUT_PREFIX, REGIME_VARIANT))
cat(sprintf("  - %sblock_abs_%s.csv\n",       XGB_OUT_PREFIX, REGIME_VARIANT))
cat(sprintf("  - %sshares_%s.csv\n",          XGB_OUT_PREFIX, REGIME_VARIANT))
cat(sprintf("  - %sregime_summary_%s.csv\n",  XGB_OUT_PREFIX, REGIME_VARIANT))
cat(sprintf("  - %sglobal_summary_%s.csv\n",  XGB_OUT_PREFIX, REGIME_VARIANT))

# --- Sanity check: A-shares sum-til-1 (samme som linje 159-187) -------------
abs_xgb_cols <- c("A_AR", "A_PPI", "A_Monetary",
                  "A_FX", "A_Oil", "A_Trade", "A_Labour")

xgb_share_sums <- xgb_shares %>%
  mutate(sum_A_xgb = rowSums(across(all_of(abs_xgb_cols)), na.rm = TRUE)) %>%
  select(Date, sum_A_xgb)

cat("=== XGBoost: absolute shares sum-to-one check ===\n")
cat(sprintf("XGBoost - max |sum - 1| = %.2e  | rows OK: %d / %d\n",
            max(abs(xgb_share_sums$sum_A_xgb - 1), na.rm = TRUE),
            sum(abs(xgb_share_sums$sum_A_xgb - 1) < tol, na.rm = TRUE),
            sum(!is.na(xgb_share_sums$sum_A_xgb))))








