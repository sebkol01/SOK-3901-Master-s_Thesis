# Replication Code — Norwegian CPI Inflation Forecasting
## Master's Thesis: Walker TVP vs XGBoost SHAP

**Authors:** Seb & Amund  
**Language:** R

---

## Overview

This repository contains the R code used to estimate the Time-Varying Parameter (Walker) model, compute and harmonise attribution contributions, compare them with XGBoost SHAP values, and produce all figures reported in the thesis.

The single entry-point script is:

```
Walker_thesis_pipeline.R
```

All parts are designed to run top-to-bottom in one R session. The MCMC backtests (Parts 2–3) take approximately 1.5 hours each and write their outputs to CSV files. All downstream parts read from those files, so they can be re-run independently without re-estimating.

---

## Required Input Files

The following files must be present in the working directory before running the script:

| File | Description |
|------|-------------|
| `master_data.csv` | Monthly macro data (KPI, PPI, oil, FX, trade, unemployment, policy rate). Must contain columns `oljepris_USD`, `import_USD`, `eksport_USD` for the USD robustness specification. |
| `xgb_v8_harmonized_shap_group_signed.csv` | XGBoost SHAP block contributions — main specification (from Python pipeline). |
| `xgb_harmonized_robustness_shap_group_signed.csv` | XGBoost SHAP block contributions — USD robustness specification. |
| `xgb_v8_harmonized_predictions.csv` | XGBoost point forecasts — used for the accuracy table and Walker vs XGBoost figure. |

---

## R Package Dependencies

```r
install.packages(c(
  "walker",       # time-varying parameter model (wraps Stan)
  "rstan",        # Stan interface
  "dplyr",        # data manipulation
  "tidyr",        # reshaping
  "readr",        # fast CSV I/O
  "ggplot2",      # figures
  "lubridate",    # date arithmetic
  "scales",       # axis formatting
  "knitr",        # table rendering
  "kableExtra"    # LaTeX table formatting
))
```

**Note on `walker` and `rstan`:** The `walker` package requires a working Stan installation. Follow the instructions at <https://mc-stan.org/rstan/> before running Parts 2–3.

---

## Script Structure

### Part 0 — Configuration
Sets all fixed hyperparameters and global constants. **Do not change these values** if you need to reproduce the original results exactly.

| Parameter | Value | Description |
|-----------|-------|-------------|
| `H` | 3 | Forecast horizon (months ahead) |
| `TEST_START` | 2020-01-01 | First expanding-window forecast origin |
| `MIN_TRAIN_N` | 50 | Minimum training observations required |
| `N_CHAINS` | 4 | MCMC chains per fit |
| `N_ITER` | 2000 | Total MCMC iterations per chain |
| `N_WARMUP` | 1000 | Warmup (burn-in) iterations |
| `SEED` | 42 | Random seed (set via `set.seed`) |
| `adapt_delta` | 0.95 | Stan sampler control parameter |
| `max_treedepth` | 12 | Stan sampler control parameter |

### Part 0b — Loop Helper Functions
Four small functions called at every forecast origin inside the MCMC loops:

- `scale_with_train_stats` — standardises features using training-window means/SDs only (enforces the real-time information constraint).
- `extract_fit_diagnostics` — extracts divergent transitions, worst R-hat, and lowest ESS.
- `predict_last_state` — computes the one-step-ahead forecast from filtered (last time-point) posterior mean coefficients.
- `extract_filtered_coefs` — saves the filtered posterior summary at the forecast origin.

Two accuracy functions used in Part 9:
- `safe_mape` — MAPE with protection against zero actuals.
- `forecast_metrics` — collects RMSE, MAE, MAPE, Bias, and N for one model.

### Part 1 — Data Loading and Feature Engineering
Reads `master_data.csv` and builds all nine lag-1 features. The feature engineering follows this order:

1. Compute YoY growth rates at time *t* (for KPI, PPI, oil, imports, exports).
2. Lag all features by 1 month — this is the lag-1 information set.
3. Set the target as `lead(kpi_yoy_raw, H = 3)`.

This ensures the forecast at each origin uses only information available at that date, matching the XGBoost specification.

### Part 2 — Walker Backtest (Main Specification)
Expanding window MCMC backtest from January 2020. At each forecast origin:

1. Split: `train = {t : t < origin}`, `new_x = {t = origin}`.
2. Standardise all nine features using training means/SDs.
3. Fit the Walker model with the formula:
   ```
   target ~ -1 + rw1(~ 1 + [9 features], beta=c(0,10), sigma=c(2,0.01))
   ```
4. Extract filtered (last time-point) posterior mean coefficients.
5. Compute the forecast: `y_hat = beta_0(T) + sum_j beta_j(T) * x_j(T)`.
6. Store real-time contributions: `c_jt = beta_jt * x_jt` (intercept contribution = `beta_0t`).

**Outputs:**
- `walker_forecasts.csv` — point forecasts and MCMC diagnostics
- `walker_filtered_coefs.csv` — filtered posterior at each origin
- `walker_x_values_at_origin.csv` — standardised feature values at each origin
- `walker_contributions.csv` — real-time feature-level contributions
- `walker_diagnostics.txt` — run log

### Part 3 — Walker Backtest (USD External Robustness)
Identical to Part 2, but replaces NOK-denominated oil price and trade with USD-denominated equivalents (`oljepris_USD`, `import_USD`, `eksport_USD`). Output files use the prefix `walker_usd_external_`.

### Part 4 — Harmonisation: Walker Main Spec
Aggregates the nine feature-level contributions to seven economic blocks:

| Block | Features |
|-------|----------|
| AR | `kpi_yoy_lag1` |
| PPI | `ppi_yoy_lag1` |
| Oil | `oil_yoy_lag1` |
| FX | `usd_nok_lag1`, `eur_nok_lag1` |
| Trade | `import_yoy_lag1`, `eksport_yoy_lag1` |
| Labour | `unemp_lag1` |
| Monetary | `rente_lag1` |

Computes signed attribution shares `S_k = C_k / Σ_m |C_m|` and absolute shares `A_k = |C_k| / Σ_m |C_m|` per block per date. Also assigns dates to one of four regimes and produces regime-conditional summaries.

**Regime cutoffs (main):**

| Regime | Period |
|--------|--------|
| COVID (2020–2021) | Before June 2021 |
| Energy Crisis | June 2021 – December 2022 |
| Disinflation | January 2023 – May 2024 |
| Normalization | June 2024 onward |

**Outputs:** `walker_harmonized_*.csv`

### Part 5 — Harmonisation: Walker USD External Spec
Same procedure as Part 4 applied to the USD robustness contributions. **Outputs:** `walker_usd_external_harmonized_*.csv`

### Part 6 — Harmonisation: XGBoost SHAP Outputs
Reads the Python-produced SHAP files and renames the long block labels (e.g., `"AR (inflation)"`) to the short names used internally (`"AR"`). Computes the same shares and regime summaries as the Walker harmonisation. **Outputs:** `xgb_main_harmonized_*.csv`, `xgb_robustness_harmonized_*.csv`

### Part 7 — SHAP vs Walker Comparison
Builds joint attribution tables for both specifications by inner-joining the XGBoost SHAP and Walker block-signed tables on `Date`. Computes:

- Sign agreement (fraction of periods where both models agree on direction).
- Pearson correlation (signed and absolute) per block.
- Mean signed and absolute attribution shares per block.
- Block importance by regime (mean absolute share × 100%).

**Outputs:** `comparison_harmonised_attributions.csv`, `comparison_routing_compare.csv`, `comparison_mean_shares_*.csv`, `comparison_block_importance_*.csv`

### Part 8 — Residual Decomposition
Computes the attribution residual `Delta_kt = C_kt^XGB - C_kt^TVP` per block per date and regime-conditional summaries of its variance and mean absolute magnitude.

**Outputs:** `residual_long.csv`, `residual_metrics_long.csv`

### Part 9 — AR(1) Baseline and Accuracy Table
Runs an expanding-window AR(1) backtest (`target ~ kpi_yoy_lag1`) and combines all model forecasts into a single accuracy table (RMSE, MAE, MAPE, Bias).

**Outputs:** `ar1_forecasts.csv`, `accuracy_table.csv`

### Part 10 — Descriptive Figures
Produces time-series plots of all macro variables from `master_data.csv`. No API calls required — all figures read from the pre-built CSV.

**Outputs:** `cpi.png`, `PPI.png`, `unemployment.png`, `oljeplot.png`, `olje_usd.png`, `usd_nok.png`, `eur_nok.png`, `policy_rate.png`, `trade.png`, `trade_usd.png`

### Part 11 — Model Output Figures
Produces all model-related figures:

- **`walker_forecast_plot.png`** — Walker vs random walk vs actual.
- **`walker_error_plot.png`** — Forecast errors over time.
- **`walker_vs_xgboost.png`** — Walker vs XGBoost vs actual (with RMSE in legend).
- **`comparison_plot_contributions.png`** — Block-level SHAP vs Walker contributions.
- **`residual_plot.png`** — Attribution residuals by block, coloured by regime.

---

## Other R Scripts (legacy/partial)

These files are kept for reference but are superseded by `Walker_thesis_pipeline.R`:

| File | Description |
|------|-------------|
| `Walker estimation.r` | Original development version (v6). Contains the pipeline in a single unsectioned file with some exploratory/incomplete code. |
| `Walker_estimation_v5.R` | Earlier version (v5) — no real-time contributions, uses a full-sample fit. |
| `walker_harmonize_from_outputs.R` | Standalone harmonisation script (used a different regime scheme). |
| `Walker_harmonize_robustness.R` | Standalone robustness harmonisation + XGBoost comparison. |
| `Walker_robustness_USD.R` | Standalone USD external backtest. |
| `acc_table.R` | Standalone accuracy table builder. |
| `MasterScript.R` | Data collection via live APIs (SSB, Norges Bank, OECD, FRED). Run this to rebuild `master_data.csv` from scratch. |
| `CompiledCode_Claude.R` | Earlier compilation attempt (more function-heavy). |

---

## Reproducibility Notes

- The random seed (`SEED = 42`) is set at the top of the script via `set.seed(42)`. Stan's MCMC sampler is additionally controlled by `rstan_options(auto_write = TRUE)` and the per-chain seeds derived from the global seed.
- Results may differ slightly across platforms or Stan versions due to floating-point arithmetic, but the economic conclusions are robust.
- All real-time standardisation is done strictly within the training window at each forecast origin — there is no look-ahead bias in the scaling step.
- The contribution formula is `c_jt = beta_jt^filtered * x_jt^scaled`. The intercept contribution equals the filtered intercept `beta_0t` directly (not multiplied by any x-value).
