# ============================================================
# Descriptive statistics for transformed model variables
# Matches the Walker/XGBoost harmonized specification
# ============================================================

library(dplyr)
library(tidyr)
library(readr)
library(knitr)

# --- Load data ---
df <- read.csv("master_data.csv", stringsAsFactors = FALSE)
df$Date <- as.Date(df$Date)
df <- df %>% arrange(Date)

# --- Forecast horizon ---
H <- 3

# --- Construct transformed variables exactly as in your Walker code ---
df_model <- df %>%
  mutate(
    # Raw YoY series
    cpi_yoy_raw     = (kpi / lag(kpi, 12) - 1) * 100,
    ppi_yoy_raw     = (ppi / lag(ppi, 12) - 1) * 100,
    oil_yoy_raw     = (oil_price_nok / lag(oil_price_nok, 12) - 1) * 100,
    import_yoy_raw  = (import / lag(import, 12) - 1) * 100,
    export_yoy_raw  = (eksport / lag(eksport, 12) - 1) * 100,
    
    # Lagged predictors
    cpi_yoy_lag1     = lag(cpi_yoy_raw, 1),
    ppi_yoy_lag1     = lag(ppi_yoy_raw, 1),
    oil_yoy_lag1     = lag(oil_yoy_raw, 1),
    usd_nok_lag1     = lag(usd_nok, 1),
    eur_nok_lag1     = lag(eur_nok, 1),
    import_yoy_lag1  = lag(import_yoy_raw, 1),
    export_yoy_lag1  = lag(export_yoy_raw, 1),
    unemp_lag1       = lag(unemployment, 1),
    rente_lag1       = lag(rente, 1),
    
    # Target
    target = lead(cpi_yoy_raw, H)
  ) %>%
  select(
    Date, target,
    cpi_yoy_lag1, ppi_yoy_lag1, oil_yoy_lag1,
    usd_nok_lag1, eur_nok_lag1,
    import_yoy_lag1, export_yoy_lag1,
    unemp_lag1, rente_lag1
  ) %>%
  filter(complete.cases(.))

# --- Check sample ---
cat("Effective sample starts:", min(df_model$Date), "\n")
cat("Effective sample ends:", max(df_model$Date), "\n")
cat("N =", nrow(df_model), "\n")

# --- Function for descriptive stats ---
desc_stats <- function(x) {
  c(
    Mean   = mean(x, na.rm = TRUE),
    Std    = sd(x, na.rm = TRUE),
    Min    = min(x, na.rm = TRUE),
    Q1     = quantile(x, 0.25, na.rm = TRUE),
    Median = median(x, na.rm = TRUE),
    Q3     = quantile(x, 0.75, na.rm = TRUE),
    Max    = max(x, na.rm = TRUE)
  )
}

# --- Variables to report ---
vars_for_table <- df_model %>%
  select(
    cpi_yoy_lag1,
    ppi_yoy_lag1,
    oil_yoy_lag1,
    usd_nok_lag1,
    eur_nok_lag1,
    rente_lag1,
    unemp_lag1,
    import_yoy_lag1,
    export_yoy_lag1
  )

# --- Compute table ---
desc_table <- t(sapply(vars_for_table, desc_stats)) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Variable")

# --- Pretty variable names for thesis ---
desc_table$Variable <- c(
  "CPI inflation, YoY (\\%)",
  "PPI inflation, YoY (\\%)",
  "Oil price, YoY (\\%)",
  "USD/NOK",
  "EUR/NOK",
  "Policy rate (\\%)",
  "Unemployment (\\%)",
  "Import growth, YoY (\\%)",
  "Export growth, YoY (\\%)"
)

# --- Round values ---
desc_table <- desc_table %>%
  mutate(across(-Variable, ~ round(.x, 2)))

# --- Print in R console ---
print(desc_table)

# --- Optional: save as csv ---
write.csv(desc_table, "descriptive_stats_transformed.csv", row.names = FALSE)

# --- Optional: produce LaTeX-ready table body ---
kable(
  desc_table,
  format = "latex",
  booktabs = TRUE,
  align = "lrrrrrrr",
  caption = "Descriptive statistics for transformed model variables, monthly data, 2001:02--2025:09",
  col.names = c("Variable", "Mean", "Std", "Min", "Q1", "Median", "Q3", "Max")
)

