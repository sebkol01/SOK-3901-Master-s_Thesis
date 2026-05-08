

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

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

# ------------------------------------------------------------
# Read Walker forecasts
# ------------------------------------------------------------

walker <- read_csv("walker_forecasts.csv", show_col_types = FALSE) %>%
  mutate(Date = as.Date(Date)) %>%
  rename(
    actual_walker = y_actual,
    walker = y_hat,
    rw_walker = y_rw
  )

# ------------------------------------------------------------
# Add AR(1) forecasts
# ------------------------------------------------------------

ar1_results <- ar1_results %>%
  mutate(Date = as.Date(Date, origin = "1970-01-01"))

walker <- walker %>%
  left_join(ar1_results, by = "Date") %>%
  rename(ar1 = y_ar1)

# ------------------------------------------------------------
# Read XGBoost forecasts
# ------------------------------------------------------------


xgb <- read_csv("~/Library/CloudStorage/OneDrive-UiTOffice365/Amund Skjalgsønn Bech's files - Data_master/Python/results_harmonized/xgb_v8_harmonized_predictions.csv",
                show_col_types = FALSE) %>%
  mutate(date = as.Date(date)) %>%
  rename(
    Date = date,
    actual_xgb = actual,
    xgboost = predicted_raw,
    rw_xgb = y_rw
  )

# ------------------------------------------------------------
# Merge Walker, XGBoost, RW and AR(1)
# ------------------------------------------------------------

combined <- walker %>%
  inner_join(
    xgb %>% select(Date, actual_xgb, xgboost, rw_xgb),
    by = "Date"
  )

# Check whether actual values match across files
combined <- combined %>%
  mutate(actual_diff = actual_walker - actual_xgb)

cat("Max absolute actual difference:",
    max(abs(combined$actual_diff), na.rm = TRUE), "\n")

cat("Common sample:",
    as.character(min(combined$Date)),
    "to",
    as.character(max(combined$Date)), "\n")

cat("N common observations:", nrow(combined), "\n")

# ------------------------------------------------------------
# Use common actual series
# ------------------------------------------------------------

combined_eval <- combined %>%
  transmute(
    Date = Date,
    actual = actual_walker,
    `Random walk` = rw_walker,
    `AR(1)` = ar1,
    Walker = walker,
    XGBoost = xgboost
  )

# ------------------------------------------------------------
# Final accuracy table
# ------------------------------------------------------------

forecast_accuracy_table <- bind_rows(
  forecast_metrics(combined_eval$actual, combined_eval$`Random walk`, "Random walk"),
  forecast_metrics(combined_eval$actual, combined_eval$`AR(1)`, "AR(1)"),
  forecast_metrics(combined_eval$actual, combined_eval$Walker, "Walker"),
  forecast_metrics(combined_eval$actual, combined_eval$XGBoost, "XGBoost")
)

print(forecast_accuracy_table)
