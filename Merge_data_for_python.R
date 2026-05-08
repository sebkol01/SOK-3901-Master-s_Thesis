# ============================================================
# merge_dataset.R
# Kjør ETTER at du har kjørt hoveddatascriptet ditt,
# slik at alle data-objektene allerede finnes i environment.
# ============================================================

library(tidyverse)
library(lubridate)
library(zoo)

# ──────────────────────────────────────────────
# 1. KPI – månedlig
# ──────────────────────────────────────────────
kpi <- kpi_data %>%
  select(Date, kpi = value)

# ──────────────────────────────────────────────
# 2. PPI – månedlig, rebasert til 2015=100
# ──────────────────────────────────────────────
ppi <- ppi_data %>%
  select(Date, ppi = value_2015_100)

# ──────────────────────────────────────────────
# 3. Arbeidsledighet (OECD) – månedlig
# ──────────────────────────────────────────────
unemp <- arbeidsledige_data %>%
  select(Date = Date, unemployment = arbledig)

# ──────────────────────────────────────────────
# 4. Oljepris NOK + USD/NOK – daglig → månedlig
# ──────────────────────────────────────────────
# olje <- oljepris %>%
#   mutate(ym = floor_date(Date, "month")) %>%
#   group_by(ym) %>%
#   summarise(
#     oil_price_nok = mean(Price_NOK, na.rm = TRUE),
#     usd_nok       = mean(usd_nok, na.rm = TRUE),
#     .groups = "drop"
#   ) %>%
#   rename(Date = ym)

# ──────────────────────────────────────────────
# 5. Import/Eksport – årlig → spredd til måneder
# ──────────────────────────────────────────────
handel <- inthandel_data %>%
  select(Date, import, eksport)

# ──────────────────────────────────────────────
# 6. Styringsrenten – allerede i environment som `ir_data`
# ──────────────────────────────────────────────
rente <- styringsrente %>%
  select(Date, rente)

# ──────────────────────────────────────────────
# 7. EUR/NOK – allerede i environment som `euro`
# ──────────────────────────────────────────────
euro <- euro %>%
  select(Date = date, eur_nok = value)

# ──────────────────────────────────────────────
# 8. Boligprisindeks – allerede i environment som `boligindex_data`
#    (forventer kolonner: dato (yearqtr), value)
# ──────────────────────────────────────────────
bolig <- boligindex_data %>%
  select(Date,  boligindex = value) %>%
  arrange(Date) %>%
  # Ekspander kvartalsdata til månedlig med lineær interpolering
  complete(Date = seq(min(Date), max(Date), by = "month"))

# ──────────────────────────────────────────────
# 9. MERGE ALT
# ──────────────────────────────────────────────
master <- kpi %>%
  left_join(ppi,    by = "Date") %>%
  left_join(unemp,  by = "Date") %>%
  left_join(olje,   by = "Date") %>%
  left_join(handel,  by = "Date") %>%
  left_join(rente,   by = "Date") %>%
  left_join(euro,    by = "Date") %>%
  left_join(bolig,   by = "Date") %>%
  arrange(Date)

# ──────────────────────────────────────────────
# 10. OVERSIKT
# ──────────────────────────────────────────────
cat("\n========================================\n")
cat("Samlet datasett:\n")
cat(sprintf("  Rader:    %d\n", nrow(master)))
cat(sprintf("  Kolonner: %d\n", ncol(master)))
cat(sprintf("  Periode:  %s → %s\n", min(master$Date), max(master$Date)))
cat("\nManglende verdier per kolonne:\n")
print(colSums(is.na(master)))

cat("\nFørste rader med komplett data:\n")
print(head(master %>% drop_na(), 10))

# ──────────────────────────────────────────────
# 11. LAGRE
# ──────────────────────────────────────────────
write_csv(master, "master_data.csv")
cat("\n✓ Lagret: master_data.csv\n")

# Versjon uten NA (klar til modellering)
master_complete <- master %>% drop_na()
write_csv(master_complete, "master_data_complete.csv")
cat(sprintf("✓ Lagret: master_data_complete.csv (%d rader)\n", nrow(master_complete)))
