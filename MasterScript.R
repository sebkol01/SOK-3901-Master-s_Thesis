
# setwd for Amund
setwd("/Users/amundbech1/Library/CloudStorage/OneDrive-UiTOffice365/Data_master/MasterOppgave")

# setwd for Seb
setwd("~/Library/CloudStorage/OneDrive-UiTOffice365/Amund Skjalgsønn Bech's files - Data_master/MasterOppgave")

# laster pakker 
library(tidyverse)
library(rjstat)
library(httr)
library(lubridate)
library(scales)
library(rsdmx)
library(zoo)
library(tseries)
library(jsonlite)
library(cowplot)
library(knitr)
library(kableExtra)
library(urca)
library(fredr)



# Henter data for KPI

url_kpi <- "https://data.ssb.no/api/v0/no/table/08981/"

kpi_query <- '{
  "query": [{"code": "Maaned","selection": {"filter": "item","values": [ "01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12" ] } } ],  "response": { "format": "json-stat2"  }}'

d.tmp1 <- POST(url_kpi , body = kpi_query, encode = "json", verbose())

kpi_data <- fromJSONstat(content(d.tmp1, "text"))

mnd_levels <- c("Januar","Februar","Mars","April","Mai","Juni",
                "Juli","August","September","Oktober","November","Desember")

kpi_data <- kpi_data %>%
  mutate(
    år_int  = as.integer(år),
    mnd_nr  = match(måned, mnd_levels),
    Date    = as.Date(sprintf("%d-%02d-01", år_int, mnd_nr)),
    value   = as.numeric(value)
  ) %>%
  arrange(Date)

kpi_data <- kpi_data %>% 
  mutate(Date = as.Date(Date))

sum(is.na(kpi_data$value))

kpi_data <- kpi_data %>% 
  na.omit()

boligindex_url <- "https://sdmx.oecd.org/public/rest/data/OECD.SDD.TPS,DSD_RHPI@DF_RHPI_ALL,1.0/COU.NOR.Q.RHPI.PC.S..._T?startPeriod=1990-Q1&endPeriod=2025-Q4&dimensionAtObservation=AllDimensions"

boligindex_data <- readSDMX(boligindex_url)

boligindex_data <- as.data.frame(boligindex_data)


boligindex_data <- boligindex_data %>% 
  select(TIME_PERIOD, obsValue) %>% 
  rename(Date = TIME_PERIOD, value = obsValue) %>% 
  mutate(Date = as.yearqtr(Date, format = "%Y-Q%q")) %>%  # hvis dato er tekst
  arrange(Date)

# endre datoformat på boligindex så det kan joines med andre datasett
boligindex_data <- boligindex_data %>% 
  mutate(Date = as.Date(as.yearqtr(Date, format = "%Y Q%q")))

url_ppi <- "https://data.ssb.no/api/v0/no/table/12462/"

ppi_query <- '{
  "query": [
    {"code": "Marked","selection": {"filter": "item","values": ["00"]}},
    {"code": "NaringUtenriks","selection": {"filter": "vs:NaringPPI1","values": ["SNN0"]}},
    {"code": "ContentsCode","selection": {"filter": "item","values": ["Indeksnivo"]}}
  ],
  "response": {"format": "json-stat2"}
}'

d.tmppi <- POST(url_ppi , body = ppi_query, encode = "json", verbose())

ppi_data <- fromJSONstat(content(d.tmppi, "text"))


ppi_data$date <- as.Date(as.yearmon(ppi_data$måned, format = "%YM%m"))


ppi_data <- ppi_data %>%
  select(date, ppi = value)

ppi_data <- ppi_data %>% 
  mutate(Date = as.Date(date))



# Rebase ppi til 2015 = 100
ppi_data <- ppi_data %>%
  mutate(
    date = as.Date(date),
    value_2015_100 = ppi / mean(ppi[year(date) == 2015], na.rm = TRUE) * 100
  )


oljepris <- read_csv("Brent_oil.csv")

head(oljepris)

oljepris <- oljepris %>% 
  mutate(
    Date = as.Date(Date, format = "%m/%d/%Y"),  # evt. "%d/%m/%Y" hvis det er det riktige
    Price = as.numeric(Price)
  ) %>%
  filter(!is.na(Date), !is.na(Price)) %>%
  arrange(Date)


# Prosess for å endre til NOK for hver dato
startP <- as.character(min(oljepris$Date, na.rm = TRUE))
endP   <- as.character(max(oljepris$Date, na.rm = TRUE))

fx_url <- paste0(
  "https://data.norges-bank.no/api/data/EXR/B.USD.NOK.SP",
  "?format=csv&locale=en&startPeriod=", startP,
  "&endPeriod=", endP
)

fx_raw <- read_delim(fx_url, delim = ";", show_col_types = FALSE)

fx <- fx_raw %>%
  transmute(
    Date = as.Date(TIME_PERIOD),
    usd_nok = as.numeric(OBS_VALUE)
  ) %>%
  arrange(Date)

fx <- fx %>% 
  mutate(Date = as.Date(Date))

oljepris <- oljepris %>%
  arrange(Date) %>%
  left_join(fx, by = "Date") %>%
  fill(usd_nok, .direction = "down") %>%      # fyll inn på helger/helligdager
  mutate(Price_NOK = Price * usd_nok)
                            



arbeidsledige_url <- "https://sdmx.oecd.org/public/rest/data/OECD.SDD.TPS,DSD_LFS@DF_IALFS_UNE_M,1.0/NOR..._Z.Y._T.Y_GE15..M?startPeriod=1989-01&endPeriod=2025-12&dimensionAtObservation=AllDimensions"

arbeidsledige_data <- readSDMX(arbeidsledige_url)

arbeidsledige_data <- as.data.frame(arbeidsledige_data)

arbeidsledige_data <- arbeidsledige_data %>%
  mutate(
    Date = ym(str_replace(TIME_PERIOD, "M", "-"))  # 2006M01 -> 2006-01 -> Date (1. i måneden)
  ) %>% 
  select(Date, obsValue) %>% 
  rename(arbledig = obsValue)

# FRED API KEY
fredr_set_key("8e19d5795cdc4a1fd2d07240b313890e")

# Load in data
eksport <- fredr(series_id = "XTEXVA01NOM667S")
import <- fredr(series_id = "XTIMVA01NOM667S")


# rename ekport
eksport <- eksport %>%
  rename(eksport = value)

# Rename import
import <- import %>%
  rename(import = value)

# Merge export and import using left_join
inthandel_data <- left_join(eksport, import, by = "date")

# Select only the date, eksport, and import columns
inthandel_data <- inthandel_data %>%
  select(date, eksport, import)

# Check if date is in date format
class(inthandel_data$date)

# Changing from dollar to NOK
fx_daily <- fx %>%
  mutate(Date = as.Date(Date)) %>%
  arrange(Date) %>%
  complete(Date = seq(min(Date), max(Date), by = "day")) %>%
  fill(usd_nok, .direction = "down")

  
inthandel_data <- inthandel_data %>%
  rename(Date = date) %>%      # bytt til "date" hvis det er den du har
  mutate(Date = as.Date(Date)) %>%
  arrange(Date) %>%
  left_join(fx_daily, by = "Date") %>%
  mutate(
    eksport = eksport * usd_nok,
    import  = import  * usd_nok
  )

# Remove all value up until 1988-06-01
inthandel_data <- inthandel_data %>%
  filter(Date > as.Date("1988-06-01"))

# Select only date, eksport, import
inthandel_data <- inthandel_data %>%
  select(Date, eksport, import)

########################################################
##### EUR NOK
euro_url <- "https://data.norges-bank.no/api/data/EXR/M.EUR.NOK.SP?format=csv&startPeriod=1980-01-01&endPeriod=2025-12-01&locale=no&bom=include"

euro <- read_delim(
  euro_url,
  delim = ";",
  locale = locale(decimal_mark = ",")
) |>
  transmute(
    date  = ym(TIME_PERIOD),
    value = OBS_VALUE
  )

euro


#### RENTER
ir_url <- "https://data.norges-bank.no/api/data/IR/M.KPRA..?format=csv&startPeriod=1983-01-01&endPeriod=2025-12-01&locale=no&bom=include"

ir_data <- read.csv(ir_url, sep = ";", dec = ",")

ir_data <- ir_data %>% 
  select(OBS_VALUE, TIME_PERIOD, Løpetid)


styringsrente <- ir_data %>%
  filter(Løpetid == "Styringsrenten") %>%
  select(TIME_PERIOD, rente = OBS_VALUE) %>%
  arrange(TIME_PERIOD)


styringsrente <- styringsrente %>% 
  rename(Date = TIME_PERIOD) %>% 
  mutate(Date = as.Date(paste0(Date, "-01"))) %>% 
  arrange(Date)

# Det var tre ulike typer renter i datasettet, det er gjort om så vi kun har styringsrente.

# --------------- PLOTS -------------------

theme_master <- function(base_size = 14){
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "plain"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_blank()
    )
}

kpiplot <- kpi_data %>% 
  ggplot(aes(x = Date, y = value)) +
  geom_line() +
  labs(
    title = "Konsumprisindeksen (KPI) over tid",
    x = "Dato",
    y = "KPI"
  ) +
  theme_master()
kpiplot


boligplot <- boligindex_data %>% 
  ggplot() +
  geom_line(aes(x = Date, y = value)) +
  labs(
    title = "Boligindex over tid",
    x = "Kvartal",
    y = "Boligprisvekst"
  ) +
  theme_master()

boligplot


# Plot PPI data
ppiplot <- ppi_data %>% 
  ggplot(aes(x = Date,  y = ppi)) +
  geom_line() +
  labs(
    title = "Produsentprisindeksen (PPI) over tid (2021 = 100)",
    x = "Dato",
    y = "PPI"
  ) +
  theme_master()

ppiplot


komboplot <- ggplot() +
  geom_line(data = kpi_data, aes(x = Date, y = value, color = "KPI")) +
  geom_line(data = boligindex_data, aes(x = Date, y = value, color = "Boligprisindex")) +
  geom_line(data = ppi_data, aes(x = Date, y = ppi, color = "PPI")) +
  scale_color_manual(values = c("KPI"="black", "Boligprisindex"="red", "PPI"="blue")) +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(breaks = pretty_breaks(n = 6)) +
  labs(
    title = "KPI, PPI og boligprisindeks over tid",
    x = NULL,
    y = "Indeksverdi"
  ) +
  theme_master()
komboplot


arbledigplot <- arbeidsledige_data %>%
  ggplot(aes(x = Date, y = arbledig)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(suffix = " %", decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Arbeidsledighet som andel av arbeidsstyrken",
    x = NULL,
    y = "Prosent"
  ) +
  theme_master()
arbledigplot

inthandelplot <- ggplot(inthandel_data, aes(x = Date)) +
  geom_line(aes(y = import, color = "Import")) +
  geom_line(aes(y = eksport, color = "Eksport")) +
  scale_color_manual(values = c("Import"="#1f77b4", "Eksport"="#d62728")) +
  scale_y_continuous(
    labels = label_number(scale = 1e-9, suffix = " mrd", decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Norsk utenrikshandel",
    x = NULL,
    y = "Milliarder kroner"
  ) +
  theme_master() +
  theme(legend.position = "top")
inthandelplot

oljeplot <- ggplot(oljepris, aes(x = Date, y = Price_NOK)) +
  geom_line(color = "#2ca02c") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(suffix = " kr", decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "Oljepris (NOK)",
    x = NULL,
    y = "Pris per fat"
  ) +
  theme_master()
oljeplot

europlot <- euro %>%
  ggplot(aes(x = date, y = value)) +
  geom_line(color = "black") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  scale_y_continuous(
    labels = label_number(decimal.mark = ","),
    breaks = pretty_breaks(n = 6)
  ) +
  labs(
    title = "EUR/NOK (spot) – månedlig gjennomsnitt",
    x = NULL,
    y = "NOK per EUR"
  ) +
  theme_master()

europlot


plot_grid(komboplot, inthandelplot, arbledigplot, oljeplot)

renteplot <- styringsrente %>% 
  ggplot() +
  geom_line(aes(x = Date, y = rente))+
  theme_master()
renteplot


##### ------------------- ADF TEST -------------------------
# Funksjon for å kjøre ADF-test og returnere resultat
run_adf <- function(data, varname) {
  # Fjern NA-verdier
  clean_data <- na.omit(data)
  
  # Sjekk om det er nok data igjen
  if(length(clean_data) < 3) {
    return(data.frame(
      Variabel = varname,
      `ADF-statistikk` = NA,
      `P-verdi` = NA,
      Status = "For få observasjoner"
    ))
  }
  
  tryCatch({
    test <- adf.test(clean_data, alternative = "stationary")
    
    data.frame(
      Variabel = varname,
      `ADF-statistikk` = round(test$statistic, 4),
      `P-verdi` = round(test$p.value, 4),
      Status = ifelse(test$p.value < 0.05, "Stasjonær", "Ikke-stasjonær")
    )
  }, error = function(e) {
    data.frame(
      Variabel = varname,
      `ADF-statistikk` = NA,
      `P-verdi` = NA,
      Status = paste("Feil:", e$message)
    )
  })
}

str(tollsatser_flat)

tollsatser_flat$sats_num <- as.numeric(str_replace(tollsatser_flat$sats, ",", "."))

##### ------------------- STRUKTURELT BRUDD -------------------------

# Zivot-Andrews på eksport
za_ex <- ur.za(na.omit(inthandel_data$eksport), model = "both", lag = 2)
summary(za_ex)
plot(za_ex)



# Zivot-Andrews på PPI YoY
ppi_yoy <- diff(log(ppi_data$ppi), lag = 12) * 100
za_ppi <- ur.za(ppi_data$ppi, lag = 1, model = "both")
summary(za_ppi)
plot(za_ppi)


# Do a cointegration test between eksport and ppi

# Legger til USD oljepris og handeldata i allerede merget fil fra "Merge_data_for_python.R"

master_data <- read_csv("master_data.csv")

oljepris <- oljepris %>% 
  rename(oljepris_USD = Price) %>% 
  select(oljepris_USD, Date)

inthandel_data <- inthandel_data %>% 
  rename(eksport_USD = eksport, import_USD = import)

master_data <- master_data %>% 
  left_join(oljepris, by = "Date") %>% 
  left_join(inthandel_data %>% select(Date, import_USD), by = "Date") %>% 
  left_join(inthandel_data %>% select(Date, eksport_USD), by = "Date")


