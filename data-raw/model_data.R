## Code to prepare `model_data` data set goes here

indicator_ts_baltic <- read.csv("data-raw/indicator_timeseries_baltic.csv")
str(indicator_ts_baltic)
pressure_ts_baltic <- read.csv("data-raw/pressure_timeseries_baltic.csv")
str(pressure_ts_baltic)

# Import North Sea pressure time series for additional tests only
pressure_ts_northsea <- read.csv("data-raw/pressure_ts_northsea.csv")
str(pressure_ts_northsea)


usethis::use_data(
  indicator_ts_baltic,
  pressure_ts_baltic,
  pressure_ts_northsea,
  overwrite = TRUE
)


