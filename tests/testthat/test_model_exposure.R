
test_that("output data frame", {
  dat <- model_exposure(
    pressure_time_series = pressure_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  )
  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(8, 15))
  expect_true(all(dat$pressure %in% names(pressure_ts_baltic)))
  expect_type(dat$pressure, "character")
  expect_type(dat$exposure, "double")
  expect_type(dat$uncertainty, "double")
  expect_type(dat$comp_magnitude, "double")
  expect_type(dat$comp_frequency, "double")
  expect_type(dat$comp_trend, "double")
  expect_type(dat$comp_direction, "character")
  expect_type(dat$comp_spatial, "double")
  expect_type(dat$uncertainty_arima, "double")
  expect_type(dat$uncertainty_gam, "double")
  expect_type(dat$mean_baseline, "double")
  expect_type(dat$mean_current, "double")
  expect_type(dat$standard_deviation_baseline, "double")
  expect_type(dat$slope_linear_model, "double")
  expect_type(dat$p_value_linear_model, "double")
})


test_that("results are correct", {
  #returning trend = good
  dat <- model_exposure(
    pressure_time_series = pressure_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  )
  expect_equal(dat$comp_magnitude, c(3, 1, 2, 2, 2, 2, 2, 2))
  expect_equal(dat$comp_trend, c(3, 3, 3, 4, 3, 4, 4, 3))
  expect_equal(dat$comp_frequency, c(5, 3, 4, 3, 5, 4, 5, 5))
  expect_equal(dat$exposure, c(3.5, 2.5, 3, 3, 3.25, 3.25, 3.5, 3.25))
  expect_equal(dat$comp_direction, c("increase", "increase", "increase", "increase",
    "decrease", "decrease",  "decrease", "increase"))
  expect_equal(dat$uncertainty, c(2, 2, 2, 2, 2, 2, 2, 2))
  expect_equal(dat$uncertainty_arima, c(2, 2, 2, 2, 2, 3, 2, 3))
  expect_equal(dat$uncertainty_gam, c(2, 2, 2, 2, 3, 2, 2, 2))

  #leaving trend = good
  dat <- model_exposure(
    pressure_time_series = pressure_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016),
    trend = "leave"
  )
  expect_equal(dat$comp_magnitude, c(3, 1, 2, 2, 2, 2, 2, 2))
  expect_equal(dat$comp_trend, c(3, 3, 3, 2, 3, 2, 2, 3))
  expect_equal(dat$comp_frequency, c(5, 3, 4, 3, 5, 4, 5, 5))
  expect_equal(dat$exposure, c(3.5, 2.5, 3, 2.5, 3.25, 2.75, 3, 3.25))
  expect_equal(dat$comp_direction, c("increase", "increase", "increase", "increase",
    "decrease", "decrease", "decrease", "increase"))

})


test_that("weaker / stronger trends", {
  # Provoke weaker / stronger trends in the North Sea datasets and test output
  pressures <- pressure_ts_northsea
  set.seed(1)
  pressures$bot_temp[42:51] <- jitter(pressure_ts_northsea$bot_temp[1:10])

  dat <- model_exposure(
    pressure_time_series = pressures,
    base_years = c(start = 1970, end = 1979),
    current_years = c(start = 2011, end = 2020)
  )
  expect_equal(dat$comp_magnitude, c(1, 5, 5))
  expect_equal(dat$comp_trend, c(4, 3, 3))
  expect_equal(dat$comp_frequency, c(1, 5, 5))
  expect_equal(dat$exposure, c(2.25, 4, 4))
  expect_equal(dat$comp_direction, c("increase", "increase", "increase"))
  expect_equal(dat$uncertainty, c(2, 1, 2))
  expect_equal(dat$uncertainty_arima, c(2, 2, 3))
  expect_equal(dat$uncertainty_gam, c(2, 1, 2))

  pressures <- pressure_ts_northsea
  pressures$bot_temp[42:51] <- mean(pressure_ts_northsea$bot_temp[1:10]) + c(1:10)*1.2

  dat <- model_exposure(
    pressure_time_series = pressures,
    base_years = c(start = 1970, end = 1979),
    current_years = c(start = 2011, end = 2020)
  )
  expect_equal(dat$comp_magnitude, c(5, 5, 5))
  expect_equal(dat$comp_trend, c(4, 3, 3))
  expect_equal(dat$comp_frequency, c(5, 5, 5))
  expect_equal(dat$exposure, c(4.25, 4, 4))
  expect_equal(dat$comp_direction, c("increase", "increase", "increase"))
  expect_equal(dat$uncertainty, c(2, 1, 2))
  expect_equal(dat$uncertainty_arima, c(2, 2, 3))
  expect_equal(dat$uncertainty_gam, c(2, 1, 2))

})


# Test the new data input validations
test_that("warning and error display", {
  test_dat <- pressure_ts_baltic[ ,1:3]
  test_dat2 <- test_dat
  test_dat2$nitrogen <- as.character(test_dat2$nitrogen)
  test_dat3 <- test_dat
  test_dat3$year <- test_dat3$year[c(20:33, 1:19)]
  test_dat4 <- pressure_ts_baltic[ ,1:4]

  # 'pressure_time_series' argument not a data.frame
  expect_error(model_exposure(
    pressure_time_series = tibble::as_tibble(test_dat),
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(model_exposure(
    pressure_time_series = as.matrix(test_dat),
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  # 'pressure_time_series' columns have wrong data types
  expect_error(model_exposure(
    pressure_time_series = test_dat2,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  # the first column (time) not sorted
  expect_warning(model_exposure(
    pressure_time_series = test_dat3,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  # no vector or df provided for baseline and current period
  suppressWarnings(expect_warning(model_exposure(
    pressure_time_series = test_dat
  )))
  # wrong assignments
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years_by_press = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years = c(start = 1984, end = 1994),
    current_years_by_press = c(start = 2010, end = 2016)
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1990, 1995), end = c(1994, 2005)),
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years = c(start = 1984, end = 1994),
    current_years = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1995, 1995), end = c(2000, 2000))
  ))
  # Missing value in general periods
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years = 1984,
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years = c(start = 1984, end = 1994),
    current_years = 2010
  ))

  # 'base_years_by_press' and 'current_years_by_press' not a data frame
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years_by_press = tibble::as_tibble(data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1990, 1995), end = c(1994, 2005))),
    current_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1995, 1995), end = c(2000, 2000))
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1990, 1995), end = c(1994, 2005)),
    current_years_by_press = as.matrix(data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1995, 1995), end = c(2000, 2000)))
  ))

  # data types in 'base_years_by_press' and 'current_years_by_press' not correct
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c("1990", "1995"), end = c(1994, 2005)),
    current_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1995, 1995), end = c(2000, 2000))
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1990, 1995), end = c(1994, 2005)),
    current_years_by_press = as.matrix(data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c("1995", "1995"), end = c(2000, 2000)))
  ))
  # not all pressure names listed or wrong pressure names
  expect_error(model_exposure(
    pressure_time_series = test_dat4,
    base_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1990, 1995), end = c(1994, 2005)),
    current_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous", "surf_temp_sum"),
      start = c(1995, 1996, 1996), end = c(2000, 2000, 2005))
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat4,
    base_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous", "surf_temp_sum"),
      start = c(1990, 1995, 1995), end = c(1994, 2005, 2005)),
    current_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1995, 1996), end = c(2000, 2000))
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years_by_press = data.frame(
      press = c("blub", "phosphorous"),
      start = c(1990, 1995), end = c(1994, 2005)),
    current_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1995, 1996), end = c(2000, 2000))
  ))
  # specified years fall outside time series range
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years = c(start = 1983, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2020)
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1980, 1995), end = c(1994, 2005)),
    current_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1995, 1996), end = c(2000, 2016))
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1984, 1995), end = c(1994, 2005)),
    current_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1995, 1996), end = c(2000, 2020))
  ))
  # starting years are higher than ending years
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years = c(start = 1994, end = 1984),
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2016, end = 2010)
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1994, 1995), end = c(1984, 2005)),
    current_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1995, 1996), end = c(2000, 2016))
  ))
  expect_error(model_exposure(
    pressure_time_series = test_dat,
    base_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1984, 1995), end = c(1994, 2005)),
    current_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(1995, 2016), end = c(2000, 2005))
  ))

})

