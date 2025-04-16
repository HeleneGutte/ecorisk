

test_that("output data frame", {

  dat <- model_sensitivity(
    indicator_time_series = indicator_ts_baltic,
    pressure_time_series = pressure_ts_baltic,
    current_years = c(start = 2010, end = 2016)
  )

  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(16, 13))
  expect_true(all(dat$pressure %in% names(pressure_ts_baltic)))

  expect_type(dat$indicator, "character")
  expect_type(dat$pressure, "character")
  expect_type(dat$type, "character")
  expect_type(dat$pathway, "character")
  expect_type(dat$sensitivity, "double")
  expect_type(dat$adaptive_capacity, "double")
  expect_type(dat$uncertainty_sen, "double")
  expect_type(dat$uncertainty_ac, "logical")
  expect_type(dat$r_sq, "double")
  expect_type(dat$p_value, "double")
  expect_type(dat$edf, "double")
  expect_type(dat$uncertainty_gam, "double")
  expect_type(dat$uncertainty_arima, "double")

})


test_that("test results", {

  # Same assessment period
  dat <- model_sensitivity(
    indicator_time_series = indicator_ts_baltic,
    pressure_time_series = pressure_ts_baltic,
    current_years = c(start = 2010, end = 2016)
  )

  expect_equal(dat$sensitivity[1:10],
    c(3, 0, 0, 0, -1, -1, 0, 0, -1, 0))
  expect_equal(round(dat$r_sq[1:10], digits = 3),
    c(0.243, 0, -0.018, 0.029, 0.101, 0.167, -0.029, -0.03, -0.022, 0.04))
  expect_equal(round(dat$p_value[1:10], digits = 3),
    c(0.017, 0.327, 0.513, 0.17, 0.04, 0.071, 0.745, 0.802, 0.79, 0.136))
  expect_equal(round(dat$edf[1:10], digits = 3),
    c(2.848, 1, 1, 1, 1, 2.193, 1, 1, 1.127, 1))
  expect_equal(dat$uncertainty_sen[1:10],
    c(2, 1, 1, 1, 2, 2, 1, 2, 2, 1))
  expect_equal(dat$uncertainty_arima[1:10],
    c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2))

  # Different assessment periods
  current <- data.frame(
    ind = rep(x = names(indicator_ts_baltic[, -1]), each = ncol(pressure_ts_baltic[, -1])),
    press = rep(x = names(pressure_ts_baltic[, -1]), times = ncol(indicator_ts_baltic[, -1])),
    start = c(rep(2010, 8), rep(2005, 8)),
    end = c(rep(2014, 8), rep(2011, 8))
  )
  dat2 <- model_sensitivity(
    indicator_time_series = indicator_ts_baltic,
    pressure_time_series = pressure_ts_baltic,
    current_years_by_ind_press = current
  )

  expect_equal(dat2$sensitivity[1:10],
    c(3, 0, 0, 0, -1, -1, 0, 0, 1, 0))
  expect_equal(round(dat2$r_sq[1:10], digits = 3),
    c(0.243, 0, -0.018, 0.029, 0.101, 0.167, -0.029, -0.03, -0.022, 0.04))
  expect_equal(round(dat2$p_value[1:10], digits = 3),
    c(0.017, 0.327, 0.513, 0.17, 0.04, 0.071, 0.745, 0.802, 0.79, 0.136))
  expect_equal(round(dat2$edf[1:10], digits = 3),
    c(2.848, 1, 1, 1, 1, 2.193, 1, 1, 1.127, 1))
  expect_equal(dat2$uncertainty_sen[1:10],
    c(2, 2, 1, 2, 2, 2, 2, 2, 2, 2))
  expect_equal(dat2$uncertainty_arima[1:10],
               c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2))

})

test_that("break function if pressures or indicators are NA", {

  fish <- data.frame(
    time = pressure_ts_baltic$year,
    fish = as.numeric(rep(NA, length(pressure_ts_baltic$year)))
  )
  suppressWarnings(dat <- model_sensitivity(
    indicator_time_series = fish,
    pressure_time_series = pressure_ts_baltic[, c(1:6)],
    current_years = c(2010, 2016)
  ))
  expect_equal(dat$sensitivity, c(NA, NA, NA, NA, NA))
  expect_equal(dat$r_sq, c(NA, NA, NA, NA, NA))
  expect_equal(dat$uncertainty_sen, c(NA, NA, NA, NA, NA))

  new_press <- data.frame(
    "year" = c(1984:2016),
    "oil" = as.numeric(rep(NA, 33)),
    "mhw" = as.numeric(rep(NA, 33))
  )
  suppressWarnings(dat <- model_sensitivity(
    indicator_time_series = indicator_ts_baltic[,c(1,2)],
    pressure_time_series = new_press,
    current_years = c(2010, 2016)
  ))
  expect_equal(dat$sensitivity, c(NA, NA))
  expect_equal(dat$r_sq, c(NA, NA))
  expect_equal(dat$uncertainty_sen, c(NA, NA))
})


test_that("NaNs in p.value will be scored as NA", {

  ind <- indicator_ts_baltic[c(28:33), ]
  ind[, 3] <- as.numeric(rep(NA))
  press <- pressure_ts_baltic[c(28:33), c(1:6)]
  press[,2] <- as.numeric(rep(NA))
  suppressWarnings(dat <- model_sensitivity(
    indicator_time_series = ind,
    pressure_time_series = press,
    current_years = c(2015, 2016)
  ))

  expect_true(is.na(dat$sensitivity[6]))
  expect_true(is.na(dat$r_sq[6]))
  expect_true(is.na(dat$edf[6]))
  expect_true(is.na(dat$p_value[6]))
})


# Test the new data input validations
test_that("warning and error display", {
  test_ind <- indicator_ts_baltic[ ,1:3]
  test_press <- pressure_ts_baltic[ ,1:3]

  current_df <- data.frame(
    ind = rep(names(test_ind)[-1], each = 2),
    press = rep(names(test_press)[-1], 2),
    start = c(2010, 2010, 2011, 2011),
    end = c(2015, 2015, 2016, 2016)
  )

  test_ind2 <- test_ind
  test_press2 <- test_press
  test_ind2$zooplankton_mean_size <- as.character(test_ind2$zooplankton_mean_size)
  test_press2$nitrogen <- as.character(test_press2$nitrogen)

  test_ind3 <- test_ind
  test_press3 <- test_press3b <-  test_press
  test_ind3$year <- test_ind3$year[c(20:33, 1:19)]
  test_press3$year <- test_press3$year[c(20:33, 1:19)]
  test_press3b$year <- 1:33

  # time series not as data.frames
  expect_error(model_sensitivity(
    indicator_time_series = tibble::as_tibble(test_ind),
    pressure_time_series = test_press,
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = as.matrix(test_press),
    current_years = c(start = 2010, end = 2016)
  ))
  # time series have wrong data types
  expect_error(model_sensitivity(
    indicator_time_series = test_ind2,
    pressure_time_series = test_press,
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press2,
    current_years = c(start = 2010, end = 2016)
  ))
  # the first column (time) in time series arguments same but not sorted
  expect_warning(model_sensitivity(
    indicator_time_series = test_ind3,
    pressure_time_series = test_press3,
    current_years = c(start = 2010, end = 2016)
  ))
  # the first column (time) in time series arguments differ
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press3b,
    current_years = c(start = 2010, end = 2016)
  ))
  # no vector or df provided for current period
  expect_warning(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press
  ))
  # wrong assignments
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years_by_ind_press = c(start = 2010, end = 2016)
  ))
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years = current_df
  ))
  # Missing value in general periods
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years = 2010
  ))
  # 'current_years_by_ind_press' not a data frame
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years_by_ind_press = as.matrix(current_df)
  ))
  # data types in 'current_years_by_ind_press' not correct
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years_by_ind_press = data.frame(
      ind = rep(names(test_ind)[-1], each = 2),
      press = rep(names(test_press)[-1], 2),
      start = as.character(c(2010, 2010, 2011, 2011)),
      end = c(2015, 2015, 2016, 2016)
    )
  ))
  # not all ind/press names listed (incorrect dimensions)
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years_by_ind_press = data.frame(
      ind = rep(names(test_ind)[-1]),
      press = rep("nitrogen", 2),
      start = as.character(c(2010, 2011)),
      end = c(2015, 2016)
    )
  ))
  # not all ind/press names listed or wrong names
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years_by_ind_press = data.frame(
      ind = rep(names(test_ind)[-1], each = 2),
      press = rep(names(test_press)[1:2], 2),
      start = as.character(c(2010, 2010, 2011, 2011)),
      end = c(2015, 2015, 2016, 2016)
    )
  ))
  # specified years fall outside time series range
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years = c(start = 2010, end = 2020)
  ))
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years_by_ind_press = data.frame(
      ind = rep(names(test_ind)[-1], each = 2),
      press = rep(names(test_press)[-1], 2),
      start = c(2010, 2010, 2011, 2011),
      end = c(2015, 2015, 2020, 2020)
    )
  ))
  # starting years are higher than ending years
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years = c(start = 2016, end = 2010)
  ))
  expect_error(model_sensitivity(
    indicator_time_series = test_ind,
    pressure_time_series = test_press,
    current_years_by_ind_press = data.frame(
      ind = rep(names(test_ind)[-1], each = 2),
      press = rep(names(test_press)[-1], 2),
      start = c(2010, 2010, 2016, 2016),
      end = c(2015, 2015, 2011, 2011)
    )
  ))

})

