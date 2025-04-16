
test_that("output data frame", {
  dat <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 2010),
    current_years = c(start = 2011, end = 2016)
  )

  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(2, 3))

  expect_type(dat$indicator, "character")
  expect_type(dat$status, "character")
  expect_type(dat$score, "double")
})


test_that("results from different threshold ranges", {
  # sd
  dat <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 2010),
    current_years = c(start = 2011, end = 2016)
  )
  expect_equal(dat$score, c(-1, -1))
  # 2sd
  dat <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    range = "2sd"
  )
  expect_equal(dat$score, c(-1, -1))

  # 95th percentile
  dat <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    range = "95percentile"
  )
  expect_equal(dat$score, c(-1, -1))

  # 75th percentile
  dat <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    range = 75
  )
  expect_equal(dat$score, c(-1, -1))
  # mean only
  dat <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    range = "mean_only"
  )
  expect_equal(dat$score, c(-1, -1))
  dat <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    range = "mean_only",
    condition = "<"
  )
  expect_equal(dat$score, c(1, 1))


})


# test different signs and conditions
test_that("different settings condition and sign", {
  dat <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    sign = "+",
    condition = "<"
  )
  expect_equal(dat$score, c(1, 1))

  dat <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    sign = "-",
    condition = "<"
  )

  expect_equal(dat$score, c(-1, -1))

  dat <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    sign = "-",
    condition = ">"
  )
  expect_equal(dat$score, c(1, 1))
})


# Test the new data input validations for time series and periods
test_that("warning and error display in ts and periods", {
  test_dat <- indicator_ts_baltic
  test_dat2 <- test_dat
  test_dat2$zooplankton_mean_size <- as.character(test_dat2$zooplankton_mean_size)
  test_dat3 <- test_dat
  test_dat3$year <- test_dat3$year[c(20:33, 1:19)]

  # 'indicator_time_series' argument not a data.frame
  expect_error(status(
    indicator_time_series = tibble::as_tibble(test_dat),
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(status(
    indicator_time_series = as.matrix(test_dat),
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  # 'indicator_time_series' columns have wrong data types
  expect_error(status(
    indicator_time_series = test_dat2,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  # the first column (time) not sorted
  expect_warning(status(
    indicator_time_series = test_dat3,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  # no vector or df provided for baseline and current period
  suppressWarnings(expect_warning(status(
    indicator_time_series = test_dat
  )))
  # wrong assignments
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years = c(start = 1984, end = 1994),
    current_years_by_ind = c(start = 2010, end = 2016)
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1990, 1995), end = c(1994, 2005)),
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years = c(start = 1984, end = 1994),
    current_years = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1995, 1995), end = c(2000, 2000))
  ))
  # Missing value in general periods
  expect_error(status(
    indicator_time_series = test_dat,
    base_years = 1984,
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years = c(start = 1984, end = 1994),
    current_years = 2010
  ))

  # 'base_years_by_ind' and 'current_years_by_ind' not a data frame
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = tibble::as_tibble(data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1990, 1995), end = c(1994, 2005))),
    current_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1995, 1995), end = c(2000, 2000))
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1990, 1995), end = c(1994, 2005)),
    current_years_by_ind = as.matrix(data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1995, 1995), end = c(2000, 2000)))
  ))

  # data types in 'base_years_by_ind' and 'current_years_by_ind' not correct
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c("1990", "1995"), end = c(1994, 2005)),
    current_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1995, 1995), end = c(2000, 2000))
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1990, 1995), end = c(1994, 2005)),
    current_years_by_ind = as.matrix(data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c("1995", "1995"), end = c(2000, 2000)))
  ))
  # not all pressure names listed or wrong pressure names
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size"),
      start = c(1990), end = c(1994)),
    current_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1995, 1996), end = c(2000, 2000))
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1990, 1995), end = c(1994, 2005)),
    current_years_by_ind = data.frame(
      ind =c("eastern_baltic_cod"),
      start = c(1996), end = c(2000))
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = data.frame(
      ind =c("blub", "eastern_baltic_cod"),
      start = c(1990, 1995), end = c(1994, 2005)),
    current_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1995, 1996), end = c(2000, 2000))
  ))
  # specified years fall outside time series range
  expect_error(status(
    indicator_time_series = test_dat,
    base_years = c(start = 1983, end = 1994),
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2020)
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1980, 1995), end = c(1994, 2005)),
    current_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1995, 1996), end = c(2000, 2016))
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1984, 1995), end = c(1994, 2005)),
    current_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1995, 1996), end = c(2000, 2020))
  ))
  # starting years are higher than ending years
  expect_error(status(
    indicator_time_series = test_dat,
    base_years = c(start = 1994, end = 1984),
    current_years = c(start = 2010, end = 2016)
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2016, end = 2010)
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1994, 1995), end = c(1984, 2005)),
    current_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1995, 1996), end = c(2000, 2016))
  ))
  expect_error(status(
    indicator_time_series = test_dat,
    base_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1984, 1995), end = c(1994, 2005)),
    current_years_by_ind = data.frame(
      ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
      start = c(1995, 2016), end = c(2000, 2005))
  ))

})


test_that("warning message occurs with method settings", {
  expect_warning(status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    range = "blub"
  ))
  expect_warning(status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    range = "sd3"
  ))
  expect_warning(status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    range = "mean_under"
  ))
  expect_warning(status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    range = "mean_only",
    condition = "="
  ))
  expect_warning(status(
    indicator_time_series = indicator_ts_baltic[ , c(1,2)],
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2011, end = 2016),
    sign = ">"
  ))

})


