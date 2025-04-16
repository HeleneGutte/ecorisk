
ex_dat <- ex_expert_sensitivity

# Test output if only the first 3 arguments provided
# (adaptive_capacities, uncertainty_sens, uncertainty_ac not provided)
test_that("output data frame with optional arguments NULL", {
  dat <- calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[ ,4:8]
  )
  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(40, 14))
  expect_true(all(dat$pressure %in% ex_dat$pressure))
  expect_true(all(dat$indicator %in% ex_dat$indicator))
  expect_type(dat$type, "character")
  expect_true(all(dat$type == "direct"))
  expect_type(dat$pathway, "character")
  expect_true(all(dat$pathway == "expert"))
  expect_type(dat$sensitivity, "double")
  expect_type(dat$adaptive_capacity, "double")
  expect_true(all(dat$adaptive_capacity == 0)) # ac is 0 if not provided
  expect_true(all(is.na(dat$uncertainty_sens))) # is NA if not provided
  expect_true(all(is.na(dat$uncertainty_ac)))  # is NA if not provided
  expect_true(all(is.na(dat$ac_original.ac_general)))  # is NA if not provided
  # Check values of sensitivity score
  expect_equal(mean(dat$sensitivity), -2.06875, tolerance = 0.01)
})


# Test output if also adaptive_capacities, uncertainty_sens, uncertainty_ac
# are provided
test_that("full output data frame", {
  dat <- calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    type = ex_dat$type,
    sensitivity_traits = ex_dat[ ,4:8],
    adaptive_capacities =ex_dat[ ,9:13],
    uncertainty_sens = ex_dat[ ,14:18],
    uncertainty_ac = ex_dat[ ,19:23]
  )
  expect_equal(dim(dat), c(40, 18))
  expect_true(all(dat$pressure %in% ex_dat$pressure))
  expect_true(all(dat$indicator %in% ex_dat$indicator))
  expect_type(dat$pathway, "character")
  expect_type(dat$sensitivity, "double")
  expect_type(dat$adaptive_capacity, "double")
  expect_type(dat$uncertainty_sens, "double")
  expect_type(dat$uncertainty_ac, "double")
  expect_true(all(names(dat) == c(
    "indicator", "pressure", "type", "pathway", "sensitivity",
    "adaptive_capacity", "uncertainty_sens",
    "uncertainty_ac", "sens_original.sens_feeding",
    "sens_original.sens_behaviour", "sens_original.sens_reproduction",
    "sens_original.sens_habitat", "sens_original.sens_general",
    "ac_original.ac_feeding", "ac_original.ac_behaviour",
    "ac_original.ac_reproduction", "ac_original.ac_habitat",
    "ac_original.ac_general"))
  )
  # Check calculated scores
  expect_equal(mean(dat$sensitivity), -2.06875, tolerance = 0.01)
  expect_equal(mean(dat$adaptive_capacity), 0.0125, tolerance = 0.01)
  expect_equal(mean(dat$uncertainty_sens), 1.5, tolerance = 0.1)
  expect_equal(mean(dat$uncertainty_ac), 1.8, tolerance = 0.1)
})


# Test number of columns in output depending on the input settings
test_that("check columns in output", {
  dat1 <- calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8]
  )
  dat2 <- calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 8]
  )
  dat3 <- calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8],
    adaptive_capacities = ex_dat[, 9:13]
  )
  dat4 <- calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8],
    adaptive_capacities = ex_dat[, 9:13],
    uncertainty_sens  = ex_dat[, 14:18],
    uncertainty_ac = ex_dat[, 19:23]
  )

  expect_equal(dim(dat1), c(40, 14))
  expect_equal(dim(dat2), c(40, 10))
  expect_equal(dim(dat3), c(40, 18))
  expect_equal(dim(dat4), c(40, 18))
})


# Test method change
test_that("results from different methods", {
  dat <- calc_sensitivity(
    indicators = ex_dat$indicator[11:30],
    pressures = ex_dat$pressure[11:30],
    type = ex_dat$type[11:30],
    sensitivity_traits = ex_dat[11:30 ,4:8],
    adaptive_capacities =ex_dat[11:30 ,9:13],
    uncertainty_sens = ex_dat[11:30 ,14:18],
    uncertainty_ac = ex_dat[11:30 ,19:23],
    method = "mean"
  )
  expect_equal(dat$sensitivity, c(-1.5, -1, -1.5, 0.5, -2.25, -3.25, -2, -3.25, 1.5, -3, -2.25, -2,
    -2.5, -0.5, -3.5, -4.25, -3, -3.75, -1.75, -4.5))
  expect_equal(dat$uncertainty_sens, c(1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1))
  expect_equal(dat$adaptive_capacity, c(0, 0, 0, 0.75, -0.75, -1, 0, 0, 1, -1, -0.75, 0,
    -0.75, 0, -1, -1, -1, -1, 1, -1))
  expect_equal(dat$uncertainty_ac, c(1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1))

  dat <- calc_sensitivity(
    indicators = ex_dat$indicator[11:30],
    pressures = ex_dat$pressure[11:30],
    type = ex_dat$type[11:30],
    sensitivity_traits = ex_dat[11:30 ,4:8],
    adaptive_capacities =ex_dat[11:30 ,9:13],
    uncertainty_sens = ex_dat[11:30 ,14:18],
    uncertainty_ac = ex_dat[11:30 ,19:23],
    method = "median"
  )
  expect_equal(dat$sensitivity, c(-1.5, -1, -1.5, 0, -2, -3, -2.5, -3, 1.5, -3.5, -2.5, -2, -3, 0,
    -4.5, -4.5, -3, -3.5, -1.5, -5))
  expect_equal(dat$uncertainty_sens, c(1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1))
  expect_equal(dat$adaptive_capacity, c(0, 0, 0, 1, -1, -1, 0, 0, 1, -1, -1, 0, -1, 0, -1, -1, -1, -1, 1, -1))
  expect_equal(dat$uncertainty_ac, c(1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1))

  dat <- calc_sensitivity(
    indicators = ex_dat$indicator[11:30],
    pressures = ex_dat$pressure[11:30],
    type = ex_dat$type[11:30],
    sensitivity_traits = ex_dat[11:30 ,4:8],
    adaptive_capacities =ex_dat[11:30 ,9:13],
    uncertainty_sens = ex_dat[11:30 ,14:18],
    uncertainty_ac = ex_dat[11:30 ,19:23],
    method = "minimum"
  )
  expect_equal(dat$sensitivity, c(-3, -2, -3, 0, -5, -4, -3, -4, 0, -5, -4, -4, -4, -2, -5, -5, -5, -5, -3, -5))
  expect_equal(dat$uncertainty_sens, c(1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1))
  expect_equal(dat$adaptive_capacity, c(0, 0, 0, 0, -1, -1, 0, 0, 1, -1, -1, 0, -1, 0, -1, -1, -1, -1, 1, -1))
  expect_equal(dat$uncertainty_ac, c(1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1))

  dat <- calc_sensitivity(
    indicators = ex_dat$indicator[11:30],
    pressures = ex_dat$pressure[11:30],
    type = ex_dat$type[11:30],
    sensitivity_traits = ex_dat[11:30 ,4:8],
    adaptive_capacities =ex_dat[11:30 ,9:13],
    uncertainty_sens = ex_dat[11:30 ,14:18],
    uncertainty_ac = ex_dat[11:30 ,19:23],
    method = "maximum"
  )
  expect_equal(dat$sensitivity, c(0, 0, 0, 2, 0, -3, 0, -3, 3, 0, 0, 0, 0, 0, 0, -3, -1, -3, -1, -3))
  expect_equal(dat$uncertainty_sens, c(1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1))
  expect_equal(dat$adaptive_capacity, c(0, 0, 0, 1, 0, -1, 0, 0, 1, -1, 0, 0, 0, 0, -1, -1, -1, -1, 1, -1))
  expect_equal(dat$uncertainty_ac, c(1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 2, 2, 1))

})


# Check output if inputs are always single variables (so no aggregation takes place)
test_that("general scores", {
  dat <- calc_sensitivity(
    indicators = ex_dat$indicator[c(1:10, 31:40)],
    pressures = ex_dat$pressure[c(1:10, 31:40)],
    type = ex_dat$type[c(1:10, 31:40)],
    sensitivity_traits = ex_dat$sens_general[c(1:10, 31:40)],
    adaptive_capacities = ex_dat$ac_general[c(1:10, 31:40)],
    uncertainty_sens = ex_dat$uncertainty_sens_general[c(1:10, 31:40)],
    uncertainty_ac = ex_dat$uncertainty_ac_general[c(1:10, 31:40)],
    method = "mean"
  )
  #sensitivity
  expect_equal(dat$sensitivity, ex_dat$sens_general[c(1:10, 31:40)])
  expect_equal(dat$sens_original.sens_general, ex_dat$sens_general[c(1:10, 31:40)])
  #adaptive capacity
  expect_equal(dat$adaptive_capacity, ex_dat$ac_general[c(1:10, 31:40)])
  expect_equal(dat$ac_original.ac_general, ex_dat$ac_general[c(1:10, 31:40)])
  #uncertainty
  expect_equal(dat$uncertainty_sens, ex_dat$uncertainty_sens_general[c(1:10, 31:40)])
  expect_equal(dat$uncertainty_ac, ex_dat$uncertainty_ac_general[c(1:10, 31:40)])

})


# Compare returned scores with the example data set 'ex_output_calc_sensitivity'
test_that("compare output with demo data", {
  dat <- calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    type = ex_dat$type,
    sensitivity_traits = ex_dat[ ,4:8],
    adaptive_capacities =ex_dat[ ,9:13],
    uncertainty_sens = ex_dat[ ,14:18],
    uncertainty_ac = ex_dat[ ,19:23]
  )

  #sensitivity
  expect_equal(dat$sensitivity, ex_output_calc_sensitivity$sensitivity)
  #adaptive capacity
  expect_equal(dat$adaptive_capacity, ex_output_calc_sensitivity$adaptive_capacity)
  #uncertainty
  expect_equal(dat$uncertainty_sens, ex_output_calc_sensitivity$uncertainty_sens)
  expect_equal(dat$uncertainty_ac, ex_output_calc_sensitivity$uncertainty_ac)

})


# Test the new data input validations
test_that("warning and error display", {

  # 'indicators' not character
  expect_error(calc_sensitivity(
    indicators = 1:length(ex_dat$indicator),
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8]
  ))
  # 'pressures' not character
  expect_error(calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = 1:length(ex_dat$pressure),
    sensitivity_traits = ex_dat[, 4:8]
  ))
  # 'indicators' and 'pressures' not same length
  expect_error(calc_sensitivity(
    indicators = ex_dat$indicator[1:10],
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8]
  ))
  # 'sensitivity_traits' not numeric
  expect_error(calc_sensitivity(
    indicators = ex_dat$indicator[1:10],
    pressures = ex_dat$pressure[1:10],
    sensitivity_traits = letters[1:10]
  ))
  # 'sensitivity_traits' data frame, but not all numeric
  expect_error(calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, c(1,4:8)]
  ))
  # 'adaptive_capacities' not NULL/numeric vector/data frame or column
  # incorrect (if > 1, not same as sensitivity_traits)
  expect_error(calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8],
    adaptive_capacities = "a"
  ))
  expect_error(calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8],
    adaptive_capacities = ex_dat[, 9:11]
  ))
  # 'uncertainty_sens' not NULL/numeric vector/data frame or column
  # incorrect (if > 1, not same as sensitivity_traits)
  expect_error(calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8],
    uncertainty_sens = "a"
  ))
  expect_error(calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8],
    uncertainty_sens = ex_dat[, 14:15]
  ))
  # 'uncertainty_ac' not NULL/numeric vector/data frame or column
  # incorrect (if > 1, not same as sensitivity_traits)
  expect_error(calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8],
    uncertainty_ac = "a"
  ))
  expect_error(calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8],
    uncertainty_ac = ex_dat[, 19:21]
  ))
  # wrong 'method'
  expect_warning(calc_sensitivity(
    indicators = ex_dat$indicator,
    pressures = ex_dat$pressure,
    sensitivity_traits = ex_dat[, 4:8],
    method = "blub"
  ))

})


