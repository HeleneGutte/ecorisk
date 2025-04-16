test_that("output data frame expert path", {
  dat <- aggregate_risk(risk_results = ex_output_risk_expert)

  expect_type(dat, "list")
  expect_s3_class(dat[[1]], "data.frame")
  expect_s3_class(dat[[2]], "data.frame")
  expect_s3_class(dat[[3]], "data.frame")
  expect_equal(dim(dat[[1]]), c(30, 5))
  expect_equal(dim(dat[[2]]), c(24, 5))
  expect_equal(dim(dat[[3]]), c(6, 4))

  expect_true(all(dat[[1]]$pressure %in% ex_output_risk_expert$pressure))
  expect_true(all(dat[[2]]$indicator %in% ex_output_risk_expert$indicator))

  expect_type(dat[[1]]$pressure, "character")
  expect_type(dat[[1]]$type, "character")
  expect_type(dat[[1]]$risk, "double")
  expect_type(dat[[1]]$pathway, "character")
  expect_type(dat[[1]]$uncertainty, "double")
  expect_type(dat[[2]]$indicator, "character")
  expect_type(dat[[2]]$type, "character")
  expect_type(dat[[2]]$risk, "double")
  expect_type(dat[[2]]$pathway, "character")
  expect_type(dat[[2]]$uncertainty, "double")
  expect_type(dat[[3]]$type, "character")
  expect_type(dat[[3]]$risk, "double")
  expect_type(dat[[3]]$pathway, "character")
  expect_type(dat[[3]]$uncertainty, "double")

})

test_that("output data frame model path", {
  dat <- aggregate_risk(risk_results = ex_output_risk_model)

  expect_type(dat, "list")
  expect_s3_class(dat[[1]], "data.frame")
  expect_s3_class(dat[[2]], "data.frame")
  expect_s3_class(dat[[3]], "data.frame")
  expect_equal(dim(dat[[1]]), c(32, 5))
  expect_equal(dim(dat[[2]]), c(8, 5))
  expect_equal(dim(dat[[3]]), c(4, 4))

  expect_true(all(dat[[1]]$pressure %in% ex_output_risk_model$pressure))
  expect_true(all(dat[[2]]$indicator %in% ex_output_risk_model$indicator))

  expect_type(dat[[1]]$pressure, "character")
  expect_type(dat[[1]]$type, "character")
  expect_type(dat[[1]]$risk, "double")
  expect_type(dat[[1]]$pathway, "character")
  expect_type(dat[[1]]$uncertainty, "double")
  expect_type(dat[[2]]$indicator, "character")
  expect_type(dat[[2]]$type, "character")
  expect_type(dat[[2]]$risk, "double")
  expect_type(dat[[2]]$pathway, "character")
  expect_type(dat[[2]]$uncertainty, "double")
  expect_type(dat[[3]]$type, "character")
  expect_type(dat[[3]]$risk, "double")
  expect_type(dat[[3]]$pathway, "character")
  expect_type(dat[[3]]$uncertainty, "double")

})


test_that("results are correct depending on the chosen method", {
  # mean
  dat <- aggregate_risk(
    risk_results = ex_output_risk_expert,
    method = "mean"
  )
  expect_equal(round(dat[[1]]$risk[1:10], digits = 2),
    c(-5.06, -6.25, -5.06, -6.25, -5.66, -5.66, -1.25, -2.14, -1.25, -2.14))
  expect_equal(round(dat[[2]]$risk[1:10], digits = 1),
    c(-5.0, -7.2, -5, -7.2, -6.1, -6.1, -3.3, -4.5, -3.3, -4.5))
  expect_equal(round(dat[[3]]$risk, digits = 2),
    c(-3.79, -3.79, -2.98, -2.98, -4.61, -4.61))

  expect_equal(round(dat[[1]]$uncertainty[1:10], digits = 2),
    c(1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 2.08, 2.08, 2.08, 2.08))
  expect_equal(round(dat[[2]]$uncertainty[1:10], digits = 2),
    c(1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62, 1.62))
  expect_equal(round(dat[[3]]$uncertainty, digits = 2),
    c(1.78, 1.78, 1.78, 1.78, 1.78, 1.78))

  # median
  dat <- aggregate_risk(
    risk_results = ex_output_risk_expert,
    method = "median"
  )
  expect_equal(round(dat[[1]]$risk[1:10], digits = 1),
    c(-6.2, -7.6, -6.2, -7.6, -7.1, -7.1, -1, -3.6, -1, -3.6))
  expect_equal(round(dat[[2]]$risk[1:10], digits = 1),
    c(-5.2, -7.0, -5.2, -7.0, -6.4, -6.4, -3.4, -5.5, -3.4, -5.5))
  expect_equal(round(dat[[3]]$risk, digits = 1),
    c(-3.8, -3.8, -2.3, -2.3, -4.4, -4.4))

  expect_equal(round(dat[[1]]$uncertainty[1:10], digits = 2),
    c(1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 2.00, 2.00, 2.00, 2.00))
  expect_equal(round(dat[[2]]$uncertainty[1:10], digits = 2),
    c(1.33, 1.33, 1.33, 1.33, 1.33, 1.33, 1.33, 1.33, 1.33, 1.33))
  expect_equal(round(dat[[3]]$uncertainty, digits = 2),
    c(1.67, 1.67, 1.67, 1.67, 1.67, 1.67))

  # minimum
  dat <- aggregate_risk(
    risk_results = ex_output_risk_expert,
    method = "minimum"
  )
  expect_equal(round(dat[[1]]$risk[1:10], digits = 2),
    c(-7.75, -9.75, -7.75, -9.75, -9.75, -9.75, -3.25, -4.25, -3.25, -4.25))
  expect_equal(round(dat[[2]]$risk[1:10], digits = 2),
    c(-7.75, -9.75, -7.75, -9.75, -9.75, -9.75, -5.50, -8.00, -5.50, -8.00))
  expect_equal(round(dat[[3]]$risk, digits = 2),
    c(-9.75, -9.75, -7.75, -7.75, -9.75, -9.75))

  expect_equal(round(dat[[1]]$uncertainty[1:10], digits = 2),
    c(1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 2.00, 2.00, 2.00, 2.00))
  expect_equal(round(dat[[2]]$uncertainty[1:10], digits = 2),
    c(1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25))
  expect_equal(round(dat[[3]]$uncertainty, digits = 2),
    c(1.25, 1.25, 1.25, 1.25, 1.25, 1.25))

  # maximum
  dat <- aggregate_risk(
    risk_results = ex_output_risk_expert,
    method = "maximum"
  )
  expect_equal(round(dat[[1]]$risk[1:10], digits = 2),
    c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.31, 2.94, 0.31, 2.94))
  expect_equal(round(dat[[2]]$risk[1:10], digits = 2),
    c(-2.06, -4.00, -2.06, -4.00, -2.06, -2.06,  0.31,  2.94,  0.31,  2.94))
  expect_equal(round(dat[[3]]$risk, digits = 2),
    c(2.94, 2.94, 0.31, 0.31, 2.94, 2.94))

  expect_equal(round(dat[[1]]$uncertainty[1:10], digits = 2),
    c(1.92, 1.92, 1.92, 1.92, 1.92, 1.92, 2.33, 2.33, 2.33, 2.33))
  expect_equal(round(dat[[2]]$uncertainty[1:10], digits = 2),
    c(2.17, 2.17, 2.17, 2.17, 2.17, 2.17, 2.17, 2.17, 2.17, 2.17))
  expect_equal(round(dat[[3]]$uncertainty, digits = 2),
    c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5))

  # sum
  dat <- aggregate_risk(
    risk_results = ex_output_risk_expert,
    method = "sum"
  )
  expect_equal(round(dat[[1]]$risk[1:10], digits = 2),
    c(-20.25, -25.00, -20.25, -25.00, -45.25, -45.25, -5.00, -8.56, -5.00, -8.56))
  expect_equal(round(dat[[2]]$risk[1:10], digits = 2),
    c(-25.19, -36.00, -25.19, -36.00, -61.19, -61.19, -16.31, -22.50, -16.31, -22.50))
  expect_equal(round(dat[[3]]$risk, digits = 2),
    c(-151.75, -151.75, -59.50, -119.00, -92.25, -184.50))

  expect_equal(round(dat[[1]]$uncertainty[1:10], digits = 2),
    c(5.67, 5.67, 5.67, 5.67, 11.33, 11.33, 8.33, 8.33, 8.33, 8.33))
  expect_equal(round(dat[[2]]$uncertainty[1:10], digits = 2),
    c(8.08, 8.08, 8.08, 8.08, 16.17, 16.17, 8.08, 8.08, 8.08, 8.08))
  expect_equal(round(dat[[3]]$uncertainty, digits = 2),
    c(71.33, 71.33, 35.67, 71.33, 35.67, 71.33))

})


test_that("warning message appears correctly", {
  expect_warning(aggregate_risk(
    risk_results = ex_output_risk_expert,
    method = "blub"
  ))

})
