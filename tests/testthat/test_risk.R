
test_that("output data frame from expert pathway", {
  dat <- risk(
    vulnerability_results = ex_output_vulnerability_expert,
    status_results = ex_expert_status
  )
  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(40,8))
  expect_true(all(dat$pressure %in% ex_output_vulnerability_expert$pressure))
  expect_true(all(dat$indicator %in% ex_output_vulnerability_expert$indicator))

  expect_type(dat$pressure, "character")
  expect_type(dat$indicator, "character")
  expect_type(dat$type, "character")
  expect_type(dat$pathway, "character")
  expect_type(dat$vulnerability, "double")
  expect_type(dat$status, "character")
  expect_type(dat$risk, "double")
  expect_type(dat$uncertainty, "double")

})


test_that("output data frame from model pathway", {
  status_res <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 2010),
    current_years = c(start = 2011, end = 2016)
  )
  dat <- risk(
    vulnerability_results = ex_output_vulnerability_model,
    status_results = status_res
  )
  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(16,8))
  expect_true(all(dat$pressure %in% ex_output_vulnerability_model$pressure))
  expect_true(all(dat$indicator %in% ex_output_vulnerability_model$indicator))

  expect_type(dat$pressure, "character")
  expect_type(dat$indicator, "character")
  expect_type(dat$type, "character")
  expect_type(dat$pathway, "character")
  expect_type(dat$vulnerability, "double")
  expect_type(dat$status, "character")
  expect_type(dat$risk, "double")
  expect_type(dat$uncertainty, "double")

})


test_that("risk scores are correct", {
  dat <- risk(
    vulnerability_results = ex_output_vulnerability_expert,
    status_results = ex_expert_status
  )
  expect_equal(dat$risk[1:10],
    c(-3.75, 0, -1.25, -3.25, 0, -4.75, 0, -2.25, -4.25, 0))

  status_res <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 2010),
    current_years = c(start = 2011, end = 2016)
  )
  dat2 <- risk(
    vulnerability_results = ex_output_vulnerability_model,
    status_results = status_res
  )

  expect_equal(dat2$risk[1:10],
    c(5, 0, 0, 0, 2, 3.25, 0, 0, -5, 0))

})


test_that("warning message occurs", {
  expect_warning(risk(
    vulnerability_results = ex_output_vulnerability_expert,
    status_results = ex_expert_status[-2, ]
  ))
})


test_that("risk special case, vulnerability_results = 0 or NA", {
  vulnerability_test <- ex_output_vulnerability_expert
  vulnerability_test$vulnerability[1] <- 0
  dat <- risk(
    vulnerability_results = vulnerability_test,
    status_results = ex_expert_status
  )
  expect_equal(dat$risk[1], 0)

  vulnerability_test <- ex_output_vulnerability_expert
  vulnerability_test$vulnerability[1] <- NA
  dat <- risk(
    vulnerability_results = vulnerability_test,
    status_results = ex_expert_status
  )

  expect_true(is.na(dat$risk[1]))
})


test_that("trimming works", {
  vulnerability_test <- ex_output_vulnerability_expert
  vulnerability_test$vulnerability[c(7, 11)] <- c(10, -10)
  dat <- risk(
    vulnerability_results = vulnerability_test,
    status_results = ex_expert_status
  )
  expect_equal(dat$risk[c(7,11)], c(10, -10))
})
