test_that("output data frame", {
  a <- ex_expert_exposure
  dat <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7:9]
  )

  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(5, 3))
  expect_true(all(dat$pressure %in% a$pressure))
  expect_type(dat$exposure, "double")
  expect_type(dat$uncertainty, "double")
})


test_that("uncertainty is NA", {
  a <- ex_expert_exposure
  dat <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5]
  )
  expect_equal(dat$uncertainty[1], NA)
})


test_that("results are correct for every method", {
  a <- ex_expert_exposure
  # without probabilities
  dat <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7:9],
    method = "median"
  )
  dat2 <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7:9],
    method = "mean"
  )
  dat3 <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7:9],
    method = "minimum"
  )
  dat4 <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7:9],
    method = "maximum"
  )

  expect_false(all(dat$exposure == dat2$exposure))
  expect_equal(dat$exposure, c(2.5, 2, 1, 2, 4.5))
  expect_equal(dat2$exposure, c(2.75, 2.25, 1.25, 2.25, 4))
  expect_equal(dat3$exposure, c(1,1,1,2,2))
  expect_equal(dat4$exposure, c(5, 4, 2, 3, 5))

  # with probabilities
  a$probabilities <- c(1, 0.8, 0.5, 0.8, 1)
  dat <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7:9],
    method = "median",
    probabilities = a$probabilities
  )
  dat2 <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7:9],
    method = "mean",
    probabilities = a$probabilities
  )
  dat3 <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7:9],
    method = "minimum",
    probabilities = a$probabilities
  )
  dat4 <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7:9],
    method = "maximum",
    probabilities = a$probabilities
  )
  expect_equal(dat$exposure, c(2.5, 1.6, 0.5, 1.6, 4.5))
  expect_equal(dat2$exposure, c(2.75, 1.8, 0.625, 1.8, 4))
  expect_equal(dat3$exposure, c(1, 0.8, 0.5, 1.6, 2))
  expect_equal(dat4$exposure, c(5, 3.2, 1, 2.4, 5))

  # Uncertainty results
  expect_equal(dat$uncertainty, c(3, 2, 3, 2, 1))
  expect_equal(round(dat2$uncertainty, digits = 4), c(2.3333, 2, 2.6667, 2.3333, 1.3333))
  expect_equal(dat3$uncertainty, c(1, 2, 2, 2, 1))
  expect_equal(dat4$uncertainty, c(3, 2, 3, 3,2))

  # Uncertainty results with single uncertainty vector
  dat <- calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7],
    method = "median",
    probabilities = a$probabilities
  )
  expect_equal(dat$uncertainty, c(1, 2, 2, 2, 1))


})


test_that("error and warning display", {
  a <- ex_expert_exposure
  expect_warning(calc_exposure(
    pressures = a$pressure,
    components = a[ ,2:5],
    uncertainty = a[ ,7:9],
    method = "blub")
  )

  # pressure is not a character
  expect_error(calc_exposure(pressures = 1:3, components = a[ ,2:5]))
  # component is not completely numeric
  expect_error(calc_exposure(pressures = a$pressure, components = a[ ,c(1,5)]))
  # probabilities is not a vector
  expect_error(calc_exposure(pressures = a$pressure, components = a[ ,c(2,5)],
    probabilities = data.frame(c(1,0,0,1,0.5), c(1,1,1,1,1))))
  # probabilities is not a numeric vector
  expect_error(calc_exposure(pressures = a$pressure, components = a[ ,c(2,5)],
    probabilities = rep(TRUE, 5)))
  # probabilities range outside of 0-1
  expect_error(calc_exposure(pressures = a$pressure, components = a[ ,c(2,5)],
    probabilities = c(2,4,0,1,0.5)))
  # uncertainity is not completely numeric
  expect_error(calc_exposure(pressures = a$pressure, components = a[ ,c(2,5)],
    uncertainty = a[ ,c(1,8,9)]))
})




