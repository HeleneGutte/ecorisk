test_that("output data frame", {
  press_names <- c("temperature", "salinity", "oxygen", "nutrient",
    "fishing")
  dat <- create_template_exposure(
    pressures = press_names ,
    mode_uncertainty = "component")

  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(5, 9))
  expect_true(all(dat$pressures %in% press_names))
  expect_type(dat$pressure, "character")


  dat <- create_template_exposure(
    pressures = press_names ,
    n_components = 2,
    mode_uncertainty = "general")
  expect_equal(dim(dat), c(5, 4))
})


test_that("warning appears", {
  press_names <- c("temperature", "salinity", "oxygen", "nutrient",
    "fishing")
  expect_warning(create_template_exposure(pressures = press_names ,
    n_components = 4, mode_uncertainty = "blub"))
})
