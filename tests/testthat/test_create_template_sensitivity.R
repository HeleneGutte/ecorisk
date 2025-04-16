
ind_names <- c("phytoplankton", "herring", "cod",
  "seabirds")
press_names <- c("temperature", "salinity", "oxygen",
  "nutrient", "fishing")

test_that("output data frame with default setting", {
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names
  )

  expect_s3_class(dat, "data.frame")
  expect_equal(dim(dat), c(20, 7))
  expect_true(all(dat$pressures %in% press_names))
  expect_true(all(dat$indicators %in% ind_names))
  expect_true(all(dat$type %in% c("direct")))
  expect_type(dat$pressure, "character")
  expect_type(dat$indicator, "character")
  expect_type(dat$type, "character")
  expect_true(all(names(dat) == c("indicator", "pressure",
    "type", "sens_general", "ac_general",
    "uncertainty_sens", "uncertainty_ac")))
})

test_that("output data frame with different settings", {
  # Testing type
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    type = c("direct", "indirect"))
  expect_equal(dim(dat), c(40, 7))

  # Testing n_sensitivity_traits
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 5)
  expect_equal(dim(dat), c(20, 11))
  expect_true(all(names(dat)[4:8] == paste0("sens_trait_", 1:5)))

  # Excluding adaptive_capacity and uncertainty
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    adaptive_capacity = FALSE,
    uncertainty = FALSE
  )
  expect_equal(dim(dat), c(20, 4))

  ### Testing adaptive_capacity (without uncertainity)
  # 1 AC column for several sensitivity traits
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 2,
    adaptive_capacity = TRUE,
    mode_adaptive_capacity = "general",
    uncertainty = FALSE
  )
  expect_equal(dim(dat), c(20, 6))
  # Trait-specific AC columns for several sensitivity traits
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 3,
    adaptive_capacity = TRUE,
    mode_adaptive_capacity = "trait",
    uncertainty = FALSE
  )
  expect_equal(dim(dat), c(20, 9))
  expect_true(all(names(dat)[7:9] == paste0("ac_trait_", 1:3)))
  # 1 AC column for 1 sensitivity traits (despite mode 'trait')
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 1,
    adaptive_capacity = TRUE,
    mode_adaptive_capacity = "trait",
    uncertainty = FALSE
  )
  expect_equal(dim(dat), c(20, 5))

  ### Testing uncertainity (without adaptive_capacity )
  # 1 sensitivity-specific UNC column for several sensitivity traits
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 2,
    adaptive_capacity = FALSE,
    uncertainty = TRUE,
    mode_uncertainty = "general"
  )
  expect_equal(dim(dat), c(20, 6))
  # Trait-specific sensitivity-UNC columns for several sensitivity traits
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 3,
    adaptive_capacity = FALSE,
    uncertainty = TRUE,
    mode_uncertainty = "trait"
  )
  expect_equal(dim(dat), c(20, 9))
  expect_true(all(names(dat)[7:9] == paste0("uncertainty_sens_trait_", 1:3)))
  # 1 sensitivity-specific UNC column for 1 sensitivity traits (despite mode 'trait')
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 1,
    adaptive_capacity = FALSE,
    uncertainty = TRUE,
    mode_uncertainty = "trait"
  )
  expect_equal(dim(dat), c(20, 5))

  ### Testing uncertainty and adaptive_capacity
  # Only 1 AC-specific UNC column if AC is set to 'general'
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 2,
    adaptive_capacity = TRUE,
    mode_adaptive_capacity = "general",
    uncertainty = TRUE,
    mode_uncertainty = "trait"
  )
  expect_equal(dim(dat), c(20, 9))
  # Only 1 AC-specific UNC, even if AC is set to 'trait'
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 2,
    adaptive_capacity = TRUE,
    mode_adaptive_capacity = "trait",
    uncertainty = TRUE,
    mode_uncertainty = "general"
  )
  expect_equal(dim(dat), c(20, 9))
  # AC-specific UNC columns of both are set to 'trait' and n_sensitivity_traits > 1
  dat <- create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 2,
    adaptive_capacity = TRUE,
    mode_adaptive_capacity = "trait",
    uncertainty = TRUE,
    mode_uncertainty = "trait"
  )
  expect_equal(dim(dat), c(20, 11))

})


test_that("error display", {

  # Wrong number of sensitivity traits
  expect_error(create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 0))
  # Sensitivity trait not numeric
  expect_error(create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = TRUE))
  # If sens_trait > and AC = TRUE, mode_adaptive_capacity checked
  expect_error(create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 5,
    mode_adaptive_capacity = "blub"))
  # If sens_trait > and unc = TRUE, mode_uncertainty checked
  expect_error(create_template_sensitivity(
    indicators = ind_names,
    pressures = press_names,
    n_sensitivity_traits = 2,
    mode_uncertainty = "blub"))

})

