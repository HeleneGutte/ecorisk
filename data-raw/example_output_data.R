
# Output from calc_exposure() ----
ex_output_calc_exposure <- calc_exposure(
  pressures = ex_expert_exposure$pressure,
  components = ex_expert_exposure[, 2:5],
  uncertainties = ex_expert_exposure[, 6:9]
)

# Output from calc_sensitivity() ----
ex_output_calc_sensitivity <- calc_sensitivity(
  indicators = ex_expert_sensitivity$indicator,
  pressures = ex_expert_sensitivity$pressure,
  type = ex_expert_sensitivity$type,
  sensitivity_traits = ex_expert_sensitivity[, 4:8],
  adaptive_capacities = ex_expert_sensitivity[, 9:13],
  uncertainty_sens = ex_expert_sensitivity[, 14:18],
  uncertainty_ac = ex_expert_sensitivity[, 19:23],
  method = "mean"
)


# Output from model_exposure ----
ex_output_model_exposure <- model_exposure(
  pressure_time_series = pressure_ts_baltic,
  base_years = c(start = 1984, end = 2014),
  current_years = c(start = 2010, end = 2016)
)

# Output from model_sensitivity ----
ex_output_model_sensitivity <- model_sensitivity(
  indicator_time_series = indicator_ts_baltic,
  pressure_time_series = pressure_ts_baltic,
  current_years = c(start = 2010, end = 2016)
)


# Output from vulnerability ----
ex_output_vulnerability_expert <- vulnerability(
  exposure_results = ex_output_calc_exposure,
  sensitivity_results = ex_output_calc_sensitivity
)

ex_output_vulnerability_model <- vulnerability(
  exposure_results = ex_output_model_exposure,
  sensitivity_results= ex_output_model_sensitivity
)

# Output from status ----
ex_output_status <- status(
  indicator_time_series = indicator_ts_baltic,
  base_years = c(start = 1984, end = 2010),
  current_years = c(start = 2011, end = 2016)
)

# Output from risk ----
ex_output_risk_expert <- risk(
  vulnerability_results = ex_output_vulnerability_expert,
  status_results = ex_expert_status
)

ex_output_risk_model <- risk(
  vulnerability_results = ex_output_vulnerability_model,
  status_results = ex_output_status
)

# Output from aggregation ----
ex_output_aggregate_risk_expert <- aggregate_risk(risk_results = ex_output_risk_expert)
ex_output_aggregate_risk_model <- aggregate_risk(risk_results = ex_output_risk_model)




# Use datasets ----
usethis::use_data(
  ex_output_calc_exposure,
  ex_output_calc_sensitivity,
  ex_output_model_exposure,
  ex_output_model_sensitivity,
  ex_output_vulnerability_expert,
  ex_output_vulnerability_model,
  ex_output_status,
  ex_output_risk_expert,
  ex_output_risk_model,
  ex_output_aggregate_risk_expert,
  ex_output_aggregate_risk_model,
  overwrite = TRUE
)

