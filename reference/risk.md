# Calculate Risk Scores Using Expert-Based or Model-Derived Vulnerability and Status Scores

The `risk` function calculates risk scores using the output from of the
[`status`](https://helenegutte.github.io/ecorisk/reference/status.md)
and the
[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)
functions. For each state indicator-pressure combination the function
adds the status score to the vulnerability score to derive the risk
score.

## Usage

``` r
risk(vulnerability_results, status_results)
```

## Arguments

- vulnerability_results:

  A data frame with the output from the `vulnerability` function.

- status_results:

  A data frame with status scores for each state indicator. The first
  column MUST contain the indicator names. The second and third column
  have to be named status and score.

## Value

a data frame containing the exposure, sensitivity, adaptive capacity,
vulnerability, and risk scores as well as their associated uncertainty
for each pressure - state indicator - type combination.

## Details

Final risk scores are in a range from -10 (severe risk for the state
indicator due to the assessed pressure) to +10 (good opportunities for
the state indicator due to the assessed pressure). The risk scores are
specific for each combination of state indicator and pressure and do NOT
take into account cumulative effects. The risk scores can be aggregated
in an additive manner with the
[`aggregate_risk`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md)
function.

## See also

[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md),
[`status`](https://helenegutte.github.io/ecorisk/reference/status.md),
[`aggregate_risk`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md)

## Examples

``` r
# Using demo output data from the vulnerability() and status() functions:
risk(
  vulnerability_results = ex_output_vulnerability_model,
  status_results = ex_output_status
)
#>                indicator      pressure            type pathway vulnerability
#> 1  zooplankton_mean_size      nitrogen direct_indirect   model          6.00
#> 2  zooplankton_mean_size   phosphorous direct_indirect   model          0.00
#> 3  zooplankton_mean_size surf_temp_sum direct_indirect   model          0.00
#> 4  zooplankton_mean_size  bot_temp_ann direct_indirect   model          0.00
#> 5  zooplankton_mean_size  surf_sal_sum direct_indirect   model          3.00
#> 6  zooplankton_mean_size   bot_sal_ann direct_indirect   model          4.25
#> 7  zooplankton_mean_size   bot_oxy_ann direct_indirect   model          0.00
#> 8  zooplankton_mean_size   fishing_cod direct_indirect   model          0.00
#> 9     eastern_baltic_cod      nitrogen direct_indirect   model         -4.00
#> 10    eastern_baltic_cod   phosphorous direct_indirect   model          0.00
#> 11    eastern_baltic_cod surf_temp_sum direct_indirect   model         -3.25
#> 12    eastern_baltic_cod  bot_temp_ann direct_indirect   model         -5.75
#> 13    eastern_baltic_cod  surf_sal_sum direct_indirect   model         -5.00
#> 14    eastern_baltic_cod   bot_sal_ann direct_indirect   model          0.00
#> 15    eastern_baltic_cod   bot_oxy_ann direct_indirect   model         -3.75
#> 16    eastern_baltic_cod   fishing_cod direct_indirect   model         -4.25
#>       status  risk uncertainty
#> 1  undesired  5.00         2.0
#> 2  undesired  0.00         1.5
#> 3  undesired  0.00         1.5
#> 4  undesired  0.00         1.5
#> 5  undesired  2.00         2.0
#> 6  undesired  3.25         2.0
#> 7  undesired  0.00         1.5
#> 8  undesired  0.00         2.0
#> 9  undesired -5.00         2.0
#> 10 undesired  0.00         1.5
#> 11 undesired -4.25         1.5
#> 12 undesired -6.75         1.5
#> 13 undesired -6.00         2.0
#> 14 undesired  0.00         2.0
#> 15 undesired -4.75         1.5
#> 16 undesired -5.25         2.0

# \donttest{
  ### Demo Expert-Based Pathway
  # - using the example scoring datasets 'ex_expert_exposure',
  #   'ex_expert_sensitivity' and 'ex_expert_status'

  # Calculate (mean) exposure score:
  exp_expert <- calc_exposure(
    pressures = ex_expert_exposure$pressure,
    components = ex_expert_exposure[ ,2:5],
    uncertainty = ex_expert_exposure[ ,6:9],
    method = "mean" # default
  )
  # Calculate (mean) sensitivity (and adaptive capacity) score:
  sens_ac_expert <- calc_sensitivity(
    indicators = ex_expert_sensitivity$indicator,
    pressures = ex_expert_sensitivity$pressure,
    type = ex_expert_sensitivity$type,
    sensitivity_traits = ex_expert_sensitivity[ ,4:8],
    adaptive_capacities = ex_expert_sensitivity[ ,9:13],
    uncertainty_sens = ex_expert_sensitivity[ ,14:18],
    uncertainty_ac = ex_expert_sensitivity[ ,19:23],
    method = "mean" # default
  )
  # Calculate (mean) vulnerability score:
  vuln_expert <- vulnerability(
    exposure_results = exp_expert,
    sensitivity_results = sens_ac_expert,
    method_vulnerability = "mean", # default
    method_uncertainty = "mean" # default
  )
  # Calculate risk score:
  risk(
    vulnerability_results = vuln_expert,
    status_results = ex_expert_status
  )
#>        indicator    pressure            type pathway vulnerability    status
#> 1  phytoplankton temperature          direct  expert       -4.7500      good
#> 2  phytoplankton    salinity          direct  expert        0.0000      good
#> 3  phytoplankton      oxygen          direct  expert       -2.2500      good
#> 4  phytoplankton    nutrient          direct  expert       -4.2500      good
#> 5  phytoplankton     fishing          direct  expert        0.0000      good
#> 6  phytoplankton temperature direct_indirect  expert       -5.7500      good
#> 7  phytoplankton    salinity direct_indirect  expert        0.0000      good
#> 8  phytoplankton      oxygen direct_indirect  expert       -3.2500      good
#> 9  phytoplankton    nutrient direct_indirect  expert       -5.2500      good
#> 10 phytoplankton     fishing direct_indirect  expert        0.0000      good
#> 11       herring temperature          direct  expert       -3.5625 undesired
#> 12       herring    salinity          direct  expert       -2.1250 undesired
#> 13       herring      oxygen          direct  expert       -2.4375 undesired
#> 14       herring    nutrient          direct  expert        1.3125 undesired
#> 15       herring     fishing          direct  expert       -4.5000 undesired
#> 16       herring temperature direct_indirect  expert       -7.0000 undesired
#> 17       herring    salinity direct_indirect  expert       -3.6875 undesired
#> 18       herring      oxygen direct_indirect  expert       -4.5000 undesired
#> 19       herring    nutrient direct_indirect  expert        3.9375 undesired
#> 20       herring     fishing direct_indirect  expert       -6.2500 undesired
#> 21           cod temperature          direct  expert       -5.0625 undesired
#> 22           cod    salinity          direct  expert       -3.1250 undesired
#> 23           cod      oxygen          direct  expert       -4.1875 undesired
#> 24           cod    nutrient          direct  expert       -1.0625 undesired
#> 25           cod     fishing          direct  expert       -6.7500 undesired
#> 26           cod temperature direct_indirect  expert       -7.5000 undesired
#> 27           cod    salinity direct_indirect  expert       -6.0000 undesired
#> 28           cod      oxygen direct_indirect  expert       -5.7500 undesired
#> 29           cod    nutrient direct_indirect  expert       -3.0000 undesired
#> 30           cod     fishing direct_indirect  expert       -8.7500 undesired
#> 31      seabirds temperature          direct  expert       -3.7500      good
#> 32      seabirds    salinity          direct  expert        0.0000      good
#> 33      seabirds      oxygen          direct  expert        0.0000      good
#> 34      seabirds    nutrient          direct  expert        0.0000      good
#> 35      seabirds     fishing          direct  expert       -8.0000      good
#> 36      seabirds temperature direct_indirect  expert       -6.7500      good
#> 37      seabirds    salinity direct_indirect  expert       -4.2500      good
#> 38      seabirds      oxygen direct_indirect  expert       -3.2500      good
#> 39      seabirds    nutrient direct_indirect  expert       -4.2500      good
#> 40      seabirds     fishing direct_indirect  expert       -9.0000      good
#>       risk uncertainty
#> 1  -3.7500    1.333333
#> 2   0.0000    2.000000
#> 3  -1.2500    2.166667
#> 4  -3.2500    2.000000
#> 5   0.0000    1.250000
#> 6  -4.7500    1.333333
#> 7   0.0000    2.000000
#> 8  -2.2500    2.166667
#> 9  -4.2500    2.000000
#> 10  0.0000    1.250000
#> 11 -4.5625    1.333333
#> 12 -3.1250    1.333333
#> 13 -3.4375    2.166667
#> 14  0.3125    2.000000
#> 15 -5.5000    1.250000
#> 16 -8.0000    1.333333
#> 17 -4.6875    1.333333
#> 18 -5.5000    2.166667
#> 19  2.9375    2.000000
#> 20 -7.2500    1.250000
#> 21 -6.0625    1.333333
#> 22 -4.1250    1.333333
#> 23 -5.1875    2.166667
#> 24 -2.0625    2.000000
#> 25 -7.7500    1.250000
#> 26 -8.5000    1.333333
#> 27 -7.0000    1.333333
#> 28 -6.7500    2.166667
#> 29 -4.0000    2.000000
#> 30 -9.7500    1.250000
#> 31 -2.7500    1.666667
#> 32  0.0000    2.333333
#> 33  0.0000    2.500000
#> 34  0.0000    2.333333
#> 35 -7.0000    1.916667
#> 36 -5.7500    1.666667
#> 37 -3.2500    2.333333
#> 38 -2.2500    2.500000
#> 39 -3.2500    2.333333
#> 40 -8.0000    1.916667


  ### Demo Model-Based Pathway
  # - using the demo time series 'pressure_ts_baltic' and 'indicator_ts_baltic'

  # Model exposure score:
  exp_model <- model_exposure(
    pressure_time_series = pressure_ts_baltic,
    base_years = c(start = 1984, end = 1994),
    current_years = c(start = 2010, end = 2016)
  )

  # Model sensitivity score:
  sens_ac_model <- model_sensitivity(
    indicator_time_series = indicator_ts_baltic,
    pressure_time_series = pressure_ts_baltic,
    current_years = c(start = 2010, end = 2016)
  )
#> Please review the model diagnostics of the GAMs applied in the time series based sensitivity scoring using the function plot_diagnostic_sensitivity(). Remove models with unacceptable diagnostics from the output table of this function.
  # Add manually adaptive capacity scores (otherwise zero):
  sens_ac_model$adaptive_capacity <- c(rep(1, 8), rep(-1, 8))

  # Calculate (mean) vulnerability score:
  vuln_model <- vulnerability(
    exposure_results = exp_model,
    sensitivity_results = sens_ac_model
  )
  # Calculate status score:
  status_model <- status(
    indicator_time_series = indicator_ts_baltic,
    base_years = c(start = 1984, end = 2010),
    current_years = c(start = 2011, end = 2016)
  )
  # Calculate risk score:
  risk(
    vulnerability_results = vuln_model,
    status_results = status_model
  )
#>                indicator      pressure            type pathway vulnerability
#> 1  zooplankton_mean_size      nitrogen direct_indirect   model          7.50
#> 2  zooplankton_mean_size   phosphorous direct_indirect   model          0.00
#> 3  zooplankton_mean_size surf_temp_sum direct_indirect   model          0.00
#> 4  zooplankton_mean_size  bot_temp_ann direct_indirect   model          0.00
#> 5  zooplankton_mean_size  surf_sal_sum direct_indirect   model          3.25
#> 6  zooplankton_mean_size   bot_sal_ann direct_indirect   model          3.25
#> 7  zooplankton_mean_size   bot_oxy_ann direct_indirect   model          0.00
#> 8  zooplankton_mean_size   fishing_cod direct_indirect   model          0.00
#> 9     eastern_baltic_cod      nitrogen direct_indirect   model         -5.50
#> 10    eastern_baltic_cod   phosphorous direct_indirect   model          0.00
#> 11    eastern_baltic_cod surf_temp_sum direct_indirect   model         -5.00
#> 12    eastern_baltic_cod  bot_temp_ann direct_indirect   model         -7.00
#> 13    eastern_baltic_cod  surf_sal_sum direct_indirect   model         -5.25
#> 14    eastern_baltic_cod   bot_sal_ann direct_indirect   model          0.00
#> 15    eastern_baltic_cod   bot_oxy_ann direct_indirect   model          3.50
#> 16    eastern_baltic_cod   fishing_cod direct_indirect   model         -5.25
#>       status  risk uncertainty
#> 1  undesired  6.50         2.0
#> 2  undesired  0.00         1.5
#> 3  undesired  0.00         1.5
#> 4  undesired  0.00         1.5
#> 5  undesired  2.25         2.0
#> 6  undesired  2.25         2.0
#> 7  undesired  0.00         1.5
#> 8  undesired  0.00         2.0
#> 9  undesired -6.50         2.0
#> 10 undesired  0.00         1.5
#> 11 undesired -6.00         1.5
#> 12 undesired -8.00         1.5
#> 13 undesired -6.25         2.0
#> 14 undesired  0.00         2.0
#> 15 undesired  2.50         1.5
#> 16 undesired -6.25         2.0

# }
```
