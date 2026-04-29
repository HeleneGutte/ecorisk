# Package index

## Expert scoring pathway

Functions for the semiquantitative expert scoring pathway.

- [`create_template_exposure()`](https://helenegutte.github.io/ecorisk/reference/create_template_exposure.md)
  : Create a Template for Expert-Based Exposure Scoring
- [`create_template_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/create_template_sensitivity.md)
  : Create a Template for Expert-Based Sensitivity and Adaptive Capacity
  Scoring
- [`calc_exposure()`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md)
  : Calculate Overall Exposure Scores from Component-Specific Expert
  Ratings
- [`calc_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/calc_sensitivity.md)
  : Calculate Overall Sensitivity and Adaptive Capacity Scores from
  Trait-Specific Expert Ratings

## Modelling pathway

Functions for the quantitative time-series based pathway.

- [`model_exposure()`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md)
  : Model Overall Exposure Scores Using Time Series Data
- [`model_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md)
  : Model Overall Sensitivity Scores Using Time Series Data
- [`plot_diagnostic_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/plot_diagnostic_sensitivity.md)
  : Produce Diagnostic Plots for the

## Joint pathway

Functions shared by both pathways for vulnerability, status, risk, and
aggregation.

- [`vulnerability()`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)
  : Calculate Vulnerability Scores Using Expert-Based or Model-Derived
  Overall Exposure and Sensitivity (Including Adaptive Capacity) Scores
- [`status()`](https://helenegutte.github.io/ecorisk/reference/status.md)
  : Compute Status Scores from Time Series Data
- [`risk()`](https://helenegutte.github.io/ecorisk/reference/risk.md) :
  Calculate Risk Scores Using Expert-Based or Model-Derived
  Vulnerability and Status Scores
- [`aggregate_risk()`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md)
  : Compute High-Complexity Multi-Risk Scores and (Eco)system Risk

## Visualisation

Plotting functions for results communication.

- [`plot_radar()`](https://helenegutte.github.io/ecorisk/reference/plot_radar.md)
  : Generate Radar Charts Displaying Pressure-Specific and Overall Risks
  for Each State Indicator
- [`plot_heatmap()`](https://helenegutte.github.io/ecorisk/reference/plot_heatmap.md)
  : Generate a Heatmap Overview of Individual Risk Scores, Aggregated
  Risk Scores, and Overall Ecosystem Risk

## Example input data

Expert scores and time series used in the tutorial.

- [`ex_expert_exposure`](https://helenegutte.github.io/ecorisk/reference/ex_expert_exposure.md)
  : Expert-based exposure scores for five pressures
- [`ex_expert_sensitivity`](https://helenegutte.github.io/ecorisk/reference/ex_expert_sensitivity.md)
  : Expert-based sensitivity and adaptive capacity scores for four
  indicators and five pressures
- [`ex_expert_status`](https://helenegutte.github.io/ecorisk/reference/ex_expert_status.md)
  : Expert-based status scores for four indicators
- [`indicator_ts_baltic`](https://helenegutte.github.io/ecorisk/reference/indicator_ts_baltic.md)
  : Baltic Sea indicator time series
- [`pressure_ts_baltic`](https://helenegutte.github.io/ecorisk/reference/pressure_ts_baltic.md)
  : Baltic Sea pressure time series
- [`pressure_ts_northsea`](https://helenegutte.github.io/ecorisk/reference/pressure_ts_northsea.md)
  : North Sea pressure time series

## Example output data

Pre-computed function outputs for illustration and testing.

- [`ex_output_calc_exposure`](https://helenegutte.github.io/ecorisk/reference/ex_output_calc_exposure.md)
  :

  Example output from the
  [`calc_exposure()`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md)
  function

- [`ex_output_calc_sensitivity`](https://helenegutte.github.io/ecorisk/reference/ex_output_calc_sensitivity.md)
  :

  Example output from the
  [`calc_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/calc_sensitivity.md)
  function

- [`ex_output_model_exposure`](https://helenegutte.github.io/ecorisk/reference/ex_output_model_exposure.md)
  :

  Example output from the
  [`model_exposure()`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md)
  function based on time series

- [`ex_output_model_sensitivity`](https://helenegutte.github.io/ecorisk/reference/ex_output_model_sensitivity.md)
  :

  Example output from the
  [`model_sensitivity()`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md)
  function based on time series

- [`ex_output_vulnerability_expert`](https://helenegutte.github.io/ecorisk/reference/ex_output_vulnerability_expert.md)
  :

  Example output from the
  [`vulnerability()`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)
  function based on expert scores

- [`ex_output_vulnerability_model`](https://helenegutte.github.io/ecorisk/reference/ex_output_vulnerability_model.md)
  :

  Example output from the
  [`vulnerability()`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)
  function based on modelled scores

- [`ex_output_status`](https://helenegutte.github.io/ecorisk/reference/ex_output_status.md)
  :

  Example output from the
  [`status()`](https://helenegutte.github.io/ecorisk/reference/status.md)
  function

- [`ex_output_risk_expert`](https://helenegutte.github.io/ecorisk/reference/ex_output_risk_expert.md)
  :

  Example output from the
  [`risk()`](https://helenegutte.github.io/ecorisk/reference/risk.md)
  function based on expert scores

- [`ex_output_risk_model`](https://helenegutte.github.io/ecorisk/reference/ex_output_risk_model.md)
  :

  Example output from the
  [`aggregate_risk()`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md)
  function based on modelled scores

- [`ex_output_aggregate_risk_expert`](https://helenegutte.github.io/ecorisk/reference/ex_output_aggregate_risk_expert.md)
  :

  Example output from the
  [`aggregate_risk()`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md)
  function based on expert scores

- [`ex_output_aggregate_risk_model`](https://helenegutte.github.io/ecorisk/reference/ex_output_aggregate_risk_model.md)
  :

  Example output from the
  [`aggregate_risk()`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md)
  function based on modelled scores
