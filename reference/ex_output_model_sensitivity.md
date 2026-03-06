# Example output from the `model_sensitivity()` function based on time series

This dataset provides example output from the
[`model_sensitivity`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md)
function, with sensitivity scores and associated uncertainties for each
indicator-pressure combination. Scores are based on time series data
from
[`pressure_ts_baltic`](https://helenegutte.github.io/ecorisk/reference/pressure_ts_baltic.md)
and
[`indicator_ts_baltic`](https://helenegutte.github.io/ecorisk/reference/indicator_ts_baltic.md).

## Usage

``` r
ex_output_model_sensitivity
```

## Format

A data frame with 16 observations and 12 variables.

- indicator:

  Names of the assessed indicator.

- pressure:

  Names of the assessed pressure.

- type:

  Type of effect (always direct + indirect for modelling pathway).

- pathway:

  Pathway used to assess sensitivity.

- sensitivity:

  Overall sensitivity score (-5 to 5).

- adaptive_capacity:

  Adaptive capacity score (default is 0).

- uncertainty_sens:

  Uncertainty score associated with sensitivity assessment (1 to 3).

- uncertainty_ac:

  Uncertainty score associated with adaptive capacity (1 to 3).

- r_sq:

  R-squared values from the GAM model, used for scoring.

- p_value:

  P-values from the GAM model, determining statistical significance.

- edf:

  Effective degrees of freedom from the GAM model, used to adjust scores
  based on non-linearity risk.

- uncertainty_gam:

  Uncertainty score for sensitivity based on predicted values from a
  GAM.

- uncertainty_arima:

  Uncertainty score for sensitivity based on predicted values from an
  ARIMA using the pressure variable as external predictor.
