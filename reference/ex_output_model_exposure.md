# Example output from the `model_exposure()` function based on time series

This dataset provides example output from the
[`model_exposure`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md)
function, with component-specific exposure scores derived from the
pressure time series in
[`pressure_ts_baltic`](https://helenegutte.github.io/ecorisk/reference/pressure_ts_baltic.md).
These scores are combined into an overall exposure score, with
associated uncertainties derived from two model types.

## Usage

``` r
ex_output_model_exposure
```

## Format

A data frame with 8 observations and 10 variables.

- pressure:

  Names of the assessed pressure.

- exposure:

  Combined exposure score (1 to 5).

- uncertainty:

  Uncertainty score from exposure modelling (1 to 3).

- comp_magnitude:

  Score for the magnitude or degree of change (1 to 5).

- comp_frequency:

  Score for the frequency or duration of change (1 to 5).

- comp_trend:

  Score for the current trend of change (1 to 5).

- comp_direction:

  Direction of the trend slope (increase or decrease).

- comp_spatial:

  Score for spatial extent of the pressure (default: 3; user-defined: 1
  to 5).

- uncertainty_arima:

  Uncertainty score based on an ARIMA model.

- uncertainty_gam:

  Uncertainty score based on a GAM model.

- mean_baseline:

  Mean of the baseline conditions, used for magnitude scoring.

- mean_current:

  Mean of the current conditions, used for magnitude and frequency
  scoring.

- standard_deviation_baseline:

  Standard deviations of the baseline conditions. Used for scoring of
  magnitude and frequency.

- slope_linear_model:

  Slope of the linear model used for scoring the future trend and to
  determine the direction.

- p_value_linear_model:

  P-value of the linear model, used to score the future trend.
