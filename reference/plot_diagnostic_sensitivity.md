# Produce Diagnostic Plots for the

The function `plot_diagnostic_sensitivity()` creates diagnostic plots
for the generalized additive models, which are used as basis for the
sensitivity scoring in
[`model_sensitivity`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md).

## Usage

``` r
plot_diagnostic_sensitivity(indicator_time_series, pressure_time_series)
```

## Arguments

- indicator_time_series:

  a data frame containing only the state indicator time series. First
  column MUST be the time column.

- pressure_time_series:

  a data frame containing only the pressure variables. First column MUST
  be the time column.

## Value

The function returns a list of ggplot objects combined with patchwork.
For each state and pressure indicator combination 4 diagnostic plots are
created: Q-Q plots, residuals vs- linear predictor, a histogram of the
residuals and response vs. fitted values.

## Examples

``` r
if (FALSE) { # \dontrun{
plot_diagnostic_sensitivity(
  indicator_time_series = indicator_ts_baltic[, c(1,2)],
  pressure_time_series = pressure_ts_baltic[, c(1,2)]
 )
 } # }
```
