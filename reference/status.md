# Compute Status Scores from Time Series Data

The `status` function assesses whether a state indicator is in a desired
or undesired status during the assessment time period. For this the
function compares the current conditions to the baseline conditions. The
user specifies whether the mean of the current conditions should be
within or outside of a specific deviation from the baseline mean.

## Usage

``` r
status(
  indicator_time_series,
  base_years = NULL,
  base_years_by_ind = NULL,
  current_years = NULL,
  current_years_by_ind = NULL,
  range = "sd",
  sign = "+",
  condition = ">"
)
```

## Arguments

- indicator_time_series:

  A data frame with time series per state indicator. The first column
  MUST be the time column.

- base_years:

  A vector with two numerics, specifying the time period for the
  baseline. The first one `start` is the starting year for all state
  indicators and the second one `end` is the end of the baseline for all
  state indicators. The default is NULL. One can specify indicator
  specific baseline periods using the `base_years_by_ind` argument. If
  `base_years` and `base_years_by_ind` are NULL, then the first 5 years
  of the time series are used as baseline period.

- base_years_by_ind:

  A data frame, specifying the baseline years for each state indicator
  individually, by setting the starting year (second column) and the end
  year (third column). The first column must contain the names of the
  state indicators used in `indicator_time_series`. The default is NULL.
  If `base_years` and `base_years_by_ind` are NULL, then the first 5
  years of the time series are used as baseline period.

- current_years:

  A vector with two numerics, specifying the time period for the
  assessment period. The first one `start` is the starting year for all
  state indicators and the second one `end` is the end of the assessment
  period for all state indicators. The default is NULL. One can specify
  indicator specific assessment periods using the `current_years_by_ind`
  argument. If `current_years` and `current_years_by_ind` are NULL, then
  the last 5 years of the time series are used as baseline period.

- current_years_by_ind:

  A data frame, specifying the assessment period years for each state
  indicator individually, by setting the starting year (second column)
  and the end year (third column). The first column must contain the
  names of the state indicators used in `indicator_time_series`. The
  default is NULL. If `current_years` and `current_years_by_ind` are
  NULL, then the last 5 years of the time series are used as assessment
  period.

- range:

  A vector specifying the allowed deviance from the baseline mean. Can
  be `sd`, `2sd`, `95percentile` or an integer between 1 and 99 to
  evaluate the nth percentile. If the current mean should be compared to
  the baseline mean without any deviance, please set
  `range = mean_only`, leave the default for the sign parameter and
  specify the condition parameter if necessary. Default is `sd`.

- sign:

  A character vector containing `+` or `-`, specifying whether the upper
  or the lower part of the deviance should be analyzed. Default is `+`.

- condition:

  A character vector containing `<` or `>` specifying whether the
  current indicator should be above (\>) or below (\<) the preset
  threshold range to be in a desired status. The default is `>`.

## Value

a data frame containing the indicator name its status and the associated
score, which will be added to the indicators vulnerability to derive the
risk.

## Details

With `range`, `sign` and `condition` one defines good status for the
state indicators. By default the function evaluates whether the current
mean is above +1 standard deviation, if yes the status will be set to
desired. If the state should be within a range of ± standard deviation
and not below that, then the arguments `sign` and `condition` must be
set to '-' and '\>', this specifies that the current mean must be higher
than the mean of the baseline period - 1 standard deviation to be
considered as good status.

## See also

[`model_exposure`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md),
[`model_sensitivity`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md),
[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md),
[`risk`](https://helenegutte.github.io/ecorisk/reference/risk.md)

## Examples

``` r
### Demo with the internal dataset 'indicator_ts_baltic'

# Define a general baseline and current assessment period:
status(
 indicator_time_series = indicator_ts_baltic,
 base_years = c(start = 1984, end = 2010),
 current_years = c(start = 2011, end = 2016)
)
#>               indicator    status score
#> 1 zooplankton_mean_size undesired    -1
#> 2    eastern_baltic_cod undesired    -1

# Define indicator-specific baseline and current assessment periods:
status(
 indicator_time_series = indicator_ts_baltic,
 base_years_by_ind = data.frame(
   ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
   start = c(1984, 1990), end = c(2010, 2010)
 ),
 current_years_by_ind = data.frame(
   ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
   start = c(2011, 2012), end = c(2016, 2016)
 )
)
#>               indicator    status score
#> 1 zooplankton_mean_size undesired    -1
#> 2    eastern_baltic_cod undesired    -1
```
