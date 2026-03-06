# Model Overall Exposure Scores Using Time Series Data

This function statistically evaluates the exposure to a pressure, based
on time series data. The scoring is based on the paper of Gaichas et
al., 2014: A risk-based approach to evaluating northeast US fish
community vulnerability to climate change. The exposure scoring is split
into four components: magnitude (or degree of change), frequency of
change, the future trend of the pressure, and spatial scale. Uncertainty
of the exposure assessment is evaluated using general additive models
(GAM) and an autoregressive integrated moving average model (ARIMA).

## Usage

``` r
model_exposure(
  pressure_time_series,
  base_years = NULL,
  base_years_by_press = NULL,
  current_years = NULL,
  current_years_by_press = NULL,
  trend = "return",
  spatial = 3
)
```

## Arguments

- pressure_time_series:

  A data frame (not a tibble object) with time series of pressures to be
  evaluated. First column MUST be the time column.

- base_years:

  A vector with two numerics, specifying the time period for the
  baseline. The first one `start` is the starting year for all pressures
  and the second one `end` is the end of the baseline for all pressures.
  The default is NULL. One can specify pressure specific baseline
  periods using the `base_years_by_press` argument. If `base_years` and
  `base_years_by_ind` are NULL, then the first 5 years of the time
  series are used as baseline period.

- base_years_by_press:

  A data frame, specifying the baseline years for each pressure
  individually, by setting the starting year (second column) and the end
  year (third column). The first column must contain the names of the
  pressure indicators used in `pressure_time_series`. The default is
  NULL. If `base_years` and `base_years_by_press` are NULL, then the
  first 5 years of the time series are used as baseline period.

- current_years:

  A vector with two numerics, specifying the time period for the
  assessment period. The first one `start` is the starting year for all
  pressures and the second one `end` is the end of the assessment period
  for all pressures. The default is NULL. One can specify pressure
  specific assessment periods using the `current_years_by_press`
  argument. If `current_years` and `current_years_by_press` are NULL,
  then the last 5 years of the time series are used as assessment
  period.

- current_years_by_press:

  A data frame, specifying the assessment period for each pressure
  individually, by setting the starting year (second column) and the end
  year (third column). The first column must contain the names of the
  pressure indicators used in `pressure_time_series`. The default is
  NULL. If `current_years` and `current_years_by_press` are NULL, then
  the last 5 years of the time series are used as assessment period.

- trend:

  a character vector specifying whether a trend returning to the
  baseline conditions should be considered as good or a trend further
  leaving the baseline conditions. Possible inputs are `return` or
  `leave`.Default is not specified is `return`, meaning a return to the
  baseline is desired.

- spatial:

  a vector with scores for the spatial scale of each pressure. The
  default is 3 for each pressure, meaning that 40 - 60% of the entire
  assessment area is affected. Scores should be on a scale from 1 - 5,
  depending on the percent of area that is affected by the pressure:

  - 1: \< 20%,

  - 2: 20 - 40%,

  - 3: 40 - 60%,

  - 4: 60 - 80%,

  - 5: \> 80%.

## Value

a data frame containing the pressure names, the aggregated exposure
score and scores for magnitude, frequency, future trend and spatial
scale of the pressures, the final uncertainty score and uncertainty
scores of the ARIMA and the GAM. If default settings are used, the
following data frame will be returned:

- `pressure`:

  Name of the assessed pressure.

- `exposure`:

  Exposure score, mean of the four assessed exposure components.

- `uncertainty`:

  Uncertainty score associated with the exposure assessment.

- `comp_magnitude`:

  Score for the magnitude of change.

- `comp_frequency`:

  Score for the frequency of a significant deviation from baseline
  conditions.

- `comp_trend`:

  Score for the future trend of the pressure.

- `comp_direction`:

  Direction of the development of the pressure in the assessment period.

- `comp_spatial`:

  Score for the spatial scale, either set by the user or automatically
  set to 3.

- `uncertainty_arima`:

  Uncertainty score based on the ARIMA model.

- `uncertainty_gam`:

  Uncertainty based on the GAM.

- `mean_baseline`:

  Mean of the baseline conditions, used for magnitude scoring.

- `mean_current`:

  Mean of the current conditions, used for magnitude and frequency
  scoring.

- `standard_deviation_baseline`:

  Standard deviations of the baseline conditions. Used for scoring of
  magnitude and frequency.

- `slope_linear_model`:

  Slope of the linear model used for scoring the future trend and to
  determine the direction.

- `p_value_linear_model`:

  P-value of the linear model, used to score the future trend.

## Details

All components are scored on a scale from 1 - 5, low impact to high
impact. The degree of change compares the mean of the current time
period to the baseline time period, the score is based on standard
deviations. The frequency evaluates in how much percent of the current
time period the mean deviates more than one standard deviation from the
baseline mean. The future trend scores if the pressure will in the
future be in desired conditions or not. Usually this means the pressure
returns to the baseline conditions. The overall exposure score is the
mean of all four components. Uncertainty of exposure is evaluated using
a general additive model and an autoregressive integrated moving average
model (ARIMA). The models are fitted using the time series except the
assessment period. The assessment period is then predicted. The function
evaluates how many of the observed data points are within the predicted
95% confidence interval. If more than 66 % are within the 95% CI the
uncertainty is 1 (low), if less than 33 % are within it, the uncertainty
is set to 3 (high). Additionally the function compares the mean size of
the predicted 95% confidence interval and compares it to the maximum
range of the observed data points to account for very large confidence
intervals, which would otherwise lead to too optimistic uncertainty
scores. The lower uncertainty score is selected as final uncertainty
score. The time periods of baseline and assessment period have to be
carefully set to reflect ongoing dynamics. Especially for oscillating
pressures time periods should be longer to assess the overall trend and
not the oscillation itself.

## See also

[`model_sensitivity`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md),
[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)

## Examples

``` r
### Example with 3 pressure time series in the demo data 'pressure_ts_baltic'
#   where the first 11 years represent the general baseline period and the last
#   7 years of the time series the current assessment period:
sub_ts <- pressure_ts_baltic[ ,c("year", "surf_temp_sum", 
  "surf_sal_sum", "bot_oxy_ann")]
model_exposure(
  pressure_time_series = sub_ts ,
  base_years = c(start = 1984, end = 1994),
  current_years = c(start = 2010, end = 2016)
)
#>        pressure exposure uncertainty comp_magnitude comp_frequency comp_trend
#> 1 surf_temp_sum     3.00           2              2              4          3
#> 2  surf_sal_sum     3.25           2              2              5          3
#> 3   bot_oxy_ann     3.50           2              2              5          4
#>   comp_direction comp_spatial uncertainty_arima uncertainty_gam mean_baseline
#> 1       increase            3                 2               2     13.360030
#> 2       decrease            3                 2               3      5.935398
#> 3       decrease            3                 2               2      4.466583
#>   mean_current standard_deviation_baseline slope_linear_model
#> 1    14.895450                    1.344546         0.06988084
#> 2     5.574401                    0.193843        -0.02150520
#> 3     4.001809                    0.253705        -0.36222659
#>   p_value_linear_model
#> 1           0.74663270
#> 2           0.92121801
#> 3           0.03754337

### Example with 2 pressure time series and pressure-specific periods
sub_ts <- pressure_ts_baltic[ ,c("year", "nitrogen", "phosphorous")]
model_exposure(
  pressure_time_series = sub_ts,
  base_years_by_press = data.frame(
    press = c("nitrogen", "phosphorous"),
    start = c(1984, 1990), end = c(1994, 2000)),
  current_years_by_press = data.frame(
      press = c("nitrogen", "phosphorous"),
      start = c(2010, 2012), end = c(2016, 2016))
)
#>      pressure exposure uncertainty comp_magnitude comp_frequency comp_trend
#> 1    nitrogen     3.50           2              3              5          3
#> 2 phosphorous     2.25           2              1              2          3
#>   comp_direction comp_spatial uncertainty_arima uncertainty_gam mean_baseline
#> 1       increase            3                 2               2    19.9366961
#> 2       decrease            3                 2               2     0.5810386
#>   mean_current standard_deviation_baseline slope_linear_model
#> 1   22.6823144                  1.27570100        0.002903849
#> 2    0.6408257                  0.07574186       -0.017404547
#>   p_value_linear_model
#> 1            0.9893508
#> 2            0.9649661
```
