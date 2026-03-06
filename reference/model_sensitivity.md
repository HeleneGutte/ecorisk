# Model Overall Sensitivity Scores Using Time Series Data

The function `model_sensitivity()` uses time series of a state indicator
and a pressure variable to assess the state indicators sensitivity
towards the pressure. The relationship between pressure and state
indicator is determined using a generalized additive model (GAM).
Uncertainty is evaluated with a GAM and an ARIMA model. Use
[`plot_diagnostic_sensitivity`](https://helenegutte.github.io/ecorisk/reference/plot_diagnostic_sensitivity.md)
to review model diagnostics of the applied GAM.

## Usage

``` r
model_sensitivity(
  indicator_time_series,
  pressure_time_series,
  current_years = NULL,
  current_years_by_ind_press = NULL
)
```

## Arguments

- indicator_time_series:

  a data frame containing only the state indicator time series. First
  column MUST be the time column.

- pressure_time_series:

  a data frame containing only the pressure variables. First column MUST
  be the time column.

- current_years:

  A vector with two numerics, specifying the time period for the
  assessment period. The first one `start` is the starting year for all
  pressure-indicator pairs and the second one `end` is the end of the
  assessment period for all pressure-indicator pairs. The default is
  NULL. One can specify pair specific assessment periods using the
  `current_years_by_ind_press` argument. If `current_years` and
  `current_years_by_ind_press` are NULL, then the last 5 years of the
  time series are used as assessment period.

- current_years_by_ind_press:

  a data frame specifying for each indicator-pressure pair the starting
  (third column) and end year (fourth column) where the current
  conditions are best reflected. The default is NULL. If `current_years`
  and `current_years_by_ind_press` are NULL, then the last 5 years of
  the time series are used as assessment period.

## Value

a data frame containing indicator, pressure, type of effect, the
sensitivity score and the associated uncertainty. Positive sensitivity
scores are associated with a positive influence of the pressure on the
indicator and vice versa. Additionally the R-squared, p-values, edf
scores and the mean confidence interval percentage, which are the basis
of the scoring, are provided. The type of effect is automatically set
to`direct_indirect` as the model cannot distinguish between direct and
indirect effects. If default settings are used, the following data frame
will be returned:

- `indicator`:

  Name of the assessed state indicator.

- `pressure`:

  Name of the assessed pressure.

- `type`:

  Type of the assessed effect.

- `pathway`:

  The pathway that has been used to derive the sensitivity scores.

- `sensitivity`:

  Sensitivity score for the assessed state indicator- pressure pair.

- `adaptive_capacity`:

  Adaptive capacity score for the assessed state indicator-pressure
  pair, is automatically set to 0 and can be changed afterwards.

- `uncertainty_sens`:

  uncertainty score associated with the sensitivity scoring.

- `uncertainty_ac`:

  uncertainty score for adaptive capacity scoring. Automatically set to
  NA, can be changed afterwards.

- `r_sq`:

  R-squared value of the GAM, used for the sensitivity scoring.

- `p_value`:

  P-value of the GAM, used to identify significant relationships.
  Unsignificant relationships get a sensitivity score of 0.

- `edf`:

  Estimated degrees of freedom, used to assess non-linearity of the
  relationship between state indicator and pressure.

- `uncertainty_gam`:

  Uncertainty score for sensitivity based on predicted values from a
  GAM.

- `uncertainty_arima`:

  Uncertainty score for sensitivity based on predicted values from an
  ARIMA using the pressure variable as external predictor.

## Details

In case the relationship of one state indicator - pressure pair is not
significant the sensitivity score is 0, and thus also vulnerability and
risk will be 0. For a significant relationship the score will be set
based on the R-squared value from 1 (R-squared \< 0.2) to 5 (R-squared
\>= 0.8). Additionally, the function evaluates the edf score of the GAM
which indicates the degree of non-linearity in the relationship. Since
highly non-linear relationships are harder to predict, the risk of
reaching an undesired state increases and the sensitivity score for
nonlinear relationships will be increased by 1 (if it was not 5
already). The direction of an effect (negative influence or positive
influence of the pressure) is evaluated with the slope of a linear model
representing the assessment period. If the slope of the linear model is
negative, the direction of effect is considered negative as well, and
vice versa for the positive effect. The function assesses uncertainty
associated with the scoring based on a general additive model and an
autoregressive integrated moving average model (ARIMA). The ARIMA model
uses the pressure variable as additional external predictor. The models
are fitted using the time series except the assessment period. The
assessment period is then predicted. The function evaluates how many of
the observed data points are within the predicted 95% confidence
interval. If more than 66 % are within the 95% CI the uncertainty is 1
(low), if less than 33 % are within it, the uncertainty is set to 3
(high). Additionally the function compares the mean size of the
predicted 95% confidence interval and compares it to the maximum range
of the observed data points to account for very large confidence
intervals, which would otherwise lead to too optimistic uncertainty
scores. The lower uncertainty score is selected as final uncertainty
score.

The function also creates columns to give the opportunity to assess
adaptive capacity and its associated uncertainty for each state
indicator-pressure pair. The scores for adaptive capacity and its
associated uncertainty must be specified before the next function
[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)
is applied (see examples). If adaptive capacity and its uncertainty are
not further specified, this will influence the further application of
the ecorisk framework.

## See also

[`model_exposure`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md),
[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)

## Examples

``` r
### Example with the 2 indicators and 8 pressure time series in the Baltic Sea demo data
#   where the last 7 years of the time series represent the current assessment period:
model_sensitivity(
  indicator_time_series = indicator_ts_baltic,
  pressure_time_series = pressure_ts_baltic,
  current_years = c(start = 2010, end = 2016)
)
#> Please review the model diagnostics of the GAMs applied in the time series based sensitivity scoring using the function plot_diagnostic_sensitivity(). Remove models with unacceptable diagnostics from the output table of this function.
#>                indicator      pressure            type pathway sensitivity
#> 1  zooplankton_mean_size      nitrogen direct_indirect   model           3
#> 2  zooplankton_mean_size   phosphorous direct_indirect   model           0
#> 3  zooplankton_mean_size surf_temp_sum direct_indirect   model           0
#> 4  zooplankton_mean_size  bot_temp_ann direct_indirect   model           0
#> 5  zooplankton_mean_size  surf_sal_sum direct_indirect   model          -1
#> 6  zooplankton_mean_size   bot_sal_ann direct_indirect   model          -1
#> 7  zooplankton_mean_size   bot_oxy_ann direct_indirect   model           0
#> 8  zooplankton_mean_size   fishing_cod direct_indirect   model           0
#> 9     eastern_baltic_cod      nitrogen direct_indirect   model          -1
#> 10    eastern_baltic_cod   phosphorous direct_indirect   model           0
#> 11    eastern_baltic_cod surf_temp_sum direct_indirect   model          -1
#> 12    eastern_baltic_cod  bot_temp_ann direct_indirect   model          -3
#> 13    eastern_baltic_cod  surf_sal_sum direct_indirect   model           3
#> 14    eastern_baltic_cod   bot_sal_ann direct_indirect   model           0
#> 15    eastern_baltic_cod   bot_oxy_ann direct_indirect   model           1
#> 16    eastern_baltic_cod   fishing_cod direct_indirect   model          -1
#>    adaptive_capacity uncertainty_sens uncertainty_ac          r_sq      p_value
#> 1                  0                2             NA  0.2431778963 0.0172601903
#> 2                  0                1             NA -0.0002984071 0.3273367464
#> 3                  0                1             NA -0.0178499546 0.5125884857
#> 4                  0                1             NA  0.0294320635 0.1703455137
#> 5                  0                2             NA  0.1005289236 0.0403973732
#> 6                  0                2             NA  0.1673089042 0.0706193795
#> 7                  0                1             NA -0.0286932790 0.7452975202
#> 8                  0                2             NA -0.0301214578 0.8015015038
#> 9                  0                2             NA -0.0218030319 0.7902502903
#> 10                 0                1             NA  0.0401034383 0.1364793690
#> 11                 0                1             NA  0.0987273652 0.0418798313
#> 12                 0                1             NA  0.3762302850 0.0009652114
#> 13                 0                2             NA  0.7424928607 0.0000000000
#> 14                 0                2             NA  0.0584726393 0.0938685078
#> 15                 0                1             NA  0.1832875957 0.0517027687
#> 16                 0                2             NA  0.0100039535 0.5556690537
#>         edf uncertainty_gam uncertainty_arima
#> 1  2.847710               2                 2
#> 2  1.000000               1                 2
#> 3  1.000000               1                 2
#> 4  1.000000               1                 2
#> 5  1.000000               2                 2
#> 6  2.192887               2                 2
#> 7  1.000000               1                 2
#> 8  1.000000               2                 2
#> 9  1.126509               2                 2
#> 10 1.000000               1                 2
#> 11 1.000000               1                 2
#> 12 2.054487               1                 2
#> 13 2.320244               2                 2
#> 14 1.000000               2                 2
#> 15 2.094211               1                 2
#> 16 1.621876               2                 2

### Example with the demo data but indicator-pressure-specific assessment periods:
sens_tbl <- model_sensitivity(
  indicator_time_series = indicator_ts_baltic,
  pressure_time_series = pressure_ts_baltic,
  current_years_by_ind_press = data.frame(
    ind = rep(names(indicator_ts_baltic)[-1], each = 8),
    press = rep(names(pressure_ts_baltic)[-1], 2),
    start = c(rep(2010, 8), rep(2008, 8)),
    end = c(rep(2016, 8), rep(2015, 8))
  )
)
#> Please review the model diagnostics of the GAMs applied in the time series based sensitivity scoring using the function plot_diagnostic_sensitivity(). Remove models with unacceptable diagnostics from the output table of this function.
# add the associated uncertainty (from 1 to 3, default is NA)
sens_tbl$adaptive_capacity <- c(0,0,1,1,1,1,-1,-1, -1,-1,1,1,1,1,1,-1)
sens_tbl$uncertainty_ac <- c(2,2,1,1,1,1,2,1, 3,3,1,1,2,2,3,1)
sens_tbl
#>                indicator      pressure            type pathway sensitivity
#> 1  zooplankton_mean_size      nitrogen direct_indirect   model           3
#> 2  zooplankton_mean_size   phosphorous direct_indirect   model           0
#> 3  zooplankton_mean_size surf_temp_sum direct_indirect   model           0
#> 4  zooplankton_mean_size  bot_temp_ann direct_indirect   model           0
#> 5  zooplankton_mean_size  surf_sal_sum direct_indirect   model          -1
#> 6  zooplankton_mean_size   bot_sal_ann direct_indirect   model          -1
#> 7  zooplankton_mean_size   bot_oxy_ann direct_indirect   model           0
#> 8  zooplankton_mean_size   fishing_cod direct_indirect   model           0
#> 9     eastern_baltic_cod      nitrogen direct_indirect   model           1
#> 10    eastern_baltic_cod   phosphorous direct_indirect   model           0
#> 11    eastern_baltic_cod surf_temp_sum direct_indirect   model          -1
#> 12    eastern_baltic_cod  bot_temp_ann direct_indirect   model          -2
#> 13    eastern_baltic_cod  surf_sal_sum direct_indirect   model           3
#> 14    eastern_baltic_cod   bot_sal_ann direct_indirect   model           0
#> 15    eastern_baltic_cod   bot_oxy_ann direct_indirect   model           1
#> 16    eastern_baltic_cod   fishing_cod direct_indirect   model          -1
#>    adaptive_capacity uncertainty_sens uncertainty_ac          r_sq      p_value
#> 1                  0                2              2  0.2431778963 0.0172601903
#> 2                  0                1              2 -0.0002984071 0.3273367464
#> 3                  1                1              1 -0.0178499546 0.5125884857
#> 4                  1                1              1  0.0294320635 0.1703455137
#> 5                  1                2              1  0.1005289236 0.0403973732
#> 6                  1                2              1  0.1673089042 0.0706193795
#> 7                 -1                1              2 -0.0286932790 0.7452975202
#> 8                 -1                2              1 -0.0301214578 0.8015015038
#> 9                 -1                2              3 -0.0218030319 0.7902502903
#> 10                -1                1              3  0.0401034383 0.1364793690
#> 11                 1                1              1  0.0987273652 0.0418798313
#> 12                 1                1              1  0.3762302850 0.0009652114
#> 13                 1                2              2  0.7424928607 0.0000000000
#> 14                 1                1              2  0.0584726393 0.0938685078
#> 15                 1                1              3  0.1832875957 0.0517027687
#> 16                -1                2              1  0.0100039535 0.5556690537
#>         edf uncertainty_gam uncertainty_arima
#> 1  2.847710               2                 2
#> 2  1.000000               1                 2
#> 3  1.000000               1                 2
#> 4  1.000000               1                 2
#> 5  1.000000               2                 2
#> 6  2.192887               2                 2
#> 7  1.000000               1                 2
#> 8  1.000000               2                 2
#> 9  1.126509               2                 2
#> 10 1.000000               1                 2
#> 11 1.000000               1                 2
#> 12 2.054487               1                 2
#> 13 2.320244               2                 2
#> 14 1.000000               1                 2
#> 15 2.094211               1                 2
#> 16 1.621876               2                 2
```
