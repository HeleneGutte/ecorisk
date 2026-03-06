# Compute High-Complexity Multi-Risk Scores and (Eco)system Risk

The function `aggregate_risk` uses the output of the risk function or
the vulnerability function. The risk or vulnerability scores are
aggregated in three ways:

1.  as multi-pressure risk per state indicator, representing the overall
    effect on one indicator

2.  as multi- state indicator risk per pressure, representing the
    overall effect each pressure has on all state indicators

3.  as ecosystem risk, the combined multi-pressure risks of all
    indicators.

## Usage

``` r
aggregate_risk(risk_results, method = "mean")
```

## Arguments

- risk_results:

  Risk score for each state indicator ~ pressure ~ type combination,
  derived from the function
  [`risk`](https://helenegutte.github.io/ecorisk/reference/risk.md) or
  [`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md).

- method:

  Character indicating the method for aggregating the pressures and
  state indicators to multiple risk scores and the ecosystem risk score.
  The use can choose between (arithmetic) `mean`, `median`, `sum`,
  `maximum` and `minimum`. Default is the arithmetic mean.

## Value

The function returns a list containing three sublists, first the
multi-state indicator risk list, containing the risks and uncertainties
aggregated per pressure, type and pathway. Second the multi-pressure
risk list, where risks and uncertainties are aggregated per state
indicator and type and pathway. The third list contains the ecosystem
risks, which aggregates the multi-pressure risk and uncertainty scores
per type and pathway.

## Details

The returned lists are required in the plotting functions
[`plot_radar`](https://helenegutte.github.io/ecorisk/reference/plot_radar.md)
and
[`plot_heatmap`](https://helenegutte.github.io/ecorisk/reference/plot_heatmap.md).
The aggregated scores are calculated for each type and pathway
combination individually and across all types and pathways. If only one
type and/or one pathway has been evaluated beforehand, the results will
be the same for the different combinations.

## See also

[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md),
[`risk`](https://helenegutte.github.io/ecorisk/reference/risk.md),
[`plot_radar`](https://helenegutte.github.io/ecorisk/reference/plot_radar.md),
[`plot_heatmap`](https://helenegutte.github.io/ecorisk/reference/plot_heatmap.md)

## Examples

``` r
### Demo with example output from the risk() function based on expert scores
# (where direct and direct/indirect effects were evaluated)

# Calculate mean risks scores per indicator/pressure/ecosystem:
mean_risk <- aggregate_risk(
  risk_results = ex_output_risk_expert,
  method = "mean" # default
)
mean_risk
#> $multi_indicator_risk
#>       pressure            type  pathway      risk uncertainty
#> 5      fishing          direct   expert -5.062500    1.416667
#> 10     fishing direct_indirect   expert -6.250000    1.416667
#> 15     fishing          direct combined -5.062500    1.416667
#> 20     fishing direct_indirect combined -6.250000    1.416667
#> 25     fishing        combined   expert -5.656250    1.416667
#> 30     fishing        combined combined -5.656250    1.416667
#> 4     nutrient          direct   expert -1.250000    2.083333
#> 9     nutrient direct_indirect   expert -2.140625    2.083333
#> 14    nutrient          direct combined -1.250000    2.083333
#> 19    nutrient direct_indirect combined -2.140625    2.083333
#> 24    nutrient        combined   expert -1.695312    2.083333
#> 29    nutrient        combined combined -1.695312    2.083333
#> 3       oxygen          direct   expert -2.468750    2.250000
#> 8       oxygen direct_indirect   expert -4.187500    2.250000
#> 13      oxygen          direct combined -2.468750    2.250000
#> 18      oxygen direct_indirect combined -4.187500    2.250000
#> 23      oxygen        combined   expert -3.328125    2.250000
#> 28      oxygen        combined combined -3.328125    2.250000
#> 2     salinity          direct   expert -1.812500    1.750000
#> 7     salinity direct_indirect   expert -3.734375    1.750000
#> 12    salinity          direct combined -1.812500    1.750000
#> 17    salinity direct_indirect combined -3.734375    1.750000
#> 22    salinity        combined   expert -2.773438    1.750000
#> 27    salinity        combined combined -2.773438    1.750000
#> 1  temperature          direct   expert -4.281250    1.416667
#> 6  temperature direct_indirect   expert -6.750000    1.416667
#> 11 temperature          direct combined -4.281250    1.416667
#> 16 temperature direct_indirect combined -6.750000    1.416667
#> 21 temperature        combined   expert -5.515625    1.416667
#> 26 temperature        combined combined -5.515625    1.416667
#> 
#> $multi_pressure_risk
#>        indicator            type  pathway     risk uncertainty
#> 3            cod          direct   expert -5.03750    1.616667
#> 7            cod direct_indirect   expert -7.20000    1.616667
#> 11           cod          direct combined -5.03750    1.616667
#> 15           cod direct_indirect combined -7.20000    1.616667
#> 19           cod        combined   expert -6.11875    1.616667
#> 23           cod        combined combined -6.11875    1.616667
#> 2        herring          direct   expert -3.26250    1.616667
#> 6        herring direct_indirect   expert -4.50000    1.616667
#> 10       herring          direct combined -3.26250    1.616667
#> 14       herring direct_indirect combined -4.50000    1.616667
#> 18       herring        combined   expert -3.88125    1.616667
#> 22       herring        combined combined -3.88125    1.616667
#> 1  phytoplankton          direct   expert -1.65000    1.750000
#> 5  phytoplankton direct_indirect   expert -2.25000    1.750000
#> 9  phytoplankton          direct combined -1.65000    1.750000
#> 13 phytoplankton direct_indirect combined -2.25000    1.750000
#> 17 phytoplankton        combined   expert -1.95000    1.750000
#> 21 phytoplankton        combined combined -1.95000    1.750000
#> 4       seabirds          direct   expert -1.95000    2.150000
#> 8       seabirds direct_indirect   expert -4.50000    2.150000
#> 12      seabirds          direct combined -1.95000    2.150000
#> 16      seabirds direct_indirect combined -4.50000    2.150000
#> 20      seabirds        combined   expert -3.22500    2.150000
#> 24      seabirds        combined combined -3.22500    2.150000
#> 
#> $ecosystem_risk
#>              type  pathway     risk uncertainty
#> 5        combined   expert -3.79375    1.783333
#> 6        combined combined -3.79375    1.783333
#> 1          direct   expert -2.97500    1.783333
#> 3          direct combined -2.97500    1.783333
#> 2 direct_indirect   expert -4.61250    1.783333
#> 4 direct_indirect combined -4.61250    1.783333
#> 
# Calculate median risks scores:
aggregate_risk(
  risk_results = ex_output_risk_expert,
  method = "median"
)
#> $multi_indicator_risk
#>       pressure            type  pathway     risk uncertainty
#> 5      fishing          direct   expert -6.25000    1.250000
#> 10     fishing direct_indirect   expert -7.62500    1.250000
#> 15     fishing          direct combined -6.25000    1.250000
#> 20     fishing direct_indirect combined -7.62500    1.250000
#> 25     fishing        combined   expert -7.12500    1.250000
#> 30     fishing        combined combined -7.12500    1.250000
#> 4     nutrient          direct   expert -1.03125    2.000000
#> 9     nutrient direct_indirect   expert -3.62500    2.000000
#> 14    nutrient          direct combined -1.03125    2.000000
#> 19    nutrient direct_indirect combined -3.62500    2.000000
#> 24    nutrient        combined   expert -2.65625    2.000000
#> 29    nutrient        combined combined -2.65625    2.000000
#> 3       oxygen          direct   expert -2.34375    2.166667
#> 8       oxygen direct_indirect   expert -3.87500    2.166667
#> 13      oxygen          direct combined -2.34375    2.166667
#> 18      oxygen direct_indirect combined -3.87500    2.166667
#> 23      oxygen        combined   expert -2.84375    2.166667
#> 28      oxygen        combined combined -2.84375    2.166667
#> 2     salinity          direct   expert -1.56250    1.666667
#> 7     salinity direct_indirect   expert -3.96875    1.666667
#> 12    salinity          direct combined -1.56250    1.666667
#> 17    salinity direct_indirect combined -3.96875    1.666667
#> 22    salinity        combined   expert -3.18750    1.666667
#> 27    salinity        combined combined -3.18750    1.666667
#> 1  temperature          direct   expert -4.15625    1.333333
#> 6  temperature direct_indirect   expert -6.87500    1.333333
#> 11 temperature          direct combined -4.15625    1.333333
#> 16 temperature direct_indirect combined -6.87500    1.333333
#> 21 temperature        combined   expert -5.25000    1.333333
#> 26 temperature        combined combined -5.25000    1.333333
#> 
#> $multi_pressure_risk
#>        indicator            type  pathway     risk uncertainty
#> 3            cod          direct   expert -5.18750    1.333333
#> 7            cod direct_indirect   expert -7.00000    1.333333
#> 11           cod          direct combined -5.18750    1.333333
#> 15           cod direct_indirect combined -7.00000    1.333333
#> 19           cod        combined   expert -6.40625    1.333333
#> 23           cod        combined combined -6.40625    1.333333
#> 2        herring          direct   expert -3.43750    1.333333
#> 6        herring direct_indirect   expert -5.50000    1.333333
#> 10       herring          direct combined -3.43750    1.333333
#> 14       herring direct_indirect combined -5.50000    1.333333
#> 18       herring        combined   expert -4.62500    1.333333
#> 22       herring        combined combined -4.62500    1.333333
#> 1  phytoplankton          direct   expert -1.25000    2.000000
#> 5  phytoplankton direct_indirect   expert -2.25000    2.000000
#> 9  phytoplankton          direct combined -1.25000    2.000000
#> 13 phytoplankton direct_indirect combined -2.25000    2.000000
#> 17 phytoplankton        combined   expert -1.75000    2.000000
#> 21 phytoplankton        combined combined -1.75000    2.000000
#> 4       seabirds          direct   expert  0.00000    2.333333
#> 8       seabirds direct_indirect   expert -3.25000    2.333333
#> 12      seabirds          direct combined  0.00000    2.333333
#> 16      seabirds direct_indirect combined -3.25000    2.333333
#> 20      seabirds        combined   expert -3.00000    2.333333
#> 24      seabirds        combined combined -3.00000    2.333333
#> 
#> $ecosystem_risk
#>              type  pathway     risk uncertainty
#> 5        combined   expert -3.81250    1.666667
#> 6        combined combined -3.81250    1.666667
#> 1          direct   expert -2.34375    1.666667
#> 3          direct combined -2.34375    1.666667
#> 2 direct_indirect   expert -4.37500    1.666667
#> 4 direct_indirect combined -4.37500    1.666667
#> 
# Calculate maximum risks scores:
aggregate_risk(
  risk_results = ex_output_risk_expert,
  method = "maximum"
)
#> $multi_indicator_risk
#>       pressure            type  pathway    risk uncertainty
#> 5      fishing          direct   expert  0.0000    1.916667
#> 10     fishing direct_indirect   expert  0.0000    1.916667
#> 15     fishing          direct combined  0.0000    1.916667
#> 20     fishing direct_indirect combined  0.0000    1.916667
#> 25     fishing        combined   expert  0.0000    1.916667
#> 30     fishing        combined combined  0.0000    1.916667
#> 4     nutrient          direct   expert  0.3125    2.333333
#> 9     nutrient direct_indirect   expert  2.9375    2.333333
#> 14    nutrient          direct combined  0.3125    2.333333
#> 19    nutrient direct_indirect combined  2.9375    2.333333
#> 24    nutrient        combined   expert  2.9375    2.333333
#> 29    nutrient        combined combined  2.9375    2.333333
#> 3       oxygen          direct   expert  0.0000    2.500000
#> 8       oxygen direct_indirect   expert -2.2500    2.500000
#> 13      oxygen          direct combined  0.0000    2.500000
#> 18      oxygen direct_indirect combined -2.2500    2.500000
#> 23      oxygen        combined   expert  0.0000    2.500000
#> 28      oxygen        combined combined  0.0000    2.500000
#> 2     salinity          direct   expert  0.0000    2.333333
#> 7     salinity direct_indirect   expert  0.0000    2.333333
#> 12    salinity          direct combined  0.0000    2.333333
#> 17    salinity direct_indirect combined  0.0000    2.333333
#> 22    salinity        combined   expert  0.0000    2.333333
#> 27    salinity        combined combined  0.0000    2.333333
#> 1  temperature          direct   expert -2.7500    1.666667
#> 6  temperature direct_indirect   expert -4.7500    1.666667
#> 11 temperature          direct combined -2.7500    1.666667
#> 16 temperature direct_indirect combined -4.7500    1.666667
#> 21 temperature        combined   expert -2.7500    1.666667
#> 26 temperature        combined combined -2.7500    1.666667
#> 
#> $multi_pressure_risk
#>        indicator            type  pathway    risk uncertainty
#> 3            cod          direct   expert -2.0625    2.166667
#> 7            cod direct_indirect   expert -4.0000    2.166667
#> 11           cod          direct combined -2.0625    2.166667
#> 15           cod direct_indirect combined -4.0000    2.166667
#> 19           cod        combined   expert -2.0625    2.166667
#> 23           cod        combined combined -2.0625    2.166667
#> 2        herring          direct   expert  0.3125    2.166667
#> 6        herring direct_indirect   expert  2.9375    2.166667
#> 10       herring          direct combined  0.3125    2.166667
#> 14       herring direct_indirect combined  2.9375    2.166667
#> 18       herring        combined   expert  2.9375    2.166667
#> 22       herring        combined combined  2.9375    2.166667
#> 1  phytoplankton          direct   expert  0.0000    2.166667
#> 5  phytoplankton direct_indirect   expert  0.0000    2.166667
#> 9  phytoplankton          direct combined  0.0000    2.166667
#> 13 phytoplankton direct_indirect combined  0.0000    2.166667
#> 17 phytoplankton        combined   expert  0.0000    2.166667
#> 21 phytoplankton        combined combined  0.0000    2.166667
#> 4       seabirds          direct   expert  0.0000    2.500000
#> 8       seabirds direct_indirect   expert -2.2500    2.500000
#> 12      seabirds          direct combined  0.0000    2.500000
#> 16      seabirds direct_indirect combined -2.2500    2.500000
#> 20      seabirds        combined   expert  0.0000    2.500000
#> 24      seabirds        combined combined  0.0000    2.500000
#> 
#> $ecosystem_risk
#>              type  pathway   risk uncertainty
#> 5        combined   expert 2.9375         2.5
#> 6        combined combined 2.9375         2.5
#> 1          direct   expert 0.3125         2.5
#> 3          direct combined 0.3125         2.5
#> 2 direct_indirect   expert 2.9375         2.5
#> 4 direct_indirect combined 2.9375         2.5
#> 


### Demo with example output from the risk() function based on modelled
#   scores (where only direct/indirect effects were evaluated)

# Calculate mean risks scores:
aggregate_risk(risk_results = ex_output_risk_model)
#> $multi_indicator_risk
#>         pressure            type  pathway   risk uncertainty
#> 7    bot_oxy_ann direct_indirect    model -2.375         1.5
#> 15   bot_oxy_ann direct_indirect combined -2.375         1.5
#> 23   bot_oxy_ann        combined    model -2.375         1.5
#> 31   bot_oxy_ann        combined combined -2.375         1.5
#> 6    bot_sal_ann direct_indirect    model  1.625         2.0
#> 14   bot_sal_ann direct_indirect combined  1.625         2.0
#> 22   bot_sal_ann        combined    model  1.625         2.0
#> 30   bot_sal_ann        combined combined  1.625         2.0
#> 4   bot_temp_ann direct_indirect    model -3.375         1.5
#> 12  bot_temp_ann direct_indirect combined -3.375         1.5
#> 20  bot_temp_ann        combined    model -3.375         1.5
#> 28  bot_temp_ann        combined combined -3.375         1.5
#> 8    fishing_cod direct_indirect    model -2.625         2.0
#> 16   fishing_cod direct_indirect combined -2.625         2.0
#> 24   fishing_cod        combined    model -2.625         2.0
#> 32   fishing_cod        combined combined -2.625         2.0
#> 1       nitrogen direct_indirect    model  0.000         2.0
#> 9       nitrogen direct_indirect combined  0.000         2.0
#> 17      nitrogen        combined    model  0.000         2.0
#> 25      nitrogen        combined combined  0.000         2.0
#> 2    phosphorous direct_indirect    model  0.000         1.5
#> 10   phosphorous direct_indirect combined  0.000         1.5
#> 18   phosphorous        combined    model  0.000         1.5
#> 26   phosphorous        combined combined  0.000         1.5
#> 5   surf_sal_sum direct_indirect    model -2.000         2.0
#> 13  surf_sal_sum direct_indirect combined -2.000         2.0
#> 21  surf_sal_sum        combined    model -2.000         2.0
#> 29  surf_sal_sum        combined combined -2.000         2.0
#> 3  surf_temp_sum direct_indirect    model -2.125         1.5
#> 11 surf_temp_sum direct_indirect combined -2.125         1.5
#> 19 surf_temp_sum        combined    model -2.125         1.5
#> 27 surf_temp_sum        combined combined -2.125         1.5
#> 
#> $multi_pressure_risk
#>               indicator            type  pathway     risk uncertainty
#> 2    eastern_baltic_cod direct_indirect    model -4.00000        1.75
#> 4    eastern_baltic_cod direct_indirect combined -4.00000        1.75
#> 6    eastern_baltic_cod        combined    model -4.00000        1.75
#> 8    eastern_baltic_cod        combined combined -4.00000        1.75
#> 1 zooplankton_mean_size direct_indirect    model  1.28125        1.75
#> 3 zooplankton_mean_size direct_indirect combined  1.28125        1.75
#> 5 zooplankton_mean_size        combined    model  1.28125        1.75
#> 7 zooplankton_mean_size        combined combined  1.28125        1.75
#> 
#> $ecosystem_risk
#>              type  pathway      risk uncertainty
#> 3        combined    model -1.359375        1.75
#> 4        combined combined -1.359375        1.75
#> 1 direct_indirect    model -1.359375        1.75
#> 2 direct_indirect combined -1.359375        1.75
#> 


### Demo with combined expert-based and model-based pathways

combined_risk <- rbind(ex_output_risk_expert, ex_output_risk_model)
aggr_risk <- aggregate_risk(risk_results = combined_risk)
aggr_risk
#> $multi_indicator_risk
#>         pressure            type  pathway      risk uncertainty
#> 17   bot_oxy_ann direct_indirect    model -2.375000    1.500000
#> 35   bot_oxy_ann direct_indirect combined -2.375000    1.500000
#> 48   bot_oxy_ann        combined    model -2.375000    1.500000
#> 61   bot_oxy_ann        combined combined -2.375000    1.500000
#> 16   bot_sal_ann direct_indirect    model  1.625000    2.000000
#> 34   bot_sal_ann direct_indirect combined  1.625000    2.000000
#> 47   bot_sal_ann        combined    model  1.625000    2.000000
#> 60   bot_sal_ann        combined combined  1.625000    2.000000
#> 14  bot_temp_ann direct_indirect    model -3.375000    1.500000
#> 32  bot_temp_ann direct_indirect combined -3.375000    1.500000
#> 45  bot_temp_ann        combined    model -3.375000    1.500000
#> 58  bot_temp_ann        combined combined -3.375000    1.500000
#> 5        fishing          direct   expert -5.062500    1.416667
#> 10       fishing direct_indirect   expert -6.250000    1.416667
#> 23       fishing          direct combined -5.062500    1.416667
#> 28       fishing direct_indirect combined -6.250000    1.416667
#> 41       fishing        combined   expert -5.656250    1.416667
#> 54       fishing        combined combined -5.656250    1.416667
#> 18   fishing_cod direct_indirect    model -2.625000    2.000000
#> 36   fishing_cod direct_indirect combined -2.625000    2.000000
#> 49   fishing_cod        combined    model -2.625000    2.000000
#> 62   fishing_cod        combined combined -2.625000    2.000000
#> 11      nitrogen direct_indirect    model  0.000000    2.000000
#> 29      nitrogen direct_indirect combined  0.000000    2.000000
#> 42      nitrogen        combined    model  0.000000    2.000000
#> 55      nitrogen        combined combined  0.000000    2.000000
#> 4       nutrient          direct   expert -1.250000    2.083333
#> 9       nutrient direct_indirect   expert -2.140625    2.083333
#> 22      nutrient          direct combined -1.250000    2.083333
#> 27      nutrient direct_indirect combined -2.140625    2.083333
#> 40      nutrient        combined   expert -1.695312    2.083333
#> 53      nutrient        combined combined -1.695312    2.083333
#> 3         oxygen          direct   expert -2.468750    2.250000
#> 8         oxygen direct_indirect   expert -4.187500    2.250000
#> 21        oxygen          direct combined -2.468750    2.250000
#> 26        oxygen direct_indirect combined -4.187500    2.250000
#> 39        oxygen        combined   expert -3.328125    2.250000
#> 52        oxygen        combined combined -3.328125    2.250000
#> 12   phosphorous direct_indirect    model  0.000000    1.500000
#> 30   phosphorous direct_indirect combined  0.000000    1.500000
#> 43   phosphorous        combined    model  0.000000    1.500000
#> 56   phosphorous        combined combined  0.000000    1.500000
#> 2       salinity          direct   expert -1.812500    1.750000
#> 7       salinity direct_indirect   expert -3.734375    1.750000
#> 20      salinity          direct combined -1.812500    1.750000
#> 25      salinity direct_indirect combined -3.734375    1.750000
#> 38      salinity        combined   expert -2.773438    1.750000
#> 51      salinity        combined combined -2.773438    1.750000
#> 15  surf_sal_sum direct_indirect    model -2.000000    2.000000
#> 33  surf_sal_sum direct_indirect combined -2.000000    2.000000
#> 46  surf_sal_sum        combined    model -2.000000    2.000000
#> 59  surf_sal_sum        combined combined -2.000000    2.000000
#> 13 surf_temp_sum direct_indirect    model -2.125000    1.500000
#> 31 surf_temp_sum direct_indirect combined -2.125000    1.500000
#> 44 surf_temp_sum        combined    model -2.125000    1.500000
#> 57 surf_temp_sum        combined combined -2.125000    1.500000
#> 1    temperature          direct   expert -4.281250    1.416667
#> 6    temperature direct_indirect   expert -6.750000    1.416667
#> 19   temperature          direct combined -4.281250    1.416667
#> 24   temperature direct_indirect combined -6.750000    1.416667
#> 37   temperature        combined   expert -5.515625    1.416667
#> 50   temperature        combined combined -5.515625    1.416667
#> 
#> $multi_pressure_risk
#>                indicator            type  pathway     risk uncertainty
#> 3                    cod          direct   expert -5.03750    1.616667
#> 7                    cod direct_indirect   expert -7.20000    1.616667
#> 13                   cod          direct combined -5.03750    1.616667
#> 17                   cod direct_indirect combined -7.20000    1.616667
#> 23                   cod        combined   expert -6.11875    1.616667
#> 29                   cod        combined combined -6.11875    1.616667
#> 10    eastern_baltic_cod direct_indirect    model -4.00000    1.750000
#> 20    eastern_baltic_cod direct_indirect combined -4.00000    1.750000
#> 26    eastern_baltic_cod        combined    model -4.00000    1.750000
#> 32    eastern_baltic_cod        combined combined -4.00000    1.750000
#> 2                herring          direct   expert -3.26250    1.616667
#> 6                herring direct_indirect   expert -4.50000    1.616667
#> 12               herring          direct combined -3.26250    1.616667
#> 16               herring direct_indirect combined -4.50000    1.616667
#> 22               herring        combined   expert -3.88125    1.616667
#> 28               herring        combined combined -3.88125    1.616667
#> 1          phytoplankton          direct   expert -1.65000    1.750000
#> 5          phytoplankton direct_indirect   expert -2.25000    1.750000
#> 11         phytoplankton          direct combined -1.65000    1.750000
#> 15         phytoplankton direct_indirect combined -2.25000    1.750000
#> 21         phytoplankton        combined   expert -1.95000    1.750000
#> 27         phytoplankton        combined combined -1.95000    1.750000
#> 4               seabirds          direct   expert -1.95000    2.150000
#> 8               seabirds direct_indirect   expert -4.50000    2.150000
#> 14              seabirds          direct combined -1.95000    2.150000
#> 18              seabirds direct_indirect combined -4.50000    2.150000
#> 24              seabirds        combined   expert -3.22500    2.150000
#> 30              seabirds        combined combined -3.22500    2.150000
#> 9  zooplankton_mean_size direct_indirect    model  1.28125    1.750000
#> 19 zooplankton_mean_size direct_indirect combined  1.28125    1.750000
#> 25 zooplankton_mean_size        combined    model  1.28125    1.750000
#> 31 zooplankton_mean_size        combined combined  1.28125    1.750000
#> 
#> $ecosystem_risk
#>              type  pathway      risk uncertainty
#> 6        combined   expert -3.793750    1.783333
#> 7        combined    model -1.359375    1.750000
#> 8        combined combined -2.982292    1.772222
#> 1          direct   expert -2.975000    1.783333
#> 4          direct combined -2.975000    1.783333
#> 2 direct_indirect   expert -4.612500    1.783333
#> 3 direct_indirect    model -1.359375    1.750000
#> 5 direct_indirect combined -3.528125    1.772222
#> 

aggr_risk$multi_indicator_risk |>
  dplyr::filter(type == "combined", pathway == "combined")
#>         pressure     type  pathway      risk uncertainty
#> 1    bot_oxy_ann combined combined -2.375000    1.500000
#> 2    bot_sal_ann combined combined  1.625000    2.000000
#> 3   bot_temp_ann combined combined -3.375000    1.500000
#> 4        fishing combined combined -5.656250    1.416667
#> 5    fishing_cod combined combined -2.625000    2.000000
#> 6       nitrogen combined combined  0.000000    2.000000
#> 7       nutrient combined combined -1.695312    2.083333
#> 8         oxygen combined combined -3.328125    2.250000
#> 9    phosphorous combined combined  0.000000    1.500000
#> 10      salinity combined combined -2.773438    1.750000
#> 11  surf_sal_sum combined combined -2.000000    2.000000
#> 12 surf_temp_sum combined combined -2.125000    1.500000
#> 13   temperature combined combined -5.515625    1.416667
aggr_risk$multi_pressure_risk |>
  dplyr::filter(type == "combined", pathway == "combined")
#>               indicator     type  pathway     risk uncertainty
#> 1                   cod combined combined -6.11875    1.616667
#> 2    eastern_baltic_cod combined combined -4.00000    1.750000
#> 3               herring combined combined -3.88125    1.616667
#> 4         phytoplankton combined combined -1.95000    1.750000
#> 5              seabirds combined combined -3.22500    2.150000
#> 6 zooplankton_mean_size combined combined  1.28125    1.750000
aggr_risk$ecosystem_risk |>
  dplyr::filter(type == "combined", pathway == "combined")
#>       type  pathway      risk uncertainty
#> 1 combined combined -2.982292    1.772222


### Demo with vulnerability scores using example output data from
#   vulnerability() based on modelled scores

aggregate_risk(risk_results = ex_output_vulnerability_model)
#> $multi_indicator_risk
#>         pressure            type  pathway   risk uncertainty
#> 7    bot_oxy_ann direct_indirect    model -1.875         1.5
#> 15   bot_oxy_ann direct_indirect combined -1.875         1.5
#> 23   bot_oxy_ann        combined    model -1.875         1.5
#> 31   bot_oxy_ann        combined combined -1.875         1.5
#> 6    bot_sal_ann direct_indirect    model  2.125         2.0
#> 14   bot_sal_ann direct_indirect combined  2.125         2.0
#> 22   bot_sal_ann        combined    model  2.125         2.0
#> 30   bot_sal_ann        combined combined  2.125         2.0
#> 4   bot_temp_ann direct_indirect    model -2.875         1.5
#> 12  bot_temp_ann direct_indirect combined -2.875         1.5
#> 20  bot_temp_ann        combined    model -2.875         1.5
#> 28  bot_temp_ann        combined combined -2.875         1.5
#> 8    fishing_cod direct_indirect    model -2.125         2.0
#> 16   fishing_cod direct_indirect combined -2.125         2.0
#> 24   fishing_cod        combined    model -2.125         2.0
#> 32   fishing_cod        combined combined -2.125         2.0
#> 1       nitrogen direct_indirect    model  1.000         2.0
#> 9       nitrogen direct_indirect combined  1.000         2.0
#> 17      nitrogen        combined    model  1.000         2.0
#> 25      nitrogen        combined combined  1.000         2.0
#> 2    phosphorous direct_indirect    model  0.000         1.5
#> 10   phosphorous direct_indirect combined  0.000         1.5
#> 18   phosphorous        combined    model  0.000         1.5
#> 26   phosphorous        combined combined  0.000         1.5
#> 5   surf_sal_sum direct_indirect    model -1.000         2.0
#> 13  surf_sal_sum direct_indirect combined -1.000         2.0
#> 21  surf_sal_sum        combined    model -1.000         2.0
#> 29  surf_sal_sum        combined combined -1.000         2.0
#> 3  surf_temp_sum direct_indirect    model -1.625         1.5
#> 11 surf_temp_sum direct_indirect combined -1.625         1.5
#> 19 surf_temp_sum        combined    model -1.625         1.5
#> 27 surf_temp_sum        combined combined -1.625         1.5
#> 
#> $multi_pressure_risk
#>               indicator            type  pathway     risk uncertainty
#> 2    eastern_baltic_cod direct_indirect    model -3.25000        1.75
#> 4    eastern_baltic_cod direct_indirect combined -3.25000        1.75
#> 6    eastern_baltic_cod        combined    model -3.25000        1.75
#> 8    eastern_baltic_cod        combined combined -3.25000        1.75
#> 1 zooplankton_mean_size direct_indirect    model  1.65625        1.75
#> 3 zooplankton_mean_size direct_indirect combined  1.65625        1.75
#> 5 zooplankton_mean_size        combined    model  1.65625        1.75
#> 7 zooplankton_mean_size        combined combined  1.65625        1.75
#> 
#> $ecosystem_risk
#>              type  pathway      risk uncertainty
#> 3        combined    model -0.796875        1.75
#> 4        combined combined -0.796875        1.75
#> 1 direct_indirect    model -0.796875        1.75
#> 2 direct_indirect combined -0.796875        1.75
#> 
```
