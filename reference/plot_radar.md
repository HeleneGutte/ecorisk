# Generate Radar Charts Displaying Pressure-Specific and Overall Risks for Each State Indicator

The `plot_radar()`function creates per indicator a ggplot object. The
plot shows the risks of all types and effect directions. The associated
uncertainty can optionally be displayed. In the middle the plot displays
the multi pressure score of a chosen effect type.

## Usage

``` r
plot_radar(
  risk_scores,
  aggregated_scores,
  type = "combined",
  pathway = "combined"
)
```

## Arguments

- risk_scores:

  output from the
  [`risk`](https://helenegutte.github.io/ecorisk/reference/risk.md)
  function.

- aggregated_scores:

  output from the
  [`aggregate_risk`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md)
  function.

- type:

  character string, type used for the multi-pressure score, can be any
  type that has been evaluated. The default is `combined`.

- pathway:

  character string specifying the multi-pressure score, should be
  plotted for each pathway individual `individual` or as a combined
  score `combined`. The default is `combined`.

## Value

a list of ggplot2 objects one for each indicator, the order depends on
the order in the risk_score data set. Each plot shows the risks for one
state indicator for each pressure and type of assessment. In the center
of the plot the multi-pressure score (either in blue or in red) and the
associated aggregated uncertainty (in black) is shown. If one indicator
has been assessed with both pathways, one plot is generated for each of
the pathways. The uncertainty of each individual risk is shown as a ring
around the risks in grey.

## See also

[`risk`](https://helenegutte.github.io/ecorisk/reference/risk.md),
[`aggregate_risk`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md)
to generate result tables/output that serve here as input

## Examples

``` r
### Demo with output data from the risk() and aggregate_risk() functions
#   based on expert scores
# The examples can run for a longer time, thus they are in dontrun{}.
# Using default settings for the indicator-specific overall risk score (coloured value)
# and associated uncertainty score (black value) (i.e., combined across both types)
if (FALSE) { # \dontrun{
p_radar <- plot_radar(
  risk_scores = ex_output_risk_expert,
  aggregated_scores = ex_output_aggregate_risk_expert
)
p_radar[[1]] # display radar chart for first indicator

# Show overall risk score based on direct effects only
p_radar_direct <- plot_radar(
  risk_scores = ex_output_risk_expert,
  aggregated_scores = ex_output_aggregate_risk_expert,
  type = "direct"
)
p_radar_direct[[1]]



### Demo with combined expert-based and model-based pathways

combined_risk <- rbind(ex_output_risk_expert, ex_output_risk_model)
aggr_risk <- aggregate_risk(risk_results = combined_risk)


# Default settings (combined type and pathway)

p_radar_comb <- plot_radar(
  risk_scores = combined_risk,
  aggregated_scores = aggr_risk
)
p_radar_comb[[1]]


# Show overall risk score based on direct/indirect effects only for both
# pathways combined

p_radar_comb_dindi <- plot_radar(
  risk_scores = ex_output_risk_expert,
  aggregated_scores = ex_output_aggregate_risk_expert,
  type = "direct_indirect"
)
p_radar_comb_dindi[[1]]
} # }
```
