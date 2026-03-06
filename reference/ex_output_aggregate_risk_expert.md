# Example output from the `aggregate_risk()` function based on expert scores

This is an expert-based example output from the
[`aggregate_risk`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md)
function applied to the
[`ex_output_risk_expert`](https://helenegutte.github.io/ecorisk/reference/ex_output_risk_expert.md)
demo data.

## Usage

``` r
ex_output_aggregate_risk_expert
```

## Format

A list of three data frames.

- multi_indicator_risk:

  A data frame with 30 rows and 5 columns, containing the
  multi-indicator risk and uncertainty of each pressure per type and
  pathway.

- multi_pressure_risk:

  A data frame with 24 rows and 5 columns, containing the multi-pressure
  risk and uncertainty on each indicator per type and pathway.

- ecosystem_risk:

  A data frame with 6 rows and 4 columns, containing the aggregated
  ecosystem risk and uncertainty per type and pathway.
