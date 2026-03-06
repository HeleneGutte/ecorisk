# Example output from the `aggregate_risk()` function based on modelled scores

This dataset provides example output from the
[`aggregate_risk`](https://helenegutte.github.io/ecorisk/reference/aggregate_risk.md)
function, applied to the
[`ex_output_risk_model`](https://helenegutte.github.io/ecorisk/reference/ex_output_risk_model.md)
demo data, following the modelling pathway.

## Usage

``` r
ex_output_aggregate_risk_model
```

## Format

A list of three data frames.

- multi_indicator_risk:

  A data frame with 32 rows and 5 variables, containing the
  multi-indicator risk and uncertainty of each pressure per type and
  pathway.

- multi_pressure_risk:

  A data frame with 8 rows and 5 variables, containing the
  multi-pressure risk and uncertainty on each indicator per type and
  pathway.

- ecosystem_risk:

  A data frame with 4 rows and 4 variables, containing the aggregated
  ecosystem risk and uncertainty per type and pathway.
