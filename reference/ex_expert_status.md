# Expert-based status scores for four indicators

This demo dataset contains the status scores of four Baltic Sea
indicators based on expert knowledge. The format is the same as the
output table from the
[`status`](https://helenegutte.github.io/ecorisk/reference/status.md)
function that evaluates the status based on time series.

## Usage

``` r
ex_expert_status
```

## Format

A data frame with 4 observations and 3 variables.

- indicator:

  Name of the assessed indicator.

- status:

  Current status of each indicator, either 'good' or 'undesired'.

- score:

  Score for each status (+1 or -1), will be combined with the
  vulnerability to derive risk.
