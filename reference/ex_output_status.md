# Example output from the `status()` function

This dataset provides example output from the
[`status`](https://helenegutte.github.io/ecorisk/reference/status.md)
function, applied to four Baltic Sea indicator time series provided in
[`indicator_ts_baltic`](https://helenegutte.github.io/ecorisk/reference/indicator_ts_baltic.md).

## Usage

``` r
ex_output_status
```

## Format

A data frame with 2 rows and 3 variables.

- indicator:

  Name of the assessed indicator.

- status:

  Qualitative description of the current status compared to the
  threshold.

- score:

  Score used in the risk function to calculate risk from vulnerability
  and status.
