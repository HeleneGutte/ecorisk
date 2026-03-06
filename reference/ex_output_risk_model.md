# Example output from the `aggregate_risk()` function based on modelled scores

This dataset provides example output from the
[`risk`](https://helenegutte.github.io/ecorisk/reference/risk.md)
function, applied to the
[`ex_output_vulnerability_model`](https://helenegutte.github.io/ecorisk/reference/ex_output_vulnerability_model.md)
and
[`ex_output_status`](https://helenegutte.github.io/ecorisk/reference/ex_output_status.md)
demo datasets, following the modelling pathway.

## Usage

``` r
ex_output_risk_model
```

## Format

A data frame with 16 observations and 8 variables.

- indicator:

  Name of the assessed indicator.

- pressure:

  Name of the assessed pressure.

- type:

  Type of effect (always direct + indirect for modelling pathway).

- pathway:

  Pathway used for the exposure and sensitivity assessment.

- vulnerability:

  Vulnerability score.

- status:

  Qualitative descriptor of the current status of the indicator.

- risk:

  Risk score.

- uncertainty:

  Uncertainty score associated with the vulnerability component scoring.
