# Example output from the `risk()` function based on expert scores

This is an expert-based example output from the
[`risk`](https://helenegutte.github.io/ecorisk/reference/risk.md)
function applied to the
[`ex_output_vulnerability_expert`](https://helenegutte.github.io/ecorisk/reference/ex_output_vulnerability_expert.md)
and
[`ex_expert_status`](https://helenegutte.github.io/ecorisk/reference/ex_expert_status.md)
demo datasets.

## Usage

``` r
ex_output_risk_expert
```

## Format

A data frame with 40 observations and 8 variables.

- indicator:

  Name of the assessed indicator.

- pressure:

  Name of the assessed pressure.

- type:

  Effect type (direct, indirect, or direct + indirect).

- pathway:

  Pathway used for the exposure and sensitivity assessment.

- vulnerability:

  Vulnerability score.

- status:

  Qualitative descriptor of the current status of the indicator.

- risk:

  Risk score.

- uncertainty:

  Uncertainty score, associated with the vulnerability component
  scoring.
