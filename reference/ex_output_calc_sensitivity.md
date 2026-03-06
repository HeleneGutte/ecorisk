# Example output from the `calc_sensitivity()` function

This is an expert-based example output from the
[`calc_sensitivity`](https://helenegutte.github.io/ecorisk/reference/calc_sensitivity.md)
function applied to the
[`ex_expert_sensitivity`](https://helenegutte.github.io/ecorisk/reference/ex_expert_sensitivity.md)
demo data.

## Usage

``` r
ex_output_calc_sensitivity
```

## Format

A data frame with 40 observations and 18 variables.

- indicator:

  Name of the assessed indicator.

- pressure:

  Name of the assessed pressure.

- type:

  Effect type (direct, indirect, or direct + indirect).

- pathway:

  Pathway with which sensitivity has been assessed (expert- or
  model-based).

- sensitivity:

  Combined sensitivity score.

- adaptive_capacity:

  Combined adaptive capacity score.

- uncertainty_sens:

  Combined score of the associated sensitivity uncertainties.

- uncertainty_ac:

  Combined score of the associated adaptive capacity uncertainties.

- sens_original.sens_feeding:

  Original sensitivity score for the feeding trait.

- sens_original.sens_behaviour:

  Original sensitivity score for the behaviour trait.

- sens_original.sens_reproduction:

  Original sensitivity score for the reproduction trait.

- sens_original.sens_habitat:

  Original sensitivity score for the habitat trait.

- sens_original.sens_general:

  Original general sensitivity score.

- ac_original.ac_feeding:

  Original adaptive capacity score for the feeding trait.

- ac_original.ac_behaviour:

  Original adaptive capacity score for the behaviour trait.

- ac_original.ac_reproduction:

  Original adaptive capacity score for the reproduction trait.

- ac_original.ac_habitat:

  Original adaptive capacity score for the habitat trait.

- ac_original.ac_general:

  Original general adaptive capacity score.
