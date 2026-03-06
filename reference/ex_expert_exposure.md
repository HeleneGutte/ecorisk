# Expert-based exposure scores for five pressures

This Baltic Sea demo dataset contains expert-assigned scores for five
environmental and anthropogenic pressures, detailing individual exposure
components and their associated uncertainties. The dataset was
initialized using the template function
[`create_template_exposure`](https://helenegutte.github.io/ecorisk/reference/create_template_exposure.md).

## Usage

``` r
ex_expert_exposure
```

## Format

A data frame with 5 observations and 9 variables (values randomly
assigned).

- pressure:

  Environmental or anthropogenic pressure.

- magnitude:

  Score for the magnitude pressure change.

- frequency:

  Score for the frequency or duration of pressure effect.

- trend:

  Score for the future trend of pressure.

- spatial:

  Score for the spatial extent of pressure change.

- uncertainty_magnitude:

  Uncertainty of magnitude score.

- uncertainty_frequency:

  Uncertainty of frequency score.

- uncertainty_trend:

  Uncertainty of trend score.

- uncertainty_spatial:

  Uncertainty of spatial extent score.

## Details

Exposure scores range from 1 (low) to 5 (high), while uncertainties
range from 1 (low) to 3 (high). This dataset can be used as input for
the function
[`calc_exposure`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md).
