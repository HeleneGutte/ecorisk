# Calculate Overall Exposure Scores from Component-Specific Expert Ratings

Calculate exposure scores from individual exposure components.
Additionally calculate associated uncertainty scores.

## Usage

``` r
calc_exposure(
  pressures,
  components,
  probabilities = NULL,
  uncertainty = NULL,
  method = "mean"
)
```

## Arguments

- pressures:

  A character vector or column of a data frame containing the names of
  the pressures.

- components:

  A numeric vector or data frame containing the numeric values per
  exposure component for each pressure. Has to be in the same order as
  the pressure vector or come from the same data frame.

- probabilities:

  Optionally: a numeric vector containing the probabilities of each
  pressure (values have to be between 0 and 1); default is `NULL`. Has
  to be in the same order as the pressure vector.

- uncertainty:

  a numeric vector or data frame containing the associated uncertainty
  per component; default is `NULL`. Has to be in the same order as the
  pressure vector.

- method:

  a character string specifying the aggregation method. Available are
  mean (default), median, maximum, and minimum.

## Value

Returns a data frame containing the pressure name, the aggregated
exposure score and associated uncertainty scores. The results serve as
input to the
[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)
function. In case no uncertainty values were provided, `NA`s will be
returned as uncertainty scores.

## Details

Often exposure components include the magnitude of change compared to
baseline conditions, the frequency of this change, a future trend and a
spatial scale. These components are scored within the ecorisk framework
for each pressure by experts on a scale from 1 (low impact) to 5 (high
impact). To express their uncertainty during the process, experts can
score the associated uncertainty generally for all components of one
pressure or for each component individually. The uncertainty is scored
in the ecorisk framework on a scale from 1 (low uncertainty) to 3 (high
uncertainty). Using exposure and sensitivity scorings vulnerability is
calculated. Guidance for the scoring process can be found here:
[`create_template_exposure`](https://helenegutte.github.io/ecorisk/reference/create_template_exposure.md)
or in the vignette or in Gutte et al., 2025.

## See also

[`create_template_exposure`](https://helenegutte.github.io/ecorisk/reference/create_template_exposure.md),
[`create_template_sensitivity`](https://helenegutte.github.io/ecorisk/reference/create_template_sensitivity.md),
[`calc_sensitivity`](https://helenegutte.github.io/ecorisk/reference/calc_sensitivity.md),
[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)

## Examples

``` r
### Example using demo data with five pressures, four components and their individual
# uncertainties (probabilities are assumed to be 1):
ex_expert_exposure
#>      pressure magnitude frequency trend spatial uncertainty_magnitude
#> 1 temperature         1         1     4       5                     1
#> 2    salinity         1         4     1       3                     2
#> 3      oxygen         1         1     1       2                     2
#> 4    nutrient         2         2     3       2                     1
#> 5     fishing         5         4     5       2                     3
#>   uncertainty_frequency uncertainty_trend uncertainty_spatial
#> 1                     1                 3                   3
#> 2                     2                 2                   2
#> 3                     2                 3                   3
#> 4                     2                 2                   3
#> 5                     1                 2                   1

calc_exposure(
  pressures = ex_expert_exposure$pressure,
  components = ex_expert_exposure[ ,2:5],
  uncertainty = ex_expert_exposure[ ,6:9]
 )
#>      pressure exposure uncertainty
#> 1 temperature     2.75        2.00
#> 2    salinity     2.25        2.00
#> 3      oxygen     1.25        2.50
#> 4    nutrient     2.25        2.00
#> 5     fishing     4.00        1.75


### Example for two hazardous risks with only two components ('magnitude' and
#   'spatial'), one general uncertainty score, and associated probabilities:
hazard <- c("heat waves", "hurricanes")

# Create scoring table using the template function:
exp_tbl <- create_template_exposure(
  pressures = hazard,
  n_components = 2,
  mode_uncertainty = "general",
  probability = TRUE
)

names(exp_tbl)[2:3] <- c("magnitude", "spatial")
# Assign component-specific scores and probabilities:
exp_tbl$magnitude <- c(5,4)
exp_tbl$spatial <- c(5,3)
exp_tbl$uncertainty <- c(2,3)
exp_tbl$probability <- c(0.8,0.3)

# Calculate exposure score:
calc_exposure(
  pressures = exp_tbl$pressure,
  components = exp_tbl[ ,c("magnitude", "spatial")],
  probabilities = exp_tbl$probability,
  uncertainty = exp_tbl$uncertainty
 )
#>     pressure exposure uncertainty
#> 1 heat waves     4.00           2
#> 2 hurricanes     1.05           3
```
