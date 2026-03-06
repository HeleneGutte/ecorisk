# Create a Template for Expert-Based Sensitivity and Adaptive Capacity Scoring

The function `crt_sensitivity` creates a template for semi-quantitative,
expert-based sensitivity and, optionally, adaptive capacity scoring. The
template allows for assessing sensitivity and adaptive capacity using
either a general score for each state indicator-pressure combination or
a trait-based approach using life history traits. The latter is
particularly useful when state indicators represent individual species
and detailed biological information is available.

## Usage

``` r
create_template_sensitivity(
  indicators,
  pressures,
  type = "direct",
  n_sensitivity_traits = 1,
  adaptive_capacity = TRUE,
  mode_adaptive_capacity = "general",
  uncertainty = TRUE,
  mode_uncertainty = "general"
)
```

## Arguments

- indicators:

  A character vector specifying the names of the state indicators to
  assess.

- pressures:

  A character vector specifying the names of the pressures to assess.

- type:

  A character vector defining the type(s) of influence, such as
  `"direct"`, `"indirect"` or `"direct_indirect"`. The default is
  `"direct"`.

- n_sensitivity_traits:

  A positive integer specifying the number of traits used to assess
  sensitivity. The default is `1`, meaning that one general sensitivity
  score per state indicator-pressure-type combination is provided.

- adaptive_capacity:

  logical; should adaptive capacity be assessed? Default is `TRUE`.

- mode_adaptive_capacity:

  A character vector specifying whether adaptive capacity should be
  assessed for each trait individually (`"trait"`) or as a general score
  (`"general"`, default). Note: This parameter is only relevant when
  `n_sensitivity_traits > 1`.

- uncertainty:

  logical; should uncertainty be assessed? Default is `TRUE`. Note:
  Uncertainty is only added for components included in the assessment.
  If adaptive capacity (`adaptive_capacity = FALSE`) is not included, no
  uncertainty scoring will be applied to it.

- mode_uncertainty:

  A character vector specifying whether uncertainty should be assessed
  for each trait individually (`"trait"`) or as a general score
  (`"general"`).

## Value

A data frame where each row represents a state indicator-pressure-type
combination, containing the specified sensitivity traits, adaptive
capacity, and uncertainty columns (if selected).

- If adaptive capacity and uncertainty are assessed, the data frame
  includes either one general column or one column per trait, depending
  on the settings.

- If using trait-based scoring, the data frame includes trait-specific
  sensitivity, adaptive capacity, and uncertainty columns, which can be
  renamed as needed.

## Details

For each state indicator-pressure combination, different types of
influence can be assessed. The type of influence describes whether the
pressure acts directly, indirectly, or as a combination of both, which
is important for identifying impact pathways and potential management
measures.

Within the *ecorisk* framework, it is recommended to use negative scores
(-1 to -5) to indicate negative impacts (low to high severity) and
positive scores (1 to 5) for positive effects of a pressure on an
indicator. If an indicator is not sensitive to a pressure, the score
should be 0. Adaptive capacity is scored from -1 (no adaptive capacity)
to 1 (high adaptive capacity).

To improve the reliability of the scoring, uncertainty should also be
assessed. Uncertainty can be scored for each trait individually or as a
general score and should be rated on a scale from 1 to 3 (low to high
uncertainty).

The returned data frame can be exported as a CSV or Excel file. Column
names can be modified as needed. The completed file can be analyzed
using the
[`calc_sensitivity`](https://helenegutte.github.io/ecorisk/reference/calc_sensitivity.md)
function.

Depending on the settings, the data frame will include:

- A single "sensitivity" column, if using general scoring.

- Multiple trait-specific sensitivity columns (e.g.,
  "sensitivity_trait_1", "sensitivity_trait_2", etc.), which can be
  renamed as necessary.

Within this data frame, trait-based and general scoring can be mixed. It
is therefore recommended to set `n_sensitivity_traits` to the maximum
number of traits to be assessed for any state indicator. The
[`calc_sensitivity`](https://helenegutte.github.io/ecorisk/reference/calc_sensitivity.md)
function automatically distinguishes between general and trait-based
scoring.

## See also

[`create_template_exposure`](https://helenegutte.github.io/ecorisk/reference/create_template_exposure.md),
[`calc_exposure`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md),
[`calc_sensitivity`](https://helenegutte.github.io/ecorisk/reference/calc_sensitivity.md)

## Examples

``` r
### Create a table for two state indicators and two pressures to evaluate direct
#   effects (default). Return a general sensitivity and adaptive capacity
#   column as well as uncertainty columns for both components:
ind <- c("seabirds", "seals")
press <- c("plastic pollution", "temperature increase")

sens_ac_tbl <- create_template_sensitivity(
  indicators = ind,
  pressures = press
)
# --> Export table and re-import after completion or fill in directly in R.

# Assign sensitivity scores from -5 (strong negative response to pressure)
# to +5 (strong positive response) (0 = no sensitivity):
sens_ac_tbl$sens_general <- c(-5,3,-4,4)

# Assign adaptive capacity scores from -1 (none) to +1 (good adaptive capacity):
sens_ac_tbl$ac_general <- c(-1,1,-1,1)

# Assign uncertainty scores from 1 (low) to 3 (high uncertainty):
sens_ac_tbl$uncertainty_sens <- c(1,2,1,1)
sens_ac_tbl$uncertainty_ac <- c(3,2,3,2)


### Create a table for four indicators and three pressures to evaluate both direct
#   and indirect effects. Return columns for five trait-specific sensitivities
#   and their respective uncertainties, but no adaptive capacity:
ind <- c("cod", "herring", "seabirds", "seals")
press <- c("fishing", "temperature increase", "salinity decrease")
sens_ac_tbl <- create_template_sensitivity(
  indicators = ind,
  pressures = press,
  type = c("direct", "direct_indirect"),
  n_sensitivity_traits = 5,
  adaptive_capacity = FALSE,
  uncertainty = TRUE,
  mode_uncertainty = "trait"
)
sens_ac_tbl
#>    indicator             pressure            type sens_trait_1 sens_trait_2
#> 1        cod              fishing          direct           NA           NA
#> 2        cod temperature increase          direct           NA           NA
#> 3        cod    salinity decrease          direct           NA           NA
#> 4        cod              fishing direct_indirect           NA           NA
#> 5        cod temperature increase direct_indirect           NA           NA
#> 6        cod    salinity decrease direct_indirect           NA           NA
#> 7    herring              fishing          direct           NA           NA
#> 8    herring temperature increase          direct           NA           NA
#> 9    herring    salinity decrease          direct           NA           NA
#> 10   herring              fishing direct_indirect           NA           NA
#> 11   herring temperature increase direct_indirect           NA           NA
#> 12   herring    salinity decrease direct_indirect           NA           NA
#> 13  seabirds              fishing          direct           NA           NA
#> 14  seabirds temperature increase          direct           NA           NA
#> 15  seabirds    salinity decrease          direct           NA           NA
#> 16  seabirds              fishing direct_indirect           NA           NA
#> 17  seabirds temperature increase direct_indirect           NA           NA
#> 18  seabirds    salinity decrease direct_indirect           NA           NA
#> 19     seals              fishing          direct           NA           NA
#> 20     seals temperature increase          direct           NA           NA
#> 21     seals    salinity decrease          direct           NA           NA
#> 22     seals              fishing direct_indirect           NA           NA
#> 23     seals temperature increase direct_indirect           NA           NA
#> 24     seals    salinity decrease direct_indirect           NA           NA
#>    sens_trait_3 sens_trait_4 sens_trait_5 uncertainty_sens_trait_1
#> 1            NA           NA           NA                       NA
#> 2            NA           NA           NA                       NA
#> 3            NA           NA           NA                       NA
#> 4            NA           NA           NA                       NA
#> 5            NA           NA           NA                       NA
#> 6            NA           NA           NA                       NA
#> 7            NA           NA           NA                       NA
#> 8            NA           NA           NA                       NA
#> 9            NA           NA           NA                       NA
#> 10           NA           NA           NA                       NA
#> 11           NA           NA           NA                       NA
#> 12           NA           NA           NA                       NA
#> 13           NA           NA           NA                       NA
#> 14           NA           NA           NA                       NA
#> 15           NA           NA           NA                       NA
#> 16           NA           NA           NA                       NA
#> 17           NA           NA           NA                       NA
#> 18           NA           NA           NA                       NA
#> 19           NA           NA           NA                       NA
#> 20           NA           NA           NA                       NA
#> 21           NA           NA           NA                       NA
#> 22           NA           NA           NA                       NA
#> 23           NA           NA           NA                       NA
#> 24           NA           NA           NA                       NA
#>    uncertainty_sens_trait_2 uncertainty_sens_trait_3 uncertainty_sens_trait_4
#> 1                        NA                       NA                       NA
#> 2                        NA                       NA                       NA
#> 3                        NA                       NA                       NA
#> 4                        NA                       NA                       NA
#> 5                        NA                       NA                       NA
#> 6                        NA                       NA                       NA
#> 7                        NA                       NA                       NA
#> 8                        NA                       NA                       NA
#> 9                        NA                       NA                       NA
#> 10                       NA                       NA                       NA
#> 11                       NA                       NA                       NA
#> 12                       NA                       NA                       NA
#> 13                       NA                       NA                       NA
#> 14                       NA                       NA                       NA
#> 15                       NA                       NA                       NA
#> 16                       NA                       NA                       NA
#> 17                       NA                       NA                       NA
#> 18                       NA                       NA                       NA
#> 19                       NA                       NA                       NA
#> 20                       NA                       NA                       NA
#> 21                       NA                       NA                       NA
#> 22                       NA                       NA                       NA
#> 23                       NA                       NA                       NA
#> 24                       NA                       NA                       NA
#>    uncertainty_sens_trait_5
#> 1                        NA
#> 2                        NA
#> 3                        NA
#> 4                        NA
#> 5                        NA
#> 6                        NA
#> 7                        NA
#> 8                        NA
#> 9                        NA
#> 10                       NA
#> 11                       NA
#> 12                       NA
#> 13                       NA
#> 14                       NA
#> 15                       NA
#> 16                       NA
#> 17                       NA
#> 18                       NA
#> 19                       NA
#> 20                       NA
#> 21                       NA
#> 22                       NA
#> 23                       NA
#> 24                       NA
# --> You might want to rename the generic trait columns with specific traits.
# --> Export table as e.g. CSV-file and re-import again after completion.

### Create a mixed table for two indicators and two pressures, where for one
#   indicator sensitivity is scored overall and for one sensitivity is scored
#   by individual traits:
ind <- c("phytoplankton", "herring")
press <- c("temperature", "salinity")

sens_ac_tbl <- create_template_sensitivity(
  indicators = ind,
  pressures = press,
  n_sensitivity_traits = 4,
  adaptive_capacity = FALSE,
  uncertainty = TRUE,
  mode_uncertainty = "general"
)
# Rename trait columns:
names(sens_ac_tbl)[4:7] <- paste0("sens_",
  c("feeding", "behaviour", "reproduction", "general"))

# Give overall sensitivity score for phytoplankton
# (keep NAs for herring):
sens_ac_tbl$sens_general[1:2] <- c(-3,0)

# Give trait-specific sensitivity scores for herring
# (keep NAs for phytoplankton):
sens_ac_tbl$sens_feeding[3:4] <- c(0,0)
sens_ac_tbl$sens_behaviour[3:4] <- c(-1,0)
sens_ac_tbl$sens_reproduction[3:4] <- c(-2,-2)

# Give overall uncertainty score:
sens_ac_tbl$uncertainty_sens <- c(1,2,1,1)
```
