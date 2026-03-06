# Create a Template for Expert-Based Exposure Scoring

The function `create_template_exposure` generates a template for a
semi-quantitative, expert-based scoring of individual exposure
components. This template is designed for use within the *ecorisk*
workflow. The user can define the number of exposure components (e.g.,
magnitude, frequency, trend, and spatial scale) to be included in the
scoring process and whether uncertainty should be assessed.

## Usage

``` r
create_template_exposure(
  pressures,
  n_components = 4,
  uncertainty = TRUE,
  mode_uncertainty = "general",
  probability = FALSE
)
```

## Arguments

- pressures:

  A character vector specifying the names of the pressures to be
  assessed.

- n_components:

  A positive integer or integer indicating the number of exposure
  components to be included in the assessment table. Within the ecorisk
  framework four components are usually scored:magnitude, frequency,
  trend and spatial scale. The default includes therefore four
  components.

- uncertainty:

  logical; should uncertainty be assessed? Default is `TRUE`.

- mode_uncertainty:

  character; if uncertainty is assessed, it can be scored for each
  exposure component and pressure individually (`"component"`) or as a
  general uncertainty score for each pressure (`"general"`, default).
  The default assigns one general uncertainty score per pressure.

- probability:

  logical; for hazardous risks that may have a lower probability of
  occurring within the assessed future time period, probabilities
  between 0 and 1 can be assigned to each pressure and considered in the
  [`calc_exposure`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md)
  function.

## Value

The function returns a data frame containing the names of the pressures
and the specified number of exposure components. If uncertainty is
assessed, the data frame also includes either one general uncertainty
column or one uncertainty column per component, depending on the
selected uncertainty assessment mode.

## Details

By default, the function creates a template with four exposure
components and a general uncertainty scoring column. Within the
*ecorisk* framework, exposure components should be scored on a scale
from 1 to 5 (low to high impact), while uncertainty should be scored
from 1 to 3 (low to high uncertainty).

The returned data frame can be exported as a CSV or Excel file.
Components can be renamed and additional components can be added if
necessary. The completed file can then be analyzed using the
[`calc_exposure`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md)
function.

The default scoring system is as follows:

- Exposure components: Scale of 1 to 5 (low to high impact).

- Uncertainty: Scale of 1 to 3 (low to high uncertainty).

## See also

[`create_template_sensitivity`](https://helenegutte.github.io/ecorisk/reference/create_template_sensitivity.md),
[`calc_exposure`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md),
[`calc_exposure`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md)

## Examples

``` r
### Create a full template for three pressures, including all four components
#   and a general uncertainty column:
press <- c("fishing", "temperature increase", "salinity decrease")

exp_tbl <- create_template_exposure(pressures = press)
# --> Export table and re-import after completion or fill in directly in R.

# Rename exposure components
names(exp_tbl)[2:5] <- c("magnitude", "frequency", "trend", "spatial")

# Assign individual scores between 1 (low) and 5 (high impact):
exp_tbl$magnitude <- c(1,5,4)
exp_tbl$frequency <- c(1,5,3)
exp_tbl$trend <- c(1,4,2)
exp_tbl$spatial <- c(1,5,5)

# Assign uncertainty scores from 1 (low) to 3 (high uncertainty):
exp_tbl$uncertainty <- c(1,1,3)


### Create a template for two more hazardous risks with only two components
#   ('magnitude' and 'spatial'), including component-specific uncertainties,
#   and a probability column:
hazard <- c("heat waves", "hurricanes")

exp_tbl <- create_template_exposure(
  pressures = hazard,
  n_components = 2,
  mode_uncertainty = "component",
  probability = TRUE
)

# Rename components and uncertainties
names(exp_tbl)[2:5] <- c("magnitude", "spatial", "uncertainty_magnitude",
                         "uncertainty_spatial")
exp_tbl$magnitude <- c(5,4)
exp_tbl$spatial <- c(5,3)
exp_tbl$uncertainty_magnitude <- c(1,1)
exp_tbl$uncertainty_spatial <- c(3,3)

# Assign probabilities of their occurrence within the assessed future period:
exp_tbl$probability <- c(0.8,0.3)
```
