# Calculate Overall Sensitivity and Adaptive Capacity Scores from Trait-Specific Expert Ratings

The function `calc_sensitivity` calculates aggregated sensitivity and
adaptive capacity scores. Additionally it prepares the scoring data for
the further usage in the `vulnerability` function. The scores for
sensitivity and adaptive capacity can be trait-based or general (one
score per state indicator and pressure combination).

## Usage

``` r
calc_sensitivity(
  indicators,
  pressures,
  type = "direct",
  sensitivity_traits,
  adaptive_capacities = NULL,
  uncertainty_sens = NULL,
  uncertainty_ac = NULL,
  method = "mean"
)
```

## Arguments

- indicators:

  A character vector or column of a data frame containing the names of
  the state indicators.

- pressures:

  A character vector or column of a data frame containing the names of
  the pressures.

- type:

  A character vector or column of a data frame specifying the effect
  type. Effects could be direct or indirect or a combination of both.
  Default is `direct`.

- sensitivity_traits:

  A data frame containing the numeric sensitivity values per species
  trait or as a general value.

- adaptive_capacities:

  A data frame (or single vector) containing the numeric values for the
  adaptive capacity. Either per trait or as one general value. When
  values are given per trait, they have to be in the same order as the
  values for the sensitivity. The default is `NULL`.

- uncertainty_sens:

  A data frame (or a vector) containing the numeric uncertainty values
  associated with the sensitivity; default is `NULL`.

- uncertainty_ac:

  A data frame (or a vector) containing the numeric uncertainty values
  associated with the adaptive capacity; default is `NULL`.

- method:

  A character string specifying the method of the aggregation of the
  traits. Available are `median`, `minimum`, `maximum` and `mean`, the
  `mean` is default.

## Value

a data frame containing the indicator, pressure and effect type, the
aggregated sensitivity and adaptive capacity as well as their associated
uncertainty scores. Additionally, the trait specific sensitivity and
adaptive capacity scores are stored and used later as input for the
`vulnerability` function.

## Details

The function calculates per state indicator and pressure combination one
aggregated sensitivity and adaptive capacity score, the aggregation
method can be determined with the parameter `method`. The assessment of
adaptive capacity is optionally, if no scores for adaptive capacity are
provided the function calculates an aggregated sensitivity score and
prepares only the sensitivity scores for the
[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)
function. Guidance for the scoring process can be found here:
[`create_template_sensitivity`](https://helenegutte.github.io/ecorisk/reference/create_template_sensitivity.md)
or in the vignette or in Gutte et al., 2025. Using exposure and
sensitivity scorings, vulnerability is calculated.

## See also

[`create_template_exposure`](https://helenegutte.github.io/ecorisk/reference/create_template_exposure.md),
[`create_template_sensitivity`](https://helenegutte.github.io/ecorisk/reference/create_template_sensitivity.md),
[`calc_exposure`](https://helenegutte.github.io/ecorisk/reference/calc_exposure.md),
[`vulnerability`](https://helenegutte.github.io/ecorisk/reference/vulnerability.md)

## Examples

``` r
### Example using demo data with four indicators and five pressures with
#   scores for direct as well as combined direct-indirect effects based on
#   the template function create_template_sensitivity(). For two
#   indicators, sensitivity, adaptive capacity, and their uncertainties are
#   provided as general scores, while for the other two, they are based on
#   individual traits.
ex_expert_sensitivity
#>        indicator    pressure            type sens_feeding sens_behaviour
#> 1  phytoplankton temperature          direct           NA             NA
#> 2  phytoplankton    salinity          direct           NA             NA
#> 3  phytoplankton      oxygen          direct           NA             NA
#> 4  phytoplankton    nutrient          direct           NA             NA
#> 5  phytoplankton     fishing          direct           NA             NA
#> 6  phytoplankton temperature direct_indirect           NA             NA
#> 7  phytoplankton    salinity direct_indirect           NA             NA
#> 8  phytoplankton      oxygen direct_indirect           NA             NA
#> 9  phytoplankton    nutrient direct_indirect           NA             NA
#> 10 phytoplankton     fishing direct_indirect           NA             NA
#> 11       herring temperature          direct            0             -1
#> 12       herring    salinity          direct            0              0
#> 13       herring      oxygen          direct            0             -1
#> 14       herring    nutrient          direct            0              0
#> 15       herring     fishing          direct            0             -4
#> 16       herring temperature direct_indirect           -3             -3
#> 17       herring    salinity direct_indirect           -2              0
#> 18       herring      oxygen direct_indirect           -3             -3
#> 19       herring    nutrient direct_indirect            1              0
#> 20       herring     fishing direct_indirect           -2             -5
#> 21           cod temperature          direct            0             -2
#> 22           cod    salinity          direct            0              0
#> 23           cod      oxygen          direct            0             -2
#> 24           cod    nutrient          direct            0              0
#> 25           cod     fishing          direct            0             -4
#> 26           cod temperature direct_indirect           -4             -3
#> 27           cod    salinity direct_indirect           -2             -1
#> 28           cod      oxygen direct_indirect           -3             -3
#> 29           cod    nutrient direct_indirect           -1             -1
#> 30           cod     fishing direct_indirect           -3             -5
#> 31      seabirds temperature          direct           NA             NA
#> 32      seabirds    salinity          direct           NA             NA
#> 33      seabirds      oxygen          direct           NA             NA
#> 34      seabirds    nutrient          direct           NA             NA
#> 35      seabirds     fishing          direct           NA             NA
#> 36      seabirds temperature direct_indirect           NA             NA
#> 37      seabirds    salinity direct_indirect           NA             NA
#> 38      seabirds      oxygen direct_indirect           NA             NA
#> 39      seabirds    nutrient direct_indirect           NA             NA
#> 40      seabirds     fishing direct_indirect           NA             NA
#>    sens_reproduction sens_habitat sens_general ac_feeding ac_behaviour
#> 1                 NA           NA           -3         NA           NA
#> 2                 NA           NA            0         NA           NA
#> 3                 NA           NA           -2         NA           NA
#> 4                 NA           NA           -3         NA           NA
#> 5                 NA           NA            0         NA           NA
#> 6                 NA           NA           -4         NA           NA
#> 7                 NA           NA            0         NA           NA
#> 8                 NA           NA           -3         NA           NA
#> 9                 NA           NA           -4         NA           NA
#> 10                NA           NA            0         NA           NA
#> 11                -2           -3           NA          0            0
#> 12                -2           -2           NA          0            0
#> 13                -2           -3           NA          0            0
#> 14                 0            2           NA          0            1
#> 15                -5            0           NA          0           -1
#> 16                -3           -4           NA         -1           -1
#> 17                -3           -3           NA          0            0
#> 18                -3           -4           NA          0            0
#> 19                 2            3           NA          1            1
#> 20                -5            0           NA         -1           -1
#> 21                -3           -4           NA          0           -1
#> 22                -4           -4           NA          0            0
#> 23                -4           -4           NA          0           -1
#> 24                 0           -2           NA          0            0
#> 25                -5           -5           NA         -1           -1
#> 26                -5           -5           NA         -1           -1
#> 27                -4           -5           NA         -1           -1
#> 28                -4           -5           NA         -1           -1
#> 29                -2           -3           NA          1            1
#> 30                -5           -5           NA         -1           -1
#> 31                NA           NA           -2         NA           NA
#> 32                NA           NA            0         NA           NA
#> 33                NA           NA            0         NA           NA
#> 34                NA           NA            0         NA           NA
#> 35                NA           NA           -3         NA           NA
#> 36                NA           NA           -4         NA           NA
#> 37                NA           NA           -2         NA           NA
#> 38                NA           NA           -2         NA           NA
#> 39                NA           NA           -2         NA           NA
#> 40                NA           NA           -5         NA           NA
#>    ac_reproduction ac_habitat ac_general uncertainty_sens_feeding
#> 1               NA         NA          1                       NA
#> 2               NA         NA          1                       NA
#> 3               NA         NA          1                       NA
#> 4               NA         NA          1                       NA
#> 5               NA         NA          0                       NA
#> 6               NA         NA          1                       NA
#> 7               NA         NA          1                       NA
#> 8               NA         NA          1                       NA
#> 9               NA         NA          1                       NA
#> 10              NA         NA          0                       NA
#> 11               0          0         NA                        1
#> 12               0          0         NA                        1
#> 13               0          0         NA                        2
#> 14               1          1         NA                        2
#> 15              -1         -1         NA                        1
#> 16              -1         -1         NA                        1
#> 17               0          0         NA                        1
#> 18               0          0         NA                        2
#> 19               1          1         NA                        2
#> 20              -1         -1         NA                        1
#> 21              -1         -1         NA                        1
#> 22               0          0         NA                        1
#> 23              -1         -1         NA                        2
#> 24               0          0         NA                        2
#> 25              -1         -1         NA                        1
#> 26              -1         -1         NA                        1
#> 27              -1         -1         NA                        1
#> 28              -1         -1         NA                        2
#> 29               1          1         NA                        2
#> 30              -1         -1         NA                        1
#> 31              NA         NA          1                       NA
#> 32              NA         NA          0                       NA
#> 33              NA         NA          0                       NA
#> 34              NA         NA          0                       NA
#> 35              NA         NA         -1                       NA
#> 36              NA         NA          0                       NA
#> 37              NA         NA          0                       NA
#> 38              NA         NA          0                       NA
#> 39              NA         NA          0                       NA
#> 40              NA         NA         -1                       NA
#>    uncertainty_sens_behaviour uncertainty_sens_reproduction
#> 1                          NA                            NA
#> 2                          NA                            NA
#> 3                          NA                            NA
#> 4                          NA                            NA
#> 5                          NA                            NA
#> 6                          NA                            NA
#> 7                          NA                            NA
#> 8                          NA                            NA
#> 9                          NA                            NA
#> 10                         NA                            NA
#> 11                          1                             1
#> 12                          1                             1
#> 13                          2                             2
#> 14                          2                             2
#> 15                          1                             1
#> 16                          1                             1
#> 17                          1                             1
#> 18                          2                             2
#> 19                          2                             2
#> 20                          1                             1
#> 21                          1                             1
#> 22                          1                             1
#> 23                          2                             2
#> 24                          2                             2
#> 25                          1                             1
#> 26                          1                             1
#> 27                          1                             1
#> 28                          2                             2
#> 29                          2                             2
#> 30                          1                             1
#> 31                         NA                            NA
#> 32                         NA                            NA
#> 33                         NA                            NA
#> 34                         NA                            NA
#> 35                         NA                            NA
#> 36                         NA                            NA
#> 37                         NA                            NA
#> 38                         NA                            NA
#> 39                         NA                            NA
#> 40                         NA                            NA
#>    uncertainty_sens_habitat uncertainty_sens_general uncertainty_ac_feeding
#> 1                        NA                        1                     NA
#> 2                        NA                        2                     NA
#> 3                        NA                        2                     NA
#> 4                        NA                        2                     NA
#> 5                        NA                        1                     NA
#> 6                        NA                        1                     NA
#> 7                        NA                        2                     NA
#> 8                        NA                        2                     NA
#> 9                        NA                        2                     NA
#> 10                       NA                        1                     NA
#> 11                        1                       NA                      1
#> 12                        1                       NA                      1
#> 13                        2                       NA                      2
#> 14                        2                       NA                      2
#> 15                        1                       NA                      1
#> 16                        1                       NA                      1
#> 17                        1                       NA                      1
#> 18                        2                       NA                      2
#> 19                        2                       NA                      2
#> 20                        1                       NA                      1
#> 21                        1                       NA                      1
#> 22                        1                       NA                      1
#> 23                        2                       NA                      2
#> 24                        2                       NA                      2
#> 25                        1                       NA                      1
#> 26                        1                       NA                      1
#> 27                        1                       NA                      1
#> 28                        2                       NA                      2
#> 29                        2                       NA                      2
#> 30                        1                       NA                      1
#> 31                       NA                        1                     NA
#> 32                       NA                        2                     NA
#> 33                       NA                        2                     NA
#> 34                       NA                        2                     NA
#> 35                       NA                        1                     NA
#> 36                       NA                        1                     NA
#> 37                       NA                        2                     NA
#> 38                       NA                        2                     NA
#> 39                       NA                        2                     NA
#> 40                       NA                        1                     NA
#>    uncertainty_ac_behaviour uncertainty_ac_reproduction uncertainty_ac_habitat
#> 1                        NA                          NA                     NA
#> 2                        NA                          NA                     NA
#> 3                        NA                          NA                     NA
#> 4                        NA                          NA                     NA
#> 5                        NA                          NA                     NA
#> 6                        NA                          NA                     NA
#> 7                        NA                          NA                     NA
#> 8                        NA                          NA                     NA
#> 9                        NA                          NA                     NA
#> 10                       NA                          NA                     NA
#> 11                        1                           1                      1
#> 12                        1                           1                      1
#> 13                        2                           2                      2
#> 14                        2                           2                      2
#> 15                        1                           1                      1
#> 16                        1                           1                      1
#> 17                        1                           1                      1
#> 18                        2                           2                      2
#> 19                        2                           2                      2
#> 20                        1                           1                      1
#> 21                        1                           1                      1
#> 22                        1                           1                      1
#> 23                        2                           2                      2
#> 24                        2                           2                      2
#> 25                        1                           1                      1
#> 26                        1                           1                      1
#> 27                        1                           1                      1
#> 28                        2                           2                      2
#> 29                        2                           2                      2
#> 30                        1                           1                      1
#> 31                       NA                          NA                     NA
#> 32                       NA                          NA                     NA
#> 33                       NA                          NA                     NA
#> 34                       NA                          NA                     NA
#> 35                       NA                          NA                     NA
#> 36                       NA                          NA                     NA
#> 37                       NA                          NA                     NA
#> 38                       NA                          NA                     NA
#> 39                       NA                          NA                     NA
#> 40                       NA                          NA                     NA
#>    uncertainty_ac_general
#> 1                       1
#> 2                       2
#> 3                       2
#> 4                       2
#> 5                       1
#> 6                       1
#> 7                       2
#> 8                       2
#> 9                       2
#> 10                      1
#> 11                     NA
#> 12                     NA
#> 13                     NA
#> 14                     NA
#> 15                     NA
#> 16                     NA
#> 17                     NA
#> 18                     NA
#> 19                     NA
#> 20                     NA
#> 21                     NA
#> 22                     NA
#> 23                     NA
#> 24                     NA
#> 25                     NA
#> 26                     NA
#> 27                     NA
#> 28                     NA
#> 29                     NA
#> 30                     NA
#> 31                      2
#> 32                      3
#> 33                      3
#> 34                      3
#> 35                      3
#> 36                      2
#> 37                      3
#> 38                      3
#> 39                      3
#> 40                      3

# Calculate only mean sensitivity scores:
calc_sensitivity(
  indicators = ex_expert_sensitivity$indicator,
  pressures = ex_expert_sensitivity$pressure,
  sensitivity_traits = ex_expert_sensitivity[ ,4:8],
  adaptive_capacities = NULL,   # (default)
  uncertainty_sens  = NULL,     # (default)
  uncertainty_ac = NULL,        # (default)
  method = "mean"               # (default)
 )
#>        indicator    pressure   type pathway sensitivity adaptive_capacity
#> 1  phytoplankton temperature direct  expert       -3.00                 0
#> 2  phytoplankton    salinity direct  expert        0.00                 0
#> 3  phytoplankton      oxygen direct  expert       -2.00                 0
#> 4  phytoplankton    nutrient direct  expert       -3.00                 0
#> 5  phytoplankton     fishing direct  expert        0.00                 0
#> 6  phytoplankton temperature direct  expert       -4.00                 0
#> 7  phytoplankton    salinity direct  expert        0.00                 0
#> 8  phytoplankton      oxygen direct  expert       -3.00                 0
#> 9  phytoplankton    nutrient direct  expert       -4.00                 0
#> 10 phytoplankton     fishing direct  expert        0.00                 0
#> 11       herring temperature direct  expert       -1.50                 0
#> 12       herring    salinity direct  expert       -1.00                 0
#> 13       herring      oxygen direct  expert       -1.50                 0
#> 14       herring    nutrient direct  expert        0.50                 0
#> 15       herring     fishing direct  expert       -2.25                 0
#> 16       herring temperature direct  expert       -3.25                 0
#> 17       herring    salinity direct  expert       -2.00                 0
#> 18       herring      oxygen direct  expert       -3.25                 0
#> 19       herring    nutrient direct  expert        1.50                 0
#> 20       herring     fishing direct  expert       -3.00                 0
#> 21           cod temperature direct  expert       -2.25                 0
#> 22           cod    salinity direct  expert       -2.00                 0
#> 23           cod      oxygen direct  expert       -2.50                 0
#> 24           cod    nutrient direct  expert       -0.50                 0
#> 25           cod     fishing direct  expert       -3.50                 0
#> 26           cod temperature direct  expert       -4.25                 0
#> 27           cod    salinity direct  expert       -3.00                 0
#> 28           cod      oxygen direct  expert       -3.75                 0
#> 29           cod    nutrient direct  expert       -1.75                 0
#> 30           cod     fishing direct  expert       -4.50                 0
#> 31      seabirds temperature direct  expert       -2.00                 0
#> 32      seabirds    salinity direct  expert        0.00                 0
#> 33      seabirds      oxygen direct  expert        0.00                 0
#> 34      seabirds    nutrient direct  expert        0.00                 0
#> 35      seabirds     fishing direct  expert       -3.00                 0
#> 36      seabirds temperature direct  expert       -4.00                 0
#> 37      seabirds    salinity direct  expert       -2.00                 0
#> 38      seabirds      oxygen direct  expert       -2.00                 0
#> 39      seabirds    nutrient direct  expert       -2.00                 0
#> 40      seabirds     fishing direct  expert       -5.00                 0
#>    uncertainty_sens uncertainty_ac sens_original.sens_feeding
#> 1                NA             NA                         NA
#> 2                NA             NA                         NA
#> 3                NA             NA                         NA
#> 4                NA             NA                         NA
#> 5                NA             NA                         NA
#> 6                NA             NA                         NA
#> 7                NA             NA                         NA
#> 8                NA             NA                         NA
#> 9                NA             NA                         NA
#> 10               NA             NA                         NA
#> 11               NA             NA                          0
#> 12               NA             NA                          0
#> 13               NA             NA                          0
#> 14               NA             NA                          0
#> 15               NA             NA                          0
#> 16               NA             NA                         -3
#> 17               NA             NA                         -2
#> 18               NA             NA                         -3
#> 19               NA             NA                          1
#> 20               NA             NA                         -2
#> 21               NA             NA                          0
#> 22               NA             NA                          0
#> 23               NA             NA                          0
#> 24               NA             NA                          0
#> 25               NA             NA                          0
#> 26               NA             NA                         -4
#> 27               NA             NA                         -2
#> 28               NA             NA                         -3
#> 29               NA             NA                         -1
#> 30               NA             NA                         -3
#> 31               NA             NA                         NA
#> 32               NA             NA                         NA
#> 33               NA             NA                         NA
#> 34               NA             NA                         NA
#> 35               NA             NA                         NA
#> 36               NA             NA                         NA
#> 37               NA             NA                         NA
#> 38               NA             NA                         NA
#> 39               NA             NA                         NA
#> 40               NA             NA                         NA
#>    sens_original.sens_behaviour sens_original.sens_reproduction
#> 1                            NA                              NA
#> 2                            NA                              NA
#> 3                            NA                              NA
#> 4                            NA                              NA
#> 5                            NA                              NA
#> 6                            NA                              NA
#> 7                            NA                              NA
#> 8                            NA                              NA
#> 9                            NA                              NA
#> 10                           NA                              NA
#> 11                           -1                              -2
#> 12                            0                              -2
#> 13                           -1                              -2
#> 14                            0                               0
#> 15                           -4                              -5
#> 16                           -3                              -3
#> 17                            0                              -3
#> 18                           -3                              -3
#> 19                            0                               2
#> 20                           -5                              -5
#> 21                           -2                              -3
#> 22                            0                              -4
#> 23                           -2                              -4
#> 24                            0                               0
#> 25                           -4                              -5
#> 26                           -3                              -5
#> 27                           -1                              -4
#> 28                           -3                              -4
#> 29                           -1                              -2
#> 30                           -5                              -5
#> 31                           NA                              NA
#> 32                           NA                              NA
#> 33                           NA                              NA
#> 34                           NA                              NA
#> 35                           NA                              NA
#> 36                           NA                              NA
#> 37                           NA                              NA
#> 38                           NA                              NA
#> 39                           NA                              NA
#> 40                           NA                              NA
#>    sens_original.sens_habitat sens_original.sens_general ac_original.ac_general
#> 1                          NA                         -3                     NA
#> 2                          NA                          0                     NA
#> 3                          NA                         -2                     NA
#> 4                          NA                         -3                     NA
#> 5                          NA                          0                     NA
#> 6                          NA                         -4                     NA
#> 7                          NA                          0                     NA
#> 8                          NA                         -3                     NA
#> 9                          NA                         -4                     NA
#> 10                         NA                          0                     NA
#> 11                         -3                         NA                     NA
#> 12                         -2                         NA                     NA
#> 13                         -3                         NA                     NA
#> 14                          2                         NA                     NA
#> 15                          0                         NA                     NA
#> 16                         -4                         NA                     NA
#> 17                         -3                         NA                     NA
#> 18                         -4                         NA                     NA
#> 19                          3                         NA                     NA
#> 20                          0                         NA                     NA
#> 21                         -4                         NA                     NA
#> 22                         -4                         NA                     NA
#> 23                         -4                         NA                     NA
#> 24                         -2                         NA                     NA
#> 25                         -5                         NA                     NA
#> 26                         -5                         NA                     NA
#> 27                         -5                         NA                     NA
#> 28                         -5                         NA                     NA
#> 29                         -3                         NA                     NA
#> 30                         -5                         NA                     NA
#> 31                         NA                         -2                     NA
#> 32                         NA                          0                     NA
#> 33                         NA                          0                     NA
#> 34                         NA                          0                     NA
#> 35                         NA                         -3                     NA
#> 36                         NA                         -4                     NA
#> 37                         NA                         -2                     NA
#> 38                         NA                         -2                     NA
#> 39                         NA                         -2                     NA
#> 40                         NA                         -5                     NA

# Calculate mean scores for sensitivity, adaptive capacity and
# associated uncertainties:
calc_sensitivity(
  indicators = ex_expert_sensitivity$indicator,
  pressures = ex_expert_sensitivity$pressure,
  type = ex_expert_sensitivity$type,
  sensitivity_traits = ex_expert_sensitivity[ ,4:8],
  adaptive_capacities = ex_expert_sensitivity[ ,9:13],
  uncertainty_sens  = ex_expert_sensitivity[ ,14:18],
  uncertainty_ac = ex_expert_sensitivity[ ,19:23]
 )
#>        indicator    pressure            type pathway sensitivity
#> 1  phytoplankton temperature          direct  expert       -3.00
#> 2  phytoplankton    salinity          direct  expert        0.00
#> 3  phytoplankton      oxygen          direct  expert       -2.00
#> 4  phytoplankton    nutrient          direct  expert       -3.00
#> 5  phytoplankton     fishing          direct  expert        0.00
#> 6  phytoplankton temperature direct_indirect  expert       -4.00
#> 7  phytoplankton    salinity direct_indirect  expert        0.00
#> 8  phytoplankton      oxygen direct_indirect  expert       -3.00
#> 9  phytoplankton    nutrient direct_indirect  expert       -4.00
#> 10 phytoplankton     fishing direct_indirect  expert        0.00
#> 11       herring temperature          direct  expert       -1.50
#> 12       herring    salinity          direct  expert       -1.00
#> 13       herring      oxygen          direct  expert       -1.50
#> 14       herring    nutrient          direct  expert        0.50
#> 15       herring     fishing          direct  expert       -2.25
#> 16       herring temperature direct_indirect  expert       -3.25
#> 17       herring    salinity direct_indirect  expert       -2.00
#> 18       herring      oxygen direct_indirect  expert       -3.25
#> 19       herring    nutrient direct_indirect  expert        1.50
#> 20       herring     fishing direct_indirect  expert       -3.00
#> 21           cod temperature          direct  expert       -2.25
#> 22           cod    salinity          direct  expert       -2.00
#> 23           cod      oxygen          direct  expert       -2.50
#> 24           cod    nutrient          direct  expert       -0.50
#> 25           cod     fishing          direct  expert       -3.50
#> 26           cod temperature direct_indirect  expert       -4.25
#> 27           cod    salinity direct_indirect  expert       -3.00
#> 28           cod      oxygen direct_indirect  expert       -3.75
#> 29           cod    nutrient direct_indirect  expert       -1.75
#> 30           cod     fishing direct_indirect  expert       -4.50
#> 31      seabirds temperature          direct  expert       -2.00
#> 32      seabirds    salinity          direct  expert        0.00
#> 33      seabirds      oxygen          direct  expert        0.00
#> 34      seabirds    nutrient          direct  expert        0.00
#> 35      seabirds     fishing          direct  expert       -3.00
#> 36      seabirds temperature direct_indirect  expert       -4.00
#> 37      seabirds    salinity direct_indirect  expert       -2.00
#> 38      seabirds      oxygen direct_indirect  expert       -2.00
#> 39      seabirds    nutrient direct_indirect  expert       -2.00
#> 40      seabirds     fishing direct_indirect  expert       -5.00
#>    adaptive_capacity uncertainty_sens uncertainty_ac sens_original.sens_feeding
#> 1               1.00                1              1                         NA
#> 2               1.00                2              2                         NA
#> 3               1.00                2              2                         NA
#> 4               1.00                2              2                         NA
#> 5               0.00                1              1                         NA
#> 6               1.00                1              1                         NA
#> 7               1.00                2              2                         NA
#> 8               1.00                2              2                         NA
#> 9               1.00                2              2                         NA
#> 10              0.00                1              1                         NA
#> 11              0.00                1              1                          0
#> 12              0.00                1              1                          0
#> 13              0.00                2              2                          0
#> 14              0.75                2              2                          0
#> 15             -0.75                1              1                          0
#> 16             -1.00                1              1                         -3
#> 17              0.00                1              1                         -2
#> 18              0.00                2              2                         -3
#> 19              1.00                2              2                          1
#> 20             -1.00                1              1                         -2
#> 21             -0.75                1              1                          0
#> 22              0.00                1              1                          0
#> 23             -0.75                2              2                          0
#> 24              0.00                2              2                          0
#> 25             -1.00                1              1                          0
#> 26             -1.00                1              1                         -4
#> 27             -1.00                1              1                         -2
#> 28             -1.00                2              2                         -3
#> 29              1.00                2              2                         -1
#> 30             -1.00                1              1                         -3
#> 31              1.00                1              2                         NA
#> 32              0.00                2              3                         NA
#> 33              0.00                2              3                         NA
#> 34              0.00                2              3                         NA
#> 35             -1.00                1              3                         NA
#> 36              0.00                1              2                         NA
#> 37              0.00                2              3                         NA
#> 38              0.00                2              3                         NA
#> 39              0.00                2              3                         NA
#> 40             -1.00                1              3                         NA
#>    sens_original.sens_behaviour sens_original.sens_reproduction
#> 1                            NA                              NA
#> 2                            NA                              NA
#> 3                            NA                              NA
#> 4                            NA                              NA
#> 5                            NA                              NA
#> 6                            NA                              NA
#> 7                            NA                              NA
#> 8                            NA                              NA
#> 9                            NA                              NA
#> 10                           NA                              NA
#> 11                           -1                              -2
#> 12                            0                              -2
#> 13                           -1                              -2
#> 14                            0                               0
#> 15                           -4                              -5
#> 16                           -3                              -3
#> 17                            0                              -3
#> 18                           -3                              -3
#> 19                            0                               2
#> 20                           -5                              -5
#> 21                           -2                              -3
#> 22                            0                              -4
#> 23                           -2                              -4
#> 24                            0                               0
#> 25                           -4                              -5
#> 26                           -3                              -5
#> 27                           -1                              -4
#> 28                           -3                              -4
#> 29                           -1                              -2
#> 30                           -5                              -5
#> 31                           NA                              NA
#> 32                           NA                              NA
#> 33                           NA                              NA
#> 34                           NA                              NA
#> 35                           NA                              NA
#> 36                           NA                              NA
#> 37                           NA                              NA
#> 38                           NA                              NA
#> 39                           NA                              NA
#> 40                           NA                              NA
#>    sens_original.sens_habitat sens_original.sens_general ac_original.ac_feeding
#> 1                          NA                         -3                     NA
#> 2                          NA                          0                     NA
#> 3                          NA                         -2                     NA
#> 4                          NA                         -3                     NA
#> 5                          NA                          0                     NA
#> 6                          NA                         -4                     NA
#> 7                          NA                          0                     NA
#> 8                          NA                         -3                     NA
#> 9                          NA                         -4                     NA
#> 10                         NA                          0                     NA
#> 11                         -3                         NA                      0
#> 12                         -2                         NA                      0
#> 13                         -3                         NA                      0
#> 14                          2                         NA                      0
#> 15                          0                         NA                      0
#> 16                         -4                         NA                     -1
#> 17                         -3                         NA                      0
#> 18                         -4                         NA                      0
#> 19                          3                         NA                      1
#> 20                          0                         NA                     -1
#> 21                         -4                         NA                      0
#> 22                         -4                         NA                      0
#> 23                         -4                         NA                      0
#> 24                         -2                         NA                      0
#> 25                         -5                         NA                     -1
#> 26                         -5                         NA                     -1
#> 27                         -5                         NA                     -1
#> 28                         -5                         NA                     -1
#> 29                         -3                         NA                      1
#> 30                         -5                         NA                     -1
#> 31                         NA                         -2                     NA
#> 32                         NA                          0                     NA
#> 33                         NA                          0                     NA
#> 34                         NA                          0                     NA
#> 35                         NA                         -3                     NA
#> 36                         NA                         -4                     NA
#> 37                         NA                         -2                     NA
#> 38                         NA                         -2                     NA
#> 39                         NA                         -2                     NA
#> 40                         NA                         -5                     NA
#>    ac_original.ac_behaviour ac_original.ac_reproduction ac_original.ac_habitat
#> 1                        NA                          NA                     NA
#> 2                        NA                          NA                     NA
#> 3                        NA                          NA                     NA
#> 4                        NA                          NA                     NA
#> 5                        NA                          NA                     NA
#> 6                        NA                          NA                     NA
#> 7                        NA                          NA                     NA
#> 8                        NA                          NA                     NA
#> 9                        NA                          NA                     NA
#> 10                       NA                          NA                     NA
#> 11                        0                           0                      0
#> 12                        0                           0                      0
#> 13                        0                           0                      0
#> 14                        1                           1                      1
#> 15                       -1                          -1                     -1
#> 16                       -1                          -1                     -1
#> 17                        0                           0                      0
#> 18                        0                           0                      0
#> 19                        1                           1                      1
#> 20                       -1                          -1                     -1
#> 21                       -1                          -1                     -1
#> 22                        0                           0                      0
#> 23                       -1                          -1                     -1
#> 24                        0                           0                      0
#> 25                       -1                          -1                     -1
#> 26                       -1                          -1                     -1
#> 27                       -1                          -1                     -1
#> 28                       -1                          -1                     -1
#> 29                        1                           1                      1
#> 30                       -1                          -1                     -1
#> 31                       NA                          NA                     NA
#> 32                       NA                          NA                     NA
#> 33                       NA                          NA                     NA
#> 34                       NA                          NA                     NA
#> 35                       NA                          NA                     NA
#> 36                       NA                          NA                     NA
#> 37                       NA                          NA                     NA
#> 38                       NA                          NA                     NA
#> 39                       NA                          NA                     NA
#> 40                       NA                          NA                     NA
#>    ac_original.ac_general
#> 1                       1
#> 2                       1
#> 3                       1
#> 4                       1
#> 5                       0
#> 6                       1
#> 7                       1
#> 8                       1
#> 9                       1
#> 10                      0
#> 11                     NA
#> 12                     NA
#> 13                     NA
#> 14                     NA
#> 15                     NA
#> 16                     NA
#> 17                     NA
#> 18                     NA
#> 19                     NA
#> 20                     NA
#> 21                     NA
#> 22                     NA
#> 23                     NA
#> 24                     NA
#> 25                     NA
#> 26                     NA
#> 27                     NA
#> 28                     NA
#> 29                     NA
#> 30                     NA
#> 31                      1
#> 32                      0
#> 33                      0
#> 34                      0
#> 35                     -1
#> 36                      0
#> 37                      0
#> 38                      0
#> 39                      0
#> 40                     -1


### Example for one indicator and three pressures to evaluate direct
#   effects where sensitivity is scored for four individual traits:
ind <- "herring"
press <- c("fishing", "temperature increase", "salinity decrease")

# Create scoring table using the template function:
sens_ac_tbl <- create_template_sensitivity(
  indicators = ind,
  pressures = press,
  type = "direct",                      # (default)
  n_sensitivity_traits = 4,
  adaptive_capacity = TRUE,             # (default)
  mode_adaptive_capacity = "general",   # (default)
  uncertainty = TRUE,                   # (default)
  mode_uncertainty = "general"          # (default)
)

# Rename trait columns:
trait_cols <- paste0("sens_",
  c("feeding", "behaviour", "reproduction", "habitat"))
names(sens_ac_tbl)[4:7] <- trait_cols
# Give trait-specific sensitivity scores:
sens_ac_tbl$sens_feeding <- c(0,0,0)
sens_ac_tbl$sens_behaviour <- c(-1,0,-4)
sens_ac_tbl$sens_reproduction <- c(-2,-2,-5)
sens_ac_tbl$sens_habitat <- c(-3,-2,0)

# Give general adaptive capacity and uncertainty scores:
sens_ac_tbl$ac_general <- c(0,0,-1)
sens_ac_tbl$uncertainty_sens <- c(1,1,1)
sens_ac_tbl$uncertainty_ac <- c(1,1,2)

sens_ac_tbl
#>   indicator             pressure   type sens_feeding sens_behaviour
#> 1   herring              fishing direct            0             -1
#> 2   herring temperature increase direct            0              0
#> 3   herring    salinity decrease direct            0             -4
#>   sens_reproduction sens_habitat ac_general uncertainty_sens uncertainty_ac
#> 1                -2           -3          0                1              1
#> 2                -2           -2          0                1              1
#> 3                -5            0         -1                1              2

# Calculate median sensitivity scores (adaptive capacities and
# uncertainties cannot be aggregated further):
calc_sensitivity(
  indicators = sens_ac_tbl$indicator,
  pressures = sens_ac_tbl$pressure,
  sensitivity_traits = sens_ac_tbl[, trait_cols],
  adaptive_capacities = sens_ac_tbl$ac_general,
  uncertainty_sens  = sens_ac_tbl$uncertainty_sens,
  uncertainty_ac = sens_ac_tbl$uncertainty_ac,
  method = "median"
)
#>   indicator             pressure   type pathway sensitivity adaptive_capacity
#> 1   herring              fishing direct  expert        -1.5                 0
#> 2   herring temperature increase direct  expert        -1.0                 0
#> 3   herring    salinity decrease direct  expert        -2.0                -1
#>   uncertainty_sens uncertainty_ac sens_original.sens_feeding
#> 1                1              1                          0
#> 2                1              1                          0
#> 3                1              2                          0
#>   sens_original.sens_behaviour sens_original.sens_reproduction
#> 1                           -1                              -2
#> 2                            0                              -2
#> 3                           -4                              -5
#>   sens_original.sens_habitat ac_original.ac_general
#> 1                         -3                      0
#> 2                         -2                      0
#> 3                          0                     -1
```
