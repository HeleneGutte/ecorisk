# Baltic Sea pressure time series

Time series of eight environmental and anthropogenic pressures
potentially affecting the zooplankton mean size or cod spawning stock
biomass in the Eastern Baltic Sea. The time series cover the period
1984–2016 (data altered from original time series). This dataset serves
as a demo input in the
[`model_exposure`](https://helenegutte.github.io/ecorisk/reference/model_exposure.md)
and
[`model_sensitivity`](https://helenegutte.github.io/ecorisk/reference/model_sensitivity.md)
functions.

## Usage

``` r
pressure_ts_baltic
```

## Format

A data frame with 33 observations and 9 variables.

- year:

  Time variable.

- nitrogen:

  Mean total nitrogen input into the Baltic Sea per year.

- phosphorous:

  Mean total phosphorus input into the Baltic Sea per year.

- surf_temp_sum:

  Mean sea surface temperature in summer in the Baltic Sea per year (in
  °C).

- bot_temp_ann:

  Mean sea bottom temperature in the Baltic Sea per year (in °C).

- surf_sal_sum:

  Mean sea surface salinity in summer in the Baltic Sea per year.

- bot_sal_ann:

  Mean sea bottom salinity in the Baltic Sea per year.

- bot_oxy_ann:

  Mean bottom oxygen concentration in the Baltic Sea per year (in
  mg/m^3).

- fishing_cod:

  Mean eastern Baltic cod fishing pressure per year.
