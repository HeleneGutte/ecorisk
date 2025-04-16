
########################### EXPERT-BASED PATHWAY ###############################

# Exposure ---------------------------------------------------------------------

#' Expert-based exposure scores for five pressures
#'
#' This Baltic Sea demo dataset contains expert-assigned scores for five
#' environmental and anthropogenic pressures, detailing individual exposure
#' components and their associated uncertainties. The dataset was initialized
#' using the template function \code{\link{create_template_exposure}}.
#'
#' Exposure scores range from 1 (low) to 5 (high), while uncertainties range
#' from 1 (low) to 3 (high). This dataset can be used as input for the function
#' \code{\link{calc_exposure}}.
#'
#' @format A data frame with 5 observations and 9 variables (values randomly assigned).
#' \describe{
#'   \item{pressure}{Environmental or anthropogenic pressure.}
#'   \item{magnitude}{Score for the magnitude pressure change.}
#'   \item{frequency}{Score for the frequency or duration of pressure effect.}
#'   \item{trend}{Score for the future trend of pressure.}
#'   \item{spatial}{Score for the spatial extent of pressure change.}
#'   \item{uncertainty_magnitude}{Uncertainty of magnitude score.}
#'   \item{uncertainty_frequency}{Uncertainty of frequency score.}
#'   \item{uncertainty_trend}{Uncertainty of trend score.}
#'   \item{uncertainty_spatial}{Uncertainty of spatial extent score.}
#' }
"ex_expert_exposure"


#' Example output from the \code{calc_exposure()} function
#'
#' This dataset provides an expert-based example output from the
#' \code{\link{calc_exposure}} function applied to \code{\link{ex_expert_exposure}}
#' demo data.
#'
#' @format A data frame with 5 observations and 3 variables.
#' \describe{
#'   \item{pressure}{Name of the assessed pressure.}
#'   \item{exposure}{Calculated combined score of all exposure components.}
#'   \item{uncertainty}{Calculated combined score of associated uncertainties.}
#' }
"ex_output_calc_exposure"


# Sensitivity and adaptive capacity --------------------------------------------

#' Expert-based sensitivity and adaptive capacity scores for four indicators and
#' five pressures
#'
#' This demo dataset includes sensitivity and adaptive capacity scores for four
#' Baltic Sea indicators and five pressures, initialized using
#' \code{\link{create_template_sensitivity}}. Depending on the indicator, a
#' general score or trait-specific scores were assigned. This dataset serves as
#' input for \code{\link{calc_sensitivity}}.
#'
#' @format A data frame with 40 observations and 23 variables.
#' \describe{
#'   \item{indicator}{Name of assessed indicator.}
#'   \item{pressure}{Name of assessed pressure.}
#'   \item{type}{Effect type (direct, indirect, or direct + indirect).}
#'   \item{sens_feeding}{Sensitivity score for the feeding trait (-5 to 5).}
#'   \item{sens_behaviour}{Sensitivity score for the behaviour trait (-5 to 5).}
#'   \item{sens_reproduction}{Sensitivity score for the reproduction trait (-5 to 5).}
#'   \item{sens_habitat}{Sensitivity score for the habitat trait (-5 to 5)}
#'   \item{sens_general}{General sensitivity score (where trait-based scoring is not possible, -5 to 5).}
#'   \item{ac_feeding}{Adaptive capacity score for the feeding trait (-1 to 1).}
#'   \item{ac_behaviour}{Adaptive capacity score for the behaviour trait (-1 to 1).}
#'   \item{ac_reproduction}{Adaptive capacity score for the reproduction trait (-1 to 1).}
#'   \item{ac_habitat}{Adaptive capacity score for the habitat trait (-1 to 1).}
#'   \item{ac_general}{General adaptive capacity score (where trait-based scoring is not possible, -1 to 1).}
#'   \item{uncertainty_sens_feeding}{Uncertainty of sensitivity score for feeding trait (1 to 3).}
#'   \item{uncertainty_sens_behaviour}{Uncertainty of sensitivity score for behaviour trait (1 to 3).}
#'   \item{uncertainty_sens_reproduction}{Uncertainty of sensitivity score for reproduction trait (1 to 3).}
#'   \item{uncertainty_sens_habitat}{Uncertainty of sensitivity score for habitat trait (1 to 3).}
#'   \item{uncertainty_sens_general}{Uncertainty of general sensitivity score (1 to 3).}
#'   \item{uncertainty_ac_feeding}{Uncertainty of adaptive capacity score for feeding trait (1 to 3).}
#'   \item{uncertainty_ac_behaviour}{Uncertainty of adaptive capacity score for behaviour trait (1 to 3).}
#'   \item{uncertainty_ac_reproduction}{Uncertainty of adaptive capacity score for reproduction trait (1 to 3).}
#'   \item{uncertainty_ac_habitat}{Uncertainty of adaptive capacity score for habitat trait (1 to 3).}
#'   \item{uncertainty_ac_general}{Uncertainty of general adaptive capacity score (1 to 3).}
#' }
"ex_expert_sensitivity"


#' Example output from the \code{calc_sensitivity()} function
#'
#' This is an expert-based example output from the \code{\link{calc_sensitivity}} function
#' applied to the \code{\link{ex_expert_sensitivity}} demo data.
#'
#' @format A data frame with 40 observations and 18 variables.
#' \describe{
#'   \item{indicator}{Name of the assessed indicator.}
#'   \item{pressure}{Name of the assessed pressure.}
#'   \item{type}{Effect type (direct, indirect, or direct + indirect).}
#'   \item{pathway}{Pathway with which sensitivity has been assessed (expert- or model-based).}
#'   \item{sensitivity}{Combined sensitivity score.}
#'   \item{adaptive_capacity}{Combined adaptive capacity score.}
#'   \item{uncertainty_sens}{Combined score of the associated sensitivity uncertainties.}
#'   \item{uncertainty_ac}{Combined score of the associated adaptive capacity uncertainties.}
#'   \item{sens_original.sens_feeding}{Original sensitivity score for the feeding trait.}
#'   \item{sens_original.sens_behaviour}{Original sensitivity score for the behaviour trait.}
#'   \item{sens_original.sens_reproduction}{Original sensitivity score for the reproduction trait.}
#'   \item{sens_original.sens_habitat}{Original sensitivity score for the habitat trait.}
#'   \item{sens_original.sens_general}{Original general sensitivity score.}
#'   \item{ac_original.ac_feeding}{Original adaptive capacity score for the feeding trait.}
#'   \item{ac_original.ac_behaviour}{Original adaptive capacity score for the behaviour trait.}
#'   \item{ac_original.ac_reproduction}{Original adaptive capacity score for the reproduction trait.}
#'   \item{ac_original.ac_habitat}{Original adaptive capacity score for the habitat trait.}
#'   \item{ac_original.ac_general}{Original general adaptive capacity score.}
#' }
"ex_output_calc_sensitivity"


# Vulnerability and status -----------------------------------------------------

#' Example output from the \code{vulnerability()} function based on expert scores
#'
#' This dataset contains an expert-based example output from the
#' \code{\link{vulnerability}} function applied to \code{\link{ex_output_calc_exposure}}
#' and \code{\link{ex_output_calc_sensitivity}} demo datasets.
#'
#' @format A data frame with 40 observations and 6 variables.
#' \describe{
#'   \item{indicator}{Name of the assessed indicator.}
#'   \item{pressure}{Name of the assessed pressure.}
#'   \item{type}{Effect type (direct, indirect, or direct + indirect).}
#'   \item{pathway}{Pathway used for the exposure and sensitivity assessment.}
#'   \item{vulnerability}{Vulnerability score.}
#'   \item{uncertainty}{Uncertainty associated with the vulnerability score.}
#' }
"ex_output_vulnerability_expert"


#' Expert-based status scores for four indicators
#'
#' This demo dataset contains the status scores of four Baltic Sea indicators based
#' on expert knowledge. The format is the same as the output table from the
#' \code{\link{status}} function that evaluates the status based on time series.
#'
#' @format A data frame with 4 observations and 3 variables.
#' \describe{
#'   \item{indicator}{Name of the assessed indicator.}
#'   \item{status}{Current status of each indicator, either 'good' or 'undesired'.}
#'   \item{score}{Score for each status (+1 or -1), will be combined with the
#'         vulnerability to derive risk.}
#' }
"ex_expert_status"


# Risk -------------------------------------------------------------------------

#' Example output from the \code{risk()} function based on expert scores
#'
#' This is an expert-based example output from the \code{\link{risk}} function
#' applied to the \code{\link{ex_output_vulnerability_expert}} and \code{\link{ex_expert_status}}
#' demo datasets.
#'
#' @format A data frame with 40 observations and 8 variables.
#' \describe{
#'   \item{indicator}{Name of the assessed indicator.}
#'   \item{pressure}{Name of the assessed pressure.}
#'   \item{type}{Effect type (direct, indirect, or direct + indirect).}
#'   \item{pathway}{Pathway used for the exposure and sensitivity assessment.}
#'   \item{vulnerability}{Vulnerability score.}
#'   \item{status}{Qualitative descriptor of the current status of the indicator.}
#'   \item{risk}{Risk score.}
#'   \item{uncertainty}{Uncertainty score, associated with the vulnerability component scoring.}
#' }
"ex_output_risk_expert"


#' Example output from the \code{aggregate_risk()} function based on expert scores
#'
#' This is an expert-based example output from the \code{\link{aggregate_risk}} function
#' applied to the \code{\link{ex_output_risk_expert}} demo data.
#'
#' @format A list of three data frames.
#' \describe{
#'   \item{multi_indicator_risk}{A data frame with 30 rows and 5 columns, containing the
#'        multi-indicator risk and uncertainty of each pressure per type and pathway.}
#'   \item{multi_pressure_risk}{A data frame with 24 rows and 5 columns, containing the
#'        multi-pressure risk and uncertainty on each indicator per type and pathway.}
#'   \item{ecosystem_risk}{A data frame with 6 rows and 4 columns, containing the
#'        aggregated ecosystem risk and uncertainty per type and pathway.}
#' }
"ex_output_aggregate_risk_expert"






############################ MODEL-BASED PATHWAY ###############################


# Baltic Sea demo time series --------------------------------------------------

#' Baltic Sea pressure time series
#'
#' Time series of eight environmental and anthropogenic pressures potentially affecting
#' the zooplankton mean size or cod spawning stock biomass in the Eastern Baltic Sea.
#' The time series cover the period 1984–2016 (data altered from original time series).
#' This dataset serves as a demo input in the \code{\link{model_exposure}} and
#' \code{\link{model_sensitivity}} functions.
#'
#' @format A data frame with 33 observations and 9 variables.
#' \describe{
#'   \item{year}{Time variable.}
#'   \item{nitrogen}{Mean total nitrogen input into the Baltic Sea per year.}
#'   \item{phosphorous}{Mean total phosphorus input into the Baltic Sea per year.}
#'   \item{surf_temp_sum}{Mean sea surface temperature in summer in the Baltic Sea per year (in °C).}
#'   \item{bot_temp_ann}{Mean sea bottom temperature in the Baltic Sea per year (in °C).}
#'   \item{surf_sal_sum}{Mean sea surface salinity in summer in the Baltic Sea per year.}
#'   \item{bot_sal_ann}{Mean sea bottom salinity in the Baltic Sea per year.}
#'   \item{bot_oxy_ann}{Mean bottom oxygen concentration in the Baltic Sea per year (in mg/m^3).}
#'   \item{fishing_cod}{Mean eastern Baltic cod fishing pressure per year.}
#' }
"pressure_ts_baltic"

#' Baltic Sea indicator time series
#'
#' Time series of two marine indicators covering the period 1984–2016 in the
#' Eastern Baltic Sea (data altered from original time series). This dataset
#' serves as a demo input in the \code{\link{model_sensitivity}} and
#' \code{\link{model_exposure}} functions.
#'
#' @format A data frame with 33 observations and 3 variables.
#' \describe{
#'   \item{year}{Time variable.}
#'   \item{zooplankton_mean_size}{Mean size of zooplankton (in wet weight micrograms).}
#'   \item{eastern_baltic_cod}{Mean spawning stock biomass of eastern Baltic cod (units unspecified).}
#' }
"indicator_ts_baltic"

#' North Sea pressure time series
#'
#' Time series of three environmental and anthropogenic pressures in the North
#' Sea covering the period 1970–2020 (data altered from original time series).
#' This dataset serves as internal test data.
#'
#' @format A data frame with 33 observations and 9 variables.
#' \describe{
#'   \item{year}{Time variable.}
#'   \item{bot_temp}{Mean sea bottom temperature in the Baltic Sea per year (in °C).}
#'   \item{bot_sal}{Mean sea bottom salinity in the Baltic Sea per year.}
#'   \item{fishing_cod}{Mean eastern Baltic cod fishing pressure per year.}
#' }
"pressure_ts_northsea"


# Modelled exposure and sensitivity --------------------------------------------

#' Example output from the \code{model_exposure()} function based on time series
#'
#' This dataset provides example output from the \code{\link{model_exposure}} function,
#' with component-specific exposure scores derived from the pressure time series
#' in \code{\link{pressure_ts_baltic}}. These scores are combined into an
#' overall exposure score, with associated uncertainties derived from two model types.
#'
#' @format A data frame with 8 observations and 10 variables.
#' \describe{
#'   \item{pressure}{Names of the assessed pressure.}
#'   \item{exposure}{Combined exposure score (1 to 5).}
#'   \item{uncertainty}{Uncertainty score from exposure modelling (1 to 3).}
#'   \item{comp_magnitude}{Score for the magnitude or degree of change (1 to 5).}
#'   \item{comp_frequency}{Score for the frequency or duration of change (1 to 5).}
#'   \item{comp_trend}{Score for the current trend of change (1 to 5).}
#'   \item{comp_direction}{Direction of the trend slope (increase or decrease).}
#'   \item{comp_spatial}{Score for spatial extent of the pressure (default: 3; user-defined: 1 to 5).}
#'   \item{uncertainty_arima}{Uncertainty score based on an ARIMA model.}
#'   \item{uncertainty_gam}{Uncertainty score based on a GAM model.}
#'   \item{mean_baseline}{Mean of the baseline conditions, used for
#'        magnitude scoring.}
#'   \item{mean_current}{Mean of the current conditions, used for magnitude
#'        and frequency scoring.}
#'   \item{standard_deviation_baseline}{Standard deviations of the baseline
#'        conditions. Used for scoring of magnitude and frequency.}
#'   \item{slope_linear_model}{Slope of the linear model used for scoring
#'        the future trend and to determine the direction.}
#'   \item{p_value_linear_model}{P-value of the linear model, used to
#'        score the future trend.}
#' }
"ex_output_model_exposure"

#' Example output from the \code{model_sensitivity()} function based on time series
#'
#' This dataset provides example output from the \code{\link{model_sensitivity}}
#' function, with sensitivity scores and associated uncertainties for each
#' indicator-pressure combination. Scores are based on time series data from
#' \code{\link{pressure_ts_baltic}} and \code{\link{indicator_ts_baltic}}.
#'
#' @format A data frame with 16 observations and 12 variables.
#' \describe{
#'   \item{indicator}{Names of the assessed indicator.}
#'   \item{pressure}{Names of the assessed pressure.}
#'   \item{type}{Type of effect (always direct + indirect for modelling pathway).}
#'   \item{pathway}{Pathway used to assess sensitivity.}
#'   \item{sensitivity}{Overall sensitivity score (-5 to 5).}
#'   \item{adaptive_capacity}{Adaptive capacity score (default is 0).}
#'   \item{uncertainty_sens}{Uncertainty score associated with sensitivity assessment (1 to 3).}
#'   \item{uncertainty_ac}{Uncertainty score associated with adaptive capacity (1 to 3).}
#'   \item{r_sq}{R-squared values from the GAM model, used for scoring.}
#'   \item{p_value}{P-values from the GAM model, determining statistical significance.}
#'   \item{edf}{Effective degrees of freedom from the GAM model, used to adjust scores based on non-linearity risk.}
#'   \item{uncertainty_gam}{Uncertainty score for sensitivity based on
#'        predicted values from a GAM.}
#'   \item{uncertainty_arima}{Uncertainty score for sensitivity based on
#'        predicted values from an ARIMA using the pressure variable as external
#'        predictor.}
#' }
"ex_output_model_sensitivity"


# Vulnerability and status -----------------------------------------------------

#' Example output from the \code{vulnerability()} function based on modelled scores
#'
#' This dataset provides example output from the \code{\link{vulnerability}} function,
#' applied to the \code{\link{ex_output_model_exposure}} and \code{\link{ex_output_model_sensitivity}}
#' demo datasets, following the modelling pathway.
#'
#' @format A data frame with 16 observations and 6 variables.
#' \describe{
#'   \item{indicator}{Names of the assessed indicator.}
#'   \item{pressure}{Names of the assessed pressure.}
#'   \item{type}{Type of effect (always direct + indirect for modelling pathway).}
#'   \item{pathway}{Pathway used for the exposure and sensitivity assessment.}
#'   \item{vulnerability}{Vulnerability score for each pressure-indicator-type combination.}
#'   \item{uncertainty}{Uncertainty associated with the vulnerability score.}
#' }
"ex_output_vulnerability_model"


#' Example output from the \code{status()} function
#'
#' This dataset provides example output from the \code{\link{status}} function,
#' applied to four Baltic Sea indicator time series provided in \code{\link{indicator_ts_baltic}}.
#'
#' @format A data frame with 2 rows and 3 variables.
#' \describe{
#'   \item{indicator}{Name of the assessed indicator.}
#'   \item{status}{Qualitative description of the current status compared to the threshold.}
#'   \item{score}{Score used in the risk function to calculate risk from vulnerability and status.}
#' }
"ex_output_status"


# Risk -------------------------------------------------------------------------

#' Example output from the \code{aggregate_risk()} function based on modelled scores
#'
#' This dataset provides example output from the \code{\link{risk}} function,
#' applied to the \code{\link{ex_output_vulnerability_model}} and \code{\link{ex_output_status}}
#' demo datasets, following the modelling pathway.
#'
#' @format A data frame with 16 observations and 8 variables.
#' \describe{
#'   \item{indicator}{Name of the assessed indicator.}
#'   \item{pressure}{Name of the assessed pressure.}
#'   \item{type}{Type of effect (always direct + indirect for modelling pathway).}
#'   \item{pathway}{Pathway used for the exposure and sensitivity assessment.}
#'   \item{vulnerability}{Vulnerability score.}
#'   \item{status}{Qualitative descriptor of the current status of the indicator.}
#'   \item{risk}{Risk score.}
#'   \item{uncertainty}{Uncertainty score associated with the vulnerability component
#'        scoring.}
#' }
"ex_output_risk_model"


#' Example output from the \code{aggregate_risk()} function based on modelled scores
#'
#' This dataset provides example output from the \code{\link{aggregate_risk}} function,
#' applied to the \code{\link{ex_output_risk_model}} demo data, following the modelling
#' pathway.
#'
#' @format A list of three data frames.
#' \describe{
#'   \item{multi_indicator_risk}{A data frame with 32 rows and 5 variables,
#'        containing the multi-indicator risk and uncertainty of each pressure per
#'        type and pathway.}
#'   \item{multi_pressure_risk}{A data frame with 8 rows and 5 variables, containing
#'        the multi-pressure risk and uncertainty on each indicator per type and
#'        pathway.}
#'   \item{ecosystem_risk}{A data frame with 4 rows and 4 variables, containing
#'        the aggregated ecosystem risk and uncertainty per type and pathway.}
#' }
"ex_output_aggregate_risk_model"






