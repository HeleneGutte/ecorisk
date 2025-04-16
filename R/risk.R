#' Calculate Risk Scores Using Expert-Based or Model-Derived Vulnerability
#' and Status Scores
#'
#' The \code{risk} function calculates risk scores using the output from of the
#' \code{\link{status}} and the \code{\link{vulnerability}} functions.
#' For each state indicator-pressure combination the function adds the status
#' score to the vulnerability score to derive the risk score.
#'
#' @param vulnerability_results A data frame with the output from the
#'        \code{vulnerability} function.
#' @param status_results A data frame with status scores for each state indicator.
#'        The first column MUST contain the indicator names.
#'        The second and third column have to be named status and score.
#'
#' @details Final risk scores are in a range from -10 (severe risk for the state
#' indicator due to the assessed pressure) to +10 (good opportunities for the
#' state indicator due to the assessed pressure). The risk scores are specific for
#' each combination of state indicator and pressure and do NOT take into account
#' cumulative effects. The risk scores can be aggregated in an additive manner
#' with the \code{\link{aggregate_risk}} function.
#'
#' @return a data frame containing the exposure, sensitivity, adaptive capacity,
#' vulnerability, and risk scores as well as their associated uncertainty for
#' each pressure - state indicator - type combination.
#'
#' @seealso \code{\link{vulnerability}}, \code{\link{status}}, \code{\link{aggregate_risk}}
#'
#' @export
#'
#' @examples
#' # Using demo output data from the vulnerability() and status() functions:
#' risk(
#'   vulnerability_results = ex_output_vulnerability_model,
#'   status_results = ex_output_status
#' )
#'
#' \donttest{
#'   ### Demo Expert-Based Pathway
#'   # - using the example scoring datasets 'ex_expert_exposure',
#'   #   'ex_expert_sensitivity' and 'ex_expert_status'
#'
#'   # Calculate (mean) exposure score:
#'   exp_expert <- calc_exposure(
#'     pressures = ex_expert_exposure$pressure,
#'     components = ex_expert_exposure[ ,2:5],
#'     uncertainties = ex_expert_exposure[ ,6:9],
#'     method = "mean" # default
#'   )
#'   # Calculate (mean) sensitivity (and adaptive capacity) score:
#'   sens_ac_expert <- calc_sensitivity(
#'     indicators = ex_expert_sensitivity$indicator,
#'     pressures = ex_expert_sensitivity$pressure,
#'     type = ex_expert_sensitivity$type,
#'     sensitivity_traits = ex_expert_sensitivity[ ,4:8],
#'     adaptive_capacities = ex_expert_sensitivity[ ,9:13],
#'     uncertainty_sens = ex_expert_sensitivity[ ,14:18],
#'     uncertainty_ac = ex_expert_sensitivity[ ,19:23],
#'     method = "mean" # default
#'   )
#'   # Calculate (mean) vulnerability score:
#'   vuln_expert <- vulnerability(
#'     exposure_results = exp_expert,
#'     sensitivity_results = sens_ac_expert,
#'     method_vulnerability = "mean", # default
#'     method_uncertainty = "mean" # default
#'   )
#'   # Calculate risk score:
#'   risk(
#'     vulnerability_results = vuln_expert,
#'     status_results = ex_expert_status
#'   )
#'
#'
#'   ### Demo Model-Based Pathway
#'   # - using the demo time series 'pressure_ts_baltic' and 'indicator_ts_baltic'
#'
#'   # Model exposure score:
#'   exp_model <- model_exposure(
#'     pressure_time_series = pressure_ts_baltic,
#'     base_years = c(start = 1984, end = 1994),
#'     current_years = c(start = 2010, end = 2016)
#'   )
#'
#'   # Model sensitivity score:
#'   sens_ac_model <- model_sensitivity(
#'     indicator_time_series = indicator_ts_baltic,
#'     pressure_time_series = pressure_ts_baltic,
#'     current_years = c(start = 2010, end = 2016)
#'   )
#'   # Add manually adaptive capacity scores (otherwise zero):
#'   sens_ac_model$adaptive_capacity <- c(rep(1, 8), rep(-1, 8))
#'
#'   # Calculate (mean) vulnerability score:
#'   vuln_model <- vulnerability(
#'     exposure_results = exp_model,
#'     sensitivity_results = sens_ac_model
#'   )
#'   # Calculate status score:
#'   status_model <- status(
#'     indicator_time_series = indicator_ts_baltic,
#'     base_years = c(start = 1984, end = 2010),
#'     current_years = c(start = 2011, end = 2016)
#'   )
#'   # Calculate risk score:
#'   risk(
#'     vulnerability_results = vuln_model,
#'     status_results = status_model
#'   )
#'
#' }

risk <- function(vulnerability_results, status_results) {

  # Shorten arguments
  vulnerability <- vulnerability_results
  status <- status_results

  # Test if all vulnerability indicators are in status table
  if (all(unique(vulnerability$indicator) %in% unique(status[, 1])) == FALSE) {
    v_names <- unique(vulnerability$indicator)
    status_names <- unique(status[, 1])
    no_match <- which(!(v_names %in% status_names))
    warning(paste(c("The following indicator names in the vulnerability table are not included in the status table:",
      v_names[c(no_match)],
      "\nPlease ensure that all indicator names of the vulnerability table are in the status table; also check for spelling",
      "errors. The affected indicators are not analysed further."),
      collapse = " "))
    skip_names <- v_names[no_match]
    vulnerability <- vulnerability[-which(vulnerability$indicator %in% skip_names), ]
  }
  r_score <- rep(0, nrow(vulnerability))
  status_out_helper <- rep(0, nrow(vulnerability))
  status_helper <- data.frame(indicators = status[, 1], score = status[, 3])


  # Calculate risk score
  for (i in 1:nrow(vulnerability)) {
    v_helper <- vulnerability$indicator[i]
    status_out_helper[i] <- status[v_helper == status[, 1], ]$status
    if (is.na(vulnerability$vulnerability[i]) == TRUE) {
      r_score[i] <- NA
    } else if (vulnerability$vulnerability[i] == 0) {
      r_score[i] <- 0
    } else {
      r_score[i] <- vulnerability$vulnerability[i] +
        status_helper[v_helper == status_helper$indicator, ]$score

      # Test that score is not > +10 / < -10
      if (r_score[i] > 10) {
        r_score[i] <- 10
      } else if (r_score[i] < -10) {
        r_score[i] <- -10
      } else {
        r_score[i] <- r_score[i]
      }
    }
  }


  output_risk <- data.frame(
    indicator = vulnerability$indicator,
    pressure = vulnerability$pressure,
    type = vulnerability$type,
    pathway = vulnerability$pathway,
    vulnerability = vulnerability$vulnerability,
    status = status_out_helper,
    risk = r_score,
    uncertainty = vulnerability$uncertainty
  )

  return(output_risk)

}

