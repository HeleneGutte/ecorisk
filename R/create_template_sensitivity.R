#' Create a Template for Expert-Based Sensitivity and Adaptive Capacity Scoring
#'
#' The function \code{crt_sensitivity} creates a template for semi-quantitative,
#' expert-based sensitivity and, optionally, adaptive capacity scoring. The
#' template allows for assessing sensitivity and adaptive capacity using either
#' a general score for each state indicator-pressure combination or a trait-based
#' approach using life history traits. The latter is particularly useful when
#' state indicators represent individual species and detailed biological information
#' is available.
#'
#' @param indicators A character vector specifying the names of the state indicators to
#'        assess.
#' @param pressures A character vector specifying the names of the pressures to
#'        assess.
#' @param type A character vector defining the type(s) of influence, such as
#'        \code{"direct"}, \code{"indirect"} or \code{"direct_indirect"}. The
#'        default is \code{"direct"}.
#' @param n_sensitivity_traits A positive integer specifying the number of traits
#'        used to assess sensitivity. The default is \code{1}, meaning that one
#'        general sensitivity score per state indicator-pressure-type combination is
#'        provided.
#' @param adaptive_capacity logical; should adaptive capacity be assessed? Default
#'        is \code{TRUE}.
#' @param mode_adaptive_capacity A character vector specifying whether adaptive
#'        capacity should be assessed for each trait individually (\code{"trait"})
#'        or as a general score (\code{"general"}, default). Note: This parameter
#'        is only relevant when \code{n_sensitivity_traits > 1}.
#' @param uncertainty logical; should uncertainty be assessed? Default is \code{TRUE}.
#'        Note: Uncertainty is only added for components included in the assessment.
#'        If adaptive capacity (\code{adaptive_capacity = FALSE}) is not included,
#'        no uncertainty scoring will be applied to it.
#' @param mode_uncertainty A character vector specifying whether uncertainty should
#'        be assessed for each trait individually (\code{"trait"}) or as a general
#'        score (\code{"general"}).
#'
#' @details
#' For each state indicator-pressure combination, different types of influence can be
#' assessed. The type of influence describes whether the pressure acts directly,
#' indirectly, or as a combination of both, which is important for identifying
#' impact pathways and potential management measures.
#'
#' Within the *ecorisk* framework, it is recommended to use negative scores
#' (-1 to -5) to indicate negative impacts (low to high severity) and positive
#' scores (1 to 5) for positive effects of a pressure on an indicator. If an
#' indicator is not sensitive to a pressure, the score should be 0. Adaptive
#' capacity is scored from -1 (no adaptive capacity) to 1 (high adaptive capacity).
#'
#' To improve the reliability of the scoring, uncertainty should also be assessed.
#' Uncertainty can be scored for each trait individually or as a general score
#' and should be rated on a scale from 1 to 3 (low to high uncertainty).
#'
#' The returned data frame can be exported as a CSV or Excel file. Column names
#' can be modified as needed. The completed file can be analyzed using the
#' \code{\link{calc_sensitivity}} function.
#'
#' Depending on the settings, the data frame will include:
#' \itemize{
#'   \item A single "sensitivity" column, if using general scoring.
#'   \item Multiple trait-specific sensitivity columns (e.g., "sensitivity_trait_1",
#'         "sensitivity_trait_2", etc.), which can be renamed as necessary.
#' }
#'
#' Within this data frame, trait-based and general scoring can be mixed. It is
#' therefore recommended to set \code{n_sensitivity_traits} to the maximum number
#' of traits to be assessed for any state indicator. The \code{\link{calc_sensitivity}}
#' function automatically distinguishes between general and trait-based scoring.
#'
#' @return
#' A data frame where each row represents a state indicator-pressure-type combination,
#' containing the specified sensitivity traits, adaptive capacity, and uncertainty
#' columns (if selected).
#' \itemize{
#'   \item If adaptive capacity and uncertainty are assessed, the data frame
#'         includes either one general column or one column per trait, depending
#'         on the settings.
#'   \item If using trait-based scoring, the data frame includes trait-specific
#'         sensitivity, adaptive capacity, and uncertainty columns, which can be
#'         renamed as needed.
#' }
#'
#' @seealso \code{\link{create_template_exposure}}, \code{\link{calc_exposure}},
#'          \code{\link{calc_sensitivity}}
#'
#' @export
#'
#' @examples
#' ### Create a table for two state indicators and two pressures to evaluate direct
#' #   effects (default). Return a general sensitivity and adaptive capacity
#' #   column as well as uncertainty columns for both components:
#' ind <- c("seabirds", "seals")
#' press <- c("plastic pollution", "temperature increase")
#'
#' sens_ac_tbl <- create_template_sensitivity(
#'   indicators = ind,
#'   pressures = press
#' )
#' # --> Export table and re-import after completion or fill in directly in R.
#'
#' # Assign sensitivity scores from -5 (strong negative response to pressure)
#' # to +5 (strong positive response) (0 = no sensitivity):
#' sens_ac_tbl$sens_general <- c(-5,3,-4,4)
#'
#' # Assign adaptive capacity scores from -1 (none) to +1 (good adaptive capacity):
#' sens_ac_tbl$ac_general <- c(-1,1,-1,1)
#'
#' # Assign uncertainty scores from 1 (low) to 3 (high uncertainty):
#' sens_ac_tbl$uncertainty_sens <- c(1,2,1,1)
#' sens_ac_tbl$uncertainty_ac <- c(3,2,3,2)
#'
#'
#' ### Create a table for four indicators and three pressures to evaluate both direct
#' #   and indirect effects. Return columns for five trait-specific sensitivities
#' #   and their respective uncertainties, but no adaptive capacity:
#' ind <- c("cod", "herring", "seabirds", "seals")
#' press <- c("fishing", "temperature increase", "salinity decrease")

#' sens_ac_tbl <- create_template_sensitivity(
#'   indicators = ind,
#'   pressures = press,
#'   type = c("direct", "direct_indirect"),
#'   n_sensitivity_traits = 5,
#'   adaptive_capacity = FALSE,
#'   uncertainty = TRUE,
#'   mode_uncertainty = "trait"
#' )
#' sens_ac_tbl
#' # --> You might want to rename the generic trait columns with specific traits.
#' # --> Export table as e.g. CSV-file and re-import again after completion.
#'
#' ### Create a mixed table for two indicators and two pressures, where for one
#' #   indicator sensitivity is scored overall and for one sensitivity is scored
#' #   by individual traits:
#' ind <- c("phytoplankton", "herring")
#' press <- c("temperature", "salinity")
#'
#' sens_ac_tbl <- create_template_sensitivity(
#'   indicators = ind,
#'   pressures = press,
#'   n_sensitivity_traits = 4,
#'   adaptive_capacity = FALSE,
#'   uncertainty = TRUE,
#'   mode_uncertainty = "general"
#' )
#' # Rename trait columns:
#' names(sens_ac_tbl)[4:7] <- paste0("sens_",
#'   c("feeding", "behaviour", "reproduction", "general"))
#'
#' # Give overall sensitivity score for phytoplankton
#' # (keep NAs for herring):
#' sens_ac_tbl$sens_general[1:2] <- c(-3,0)
#'
#' # Give trait-specific sensitivity scores for herring
#' # (keep NAs for phytoplankton):
#' sens_ac_tbl$sens_feeding[3:4] <- c(0,0)
#' sens_ac_tbl$sens_behaviour[3:4] <- c(-1,0)
#' sens_ac_tbl$sens_reproduction[3:4] <- c(-2,-2)
#'
#' # Give overall uncertainty score:
#' sens_ac_tbl$uncertainty_sens <- c(1,2,1,1)

create_template_sensitivity <- function(indicators, pressures,
  type = "direct", n_sensitivity_traits = 1,
  adaptive_capacity = TRUE, mode_adaptive_capacity = "general",
  uncertainty = TRUE, mode_uncertainty = "general") {

  tbl <- data.frame(
    indicator = rep(indicators, each = (length(pressures) * length(type))),
    pressure = rep(pressures, times = (length(indicators) * length(type))),
    type = rep(type, each = length(pressures))
  )

  # Add sensitivity
  if (!is.numeric(n_sensitivity_traits) | n_sensitivity_traits < 1 ) {
    stop("The number of sensitivity traits must be at least 1.")
  }
  if (n_sensitivity_traits == 1) {
    tbl$sens_general <- rep(NA, nrow(tbl))
  } else {
    for (i in 1:n_sensitivity_traits) {
      tbl[, i + 3] <- rep(NA, nrow(tbl))
      names(tbl)[i + 3] <- paste0("sens_trait_", i)
    }
  }

  # Add adaptive capacity
  if (adaptive_capacity == TRUE) {
    if (n_sensitivity_traits == 1) {
      # If only 1 sensitivity component selected, only 1 adaptive capacity column
      # can be returned (independent of mode_adaptive_capacity setting)
      tbl$ac_general <- rep(NA, nrow(tbl))
    } else {
      if (mode_adaptive_capacity == "general") {
        tbl$ac_general <- rep(NA, nrow(tbl))
      } else if (mode_adaptive_capacity == "trait") {
        nr_col <- ncol(tbl)
        for (i in 1:n_sensitivity_traits) {
          tbl[, i + nr_col] <- rep(NA, nrow(tbl))
          names(tbl)[i + nr_col] <- paste0("ac_trait_", i)
        }
      } else {
        stop("Mode of adaptive capacity not recognized. Argument 'mode_adaptive_capacity' must either be 'general' or 'trait'.")
      }
    }
  }

  # Add uncertainty
  if (uncertainty == TRUE) {
    if (n_sensitivity_traits == 1) {
      # If only 1 sensitivity component selected, only 1 uncertainty column
      # for sensitivity (and 1 for adaptive capacity) can be returned
      tbl$uncertainty_sens <- rep(NA, nrow(tbl))
      if (adaptive_capacity == TRUE) {
        tbl$uncertainty_ac <- rep(NA, nrow(tbl))
      }
    } else {
      if (mode_uncertainty == "general") {
        tbl$uncertainty_sens <- rep(NA, nrow(tbl))
        if (adaptive_capacity == TRUE) {
          tbl$uncertainty_ac <- rep(NA, nrow(tbl))
        }
      } else if (mode_uncertainty == "trait") {
        nr_col <- ncol(tbl)
        for (i in 1:n_sensitivity_traits) {
          tbl[, i + nr_col] <- rep(NA, nrow(tbl))
          names(tbl)[i + nr_col] <- paste0("uncertainty_sens_trait_", i)
        }
        if (adaptive_capacity == TRUE) {
          # Trait-specific uncertainty columns can only be returned if adaptive
          # capacity is set to 'trait', otherwise return general column
          if (mode_adaptive_capacity == "trait") {
            nr_col <- ncol(tbl)
            for (i in 1:n_sensitivity_traits) {
              tbl[, i + nr_col] <- rep(NA, nrow(tbl))
              names(tbl)[i + nr_col] <- paste0("uncertainty_ac_trait_", i)
            }
          } else {
            tbl$uncertainty_ac <- rep(NA, nrow(tbl))
          }
        }
      } else {
        stop("Mode of uncertainty not recognized. Argument 'mode_uncertainty' must either be 'general' or 'component'.")
      }
    }
  }

  return(tbl)
}
