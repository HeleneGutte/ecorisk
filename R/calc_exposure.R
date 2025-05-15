#' Calculate Overall Exposure Scores from Component-Specific Expert Ratings
#'
#' Calculate exposure scores from individual exposure components. Additionally
#' calculate associated uncertainty scores.
#'
#' @param pressures A character vector or column of a data frame containing the
#'        names of the pressures.
#' @param components A numeric vector or data frame containing the numeric values
#'        per exposure component for each pressure. Has to be in the same order
#'        as the pressure vector or come from the same data frame.
#' @param probabilities Optionally: a numeric vector containing the probabilities of
#'        each pressure (values have to be between 0 and 1); default is \code{NULL}.
#'        Has to be in the same order as the pressure vector.
#' @param uncertainty a numeric vector or  data frame containing the associated
#'        uncertainty per component; default is \code{NULL}.
#'        Has to be in the same order as the pressure vector.
#' @param method a character string specifying the aggregation method. Available are
#'        mean (default), median, maximum, and minimum.
#'
#' @return
#' Returns a data frame containing the pressure name, the aggregated exposure score and
#' associated uncertainty scores. The results serve as input to the \code{\link{vulnerability}}
#' function. In case no uncertainty values were provided, \code{NA}s will be returned as uncertainty
#' scores.
#'
#' @details  Often exposure components include the magnitude of change compared
#' to baseline conditions, the frequency of this change, a future trend and a spatial scale.
#' These components are scored within the ecorisk framework for each pressure by
#' experts on a scale from 1 (low impact) to 5 (high impact). To express their
#' uncertainty during the process, experts can score the associated uncertainty
#' generally for all components of one pressure or for each component individually.
#' The uncertainty is scored in the ecorisk framework on a scale from 1 (low uncertainty)
#' to 3 (high uncertainty).
#' Using exposure and sensitivity scorings vulnerability is calculated.
#' Guidance for the scoring process can be found here: \code{\link{create_template_exposure}}
#' or in the vignette or in Gutte et al., 2025.
#'
#' @seealso \code{\link{create_template_exposure}}, \code{\link{create_template_sensitivity}},
#'          \code{\link{calc_sensitivity}}, \code{\link{vulnerability}}
#'
#' @export
#'
#'
#' @examples
#' ### Example using demo data with five pressures, four components and their individual
#' # uncertainties (probabilities are assumed to be 1):
#' ex_expert_exposure
#'
#' calc_exposure(
#'   pressures = ex_expert_exposure$pressure,
#'   components = ex_expert_exposure[ ,2:5],
#'   uncertainty = ex_expert_exposure[ ,6:9]
#'  )
#'
#'
#' ### Example for two hazardous risks with only two components ('magnitude' and
#' #   'spatial'), one general uncertainty score, and associated probabilities:
#' hazard <- c("heat waves", "hurricanes")
#'
#' # Create scoring table using the template function:
#' exp_tbl <- create_template_exposure(
#'   pressures = hazard,
#'   n_components = 2,
#'   mode_uncertainty = "general",
#'   probability = TRUE
#' )
#'
#' names(exp_tbl)[2:3] <- c("magnitude", "spatial")
#' # Assign component-specific scores and probabilities:
#' exp_tbl$magnitude <- c(5,4)
#' exp_tbl$spatial <- c(5,3)
#' exp_tbl$uncertainty <- c(2,3)
#' exp_tbl$probability <- c(0.8,0.3)
#'
#' # Calculate exposure score:
#' calc_exposure(
#'   pressures = exp_tbl$pressure,
#'   components = exp_tbl[ ,c("magnitude", "spatial")],
#'   probabilities = exp_tbl$probability,
#'   uncertainty = exp_tbl$uncertainty
#'  )

calc_exposure <- function(pressures, components, probabilities = NULL,
  uncertainty = NULL, method = "mean") {

  if (!is.character(pressures)) {
    stop("The 'pressures' argument must be a character vector of pressure names.")
  }
  if (all(sapply(components, is.numeric)) == FALSE) {
    stop("The 'components' argument must be a numeric vector or data frame with numeric variables.")
  }
  if (!is.null(probabilities)) {
    if (is.null(dim(probabilities)) == FALSE) {
      stop("The 'probabilities' argument must be NULL or a numeric vector.")
    }
    if (!is.numeric(probabilities) | (min(probabilities) < 0 | max(probabilities) > 1)) {
      stop("The 'probabilities' argument must be NULL or a numeric vector with values between 0 and 1.")
    }
  }
  if (!is.null(uncertainty) & all(sapply(uncertainty, is.numeric)) == FALSE) {
    stop("The 'uncertainty' argument must be NULL or a numeric vector or data frame with numeric variables.")
  }
  if (!(method %in% c("mean", "median", "maximum", "minimum"))) {
    warning("method not recognised. Will use method = arithmetic mean instead.")
  }

  score <- rep(0, length(pressures))
  unc <- rep(0, length(pressures))

  # Aggregate exposure score without probabilities
  if (is.null(probabilities) == TRUE) {
    if (is.null(ncol(components))) {
      score <- components
    } else {
      if (method == "median") {
        for (i in 1:length(pressures)) {
          score[i] <- stats::median(as.numeric(components[i, ]),
            na.rm = TRUE)
        }
      }
      if (method == "maximum") {
        for (i in 1:length(pressures)) {
          score[i] <- max(as.numeric(components[i, ]), na.rm = TRUE)
        }
      }
      if (method == "minimum") {
        for (i in 1:length(pressures)) {
          score[i] <- min(as.numeric(components[i, ]), na.rm = TRUE)
        }
      }
      if (!(method %in% c("median", "maximum", "minimum"))) {
        for (i in 1:length(pressures)) {
          score[i] <- mean(as.numeric(components[i, ]), na.rm = TRUE)
        }
      }
    }
  } else { # Aggregate exposure score including probabilities
    if (is.null(ncol(components))) {
      score <- components * probabilities
    } else {
      if (method == "median") {
        for (i in 1:length(pressures)) {
          score[i] <- stats::median(as.numeric(components[i, ]),
            na.rm = TRUE) * probabilities[i]
        }
      }
      if (method == "maximum") {
        for (i in 1:length(pressures)) {
          score[i] <- max(as.numeric(components[i, ]), na.rm = TRUE) *
            probabilities[i]
        }
      }
      if (method == "minimum") {
        for (i in 1:length(pressures)) {
          score[i] <- min(as.numeric(components[i, ]), na.rm = TRUE) *
            probabilities[i]
        }
      }
      if (!(method %in% c("median", "maximum", "minimum"))) {
        for (i in 1:length(pressures)) {
          score[i] <- mean(as.numeric(components[i, ]), na.rm = TRUE) *
            probabilities[i]
        }
      }
    }
  }

  # uncertainty
  if (is.null(ncol(uncertainty)) == TRUE) {
    if (is.null(uncertainty) == TRUE) {
      unc <- rep(NA, length(pressures))
    } else {
      unc <- uncertainty
    }

  } else {
    if (method == "median") {
      for (i in 1:length(pressures)) {
        unc[i] <- stats::median(as.numeric(uncertainty[i, ]), na.rm = T)
      }
    }
    if (method == "maximum") {
      for (i in 1:length(pressures)) {
        unc[i] <- max(as.numeric(uncertainty[i, ]), na.rm = T)
      }
    }
    if (method == "minimum") {
      for (i in 1:length(pressures)) {
        unc[i] <- min(as.numeric(uncertainty[i, ]), na.rm = T)
      }
    }
    if (!(method %in% c("median", "maximum", "minimum"))) {
      for (i in 1:length(pressures)) {
        unc[i] <- mean(as.numeric(uncertainty[i, ]), na.rm = T)
      }
    }
  }

  output_calc_exp <- data.frame(
    pressure = pressures,
    exposure = score,
    uncertainty = unc
  )

  return(output_calc_exp)
}
