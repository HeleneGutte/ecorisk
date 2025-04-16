#' Calculate Overall Sensitivity and Adaptive Capacity Scores from Trait-Specific Expert Ratings
#'
#' The function \code{calc_sensitivity} calculates aggregated sensitivity and adaptive
#' capacity scores. Additionally it prepares the scoring data for the further
#' usage in the \code{vulnerability} function. The scores for sensitivity and
#' adaptive capacity can be trait-based or general (one score per state indicator
#' and pressure combination).
#'
#' @param indicators A character vector or column of a data frame containing the
#'        names of the state indicators.
#' @param pressures A character vector or column of a data frame containing
#'        the names of the pressures.
#' @param type A character vector or column of a data frame specifying the effect type.
#'        Effects could be direct or indirect or a combination of both. Default is
#'        \code{direct}.
#' @param sensitivity_traits A data frame containing the numeric sensitivity
#'        values per species trait or as a general value.
#' @param adaptive_capacities A data frame (or single vector) containing the
#'        numeric values for the adaptive capacity. Either per trait or as one
#'        general value. When values are given per trait, they have to be in the
#'        same order as the values for the sensitivity. The default is \code{NULL}.
#' @param uncertainty_sens  A data frame (or a vector) containing the numeric uncertainty
#'        values associated with the sensitivity; default is \code{NULL}.
#' @param uncertainty_ac A data frame (or a vector) containing the numeric uncertainty
#'        values associated with the adaptive capacity; default is \code{NULL}.
#' @param method A character string specifying the method of the aggregation of the traits.
#'        Available are \code{median}, \code{minimum}, \code{maximum} and \code{mean},
#'        the \code{mean} is default.
#'
#' @details The function calculates per state indicator and
#' pressure combination one aggregated sensitivity and adaptive capacity score,
#' the aggregation method can be determined with the parameter \code{method}.
#' The assessment of adaptive capacity is optionally, if no scores for adaptive
#' capacity are provided the function calculates an aggregated sensitivity score
#' and prepares only the sensitivity scores for the \code{\link{vulnerability}} function.
#' Guidance for the scoring process can be found here: \code{\link{create_template_sensitivity}}
#' or in the vignette or in Gutte et al., 2025.
#' Using exposure and sensitivity scorings, vulnerability is calculated.
#'
#' @return a data frame containing the indicator, pressure and effect type,
#'  the aggregated sensitivity and adaptive capacity as well as their associated uncertainty scores.
#'  Additionally, the trait specific sensitivity and adaptive capacity scores are stored and used later
#'  as input for the \code{vulnerability} function.
#'
#' @seealso \code{\link{create_template_exposure}}, \code{\link{create_template_sensitivity}},
#'          \code{\link{calc_exposure}}, \code{\link{vulnerability}}
#'
#' @export
#'
#' @examples
#' ### Example using demo data with four indicators and five pressures with
#' #   scores for direct as well as combined direct-indirect effects based on
#' #   the template function create_template_sensitivity(). For two
#' #   indicators, sensitivity, adaptive capacity, and their uncertainties are
#' #   provided as general scores, while for the other two, they are based on
#' #   individual traits.
#' ex_expert_sensitivity
#'
#' # Calculate only mean sensitivity scores:
#' calc_sensitivity(
#'   indicators = ex_expert_sensitivity$indicator,
#'   pressures = ex_expert_sensitivity$pressure,
#'   sensitivity_traits = ex_expert_sensitivity[ ,4:8],
#'   adaptive_capacities = NULL,   # (default)
#'   uncertainty_sens  = NULL,     # (default)
#'   uncertainty_ac = NULL,        # (default)
#'   method = "mean"               # (default)
#'  )
#'
#' # Calculate mean scores for sensitivity, adaptive capacity and
#' # associated uncertainties:
#' calc_sensitivity(
#'   indicators = ex_expert_sensitivity$indicator,
#'   pressures = ex_expert_sensitivity$pressure,
#'   type = ex_expert_sensitivity$type,
#'   sensitivity_traits = ex_expert_sensitivity[ ,4:8],
#'   adaptive_capacities = ex_expert_sensitivity[ ,9:13],
#'   uncertainty_sens  = ex_expert_sensitivity[ ,14:18],
#'   uncertainty_ac = ex_expert_sensitivity[ ,19:23]
#'  )
#'
#'
#' ### Example for one indicator and three pressures to evaluate direct
#' #   effects where sensitivity is scored for four individual traits:
#' ind <- "herring"
#' press <- c("fishing", "temperature increase", "salinity decrease")
#'
#' # Create scoring table using the template function:
#' sens_ac_tbl <- create_template_sensitivity(
#'   indicators = ind,
#'   pressures = press,
#'   type = "direct",                      # (default)
#'   n_sensitivity_traits = 4,
#'   adaptive_capacity = TRUE,             # (default)
#'   mode_adaptive_capacity = "general",   # (default)
#'   uncertainty = TRUE,                   # (default)
#'   mode_uncertainty = "general"          # (default)
#' )
#'
#' # Rename trait columns:
#' trait_cols <- paste0("sens_",
#'   c("feeding", "behaviour", "reproduction", "habitat"))
#' names(sens_ac_tbl)[4:7] <- trait_cols

#' # Give trait-specific sensitivity scores:
#' sens_ac_tbl$sens_feeding <- c(0,0,0)
#' sens_ac_tbl$sens_behaviour <- c(-1,0,-4)
#' sens_ac_tbl$sens_reproduction <- c(-2,-2,-5)
#' sens_ac_tbl$sens_habitat <- c(-3,-2,0)
#'
#' # Give general adaptive capacity and uncertainty scores:
#' sens_ac_tbl$ac_general <- c(0,0,-1)
#' sens_ac_tbl$uncertainty_sens <- c(1,1,1)
#' sens_ac_tbl$uncertainty_ac <- c(1,1,2)
#'
#' sens_ac_tbl
#'
#' # Calculate median sensitivity scores (adaptive capacities and
#' # uncertainties cannot be aggregated further):
#' calc_sensitivity(
#'   indicators = sens_ac_tbl$indicator,
#'   pressures = sens_ac_tbl$pressure,
#'   sensitivity_traits = sens_ac_tbl[, trait_cols],
#'   adaptive_capacities = sens_ac_tbl$ac_general,
#'   uncertainty_sens  = sens_ac_tbl$uncertainty_sens,
#'   uncertainty_ac = sens_ac_tbl$uncertainty_ac,
#'   method = "median"
#' )

calc_sensitivity <- function(indicators, pressures, type = "direct",
  sensitivity_traits, adaptive_capacities = NULL, uncertainty_sens = NULL,
  uncertainty_ac = NULL, method = "mean") {

  check_dim_sens <- as.data.frame(sensitivity_traits)
  if (!is.null(adaptive_capacities)) check_dim_ac <- as.data.frame(adaptive_capacities)
  if (!is.null(uncertainty_sens)) check_dim_unc_sens <- as.data.frame(uncertainty_sens)
  if (!is.null(uncertainty_ac)) check_dim_unc_ac <- as.data.frame(uncertainty_ac)


  if (!is.character(indicators)) {
    stop("The 'indicators' argument must be a character vector containing indicator names.")
  }
  if (!is.character(pressures)) {
    stop("The 'pressures' argument must be a character vector containing pressure names.")
  }
  if (length(indicators) != length(pressures)) {
    stop("'indicators' and 'pressures' must be vectors of the same length (see example output from create_template_sensitivity()).")
  }

  if (all(sapply(sensitivity_traits, is.numeric)) == FALSE) {
    stop("The 'sensitivity_traits' argument must be either a numeric vector or a data frame where all columns are numeric.")
  }
  if (!is.null(adaptive_capacities)) {
    if (all(sapply(adaptive_capacities, is.numeric)) == FALSE) {
      stop("The 'adaptive_capacities' argument must be NULL, a numeric vector, or a data frame where all columns are numeric.")
    }
    if (ncol(check_dim_ac) > 1 & (ncol(check_dim_ac) != ncol(check_dim_sens))) {
      stop("If 'adaptive_capacities' is a data frame, it must contain either a single column or the same number of columns as 'sensitivity_traits'.")
    }
  }
  if (!is.null(uncertainty_sens)) {
    if (all(sapply(uncertainty_sens, is.numeric)) == FALSE) {
      stop("The 'uncertainty_sens' argument must be NULL, a numeric vector, or a data frame where all columns are numeric.")
    }
    if (ncol(check_dim_unc_sens) > 1 & (ncol(check_dim_unc_sens) != ncol(check_dim_sens))) {
      stop("If 'uncertainty_sens' is a data frame, it must contain either a single column or the same number of columns as 'sensitivity_traits'.")
    }
  }
  if (!is.null(uncertainty_ac)) {
    if (all(sapply(uncertainty_ac, is.numeric)) == FALSE) {
      stop("The 'uncertainty_ac' argument must be NULL, a numeric vector, or a data frame where all columns are numeric.")
    }
    if (ncol(check_dim_unc_ac) > 1 & (ncol(check_dim_unc_ac) != ncol(check_dim_ac))) {
      stop("If 'uncertainty_ac' is a data frame, it must contain either a single column or the same number of columns as 'adaptive_capacities'.")
    }
  }
  if (!(method %in% c("mean", "median", "maximum", "minimum"))) {
    warning("The specified method is not recognized. Defaulting to 'mean'.")
  }


  sens_aggr <- rep(0, length(indicators))
  ac_aggr <- rep(0, length(indicators))
  unc_sens_aggr <- rep(0, length(indicators))
  unc_ac_aggr <- rep(0, length(indicators))
  pathway <- rep("expert", length(indicators))

  # Abbreviate arguments
  traits <- sensitivity_traits
  ac <- adaptive_capacities


  # In case of a general sensitivity score, duplicate it to ensure
  # functioning of the for-loop
  if (is.null(ncol(traits)) == TRUE) {
    traits <- as.data.frame(replicate(2, traits))
    general_sens_score <- TRUE
  } else {
    general_sens_score <- FALSE
  }

  for (i in 1:length(indicators)) {

    if (is.null(ncol(ac)) == TRUE) {
      if (is.null(ac) == TRUE) {
        # if ac score is not provided by the user, set ac to 0, that
        # way it has no effect on the resulting sensitivity score
        n_traits <- ncol(traits)
        ac <- as.data.frame(replicate(n_traits, rep(0, length(indicators))))
      } else {
        # in case of a general ac score, duplicate the ac score until
        # it matches number of traits
        n_traits <- ncol(traits)
        ac <- as.data.frame(replicate(n_traits, ac))
      }
    }

    # Aggregate sensitivity and adaptive capacity
    if (method == "median") {
      sens_aggr[i] <- stats::median(as.numeric(traits[i, ]), na.rm = TRUE)
      ac_aggr[i] <- stats::median(as.numeric(ac[i, ]), na.rm = TRUE)
    }
    if (method == "maximum") {
      sens_aggr[i] <- max(as.numeric(traits[i, ]), na.rm = TRUE)
      ac_aggr[i] <- max(as.numeric(ac[i, ]), na.rm = TRUE)
    }
    if (method == "minimum") {
      sens_aggr[i] <- min(as.numeric(traits[i, ]), na.rm = TRUE)
      ac_aggr[i] <- min(as.numeric(ac[i, ]), na.rm = TRUE)
    }
    if (!(method %in% c("median", "maximum", "minimum"))) {
      sens_aggr[i] <- mean(as.numeric(traits[i, ]), na.rm = TRUE)
      ac_aggr[i] <- mean(as.numeric(ac[i, ]), na.rm = TRUE)
    }


    # Uncertainty
    if (is.null(ncol(uncertainty_sens )) == TRUE) {
      if (is.null(uncertainty_sens ) == TRUE) {
        unc_sens_aggr <- rep(NA, length(indicators))
      } else {
        unc_sens_aggr <- uncertainty_sens
      }
    } else {
      if (method == "median") {
        unc_sens_aggr[i] <- stats::median(as.numeric(uncertainty_sens [i, ]),
          na.rm = TRUE)
      }
      if (method == "maximum") {
        unc_sens_aggr[i] <- max(as.numeric(uncertainty_sens [i, ]),
          na.rm = TRUE)
      }
      if (method == "minimum") {
        unc_sens_aggr[i] <- min(as.numeric(uncertainty_sens [i, ]),
          na.rm = TRUE)
      }
      if (!(method %in% c("median", "maximum", "minimum"))) {
        unc_sens_aggr[i] <- mean(as.numeric(uncertainty_sens [i, ]),
          na.rm = TRUE)
      }
    }

    if (is.null(ncol(uncertainty_ac)) == TRUE) {
      if (is.null(uncertainty_ac) == TRUE) {
        unc_ac_aggr <- rep(NA, length(indicators))
      } else {
        unc_ac_aggr <- uncertainty_ac
      }

    } else {
      if (method == "median") {
        unc_ac_aggr[i] <- stats::median(as.numeric(uncertainty_ac[i, ]),
          na.rm = TRUE)
      }
      if (method == "maximum") {
        unc_ac_aggr[i] <- max(as.numeric(uncertainty_ac[i, ]),
          na.rm = TRUE)
      }
      if (method == "minimum") {
        unc_ac_aggr[i] <- min(as.numeric(uncertainty_ac[i, ]),
          na.rm = TRUE)
      }
      if (!(method %in% c("median", "maximum", "minimum"))) {
        unc_ac_aggr[i] <- mean(as.numeric(uncertainty_ac[i, ]),
          na.rm = TRUE)
      }
    }

  }

  if (general_sens_score == TRUE) {
    names(traits)[names(traits) == "V1"] <- "sens_original.sens_general"
    names(ac)[names(ac) == "V1"] <- "ac_original.ac_general"
    output_calc_sens <- data.frame(
      indicator = indicators,
      pressure = pressures,
      type = type,
      pathway = pathway,
      sensitivity = sens_aggr,
      adaptive_capacity = ac_aggr,
      uncertainty_sens  = unc_sens_aggr,
      uncertainty_ac = unc_ac_aggr,
      sens_original  = traits[1],
      ac_original = ac[1]
    )
  } else {
    # If adaptive_capacities == NULL, return NA as original input, else
    # return only first column and rename
    if (is.null(adaptive_capacities)) {
      ac <- data.frame(ac_original.ac_general = rep(NA, length(indicators)))
    } else {
      if (ncol(check_dim_ac) == 1) {
        names(ac)[names(ac) == "V1"] <- "ac_original.ac_general"
        ac <- ac[1]
      }
    }
    output_calc_sens <- data.frame(
      indicator = indicators,
      pressure = pressures,
      type = type,
      pathway = pathway,
      sensitivity = sens_aggr,
      adaptive_capacity = ac_aggr,
      uncertainty_sens  = unc_sens_aggr,
      uncertainty_ac = unc_ac_aggr,
      sens_original  = traits,
      ac_original = ac
    )
  }

  return(output_calc_sens)
}

