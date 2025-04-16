#' Compute Status Scores from Time Series Data
#'
#' The `status` function assesses whether a state indicator is in a desired
#' or undesired status during the assessment time period. For this the function
#' compares the current conditions to the baseline conditions. The user specifies
#' whether the mean of the current conditions should be within or outside of a
#' specific deviation from the baseline mean.
#'
#' @param indicator_time_series A data frame with time series per state indicator.
#'        The first column MUST be the time column.
#' @param base_years A vector with two numerics, specifying the time period for the baseline.
#'        The first one \code{start} is the starting year for all state indicators
#'        and the second one \code{end} is the end of the baseline for all state
#'        indicators. The default is NULL. One can specify indicator specific
#'        baseline periods using the \code{base_years_by_ind} argument.
#'        If \code{base_years} and \code{base_years_by_ind} are NULL, then the
#'        first 5 years of the time series are used as baseline period.
#' @param base_years_by_ind A data frame, specifying the baseline years for each state
#'        indicator individually, by setting the starting year (second column)
#'        and the end year (third column). The first column must contain the
#'        names of the state indicators used in \code{indicator_time_series}.
#'        The default is NULL. If \code{base_years} and \code{base_years_by_ind}
#'        are NULL, then the first 5 years of the time series are used as baseline period.
#' @param current_years A vector with two numerics, specifying the time period for
#'        the assessment period. The first one \code{start} is the starting year
#'        for all state indicators and the second one \code{end} is the end of the
#'        assessment period for all state indicators. The default is NULL.
#'        One can specify indicator specific assessment periods using the
#'        \code{current_years_by_ind} argument.
#'        If \code{current_years} and \code{current_years_by_ind} are NULL, then the
#'        last 5 years of the time series are used as baseline period.
#' @param current_years_by_ind A data frame, specifying the assessment period years
#'        for each state indicator individually, by setting the starting year (second column)
#'        and the end year (third column). The first column must contain the
#'        names of the state indicators used in \code{indicator_time_series}.
#'        The default is NULL. If \code{current_years} and \code{current_years_by_ind}
#'        are NULL, then the last 5 years of the time series are used as assessment period.
#' @param range A vector specifying the allowed deviance from the baseline mean.
#'        Can be \code{sd}, \code{2sd}, \code{95percentile} or an integer between 1 and 99
#'        to evaluate the nth percentile. If the current mean should be compared to the
#'        baseline mean without any deviance, please set \code{range = mean_only},
#'        leave the default for the sign parameter and specify the condition
#'        parameter if necessary. Default is \code{sd}.
#' @param sign A character vector containing \code{+} or  \code{-}, specifying
#'        whether the upper or the lower part of the deviance should be analyzed.
#'        Default is \code{+}.
#' @param condition A character vector containing \code{<} or \code{>} specifying
#'        whether the current indicator should be above (>) or below (<) the preset
#'        threshold range to be in a desired status. The default is \code{>}.
#'
#' @details With \code{range}, \code{sign} and \code{condition} one defines good status
#' for the state indicators. By default the function evaluates whether the current
#' mean is above +1 standard deviation, if yes the status will be set to desired.
#' If the state should be within a range of ± standard deviation and not below that,
#' then the arguments \code{sign} and \code{condition} must be set to '-' and '>', this specifies
#' that the current mean must be higher than the mean of the baseline period - 1 standard
#' deviation to be considered as good status.
#'
#' @return a data frame containing the indicator name its status and the associated score,
#'  which will be added to the indicators vulnerability to derive the risk.
#'
#' @seealso \code{\link{model_exposure}}, \code{\link{model_sensitivity}},
#'          \code{\link{vulnerability}}, \code{\link{risk}}
#'
#' @export
#'
#' @examples
#'### Demo with the internal dataset 'indicator_ts_baltic'
#'
#' # Define a general baseline and current assessment period:
#' status(
#'  indicator_time_series = indicator_ts_baltic,
#'  base_years = c(start = 1984, end = 2010),
#'  current_years = c(start = 2011, end = 2016)
#')
#'
#' # Define indicator-specific baseline and current assessment periods:
#' status(
#'  indicator_time_series = indicator_ts_baltic,
#'  base_years_by_ind = data.frame(
#'    ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
#'    start = c(1984, 1990), end = c(2010, 2010)
#'  ),
#'  current_years_by_ind = data.frame(
#'    ind =c("zooplankton_mean_size", "eastern_baltic_cod"),
#'    start = c(2011, 2012), end = c(2016, 2016)
#'  )
#')

status <- function(indicator_time_series,
  base_years = NULL, base_years_by_ind = NULL,
  current_years = NULL, current_years_by_ind = NULL,
  range = "sd", sign = "+", condition = ">") {

### Data input validations
  # Time series
  if (any(class(indicator_time_series) != "data.frame")) {
    stop("The 'indicator_time_series' argument must be a data frame (not a
         matrix or tibble) in which all columns are numeric. The first column
         should represent the time component (e.g., years).")
  }
  if (all(sapply(indicator_time_series, is.numeric)) == FALSE) {
    stop("The 'indicator_time_series' argument must be a data frame where all
         columns are numeric. The first column should represent the time
         component (e.g., years).")
  }
  if (all(indicator_time_series[, 1] == sort(indicator_time_series[, 1])) == FALSE) {
    warning("The first column in 'indicator_time_series' appears to be out of
            order. Are you sure it represents the time component (e.g., years)?")
  }

  # Baseline period
  if (is.null(base_years) & is.null(base_years_by_ind)) {
    warning(paste0("You did not provide the first and last year of the baseline period. ",
      "You can either define this generally ('base_years = c(firstyear, lastyear)') or specify it ",
      "for each indicator using 'base_years_by_ind'. The default will ",
      "be the first 5 years of the time series. Refer to the R documentation for further guidance."))
    yrs <- indicator_time_series[ ,1]
    base_years <- c(min(yrs), min(yrs)+4)
  }

  if (!is.null(base_years)) {
    if (!is.null(dim(base_years))) {
      stop("'base_years' must be a numeric vector specifying the start and end
           year of the baseline period, which will be applied to all indicators.")
    } else {
      if (length(base_years) != 2) {
        stop("'base_years' must be a numeric vector specifying the start and end
             year of the baseline period, which will be applied to all indicators.")
      }
    }
    if (all(base_years %in% indicator_time_series[ ,1]) == FALSE ) {
      stop("The years defined in 'base_years' fall outside the available time
           series range.")
    }
    if (base_years[1] > base_years[2]) {
      stop("The start year cannot be later than the end year in 'base_years'.")
    }
  }

  if (!is.null(base_years_by_ind)) {
    by_format_txt <- paste0("The 'base_years_by_ind' argument must be a data frame with three columns: ",
      "'ind' (a character vector containing the indicator names) and two numeric columns, ",
      "'start' and 'end' (the first and last years of the baseline period). See the examples ",
      "section in the R documentation for further guidance.")
    if (any(class(base_years_by_ind) != "data.frame")) {
      stop(by_format_txt)
    } else {
      if (is.character(base_years_by_ind[, 1]) == FALSE |
          all(sapply(base_years_by_ind[, 2:3], is.numeric)) == FALSE |
          !identical(names(base_years_by_ind), c("ind", "start", "end")) ) {
        stop(by_format_txt)
      }
      if (nrow(base_years_by_ind) != ncol(indicator_time_series[ ,-1]) |
          all(names(indicator_time_series)[-1] %in% base_years_by_ind$ind) == FALSE ) {
        stop("The data frame provided in 'base_years_by_ind' must include all
             indicators listed in 'indicator_time_series'.")
      }
      if ( all(base_years_by_ind$start %in% indicator_time_series[ ,1]) == FALSE |
          all(base_years_by_ind$end %in% indicator_time_series[ ,1]) == FALSE ) {
        stop("Some of the years defined in 'base_years_by_ind' fall outside the
             available time series range.")
      }
      if (any(base_years_by_ind$start > base_years_by_ind$end)) {
        stop("The start year cannot be later than the end year in 'base_years_by_ind'.")
      }
    }
  }

  # Current (assessment) period
  if (is.null(current_years) & is.null(current_years_by_ind)) {
    warning(paste0("You did not provide the first and last year of the current assessment period. ",
      "You can either define this generally ('current_years = c(firstyear, lastyear)') or specify it ",
      "for each indicator using 'current_years_by_ind'. The default will ",
      "be the last 5 years of the time series. Refer to the R documentation for further guidance."))
    yrs <- indicator_time_series[ ,1]
    current_years <- c(max(yrs)-4, max(yrs))
  }

  if (!is.null(current_years)) {
    if (!is.null(dim(current_years))) {
      stop("'current_years' must be a numeric vector specifying the start and
           end year of the current assessment period, which will be applied to
           all indicators.")
    } else {
      if (length(current_years) != 2) {
        stop("'current_years' must be a numeric vector specifying the start and
             end year of the current assessment period, which will be applied
             to all indicators.")
      }
    }
    if (all(current_years %in% indicator_time_series[ ,1]) == FALSE ) {
      stop("The years defined in 'current_years' fall outside the available
           time series range.")
    }
    if (current_years[1] > current_years[2]) {
      stop("The start year cannot be later than the end year in 'current_years'.")
    }
  }

  if (!is.null(current_years_by_ind)) {
    cy_format_txt <- paste0("The 'current_years_by_ind' argument must be a data frame with three columns: ",
      "'ind' (a character vector containing the indicator names) and two numeric columns, ",
      "'start' and 'end' (the first and last years of the current assessment period). See the examples ",
      "section in the R documentation for further guidance.")
    if (any(class(current_years_by_ind) != "data.frame")) {
      stop(cy_format_txt)
    } else {
      if (is.character(current_years_by_ind[, 1]) == FALSE |
          all(sapply(current_years_by_ind[, 2:3], is.numeric)) == FALSE |
          !identical(names(current_years_by_ind), c("ind", "start", "end")) ) {
        stop(cy_format_txt)
      }
      if (nrow(current_years_by_ind) != ncol(indicator_time_series[ ,-1]) |
          all(names(indicator_time_series)[-1] %in% current_years_by_ind$ind) == FALSE ) {
        stop("The data frame provided in 'current_years_by_ind' must include all
             indicators listed in 'indicator_time_series'.")
      }
      if ( all(current_years_by_ind$start %in% indicator_time_series[ ,1]) == FALSE |
          all(current_years_by_ind$end %in% indicator_time_series[ ,1]) == FALSE ) {
        stop("Some of the years defined in 'current_years_by_ind' fall outside
             the available time series range.")
      }
      if (any(current_years_by_ind$start > current_years_by_ind$end)) {
        stop("The start year cannot be later than the end year in 'current_years_by_ind'.")
      }
    }
  }

  # Method settings
  if (all(sign %in% c("+", "-")) == FALSE) {
    warning("'sign' can only be set to '+' or '-'. The default value '+' will be used.")
  }
  if (all(condition %in% c(">", "<")) == FALSE) {
    warning("'condition' can only be set to '<' or '>'. The default value '>' will be used.")
  }


  # ---
  ind_ts <- indicator_time_series

# Calculate for each indicator mean and assigned range in baseline and current
  baseline <- data.frame(
    indicator = names(ind_ts[-1]),
    mean = rep(0, ncol(ind_ts) - 1),
    dev = rep(0, ncol(ind_ts) - 1)
  )

  current <- data.frame(
    indicator = names(ind_ts[-1]),
    mean = rep(0, ncol(ind_ts) - 1),
    dev = rep(0, ncol(ind_ts) - 1)
  )

  # Check input length of range
  if (length(range) == 1) {
    range <- rep(range, ncol(ind_ts) - 1)
  }


  for (i in 1:(ncol(ind_ts) - 1)) {

    ### Baseline period
    # Create matching baseline years vector
    helper_ind <- names(ind_ts[i + 1])
    if (!is.null(base_years)) {
      helper_by <- base_years[1]:base_years[2]
    } else {
      helper_by <- base_years_by_ind$start[base_years_by_ind$ind ==
          helper_ind]:base_years_by_ind$end[base_years_by_ind$ind == helper_ind]
    }

    baseline$mean[i] <- mean(as.numeric(ind_ts[ind_ts[, 1] %in%
        helper_by, i + 1]), na.rm = TRUE)

    if (typeof(range[i]) == "character") {
      if (!(range[i] %in% c("sd", "2sd", "95percentile", "mean_only"))) {
        warning("Deviance range not recognised. Will use range = 'sd' as default.")
        range <- rep("sd", ncol(ind_ts) - 1)
      }

      if (range[i] == "sd") {
        baseline$dev[i] <- stats::sd(as.numeric(ind_ts[ind_ts[,
          1] %in% helper_by, i + 1]), na.rm = TRUE)
      } else if (range[i] == "2sd") {
        baseline$dev[i] <- 2 * (stats::sd(as.numeric(ind_ts[ind_ts[,
          1] %in% helper_by, i + 1]), na.rm = TRUE))
      } else if (range[i] == "95percentile") {
        baseline$dev[i] <- stats::quantile(as.numeric(ind_ts[ind_ts[,
          1] %in% helper_by, i + 1]), na.rm = TRUE, probs = 0.95)
      } else {
        # mean only
        baseline$dev[i] <- 0
      }
    } else if (typeof(range[i]) == "double") {
      baseline$dev[i] <- stats::quantile(as.numeric(ind_ts[ind_ts[,
        1] %in% helper_by, i + 1]), na.rm = TRUE, probs = (range[i]/100))
    }

    ### Current period
    # Create matching current years vector
    helper_ind <- names(ind_ts[i + 1])
    if (!is.null(current_years)) {
      helper_cy <- current_years[1]:current_years[2]
    } else {
      helper_cy <- current_years_by_ind$start[current_years_by_ind$ind ==
          helper_ind]:current_years_by_ind$end[current_years_by_ind$ind == helper_ind]
    }

    current$mean[i] <- mean(as.numeric(ind_ts[ind_ts[, 1] %in%
        helper_cy, i + 1]), na.rm = TRUE)

  }


  # Test if a specified condition and sign is set for each indicator
  # or if general settings is selected
  if (length(condition) == 1) {
    condition <- rep(condition, ncol(ind_ts) - 1)
  }
  if (length(sign) == 1) {
    sign <- rep(sign, ncol(ind_ts) - 1)
  }

  # Now set status score for each indicator: can be either
  #   good -> risk = vulnerability decreased by one
  #   undesired -> risk = vulnerability increased by one

  status_helper <- data.frame(
    indicator = current$indicator,
    status = rep(NA, nrow(current)),
    score = rep(0, nrow(current))
  )

  for (i in 1:nrow(current)) {
    if (range[i] == "mean_only") {
      if (condition[i] == ">") {
        # above baseline mean
        if (current$mean[i] > baseline$mean[i]) {
          status_helper$score[i] <- 1
          status_helper$status[i] <- "good"
        } else {
          status_helper$score[i] <- -1
          status_helper$status[i] <- "undesired"
        }
      } else {
        # below baseline mean
        if (current$mean[i] < baseline$mean[i]) {
          status_helper$score[i] <- 1
          status_helper$status[i] <- "good"
        } else {
          status_helper$score[i] <- -1
          status_helper$status[i] <- "undesired"
        }
      }
    } else {
      if (sign[i] == "+") {
        # above baseline mean above deviance
        if (condition[i] == ">") {
          if (current$mean[i] > (baseline$mean[i] + baseline$dev[i])) {
            status_helper$score[i] <- 1
            status_helper$status[i] <- "good"
          } else {
            status_helper$score[i] <- -1
            status_helper$status[i] <- "undesired"
          }
        } else {
          # below deviance
          if (current$mean[i] < (baseline$mean[i] + baseline$dev[i])) {
            status_helper$score[i] <- 1
            status_helper$status[i] <- "good"
          } else {
            status_helper$score[i] <- -1
            status_helper$status[i] <- "undesired"
          }
        }
      } else {
        # below baseline mean above deviance
        if (condition[i] == ">") {
          if (current$mean[i] > (baseline$mean[i] - baseline$dev[i])) {
            status_helper$score[i] <- 1
            status_helper$status[i] <- "good"
          } else {
            status_helper$score[i] <- -1
            status_helper$status[i] <- "undesired"
          }
        } else {
          # below deviance
          if (current$mean[i] < (baseline$mean[i] - baseline$dev[i])) {
            status_helper$score[i] <- 1
            status_helper$status[i] <- "good"
          } else {
            status_helper$score[i] <- -1
            status_helper$status[i] <- "undesired"
          }
        }
      }
    }
  }  # end status scoring loop

  output_status <- status_helper

  return(output_status)
}
