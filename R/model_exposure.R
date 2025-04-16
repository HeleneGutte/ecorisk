#' Model Overall Exposure Scores Using Time Series Data
#'
#' This function statistically evaluates the exposure to a pressure, based on time
#' series data. The scoring is based on the paper of Gaichas et al., 2014:
#' A risk-based approach to evaluating northeast US fish community vulnerability
#' to climate change. The exposure scoring is split into four components:
#' magnitude (or degree of change), frequency of change, the future trend of
#' the pressure, and spatial scale. Uncertainty of the exposure assessment is
#' evaluated using general additive models (GAM) and an autoregressive integrated
#' moving average model (ARIMA).
#'
#' @param pressure_time_series A data frame (not a tibble object) with time series
#'        of pressures to be evaluated. First column MUST be the time column.
#' @param base_years A vector with two numerics, specifying the time period for
#'        the baseline. The first one \code{start} is the starting year for all
#'        pressures and the second one \code{end} is the end of the baseline for
#'        all pressures. The default is NULL. One can specify pressure specific
#'        baseline periods using the \code{base_years_by_press} argument.
#'        If \code{base_years} and \code{base_years_by_ind} are NULL, then the
#'        first 5 years of the time series are used as baseline period.
#' @param base_years_by_press A data frame, specifying the baseline years for each
#'        pressure individually, by setting the starting year (second column)
#'        and the end year (third column). The first column must contain the
#'        names of the pressure indicators used in \code{pressure_time_series}.
#'        The default is NULL. If \code{base_years} and \code{base_years_by_press}
#'        are NULL, then the first 5 years of the time series are used as baseline period.
#' @param current_years A vector with two numerics, specifying the time period for
#'        the assessment period. The first one \code{start} is the starting year
#'        for all pressures and the second one \code{end} is the end of the
#'        assessment period for all pressures. The default is NULL.
#'        One can specify pressure specific assessment periods using the
#'        \code{current_years_by_press} argument.
#'        If \code{current_years} and \code{current_years_by_press} are NULL, then the
#'        last 5 years of the time series are used as assessment period.
#' @param current_years_by_press A data frame, specifying the assessment period
#'        for each pressure individually, by setting the starting year (second column)
#'        and the end year (third column). The first column must contain the
#'        names of the pressure indicators used in \code{pressure_time_series}.
#'        The default is NULL. If \code{current_years} and \code{current_years_by_press}
#'        are NULL, then the last 5 years of the time series are used as assessment period.
#' @param trend a character vector specifying whether a trend returning to the
#'        baseline conditions should be considered as good or a trend further
#'        leaving the baseline conditions. Possible inputs are \code{return} or
#'        \code{leave}.Default is not specified is \code{return}, meaning a
#'        return to the baseline is desired.
#' @param spatial a vector with scores for the spatial scale of each pressure.
#'        The default is 3 for each pressure, meaning that 40 - 60% of the entire
#'        assessment area is affected. Scores should be on a scale from 1 - 5,
#'        depending on the percent of area that is affected by the pressure:
#'        \itemize{
#'          \item 1: < 20%,
#'          \item 2: 20 - 40%,
#'          \item 3: 40 - 60%,
#'          \item 4: 60 - 80%,
#'          \item 5: > 80%.
#'          }
#'
#' @details All components are scored on a scale from 1 - 5, low impact to high
#' impact. The degree of change compares the mean of the current time period to the
#' baseline time period, the score is based on standard deviations.
#' The frequency evaluates in how much percent of the current time period the
#' mean deviates more than one standard deviation from the baseline mean.
#' The future trend scores if the pressure will in the future be in desired
#' conditions or not. Usually this means the pressure returns to the baseline
#' conditions. The overall exposure score is the mean of all four components.
#' Uncertainty of exposure is evaluated using a general additive model and
#' an autoregressive integrated moving average model (ARIMA). The models are
#' fitted using the time series except the assessment period. The assessment
#' period is then predicted. The function evaluates how many of the observed
#' data points are within the predicted 95% confidence interval.
#' If more than 66 % are within the 95% CI the uncertainty is 1 (low), if less
#' than 33 % are within it, the uncertainty is set to 3 (high).
#' Additionally the function compares the mean size of the predicted 95%
#' confidence interval and compares it to the maximum range of
#' the observed data points to account for very large confidence intervals, which
#' would otherwise lead to too optimistic uncertainty scores. The lower uncertainty
#' score is selected as final uncertainty score.
#' The time periods of baseline and assessment period have to be carefully set to
#' reflect ongoing dynamics. Especially for oscillating pressures time periods
#' should be longer to assess the overall trend and not the oscillation itself.
#'
#' @return a data frame containing the pressure names, the aggregated exposure score
#' and scores for magnitude, frequency, future trend and spatial scale of the
#' pressures, the final uncertainty score and uncertainty scores of the ARIMA and the GAM.
#' If default settings are used, the following data frame will be returned:
#'
#' \describe{
#'   \item{\code{pressure}}{Name of the assessed pressure.}
#'   \item{\code{exposure}}{Exposure score, mean of the four assessed exposure
#'        components.}
#'   \item{\code{uncertainty}}{Uncertainty score associated with the exposure
#'        assessment.}
#'   \item{\code{comp_magnitude}}{Score for the magnitude of change.}
#'   \item{\code{comp_frequency}}{Score for the frequency of a significant deviation
#'        from baseline conditions.}
#'   \item{\code{comp_trend}}{Score for the future trend of the pressure.}
#'   \item{\code{comp_direction}}{Direction of the development of the pressure in
#'        the assessment period.}
#'   \item{\code{comp_spatial}}{Score for the spatial scale, either set by the user
#'        or automatically set to 3.}
#'   \item{\code{uncertainty_arima}}{Uncertainty score based on the ARIMA model.}
#'   \item{\code{uncertainty_gam}}{Uncertainty based on the GAM.}
#'   \item{\code{mean_baseline}}{Mean of the baseline conditions, used for
#'        magnitude scoring.}
#'   \item{\code{mean_current}}{Mean of the current conditions, used for magnitude
#'        and frequency scoring.}
#'   \item{\code{standard_deviation_baseline}}{Standard deviations of the baseline
#'        conditions. Used for scoring of magnitude and frequency.}
#'   \item{\code{slope_linear_model}}{Slope of the linear model used for scoring
#'        the future trend and to determine the direction.}
#'   \item{\code{p_value_linear_model}}{P-value of the linear model, used to
#'        score the future trend.}
#' }
#'
#' @seealso \code{\link{model_sensitivity}}, \code{\link{vulnerability}}
#'
#' @export
#'
#' @examples
#' ### Example with 8 pressure time series in the demo data 'pressure_ts_baltic'
#' #   where the first 11 years represent the general baseline period and the last
#' #   7 years of the time series the current assessment period:
#' model_exposure(
#'   pressure_time_series = pressure_ts_baltic,
#'   base_years = c(start = 1984, end = 1994),
#'   current_years = c(start = 2010, end = 2016)
#' )
#'
#' ### Example with 2 pressure time series and pressure-specific periods
#' sub_ts <- pressure_ts_baltic[ ,c("year", "nitrogen", "phosphorous")]
#' model_exposure(
#'   pressure_time_series = sub_ts,
#'   base_years_by_press = data.frame(
#'     press = c("nitrogen", "phosphorous"),
#'     start = c(1984, 1990), end = c(1994, 2000)),
#'   current_years_by_press = data.frame(
#'       press = c("nitrogen", "phosphorous"),
#'       start = c(2010, 2012), end = c(2016, 2016))
#' )

model_exposure <- function(pressure_time_series,
  base_years = NULL, base_years_by_press = NULL,
  current_years = NULL, current_years_by_press = NULL,
  trend = "return", spatial = 3) {

  ### Data input validations
  # Time series
  if (any(class(pressure_time_series) != "data.frame")) {
    stop("The 'pressure_time_series' argument must be a data frame
         (not a matrix or tibble) in which all columns are numeric.
         The first column should represent the time component (e.g., years).")
  }
  if (all(sapply(pressure_time_series, is.numeric)) == FALSE) {
    stop("The 'pressure_time_series' argument must be a data frame where all
         columns are numeric. The first column should represent the time
         component (e.g., years).")
  }
  if (all(pressure_time_series[, 1] == sort(pressure_time_series[, 1])) == FALSE) {
    warning("The first column in 'pressure_time_series' appears to be out of order.
            Are you sure it represents the time component (e.g., years)?")
  }

  # Baseline period
  if (is.null(base_years) & is.null(base_years_by_press)) {
    warning(paste0("You did not provide the first and last year of the baseline period. ",
      "You can either define this generally ('base_years = c(firstyear, lastyear)') or specify it ",
      "for each pressure using 'base_years_by_press'. The default will ",
      "be the first 5 years of the time series. Refer to the R documentation for further guidance."))
    yrs <- pressure_time_series[ ,1]
    base_years <- c(min(yrs), min(yrs)+4)
  }

  if (!is.null(base_years)) {
    if (!is.null(dim(base_years))) {
      stop("'base_years' must be a numeric vector specifying the start and end
           year of the baseline period, which will be applied to all pressures.")
    } else {
      if (length(base_years) != 2) {
        stop("'base_years' must be a numeric vector specifying the start and end
             year of the baseline period, which will be applied to all pressures.")
      }
    }
    if (all(base_years %in% pressure_time_series[ ,1]) == FALSE ) {
      stop("The years defined in 'base_years' fall outside the available time series range.")
    }
    if (base_years[1] > base_years[2]) {
      stop("The start year cannot be later than the end year in 'base_years'.")
    }
  }

  if (!is.null(base_years_by_press)) {
    by_format_txt <- paste0("The 'base_years_by_press' argument must be a data frame with three columns: ",
      "'press' (a character vector containing the pressure names) and two numeric columns, ",
      "'start' and 'end' (the first and last years of the baseline period). See the examples ",
      "section in the R documentation for further guidance.")
    if (any(class(base_years_by_press) != "data.frame")) {
      stop(by_format_txt)
    } else {
      if (is.character(base_years_by_press[, 1]) == FALSE |
          all(sapply(base_years_by_press[, 2:3], is.numeric)) == FALSE |
          !identical(names(base_years_by_press), c("press", "start", "end")) ) {
        stop(by_format_txt)
      }
      if (nrow(base_years_by_press) != ncol(pressure_time_series[ ,-1]) |
          all(names(pressure_time_series)[-1] %in% base_years_by_press$press) == FALSE ) {
        stop("The data frame provided in 'base_years_by_press' must include all
             pressures listed in 'pressure_time_series'.")
      }
      if ( all(base_years_by_press$start %in% pressure_time_series[ ,1]) == FALSE |
          all(base_years_by_press$end %in% pressure_time_series[ ,1]) == FALSE ) {
        stop("Some of the years defined in 'base_years_by_press' fall outside
             the available time series range.")
      }
      if (any(base_years_by_press$start > base_years_by_press$end)) {
        stop("The start year cannot be later than the end year in 'base_years_by_press'.")
      }
    }
  }

  # Current (assessment) period
  if (is.null(current_years) & is.null(current_years_by_press)) {
    warning(paste0("You did not provide the first and last year of the current assessment period. ",
      "You can either define this generally ('current_years = c(firstyear, lastyear)') or specify it ",
      "for each pressure using 'current_years_by_press'. The default will ",
      "be the last 5 years of the time series. Refer to the R documentation for further guidance."))
    yrs <- pressure_time_series[ ,1]
    current_years <- c(max(yrs)-4, max(yrs))
  }

  if (!is.null(current_years)) {
    if (!is.null(dim(current_years))) {
      stop("'current_years' must be a numeric vector specifying the start and end year of the current assessment period, which will be applied to all pressures.")
    } else {
      if (length(current_years) != 2) {
        stop("'current_years' must be a numeric vector specifying the start and end year of the current assessment period, which will be applied to all pressures.")
      }
    }
    if (all(current_years %in% pressure_time_series[ ,1]) == FALSE ) {
      stop("The years defined in 'current_years' fall outside the available time series range.")
    }
    if (current_years[1] > current_years[2]) {
      stop("The start year cannot be later than the end year in 'current_years'.")
    }
  }

  if (!is.null(current_years_by_press)) {
    cy_format_txt <- paste0("The 'current_years_by_press' argument must be a data frame with three columns: ",
      "'press' (a character vector containing the pressure names) and two numeric columns, ",
      "'start' and 'end' (the first and last years of the current assessment period). See the examples ",
      "section in the R documentation for further guidance.")
    if (any(class(current_years_by_press) != "data.frame")) {
      stop(cy_format_txt)
    } else {
      if (is.character(current_years_by_press[, 1]) == FALSE |
          all(sapply(current_years_by_press[, 2:3], is.numeric)) == FALSE |
          !identical(names(current_years_by_press), c("press", "start", "end")) ) {
        stop(cy_format_txt)
      }
      if (nrow(current_years_by_press) != ncol(pressure_time_series[ ,-1]) |
          all(names(pressure_time_series)[-1] %in% current_years_by_press$press) == FALSE ) {
        stop("The data frame provided in 'current_years_by_press' must include all pressures listed in 'pressure_time_series'.")
      }
      if ( all(current_years_by_press$start %in% pressure_time_series[ ,1]) == FALSE |
          all(current_years_by_press$end %in% pressure_time_series[ ,1]) == FALSE ) {
        stop("Some of the years defined in 'current_years_by_press' fall outside the available time series range.")
      }
      if (any(current_years_by_press$start > current_years_by_press$end)) {
        stop("The start year cannot be later than the end year in 'current_years_by_press'.")
      }
    }
  }


  # ---
  press_ts <- pressure_time_series

  baseline <- data.frame(
    pressures = names(press_ts[-1]),
    mean = rep(0, ncol(press_ts) - 1),
    sd = rep(0, ncol(press_ts) - 1)
  )

  current <- data.frame(
    pressures = names(press_ts[-1]),
    mean = rep(0, ncol(press_ts) - 1),
    sd = rep(0, ncol(press_ts) - 1)
  )


  for (i in 1:(ncol(press_ts) - 1)) {

    ### Baseline period
    # Create matching baseline years vector
    helper_pres <- names(press_ts[i + 1])
    if (!is.null(base_years)) {
      helper_by <- base_years[1]:base_years[2]
    } else {
      helper_by <- base_years_by_press$start[base_years_by_press$press ==
          helper_pres]:base_years_by_press$end[base_years_by_press$press == helper_pres]
    }

    baseline$mean[i] <- mean(as.numeric(press_ts[press_ts[, 1] %in%
        helper_by, i+1]), na.rm = TRUE)
    baseline$sd[i] <- stats::sd(as.numeric(press_ts[press_ts[, 1] %in%
        helper_by, i+1]), na.rm = TRUE)

    # Extract 95% confidence interval of both periods to evaluate whether
    # one of the means falls with in the CI of the other period:
    it <- 10000
    # empty vector:
    sm <- numeric(it)
    # loop:
    for (l in 1:it) {
      xs <- sample(as.numeric(press_ts[press_ts[, 1] %in% helper_by, i+1]),
        replace = T)
      sm[l] <- mean(xs, na.rm = TRUE)
    }
    cis <- stats::quantile(x = sm, probs = c(0.025, 0.975), na.rm = TRUE)
    baseline$ci_low[i] <- unname(cis[1])
    baseline$ci_up[i] <- unname(cis[2])


    ### Current period
    # Create matching current years vector
    helper_pres <- names(press_ts[i + 1])
    if (!is.null(current_years)) {
      helper_cy <- current_years[1]:current_years[2]
    } else {
      helper_cy <- current_years_by_press$start[current_years_by_press$press ==
          helper_pres]:current_years_by_press$end[current_years_by_press$press == helper_pres]
    }

    current$mean[i] <- mean(as.numeric(press_ts[press_ts[ ,1] %in%
        helper_cy, i+1]), na.rm = TRUE)
    current$sd[i] <- stats::sd(as.numeric(press_ts[press_ts[ ,1] %in%
        helper_cy, i+1]), na.rm = TRUE)

    # Extract 95% confidence interval of both periods to evaluate whether
    # one of the means falls with in the CI of the other period:
    it <- 10000
    # empty vector:
    sm <- numeric(it)
    # loop:
    for (l in 1:it) {
      xs <- sample(as.numeric(press_ts[press_ts[, 1] %in% helper_cy, i+1]),
        replace = T)
      sm[l] <- mean(xs, na.rm = TRUE)
    }
    cis <- stats::quantile(x = sm, probs = c(0.025, 0.975), na.rm = TRUE)
    current$ci_low[i] <- unname(cis[1])
    current$ci_up[i] <- unname(cis[2])
  }

  if (length(trend) == 1) {
    trend <- rep(trend, ncol(press_ts) - 1)
  }
  if (length(spatial) == 1) {
    spatial <- rep(spatial, ncol(press_ts) - 1)
  }


  ### Scoring ----
  scores <- data.frame(
    pressure = current$pressures,
    exposure = rep(0, nrow(current)),
    uncertainty = rep(NA, nrow(current)),
    comp_magnitude = rep(0, nrow(current)),
    comp_frequency = rep(0, nrow(current)),
    comp_trend = rep(0, nrow(current)),
    comp_direction = rep(0, nrow(current)),
    comp_spatial = rep(0, nrow(current)),
    uncertainty_arima = rep(NA, nrow(current)),
    uncertainty_gam = rep(NA, nrow(current)),
    mean_baseline = rep(NA, nrow(current)),
    mean_current = rep(NA, nrow(current)),
    standard_deviation_baseline = rep(NA, nrow(current)),
    slope_linear_model = rep(NA, nrow(current)),
    p_value_linear_model = rep(NA, nrow(current))
  )


  for (i in 1:nrow(baseline)) {
    while (is.na(current$sd[i]) == TRUE & i < nrow(baseline)) {
      i <- i + 1
    }
    if (is.na(current$sd[i]) == TRUE & i == nrow(baseline)) {
      break
    }
    # Magnitude (degree) ---- above baseline
    if (current$mean[i] > baseline$mean[i]) {
      if (current$mean[i] > (baseline$mean[i] + 4 * baseline$sd[i])) {
        scores$comp_magnitude[i] <- 5
      } else if (current$mean[i] > (baseline$mean[i] + 3 * baseline$sd[i])) {
        scores$comp_magnitude[i] <- 4
      } else if (current$mean[i] > (baseline$mean[i] + 2 * baseline$sd[i])) {
        scores$comp_magnitude[i] <- 3
      } else if (current$mean[i] > (baseline$mean[i] + 1 * baseline$sd[i])) {
        scores$comp_magnitude[i] <- 2
      } else {
        scores$comp_magnitude[i] <- 1
      }
    } else {
      # below baseline
      if (current$mean[i] < (baseline$mean[i] - 4 * baseline$sd[i])) {
        scores$comp_magnitude[i] <- 5
      } else if (current$mean[i] < (baseline$mean[i] - 3 * baseline$sd[i])) {
        scores$comp_magnitude[i] <- 4
      } else if (current$mean[i] < (baseline$mean[i] - 2 * baseline$sd[i])) {
        scores$comp_magnitude[i] <- 3
      } else if (current$mean[i] < (baseline$mean[i] - 1 * baseline$sd[i])) {
        scores$comp_magnitude[i] <- 2
      } else {
        scores$comp_magnitude[i] <- 1
      }
    }

    # Evaluate whether means overlap with the CI of the other period
    if (baseline$mean[i] >= current$ci_low[i] &
        baseline$mean[i] <= current$ci_up[i]) {
      if (scores$comp_magnitude[i] == 1) {
        scores$comp_magnitude[i] <- 1
      } else {
        scores$comp_magnitude[i] <- scores$comp_magnitude[i] - 1
      }
    } else if (current$mean[i] >= baseline$ci_low[i] &
        current$mean[i] <= baseline$ci_up[i]) {
      if (scores$comp_magnitude[i] == 1) {
        scores$comp_magnitude[i] <- 1
      } else {
        scores$comp_magnitude[i] <- scores$comp_magnitude[i] - 1
      }
    }


    # Frequency (duration)  ----
    helper_pres <- names(press_ts[i + 1])
    if (!is.null(current_years)) {
      helper_cy <- current_years[1]:current_years[2]
    } else {
      helper_cy <- current_years_by_press$start[current_years_by_press$press ==
          helper_pres]:current_years_by_press$end[current_years_by_press$press == helper_pres]
    }

    curr_cond <- data.frame(
      pressure = press_ts[press_ts[, 1] %in% helper_cy, i + 1],
      year = press_ts[press_ts[, 1] %in% helper_cy, 1]
    )
    if (current$mean[i] > baseline$mean[i]) {
      temp <- sum(curr_cond[, 1] > baseline$mean[i] + baseline$sd[i], na.rm = TRUE)
    } else {
      temp <- sum(curr_cond[, 1] < baseline$mean[i] - baseline$sd[i], na.rm = TRUE)
    }

    length_current <- nrow(curr_cond)
    if (temp <= (0.2 * length_current)) {
      scores$comp_frequency[i] <- 1
    } else if (temp <= (0.4 * length_current)) {
      scores$comp_frequency[i] <- 2
    } else if (temp <= (0.6 * length_current)) {
      scores$comp_frequency[i] <- 3
    } else if (temp <= (0.8 * length_current)) {
      scores$comp_frequency[i] <- 4
    } else {
      scores$comp_frequency[i] <- 5
    }

    # Trend analysis with linear model ----
    fm <- stats::lm(scale(curr_cond[, 1]) ~ curr_cond$year)  #standardize data with z-transformation
    sm_fm <- summary(fm)
    # Extract p-value and slope
    p_val <- sm_fm$coefficients[8]
    slope <- sm_fm$coefficients[2]
    intercept <- sm_fm$coefficients[1]

    if (slope > 0) {
      scores$comp_direction[i] <- "increase"
    } else {
      scores$comp_direction[i] <- "decrease"
    }

    if (trend[i] == "return") {
      if (current$mean[i] < baseline$mean[i]) {
        # current below baseline test if returning to baseline:
        if (p_val > 0.05) {
          scores$comp_trend[i] <- 3
        } else if (slope > 1) {
          scores$comp_trend[i] <- 1
        } else if (slope > 0) {
          scores$comp_trend[i] <- 2
        } else if (slope > -1) {
          scores$comp_trend[i] <- 4
        } else {
          scores$comp_trend[i] <- 5
        }
      } else {
        # current above baseline test if returning to baseline:
        if (p_val > 0.05) {
          scores$comp_trend[i] <- 3
        } else if (slope > 1) {
          scores$comp_trend[i] <- 5
        } else if (slope > 0) {
          scores$comp_trend[i] <- 4
        } else if (slope > -1) {
          scores$comp_trend[i] <- 2
        } else {
          scores$comp_trend[i] <- 1
        }
      }
    } else {
      # Leaving from baseline below baseline test if returning to baseline:
      if (current$mean[i] < baseline$mean[i]) {
        if (p_val > 0.05) {
          scores$comp_trend[i] <- 3
        } else if (slope > 1) {
          scores$comp_trend[i] <- 5
        } else if (slope > 0) {
          scores$comp_trend[i] <- 4
        } else if (slope > -1) {
          scores$comp_trend[i] <- 2
        } else {
          scores$comp_trend[i] <- 1
        }
      } else {
        # above baseline test if returning to baseline:
        if (p_val > 0.05) {
          scores$comp_trend[i] <- 3
        } else if (slope > 1) {
          scores$comp_trend[i] <- 1
        } else if (slope > 0) {
          scores$comp_trend[i] <- 2
        } else if (slope > -1) {
          scores$comp_trend[i] <- 4
        } else {
          scores$comp_trend[i] <- 5
        }
      }
    }

    # Spatial scale
    scores$comp_spatial[i] <- spatial[i]

    # Uncertainty ----

    # GAM for the entire time series except the current conditions
    not_curr_cond <- data.frame(
      pressure = press_ts[!(press_ts[, 1] %in% helper_cy), i + 1],
      year = press_ts[!(press_ts[, 1] %in% helper_cy), 1]
    )
    temp_gam <- mgcv::gam(pressure ~ 1 + s(year, k = 4),
      data = not_curr_cond, family = stats::gaussian())
    # ARIMA model for the entire time series except the current conditions
    auto_arima <- forecast::auto.arima(y = stats::ts(not_curr_cond$pressure,
      start = min(not_curr_cond$year), end = max(not_curr_cond$year),
      frequency = 1))
    p <- length(auto_arima$model$phi)
    d <- length(auto_arima$model$Delta)
    q <- length(auto_arima$model$theta)
    arima_model <- stats::arima(x = stats::ts(not_curr_cond$pressure,
      start = min(not_curr_cond$year), end = max(not_curr_cond$year),
      frequency = 1), order = c(p, d, q))

    # Calculate ARIMA uncertainty predict current conditions
    pred_arima <- stats::predict(arima_model, n.ahead = length(helper_cy),
      se.fit = TRUE)
    ci_low_arima <- pred_arima$pred - 1.96 * pred_arima$se
    ci_up_arima <- pred_arima$pred + 1.96 * pred_arima$se
    # Count how many of the observed data points are within the
    # ± 95% CI range of the predicted values
    counter_arima <- sum(curr_cond[, 1] > ci_low_arima & curr_cond[,1] <
                           ci_up_arima, na.rm = TRUE)

    # Calculate mean 95 % confidence interval and compare to maximum
    # range of supplied values
    max_range_press <- max(press_ts[, i + 1], na.rm = TRUE) -
      min(press_ts[,i + 1], na.rm = TRUE)
    mean_confid_arima <- mean(ci_up_arima - ci_low_arima,
      na.rm = TRUE)/max_range_press

    if (mean_confid_arima < (1/3) & counter_arima >= (2/3 * length_current)) {
      scores$uncertainty_arima[i] <- 1
    } else if (mean_confid_arima < (1/3) & counter_arima >= (1/3 * length_current)) {
      scores$uncertainty_arima[i] <- 2
    } else if (mean_confid_arima < (1/3) & counter_arima < (1/3 * length_current)) {
      scores$uncertainty_arima[i] <- 2
    } else if (mean_confid_arima < (2/3) & counter_arima >= (2/3 * length_current)) {
      scores$uncertainty_arima[i] <- 2
    } else if (mean_confid_arima < (2/3) & counter_arima >= (1/3 * length_current)) {
      scores$uncertainty_arima[i] <- 2
    } else if (mean_confid_arima < (2/3) & counter_arima < (1/3 * length_current)) {
      scores$uncertainty_arima[i] <- 3
    } else if (mean_confid_arima >= (2/3) & counter_arima >= (2/3 * length_current)) {
      scores$uncertainty_arima[i] <- 2
    } else if (mean_confid_arima >= (2/3) & counter_arima >= (1/3 * length_current)) {
      scores$uncertainty_arima[i] <- 3
    } else if (mean_confid_arima >= (2/3) & counter_arima < (1/3 * length_current)) {
      scores$uncertainty_arima[i] <- 3
    }

    # Calculate GAM uncertainty predict current conditions
    press_seq <- press_ts[press_ts[, 1] %in% helper_cy, 1]
    new_data <- data.frame(press = NA, year = press_seq)
    pred_gam <- mgcv::predict.gam(temp_gam, type = "response", se.fit = TRUE,
      newdata = new_data)
    ci_low_gam <- pred_gam$fit - 1.96 * pred_gam$se.fit
    ci_up_gam <- pred_gam$fit + 1.96 * pred_gam$se.fit
    # Count how many of the observed data points are within the
    # ± 95% CI range of the predicted values
    counter_gam <- sum(curr_cond[, 1] > ci_low_gam & curr_cond[, 1] <
        ci_up_gam, na.rm = TRUE)

    # Calculate mean 95 % confidence interval and compare to maximum
    # range of supplied values
    max_range_press <- max(press_ts[, i + 1], na.rm = TRUE) -
      min(press_ts[,i + 1], na.rm = TRUE)
    mean_confid_gam <- mean(ci_up_gam - ci_low_gam,
      na.rm = TRUE)/max_range_press

    if (mean_confid_gam < (1/3) & counter_gam >= (2/3 * length_current)) {
      scores$uncertainty_gam[i] <- 1
    } else if (mean_confid_gam < (1/3) & counter_gam >= (1/3 * length_current)) {
      scores$uncertainty_gam[i] <- 2
    } else if (mean_confid_gam < (1/3) & counter_gam < (1/3 * length_current)) {
      scores$uncertainty_gam[i] <- 2
    } else if (mean_confid_gam < (2/3) & counter_gam >= (2/3 * length_current)) {
      scores$uncertainty_gam[i] <- 2
    } else if (mean_confid_gam < (2/3) & counter_gam >= (1/3 * length_current)) {
      scores$uncertainty_gam[i] <- 2
    } else if (mean_confid_gam < (2/3) & counter_gam < (1/3 * length_current)) {
      scores$uncertainty_gam[i] <- 3
    } else if (mean_confid_gam >= (2/3) & counter_gam >= (2/3 * length_current)) {
      scores$uncertainty_gam[i] <- 2
    } else if (mean_confid_gam >= (2/3) & counter_gam >= (1/3 * length_current)) {
      scores$uncertainty_gam[i] <- 3
    } else if (mean_confid_gam >= (2/3) & counter_gam < (1/3 * length_current)) {
      scores$uncertainty_gam[i] <- 3
    }


    if (scores$uncertainty_arima[i] < scores$uncertainty_gam[i]) {
      scores$uncertainty[i] <- scores$uncertainty_arima[i]
    } else {
      scores$uncertainty[i] <- scores$uncertainty_gam[i]
    }

    # extract diagnostics
    scores$mean_baseline[i] <- baseline$mean[i]
    scores$mean_current[i] <- current$mean[i]
    scores$standard_deviation_baseline[i] <- baseline$sd[i]
    scores$slope_linear_model[i] <- slope
    scores$p_value_linear_model[i] <- p_val

  }  # End of scoring loop

  output_model_exp <- scores
  exposure <- ecorisk::calc_exposure(
    pressure = scores$pressure,
    components = scores[ ,c(4:6,8)]
  )
  output_model_exp$exposure <- exposure$exposure

  return(output_model_exp)
}

