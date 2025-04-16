#' Model Overall Sensitivity Scores Using Time Series Data
#'
#' The function `model_sensitivity()` uses time series of a state indicator and
#' a pressure variable to assess the state indicators sensitivity towards the
#' pressure. The relationship between pressure and state indicator is determined
#' using a generalized additive model (GAM). Uncertainty is evaluated with a GAM
#' and an ARIMA model.
#'
#' @param indicator_time_series a data frame containing only the state indicator
#'        time series. First column MUST be the time column.
#' @param pressure_time_series a data frame containing only the pressure variables.
#'        First column MUST be the time column.
#' @param current_years A vector with two numerics, specifying the time period for
#'        the assessment period. The first one \code{start} is the starting year
#'        for all pressure-indicator pairs and the second one \code{end} is the end of the
#'        assessment period for all pressure-indicator pairs. The default is NULL.
#'        One can specify pair specific assessment periods using the
#'        \code{current_years_by_ind_press} argument.
#'        If \code{current_years} and \code{current_years_by_ind_press} are NULL,
#'        then the last 5 years of the time series are used as assessment period.
#' @param current_years_by_ind_press a data frame specifying for each indicator-pressure
#'        pair the starting (third column) and end year (fourth column) where the
#'        current conditions are best reflected. The default is NULL.
#'        If \code{current_years} and \code{current_years_by_ind_press} are NULL,
#'        then the last 5 years of the time series are used as assessment period.
#'
#' @details In case the relationship of one state indicator - pressure pair is not
#' significant the sensitivity score is 0, and thus also vulnerability and risk
#' will be 0. For a significant relationship the score will be set based on the
#' R-squared value from 1 (R-squared < 0.2) to 5 (R-squared ≥ 0.8). Additionally,
#' the function evaluates the edf score of the GAM which indicates the degree of
#' non-linearity in the relationship. Since highly non-linear relationships are
#' harder to predict, the risk of reaching an undesired state increases and the
#' sensitivity score for nonlinear relationships will be increased by 1
#' (if it was not 5 already). The direction of an effect (negative influence or
#' positive influence of the pressure) is evaluated with the slope of a linear
#' model representing the assessment period. If the slope of the linear model is
#' negative, the direction of effect is considered negative as well, and vice
#' versa for the positive effect.
#' The function assesses uncertainty associated with the scoring based on a
#' general additive model and an autoregressive integrated moving average model
#' (ARIMA). The ARIMA model uses the pressure variable as additional external predictor.
#' The models are fitted using the time series except the assessment period.
#' The assessment period is then predicted. The function evaluates how many of
#' the observed data points are within the predicted 95% confidence interval.
#' If more than 66 % are within the 95% CI the uncertainty is 1 (low), if less
#' than 33 % are within it, the uncertainty is set to 3 (high).
#' Additionally the function compares the mean size of the predicted 95%
#' confidence interval and compares it to the maximum range of
#' the observed data points to account for very large confidence intervals, which
#' would otherwise lead to too optimistic uncertainty scores. The lower uncertainty
#' score is selected as final uncertainty score.
#'
#' The function also creates columns to give the opportunity to assess adaptive
#' capacity and its associated uncertainty for each state indicator-pressure pair.
#' The scores for adaptive capacity and its associated uncertainty must be specified
#' before the next function \code{\link{vulnerability}} is applied (see examples).
#' If adaptive capacity and its uncertainty are not further specified, this will
#' influence the further application of the ecorisk framework.
#'
#' @return a data frame containing indicator, pressure, type of effect, the
#' sensitivity score and the associated uncertainty. Positive sensitivity scores
#' are associated with a positive influence of the pressure on the indicator and vice versa.
#' Additionally the R-squared, p-values, edf scores and the mean confidence
#' interval percentage, which are the basis of the scoring, are provided.
#' The type of effect is automatically set to\code{direct_indirect} as the model
#' cannot distinguish between direct and indirect effects.
#' If default settings are used, the following data frame will be returned:
#' \describe{
#'   \item{\code{indicator}}{Name of the assessed state indicator.}
#'   \item{\code{pressure}}{Name of the assessed pressure.}
#'   \item{\code{type}}{Type of the assessed effect.}
#'   \item{\code{pathway}}{The pathway that has been used to derive the sensitivity scores.}
#'   \item{\code{sensitivity}}{Sensitivity score for the assessed state indicator-
#'   pressure pair.}
#'   \item{\code{adaptive_capacity}}{Adaptive capacity score for the assessed
#'        state indicator-pressure pair, is automatically set to 0 and can be
#'        changed afterwards.}
#'   \item{\code{uncertainty_sens}}{uncertainty score associated with the sensitivity
#'        scoring.}
#'   \item{\code{uncertainty_ac}}{uncertainty score for adaptive capacity scoring.
#'        Automatically set to NA, can be changed afterwards.}
#'   \item{\code{r_sq}}{R-squared value of the GAM, used for the sensitivity scoring.}
#'   \item{\code{p_value}}{P-value of the GAM, used to identify significant
#'        relationships. Unsignificant relationships get a sensitivity score of 0.}
#'   \item{\code{edf}}{Estimated degrees of freedom, used to assess non-linearity
#'        of the relationship between state indicator and pressure. }
#'   \item{\code{uncertainty_gam}}{Uncertainty score for sensitivity based on
#'        predicted values from a GAM.}
#'   \item{\code{uncertainty_arima}}{Uncertainty score for sensitivity based on
#'        predicted values from an ARIMA using the pressure variable as external
#'        predictor.}
#' }
#'
#' @seealso \code{\link{model_exposure}}, \code{\link{vulnerability}}
#'
#' @export
#'
#' @examples
#' ### Example with the 2 indicators and 8 pressure time series in the Baltic Sea demo data
#' #   where the last 7 years of the time series represent the current assessment period:
#' model_sensitivity(
#'   indicator_time_series = indicator_ts_baltic,
#'   pressure_time_series = pressure_ts_baltic,
#'   current_years = c(start = 2010, end = 2016)
#' )
#'
#' ### Example with the demo data but indicator-pressure-specific assessment periods:
#' sens_tbl <- model_sensitivity(
#'   indicator_time_series = indicator_ts_baltic,
#'   pressure_time_series = pressure_ts_baltic,
#'   current_years_by_ind_press = data.frame(
#'     ind = rep(names(indicator_ts_baltic)[-1], each = 8),
#'     press = rep(names(pressure_ts_baltic)[-1], 2),
#'     start = c(rep(2010, 8), rep(2008, 8)),
#'     end = c(rep(2016, 8), rep(2015, 8))
#'   )
#' )
# # Adjust the adaptive capacity scores (score from -1 to 1, default is 0) and
#' # add the associated uncertainty (from 1 to 3, default is NA)
#' sens_tbl$adaptive_capacity <- c(0,0,1,1,1,1,-1,-1, -1,-1,1,1,1,1,1,-1)
#' sens_tbl$uncertainty_ac <- c(2,2,1,1,1,1,2,1, 3,3,1,1,2,2,3,1)
#' sens_tbl

model_sensitivity <- function(indicator_time_series, pressure_time_series,
  current_years = NULL, current_years_by_ind_press = NULL) {

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
  if (any(class(pressure_time_series) != "data.frame")) {
    stop("The 'pressure_time_series' argument must be a data frame (not a
         matrix or tibble) in which all columns are numeric. The first column
         should represent the time component (e.g., years).")
  }
  if (all(sapply(pressure_time_series, is.numeric)) == FALSE) {
    stop("The 'pressure_time_series' argument must be a data frame where all
         columns are numeric. The first column should represent the time
         component (e.g., years).")
  }
  if (all(indicator_time_series[ ,1] == pressure_time_series[ ,1]) == FALSE) {
    stop("The time columns in 'indicator_time_series' and 'pressure_time_series'
         are not identical.")
  }
  if (all(indicator_time_series[, 1] == sort(indicator_time_series[, 1])) == FALSE) {
    warning("The first columns in 'indicator_time_series' and 'pressure_time_series'
            appear to be out of order. Are you sure they represent the time
            component (e.g., years)?")
  }

  # Current assessment period
  if (is.null(current_years) & is.null(current_years_by_ind_press)) {
    warning(paste0("You did not provide the first and last year of the current assessment period. ",
      "You can either define this generally ('current_years = c(firstyear, lastyear)') or specify it ",
      "for each indicator-pressure combination using 'current_years_by_ind_press'. The default will ",
      "be the last 5 years of the time series. Refer to the R documentation for further guidance."))
    yrs <- indicator_time_series[ ,1]
    current_years <- c(max(yrs)-4, max(yrs))
  }

  if (!is.null(current_years)) {
    if (!is.null(dim(current_years))) {
      stop("'current_years' must be a numeric vector specifying the start and
           end year of the current assessment period, which will be applied to
           all pressures.")
    } else {
      if (length(current_years) != 2) {
        stop("'current_years' must be a numeric vector specifying the start and
             end year of the current assessment period, which will be applied to
             all pressures.")
      }
    }
    if (all(current_years %in% indicator_time_series[ ,1]) == FALSE ) {
      stop("The years defined in 'current_years' fall outside the available time
           series range.")
    }
    if (current_years[1] > current_years[2]) {
      stop("The start year cannot be later than the end year in 'current_years'.")
    }
  }

  if (!is.null(current_years_by_ind_press)) {
    cy_format_txt <- paste0("The 'current_years_by_ind_press' argument must be a data frame with four ",
      "columns: 'ind' (a character vector containing the indicator names), 'press' (a character vector ",
      "containing the pressure names), and two numeric columns, 'start' and 'end' (the first and ",
      "last years of the evaluation period). See the examples section in the R documentation for ",
      "further guidance.")
    if (any(class(current_years_by_ind_press) != "data.frame")) {
      stop(cy_format_txt)
    } else {
      if (all(sapply(current_years_by_ind_press[, 1:2], is.character)) == FALSE|
          all(sapply(current_years_by_ind_press[, 3:4], is.numeric)) == FALSE |
          !identical(names(current_years_by_ind_press), c("ind", "press", "start", "end")) ) {
        stop(cy_format_txt)
      }
      n_combo <- ncol(indicator_time_series[ ,-1])*ncol(pressure_time_series[ ,-1])
      if (nrow(current_years_by_ind_press) != n_combo |
          all(names(indicator_time_series)[-1] %in% current_years_by_ind_press$ind) == FALSE |
          all(names(pressure_time_series)[-1] %in% current_years_by_ind_press$press) == FALSE) {
        stop("The data frame provided in 'current_years_by_ind_press' must include
             all indicator-pressure combinations.")
      }
      if ( all(current_years_by_ind_press$start %in% pressure_time_series[ ,1]) == FALSE |
          all(current_years_by_ind_press$end %in% pressure_time_series[ ,1]) == FALSE ) {
        stop("Some of the years defined in 'current_years_by_ind_press' fall
             outside the available time series range.")
      }
      if (any(current_years_by_ind_press$start > current_years_by_ind_press$end)) {
        stop("The start year cannot be later than the end year in
             'current_years_by_ind_press'.")
      }
    }
  }


  # ---
  # Create a data frame where all indicators are combined with each pressure
  time <- indicator_time_series[, 1]
  indicator <- tibble::as_tibble(indicator_time_series)
  pressure <- tibble::as_tibble(pressure_time_series)

  ind <- rep(x = names(indicator)[-1], each = ncol(pressure) - 1)
  press <- rep(x = names(pressure)[-1], times = ncol(indicator) - 1)

  dat <- data.frame(ind = ind, press = press)

  # Add the data (as list-column)
  for (i in 1:nrow(dat)) {
    dat$press_data[i] <- pressure[, names(pressure) == dat$press[i]]
    dat$ind_data[i] <- indicator[, names(indicator) == dat$ind[i]]
  }

  # Create sensitivity scoring table
  scores <- 1:5
  length_criteria <- 1/length(scores)
  score_table <- data.frame(score = scores, criteria = rep(0, length(scores)))

  for (j in 1:length(scores)) {
    score_table$criteria[j] <- (((j - 1) * length_criteria):(j * length_criteria))
  }
  score_table[nrow(score_table) + 1, ] <- c(NA, 1)

  # Set up table to save GAM output and scoring
  gam <- dat[, 1:2]
  gam$type <- "direct_indirect"
  gam$sen_score <- rep(NA, nrow(gam))
  gam$r_sq <- rep(NA, nrow(gam))
  gam$edf <- rep(NA, nrow(gam))
  gam$p_value <- rep(NA, nrow(gam))
  gam$unc <- rep(NA, nrow(gam))
  gam$mean_confid <- rep(NA, nrow(gam))
  gam$pathway <- rep("model", nrow(gam))
  gam$unc_gam <- rep(NA, nrow(gam))
  gam$unc_arima <- rep(NA, nrow(gam))


  # Apply GAMs and evaluate output to derive scores
  for (l in 1:nrow(dat)) {
    # get input for gam
    ind_list <- dat$ind_data[[l]]
    press_list <- dat$press_data[[l]]
    if (all(is.na(press_list)) == TRUE & l < nrow(dat)) {
      warning(paste(c("The pressure", dat$press[l],
        "only consists of NAs, making it impossible to evaluate further."),
        collapse = " "))
      next
    }
    if (all(is.na(press_list)) == TRUE & l == nrow(dat)) {
      warning(paste(c("The pressure", dat$press[l],
        "only consists of NAs, making it impossible to evaluate further."),
        collapse = " "))
      break
    }
    if (all(is.na(ind_list)) == TRUE & l < nrow(dat)) {
      warning(paste(c("The indicator", dat$ind[l],
        "only consists of NAs, making it impossible to evaluate further."),
        collapse = " "))
      next
    }
    if (all(is.na(ind_list)) == TRUE & l == nrow(dat)) {
      warning(paste(c("The indicator", dat$ind[l],
        "only consists of NAs, , making it impossible to evaluate further."),
        collapse = " "))
      break
    }
    dat_temp <- data.frame(ind = ind_list, press = press_list)
    # Apply gam
    temp_gam <- mgcv::gam(ind ~ 1 + s(press, k = 4), family = stats::gaussian(),
      data = dat_temp)
    # Extract necessary output of GAM
    smy_gam <- mgcv::summary.gam(temp_gam)
    r_sq <- smy_gam$r.sq
    edf <- smy_gam$edf
    p_value <- smy_gam$s.table[[4]]

    # Calculate first derivative
    derivs <- function(x, y, z) {
      as.vector((x - y)/(2 * z))
    }

    # Extract starting and end point of current period
    if (!is.null(current_years)) {
      begin <- current_years[1]
      end <- current_years[2]
    } else {
      begin <- current_years_by_ind_press[current_years_by_ind_press$ind == dat$ind[l] &
          current_years_by_ind_press$press == dat$press[l], "start"]
      end <- current_years_by_ind_press[current_years_by_ind_press$ind == dat$ind[l] &
          current_years_by_ind_press$press == dat$press[l], "end"]
    }

    # Create linear model to evaluate direction and current strength of
    # ind - press relationship
    lm <- lm(scale(ind) ~ scale(press),
      data = dat_temp[c(which(time == begin):which(time == end)), ])
    smy_lm <- summary(lm)

    # Evaluate output and set sensitivity score
    if (is.nan(p_value) == TRUE) {
      gam$sen_score[l] <- NA
      gam$r_sq[l] <- NA
      gam$edf[l] <- NA
      gam$p_value[l] <- NA
    } else if (p_value < 0.05) {
      if (smy_lm$coefficients[2, 1] >= -0.5 & smy_lm$coefficients[2, 1] <= 0.5) {
        # only weak slope, thus weaker sensitivity score
        if (r_sq <= 0.4) {
          gam$sen_score[l] <- min(score_table$score, na.rm = TRUE)
          gam$r_sq[l] <- r_sq
          gam$edf[l] <- edf
          gam$p_value[l] <- p_value
        } else if (r_sq <= 0.8) {
          gam$sen_score[l] <- score_table$score[2]
          gam$r_sq[l] <- r_sq
          gam$edf[l] <- edf
          gam$p_value[l] <- p_value
        } else {
          gam$sen_score[l] <- score_table$score[3]
          gam$r_sq[l] <- r_sq
          gam$edf[l] <- edf
          gam$p_value[l] <- p_value
        }
      } else {
        # stronger slope values, thus 'normal' criteria are applied
        if (r_sq < 0) {
          gam$sen_score[l] <- min(score_table$score, na.rm = TRUE)
          gam$r_sq[l] <- r_sq
          gam$edf[l] <- edf
          gam$p_value[l] <- p_value
        } else {
          for (m in 1:nrow(score_table)) {
            if (r_sq >= score_table$criteria[m] &
                r_sq < score_table$criteria[m + 1]) {
              gam$sen_score[l] <- score_table$score[m]
              gam$r_sq[l] <- r_sq
              gam$edf[l] <- edf
              gam$p_value[l] <- p_value
            }
          }
        }
      }
    } else {
      gam$sen_score[l] <- 0
      gam$r_sq[l] <- r_sq
      gam$edf[l] <- edf
      gam$p_value[l] <- p_value
    }  #end if statement sen score
    if (edf > 1.1 & gam$sen_score[l] < 5 & is.nan(p_value) == FALSE) {
      gam$sen_score[l] <- gam$sen_score[l] + 1
    }

    # Direction of effect in last years in the linear model
    if (smy_lm$coefficients[2, 1] < 0) {
      gam$sen_score[l] <- gam$sen_score[l] * -1
    }


    # Uncertainty ----
    not_curr_cond <- dat_temp[-c(which(time == begin):which(time == end)), ]
    not_curr_cond_time <- time[-c(which(time == begin):which(time == end))]
    length_current <- length(begin:end)

    # GAM for the entire time series except the current conditions
    temp_gam_unc <- mgcv::gam(ind ~ 1 + s(press, k = 4), family = stats::gaussian(),
      data = not_curr_cond)

    press_seq <- dat_temp$press[c(which(time == begin):which(time == end))]
    new_data <- data.frame(ind = NA, press = press_seq)
    pred_gam <- mgcv::predict.gam(temp_gam_unc, type = "response",
      se.fit = TRUE, newdata = new_data)
    ci_low_gam <- pred_gam$fit - 1.96 * pred_gam$se.fit
    ci_up_gam <- pred_gam$fit + 1.96 * pred_gam$se.fit

    counter_gam <- sum(dat_temp[c(which(time == begin):which(time == end)),]$ind >
                     ci_low_gam & dat_temp[c(which(time == begin):which(time ==end)), ]$ind <
                     ci_up_gam, na.rm = TRUE)

    max_range_ind <- max(dat_temp$ind, na.rm = TRUE) - min(dat_temp$ind,
      na.rm = TRUE)
    mean_confid_gam <- mean(ci_up_gam - ci_low_gam, na.rm = TRUE)/max_range_ind

    # Evaluate mean confidence interval percentage and prediction
    # accuracy (how many data points are within the 95% CI)
    if (mean_confid_gam < (1/3) & counter_gam >= (2/3 * length_current)) {
      gam$unc_gam[l] <- 1
    } else if (mean_confid_gam < (1/3) & counter_gam >= (1/3 * length_current)) {
      gam$unc_gam[l] <- 2
    } else if (mean_confid_gam < (1/3) & counter_gam < (1/3 * length_current)) {
      gam$unc_gam[l] <- 2
    } else if (mean_confid_gam < (2/3) & counter_gam >= (2/3 * length_current)) {
      gam$unc_gam[l] <- 2
    } else if (mean_confid_gam < (2/3) & counter_gam >= (1/3 * length_current)) {
      gam$unc_gam[l] <- 2
    } else if (mean_confid_gam < (2/3) & counter_gam < (1/3 * length_current)) {
      gam$unc_gam[l] <- 3
    } else if (mean_confid_gam >= (2/3) & counter_gam >= (2/3 * length_current)) {
      gam$unc_gam[l] <- 2
    } else if (mean_confid_gam >= (2/3) & counter_gam >= (1/3 * length_current)) {
      gam$unc_gam[l] <- 3
    } else if (mean_confid_gam >= (2/3) & counter_gam < (1/3 * length_current)) {
      gam$unc_gam[l] <- 3
    }

    # Calculate ARIMA without current conditions, predict current conditions
    # and count how many data points are actually in the predicted 95 % CIs
    # account for large CIs
    auto_arima <- forecast::auto.arima(y = not_curr_cond$ind, xreg = not_curr_cond$press)
    p <- length(auto_arima$model$phi)
    d <- length(auto_arima$model$Delta)
    q <- length(auto_arima$model$theta)
    arima_model <- stats::arima(x = not_curr_cond$ind, xreg = not_curr_cond$press, order = c(p, d, q))

    pred_arima <- stats::predict(arima_model, n.ahead = length_current,
                                 se.fit = TRUE, newxreg = press_seq)
    ci_low_arima <- pred_arima$pred - 1.96 * pred_arima$se
    ci_up_arima <- pred_arima$pred + 1.96 * pred_arima$se

    counter_arima <- sum(dat_temp[c(which(time == begin):which(time == end)),]$ind >
                           ci_low_arima & dat_temp[c(which(time == begin):which(time ==end)), ]$ind <
                           ci_up_arima, na.rm = TRUE)
    mean_confid_arima <- mean(ci_up_arima - ci_low_arima,
                              na.rm = TRUE)/max_range_ind

    if (mean_confid_arima < (1/3) & counter_arima >= (2/3 * length_current)) {
      gam$unc_arima[l] <- 1
    } else if (mean_confid_arima < (1/3) & counter_arima >= (1/3 * length_current)) {
      gam$unc_arima[l] <- 2
    } else if (mean_confid_arima < (1/3) & counter_arima < (1/3 * length_current)) {
      gam$unc_arima[l] <- 2
    } else if (mean_confid_arima < (2/3) & counter_arima >= (2/3 * length_current)) {
      gam$unc_arima[l] <- 2
    } else if (mean_confid_arima < (2/3) & counter_arima >= (1/3 * length_current)) {
      gam$unc_arima[l] <- 2
    } else if (mean_confid_arima < (2/3) & counter_arima < (1/3 * length_current)) {
      gam$unc_arima[l] <- 3
    } else if (mean_confid_arima >= (2/3) & counter_arima >= (2/3 * length_current)) {
      gam$unc_arima[l] <- 2
    } else if (mean_confid_arima >= (2/3) & counter_arima >= (1/3 * length_current)) {
      gam$unc_arima[l] <- 3
    } else if (mean_confid_arima >= (2/3) & counter_arima < (1/3 * length_current)) {
      gam$unc_arima[l] <- 3
    }

    if (gam$unc_arima[l] < gam$unc_gam[l]) {
      gam$unc[l] <- gam$unc_arima[l]
    } else {
      gam$unc[l] <- gam$unc_gam[l]
    }
  }  # End of Scoring for-loop




  # Create output data frame (set adaptive_capacities as default to 0 to not
  # influence the vulnerability calculation, uncertainty_ac will be NA as default)
  output_model_sens <- data.frame(
    indicator = gam$ind,
    pressure = gam$press,
    type = gam$type,
    pathway = gam$pathway,
    sensitivity = gam$sen_score,
    adaptive_capacity = rep(0, nrow(gam)),
    uncertainty_sens  = gam$unc,
    uncertainty_ac = rep(NA, nrow(gam)),
    r_sq = gam$r_sq,
    p_value = gam$p_value,
    edf = gam$edf,
    uncertainty_gam = gam$unc_gam,
    uncertainty_arima = gam$unc_arima
  )

  return(output_model_sens)

}
