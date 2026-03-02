#' Produce Diagnostic Plots for the 
#' 
#' The function `plot_diagnostic_sensitivity()` creates diagnostic plots for the
#' generalized additive models, which are used as basis for the sensitivity scoring 
#' in \code{\link{model_sensitivity}}. 
#' 
#' @param indicator_time_series a data frame containing only the state indicator
#'        time series. First column MUST be the time column.
#' @param pressure_time_series a data frame containing only the pressure variables.
#'        First column MUST be the time column.
#'
#' @returns The function returns a list of ggplot objects combined with patchwork.
#'  For each state and pressure indicator combination 4 diagnostic plots are created:
#'  Q-Q plots, residuals vs- linear predictor, a histogram of the residuals and 
#'  response vs. fitted values. 
#'  
#' @export
#'
#' @examples
#' \dontrun{
#' plot_diagnostic_sensitivity(
#'   indicator_time_series = indicator_ts_baltic[, c(1,2)],
#'   pressure_time_series = pressure_ts_baltic[, c(1,2)]
#'  )
#'  }

plot_diagnostic_sensitivity <- function(indicator_time_series, pressure_time_series){
  
  ### check if suggested packages are installed: 
  if(!requireNamespace(c("mgcViz"), quietly = TRUE)){
    stop(
      "Package \"mgcViz\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if(!requireNamespace(c("patchwork"), quietly = TRUE)){
    stop(
      "Package \"patchwork\" must be installed to use this function.",
      call. = FALSE
    )
  }
  
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
  
  # Create a data frame where all indicators are combined with each pressure
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
  
  
  diag_plot_list <- vector(mode = "list", length = nrow(dat))
  # apply GAMs
  for(j in 1:nrow(dat)){
    # get input for gam
    ind_list <- dat$ind_data[[j]]
    press_list <- dat$press_data[[j]]
    if (all(is.na(press_list)) == TRUE & j < nrow(dat)) {
      warning(paste(c("The pressure", dat$press[j],
                      "only consists of NAs, making it impossible to evaluate further."),
                    collapse = " "))
      next
    }
    if (all(is.na(press_list)) == TRUE & j == nrow(dat)) {
      warning(paste(c("The pressure", dat$press[j],
                      "only consists of NAs, making it impossible to evaluate further."),
                    collapse = " "))
      break
    }
    if (all(is.na(ind_list)) == TRUE & j < nrow(dat)) {
      warning(paste(c("The indicator", dat$ind[j],
                      "only consists of NAs, making it impossible to evaluate further."),
                    collapse = " "))
      next
    }
    if (all(is.na(ind_list)) == TRUE & j == nrow(dat)) {
      warning(paste(c("The indicator", dat$ind[j],
                      "only consists of NAs, , making it impossible to evaluate further."),
                    collapse = " "))
      break
    }
    dat_temp <- data.frame(ind = ind_list, press = press_list)
    # Apply gam
    temp_gam <- mgcv::gam(ind ~ 1 + s(press, k = 4), family = stats::gaussian(),
                          data = dat_temp)
    
    temp_gam_viz <- mgcViz::getViz(temp_gam)
    temp_diag_plot <- suppressMessages(mgcViz::check.gamViz(temp_gam_viz)) # does not work yet

    diag_plot_list[[j]] <- patchwork::wrap_plots(temp_diag_plot) + 
      patchwork::plot_annotation(title = paste0(dat$ind[j], " ~ ", dat$press[j]))
  }
  
  # add message what to do with the plots and how to proceed, when not satisfied 
  # with the model diagnostics
  
  message("Please review the model diagnostics of the GAMs applied in the time series based sensitivity scoring. Remove models with unacceptable diagnostics from the output table of the model_sensitivity() function.")
  
  return(diag_plot_list)
  
}