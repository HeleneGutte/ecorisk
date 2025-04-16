#' Compute High-Complexity Multi-Risk Scores and (Eco)system Risk
#'
#' The function \code{aggregate_risk} uses the output of the risk function or
#' the vulnerability function. The risk or vulnerability scores are aggregated
#' in three ways:
#' 1. as multi-pressure risk per state indicator, representing the overall
#' effect on one indicator
#' 2. as multi- state indicator risk per pressure, representing the overall
#' effect each pressure has on all state indicators
#' 3. as ecosystem risk, the combined multi-pressure risks of all indicators.
#'
#' @param risk_results Risk score for each state indicator ~ pressure ~
#'        type combination, derived from the function \code{\link{risk}}
#'        or \code{\link{vulnerability}}.
#' @param method Character indicating the method for aggregating the pressures and
#'        state indicators to multiple risk scores and the ecosystem risk score.
#'        The use can choose between (arithmetic) `mean`, `median`, `sum`, `maximum`
#'        and `minimum`. Default is the arithmetic mean.
#'
#' @details
#' The returned lists are required in the plotting functions \code{\link{plot_radar}}
#' and \code{\link{plot_heatmap}}.
#' The aggregated scores are calculated for each type and pathway combination
#' individually and across all types and pathways. If only one type and/or one
#' pathway has been evaluated beforehand, the results will be the same for the
#' different combinations.
#'
#' @return
#' The function returns a list containing three sublists, first the multi-state
#' indicator risk list, containing the risks and uncertainties aggregated per
#' pressure, type and pathway.
#' Second the multi-pressure risk list, where risks and uncertainties are aggregated
#' per state indicator and type and pathway. The third list contains the ecosystem
#' risks, which aggregates the multi-pressure risk and uncertainty scores per
#' type and pathway.
#'
#'
#'
#' @seealso \code{\link{vulnerability}}, \code{\link{risk}},
#'          \code{\link{plot_radar}}, \code{\link{plot_heatmap}}
#'
#' @export
#'
#' @examples
#' ### Demo with example output from the risk() function based on expert scores
#' # (where direct and direct/indirect effects were evaluated)
#'
#' # Calculate mean risks scores per indicator/pressure/ecosystem:
#' mean_risk <- aggregate_risk(
#'   risk_results = ex_output_risk_expert,
#'   method = "mean" # default
#' )
#' mean_risk
#' # Calculate median risks scores:
#' aggregate_risk(
#'   risk_results = ex_output_risk_expert,
#'   method = "median"
#' )
#' # Calculate maximum risks scores:
#' aggregate_risk(
#'   risk_results = ex_output_risk_expert,
#'   method = "maximum"
#' )
#'
#'
#' ### Demo with example output from the risk() function based on modelled
#' #   scores (where only direct/indirect effects were evaluated)
#'
#' # Calculate mean risks scores:
#' aggregate_risk(risk_results = ex_output_risk_model)
#'
#'
#' ### Demo with combined expert-based and model-based pathways
#'
#' combined_risk <- rbind(ex_output_risk_expert, ex_output_risk_model)
#' aggr_risk <- aggregate_risk(risk_results = combined_risk)
#' aggr_risk
#'
# Look at individual results, e.g., combined pathways and combined types
#  (direct and direct/indirect):
#' aggr_risk$multi_indicator_risk |>
#'   dplyr::filter(type == "combined", pathway == "combined")
#' aggr_risk$multi_pressure_risk |>
#'   dplyr::filter(type == "combined", pathway == "combined")
#' aggr_risk$ecosystem_risk |>
#'   dplyr::filter(type == "combined", pathway == "combined")
#'
#'
#' ### Demo with vulnerability scores using example output data from
#' #   vulnerability() based on modelled scores
#'
#' aggregate_risk(risk_results = ex_output_vulnerability_model)

aggregate_risk <- function(risk_results, method = "mean") {

  if (!(method %in% c("mean", "median", "maximum", "minimum", "sum"))) {
    warning("Method not recognised. Will use method = arithmetic mean instead.")
  }

  # Extract different types and pathways, so that cumulative effects
  # are calculated for each type individually
  types <- unique(risk_results$type)
  pathways <- unique(risk_results$pathway)
  ind <- unique(risk_results$indicator)

  # Check if it is a vulnerability or a risk output and rename in
  # case it is vulnerability
  if ("vulnerability" %in% names(risk_results) == TRUE &
      !("risk" %in% names(risk_results) == TRUE)) {
    names(risk_results)[names(risk_results) == "vulnerability"] <- "risk"
  }


  # Calculate aggregated risk per indicator, pathway and type (multi-
  # pressure risk on each indicator)
  for (h in 1:length(pathways)) {
    for (i in 1:length(types)) {
      dat <- risk_results[risk_results$type == types[i] &
          risk_results$pathway ==  pathways[h], ]
      if (nrow(dat) == 0) {
        next
      }
      # get inds in this subset
      ind_helper <- unique(dat$indicator)
      risk <- rep(0, length(ind_helper))
      indicators <- rep(NA, length(ind_helper))
      type <- rep(NA, length(ind_helper))
      pathway <- rep(NA, length(ind_helper))
      uncertainty <- rep(NA, length(ind_helper))
      # calculate for each ind the aggregated risk
      for (j in 1:length(ind_helper)) {
        # choose method
        if (method == "median") {
          risk[j] <- stats::median(dat[dat$indicator == ind_helper[j], ]$risk,
            na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- stats::median(
              dat[dat$indicator == ind_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "maximum") {
          risk[j] <- max(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- max(
              dat[dat$indicator == ind_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "minimum") {
          risk[j] <- min(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- min(
              dat[dat$indicator == ind_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "sum") {
          risk[j] <- sum(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- sum(
              dat[dat$indicator == ind_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
          risk[j] <- mean(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- mean(
              dat[dat$indicator == ind_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        indicators[j] <- ind_helper[j]
        type[j] <- types[i]
        pathway[j] <- pathways[h]
        multi_press_helper <- data.frame(
          indicator = indicators,
          type = type,
          pathway = pathway,
          risk = risk,
          uncertainty = uncertainty
        )
      }
      if (i == 1) {
        multi_press <- multi_press_helper
      } else {
        multi_press <- rbind(multi_press, multi_press_helper)
      }
    }

  }  # End multi-pressure risk per indicators, type and pathway

  # Calculate multi-pressure risk per indicators and pathway and
  # across everything
  for (i in 1:length(pathways)) {
    dat <- risk_results[risk_results$pathway == pathways[i], ]
    # get inds in this subset
    ind_helper <- unique(dat$indicator)
    risk <- rep(0, length(ind_helper))
    indicators <- rep(NA, length(ind_helper))
    type <- rep(NA, length(ind_helper))
    pathway <- rep(NA, length(ind_helper))
    uncertainty <- rep(NA, length(ind_helper))
    # calculate for each ind the aggregated risk
    for (j in 1:length(ind_helper)) {
      # choose method
      if (method == "median") {
        risk[j] <- stats::median(dat[dat$indicator == ind_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- stats::median(
            dat[dat$indicator == ind_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "maximum") {
        risk[j] <- max(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- max(
            dat[dat$indicator == ind_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "minimum") {
        risk[j] <- min(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- min(
            dat[dat$indicator == ind_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "sum") {
        risk[j] <- sum(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- sum(
            dat[dat$indicator == ind_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
        risk[j] <- mean(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- mean(
            dat[dat$indicator == ind_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      indicators[j] <- ind_helper[j]
      type[j] <- "combined"
      pathway[j] <- pathways[i]
      multi_press_helper <- data.frame(
        indicator = indicators,
        type = type,
        pathway = pathway,
        risk = risk,
        uncertainty = uncertainty
      )
    }
    if (i == 1) {
      multi_press_2 <- multi_press_helper
    } else {
      multi_press_2 <- rbind(multi_press_2, multi_press_helper)
    }
    if (i == max(length(pathways))) {
      dat <- risk_results
      # get inds in this subset
      ind_helper <- unique(dat$indicator)
      risk <- rep(0, length(ind_helper))
      indicators <- rep(NA, length(ind_helper))
      type <- rep(NA, length(ind_helper))
      pathway <- rep(NA, length(ind_helper))
      uncertainty <- rep(NA, length(ind_helper))
      # calculate for each ind the aggregated risk
      for (j in 1:length(ind_helper)) {
        # choose method
        if (method == "median") {
          risk[j] <- stats::median(dat[dat$indicator == ind_helper[j], ]$risk,
            na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- stats::median(
              dat[dat$indicator == ind_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "maximum") {
          risk[j] <- max(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- max(
              dat[dat$indicator == ind_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "minimum") {
          risk[j] <- min(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- min(
              dat[dat$indicator == ind_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "sum") {
          risk[j] <- sum(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) ==
              FALSE) {
            uncertainty[j] <- sum(
              dat[dat$indicator == ind_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
          risk[j] <- mean(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- mean(
              dat[dat$indicator == ind_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        indicators[j] <- ind_helper[j]
        type[j] <- "combined"
        pathway[j] <- "combined"
        multi_press_helper <- data.frame(
          indicator = indicators,
          type = type,
          pathway = pathway,
          risk = risk,
          uncertainty = uncertainty
        )
      }
      multi_press_2 <- rbind(multi_press_2, multi_press_helper)

    }
  }
  # Calculate per type, regardless off pathway
  for (i in 1:length(types)) {
    dat <- risk_results[risk_results$type == types[i], ]
    # get inds in this subset
    ind_helper <- unique(dat$indicator)
    risk <- rep(0, length(ind_helper))
    indicators <- rep(NA, length(ind_helper))
    type <- rep(NA, length(ind_helper))
    pathway <- rep(NA, length(ind_helper))
    uncertainty <- rep(NA, length(ind_helper))
    # calculate for each ind the aggregated risk
    for (j in 1:length(ind_helper)) {
      # choose method
      if (method == "median") {
        risk[j] <- stats::median(dat[dat$indicator == ind_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- stats::median(
            dat[dat$indicator == ind_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "maximum") {
        risk[j] <- max(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- max(
            dat[dat$indicator == ind_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "minimum") {
        risk[j] <- min(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- min(
            dat[dat$indicator == ind_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "sum") {
        risk[j] <- sum(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- sum(
            dat[dat$indicator == ind_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
        risk[j] <- mean(dat[dat$indicator == ind_helper[j], ]$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- mean(
            dat[dat$indicator == ind_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      indicators[j] <- ind_helper[j]
      type[j] <- types[i]
      pathway[j] <- "combined"
      multi_press_type_helper <- data.frame(
        indicator = indicators,
        type = type,
        pathway = pathway,
        risk = risk,
        uncertainty = uncertainty
      )
    }
    if (i == 1) {
      multi_press_type <- multi_press_type_helper
    } else {
      multi_press_type <- rbind(multi_press_type, multi_press_type_helper)
    }
  }

  multi_press_all <- rbind(multi_press, multi_press_type, multi_press_2)
  multi_press_all <- multi_press_all[order(multi_press_all$indicator), ]
  # End multi-pressure risks



  # Calculate aggregated risk per pressure, pathway and type (multi
  # indicator risk of each pressure)
  for (h in 1:length(pathways)) {
    for (i in 1:length(types)) {
      dat <- risk_results[risk_results$type == types[i] &
          risk_results$pathway ==  pathways[h], ]
      if (nrow(dat) == 0) {
        next
      }
      # get press in this subset
      press_helper <- unique(dat$pressure)
      risk <- rep(0, length(press_helper))
      pressures <- rep(NA, length(press_helper))
      type <- rep(NA, length(press_helper))
      pathway <- rep(NA, length(press_helper))
      uncertainty <- rep(NA, length(press_helper))

      # calculate for each ind the aggregated risk
      for (j in 1:length(press_helper)) {
        # choose method
        if (method == "median") {
          risk[j] <- stats::median(dat[dat$pressure == press_helper[j], ]$risk,
            na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- stats::median(
              dat[dat$pressure == press_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "maximum") {
          risk[j] <- max(dat[dat$pressure == press_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- max(
              dat[dat$pressure == press_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "minimum") {
          risk[j] <- min(dat[dat$pressure == press_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- min(
              dat[dat$pressure == press_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "sum") {
          risk[j] <- sum(dat[dat$pressure == press_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- sum(
              dat[dat$pressure == press_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
          risk[j] <- mean(dat[dat$pressure == press_helper[j], ]$risk, na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- mean(
              dat[dat$pressure == press_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        pressures[j] <- press_helper[j]
        type[j] <- types[i]
        pathway[j] <- pathways[h]
        multi_ind_helper <- data.frame(
          pressure = pressures,
          type = type,
          pathway = pathway,
          risk = risk,
          uncertainty = uncertainty
        )
      }
      if (i == 1) {
        multi_ind <- multi_ind_helper
      } else {
        multi_ind <- rbind(multi_ind, multi_ind_helper)
      }
    }
  }  # End multi-indicator risks per pressure, pathway and type

  # Calculate multi-indicator risk per pressure and pathway and
  # across everything
  for (i in 1:length(pathways)) {
    dat <- risk_results[risk_results$pathway == pathways[i], ]
    # get press in this subset
    press_helper <- unique(dat$pressure)
    risk <- rep(0, length(press_helper))
    pressures <- rep(NA, length(press_helper))
    type <- rep(NA, length(press_helper))
    pathway <- rep(NA, length(press_helper))
    uncertainty <- rep(NA, length(press_helper))

    # calculate for each ind the aggregated risk
    for (j in 1:length(press_helper)) {
      # choose method
      if (method == "median") {
        risk[j] <- stats::median(dat[dat$pressure == press_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- stats::median(
            dat[dat$pressure == press_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "maximum") {
        risk[j] <- max(dat[dat$pressure == press_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- max(
            dat[dat$pressure == press_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "minimum") {
        risk[j] <- min(dat[dat$pressure == press_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- min(
            dat[dat$pressure == press_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "sum") {
        risk[j] <- sum(dat[dat$pressure == press_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- sum(
            dat[dat$pressure == press_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
        risk[j] <- mean(dat[dat$pressure == press_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- mean(
            dat[dat$pressure == press_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      pressures[j] <- press_helper[j]
      type[j] <- "combined"
      pathway[j] <- pathways[i]
      multi_ind_helper <- data.frame(
        pressure = pressures,
        type = type,
        pathway = pathway,
        risk = risk,
        uncertainty = uncertainty
      )
    }
    if (i == 1) {
      multi_ind_2 <- multi_ind_helper
    } else {
      multi_ind_2 <- rbind(multi_ind_2, multi_ind_helper)
    }
    if (i == max(length(pathways))) {
      dat <- risk_results
      # get inds in this subset
      press_helper <- unique(dat$pressure)
      risk <- rep(0, length(press_helper))
      pressures <- rep(NA, length(press_helper))
      type <- rep(NA, length(press_helper))
      pathway <- rep(NA, length(press_helper))
      uncertainty <- rep(NA, length(press_helper))
      # calculate for each ind the aggregated risk
      for (j in 1:length(press_helper)) {
        # choose method
        if (method == "median") {
          risk[j] <- stats::median(dat[dat$pressure == press_helper[j], ]$risk,
            na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- stats::median(
              dat[dat$pressure == press_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "maximum") {
          risk[j] <- max(dat[dat$pressure == press_helper[j], ]$risk,
            na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- max(
              dat[dat$pressure == press_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "minimum") {
          risk[j] <- min(dat[dat$pressure == press_helper[j], ]$risk,
            na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- min(
              dat[dat$pressure == press_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (method == "sum") {
          risk[j] <- sum(dat[dat$pressure == press_helper[j], ]$risk,
            na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- sum(
              dat[dat$pressure == press_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
          risk[j] <- mean(dat[dat$pressure == press_helper[j], ]$risk,
            na.rm = TRUE)
          if (all(is.na(dat$uncertainty)) == FALSE) {
            uncertainty[j] <- mean(
              dat[dat$pressure == press_helper[j], ]$uncertainty,
              na.rm = TRUE
            )
          }
        }
        pressures[j] <- press_helper[j]
        type[j] <- "combined"
        pathway[j] <- "combined"
        multi_ind_helper <- data.frame(
          pressure = pressures,
          type = type,
          pathway = pathway,
          risk = risk,
          uncertainty = uncertainty
        )
      }
      multi_ind_2 <- rbind(multi_ind_2, multi_ind_helper)
    }
  }
  # Calculate multi-indicator risk per type, regardless of pathway
  for (i in 1:length(types)) {
    dat <- risk_results[risk_results$type == types[i], ]
    # get press in this subset
    press_helper <- unique(dat$pressure)
    risk <- rep(0, length(press_helper))
    pressures <- rep(NA, length(press_helper))
    type <- rep(NA, length(press_helper))
    pathway <- rep(NA, length(press_helper))
    uncertainty <- rep(NA, length(press_helper))

    # Calculate for each ind the aggregated risk
    for (j in 1:length(press_helper)) {
      # choose method
      if (method == "median") {
        risk[j] <- stats::median(dat[dat$pressure == press_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- stats::median(
            dat[dat$pressure == press_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "maximum") {
        risk[j] <- max(dat[dat$pressure == press_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- max(
            dat[dat$pressure == press_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "minimum") {
        risk[j] <- min(dat[dat$pressure == press_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- min(
            dat[dat$pressure == press_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (method == "sum") {
        risk[j] <- sum(dat[dat$pressure == press_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- sum(
            dat[dat$pressure == press_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
        risk[j] <- mean(dat[dat$pressure == press_helper[j], ]$risk,
          na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          uncertainty[j] <- mean(
            dat[dat$pressure == press_helper[j], ]$uncertainty,
            na.rm = TRUE
          )
        }
      }
      pressures[j] <- press_helper[j]
      type[j] <- types[i]
      pathway[j] <- "combined"
      multi_ind_type_helper <- data.frame(
        pressure = pressures,
        type = type,
        pathway = pathway,
        risk = risk,
        uncertainty = uncertainty
      )
    }
    if (i == 1) {
      multi_ind_type <- multi_ind_type_helper
    } else {
      multi_ind_type <- rbind(multi_ind_type, multi_ind_type_helper)
    }
  }

  multi_ind_all <- rbind(multi_ind, multi_ind_type, multi_ind_2)
  multi_ind_all <- multi_ind_all[order(multi_ind_all$pressure), ]



  # Now calculate ecosystem-wide risk for each type and pathway
  for (h in 1:length(pathways)) {
    for (i in 1:length(types)) {
      dat <- multi_press_all[multi_press_all$type == types[i] &
          multi_press_all$pathway ==  pathways[h], ]
      if (nrow(dat) == 0) {
        next
      }
      if (method == "median") {
        risk <- stats::median(dat$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) ==
            FALSE) {
          uncertainty <- stats::median(dat$uncertainty, na.rm = TRUE)
        }
      }
      if (method == "maximum") {
        risk <- max(dat$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) ==
            FALSE) {
          uncertainty <- max(dat$uncertainty, na.rm = TRUE)
        }
      }
      if (method == "minimum") {
        risk <- min(dat$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) ==
            FALSE) {
          uncertainty <- min(dat$uncertainty, na.rm = TRUE)
        }
      }
      if (method == "sum") {
        risk <- sum(dat$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) ==
            FALSE) {
          uncertainty <- sum(dat$uncertainty, na.rm = TRUE)
        }
      }
      if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
        risk <- mean(dat$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) ==
            FALSE) {
          uncertainty <- mean(dat$uncertainty, na.rm = TRUE)
        }
      }
      ecosystem_risk_helper <- data.frame(
        type = types[i], pathway = pathways[h], risk = risk,
        uncertainty = uncertainty
      )

      if (i == 1) {
        ecosystem_risk <- ecosystem_risk_helper
      } else {
        ecosystem_risk <- rbind(ecosystem_risk, ecosystem_risk_helper)
      }
    }
  }

  # Risk per type regardless of pathway
  ecosystem_risk_type <- data.frame(
    type = rep(NA, times = length(types)),
    pathway = rep("combined"),
    risk = rep(0, times = length(types)),
    uncertainty = rep(NA, times = length(types))
  )
  for (i in 1:length(types)) {
    dat <- multi_press_all[multi_press_all$type == types[i], ]
    ecosystem_risk_type$type[i] <- types[i]
    if (method == "median") {
      ecosystem_risk_type$risk[i] <- stats::median(dat$risk, na.rm = TRUE)
      if (all(is.na(dat$uncertainty)) == FALSE) {
        ecosystem_risk_type$uncertainty[i] <- stats::median(dat$uncertainty,
          na.rm = TRUE)
      }
    }
    if (method == "maximum") {
      ecosystem_risk_type$risk[i] <- max(dat$risk, na.rm = TRUE)
      if (all(is.na(dat$uncertainty)) == FALSE) {
        ecosystem_risk_type$uncertainty[i] <- max(dat$uncertainty, na.rm = TRUE)
      }
    }
    if (method == "minimum") {
      ecosystem_risk_type$risk[i] <- min(dat$risk, na.rm = TRUE)
      if (all(is.na(dat$uncertainty)) == FALSE) {
        ecosystem_risk_type$uncertainty[i] <- min(dat$uncertainty, na.rm = TRUE)
      }
    }
    if (method == "sum") {
      ecosystem_risk_type$risk[i] <- sum(dat$risk, na.rm = TRUE)
      if (all(is.na(dat$uncertainty)) == FALSE) {
        ecosystem_risk_type$uncertainty[i] <- sum(dat$uncertainty, na.rm = TRUE)
      }
    }
    if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
      ecosystem_risk_type$risk[i] <- mean(dat$risk, na.rm = TRUE)
      if (all(is.na(dat$uncertainty)) == FALSE) {
        ecosystem_risk_type$uncertainty[i] <- mean(dat$uncertainty, na.rm = TRUE)
      }
    }

  }


  # Now for each pathway and across pathways
  for (i in 1:length(pathways)) {
    dat <- multi_press_all[multi_press_all$pathway == pathways[i] &
        multi_press_all$type == "combined", ]
    temp <- data.frame(type = "combined", pathway = pathways[i], risk = NA,
      uncertainty = NA)
    if (method == "median") {
      temp$risk <- stats::median(dat$risk, na.rm = TRUE)
      if (all(is.na(dat$uncertainty)) == FALSE) {
        temp$uncertainty <- stats::median(dat$uncertainty, na.rm = TRUE)
      }
    }
    if (method == "maximum") {
      temp$risk <- max(dat$risk, na.rm = TRUE)
      if (all(is.na(dat$uncertainty)) == FALSE) {
        temp$uncertainty <- max(dat$uncertainty, na.rm = TRUE)
      }
    }
    if (method == "minimum") {
      temp$risk <- min(dat$risk, na.rm = TRUE)
      if (all(is.na(dat$uncertainty)) == FALSE) {
        temp$uncertainty <- min(dat$uncertainty, na.rm = TRUE)
      }
    }
    if (method == "sum") {
      temp$risk <- sum(dat$risk, na.rm = TRUE)
      if (all(is.na(dat$uncertainty)) == FALSE) {
        temp$uncertainty <- sum(dat$uncertainty, na.rm = TRUE)
      }
    }
    if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
      temp$risk <- mean(dat$risk, na.rm = TRUE)
      if (all(is.na(dat$uncertainty)) == FALSE) {
        temp$uncertainty <- mean(dat$uncertainty, na.rm = TRUE)
      }

    }
    if (i == 1) {
      ecosystem_risk_2 <- temp
    } else {
      ecosystem_risk_2 <- rbind(ecosystem_risk_2, temp)
    }
    if (i == max(length(pathways))) {
      dat <- multi_press_all[multi_press_all$pathway == "combined" &
          multi_press_all$type ==  "combined", ]
      temp <- data.frame(type = "combined", pathway = "combined", risk = NA,
        uncertainty = NA)
      if (method == "median") {
        temp$risk <- stats::median(dat$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          temp$uncertainty <- stats::median(dat$uncertainty, na.rm = TRUE)
        }
      }
      if (method == "maximum") {
        temp$risk <- max(dat$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          temp$uncertainty <- max(dat$uncertainty, na.rm = TRUE)
        }
      }
      if (method == "minimum") {
        temp$risk <- min(dat$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          temp$uncertainty <- min(dat$uncertainty, na.rm = TRUE)
        }
      }
      if (method == "sum") {
        temp$risk <- sum(dat$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          temp$uncertainty <- sum(dat$uncertainty, na.rm = TRUE)
        }
      }
      if (!(method %in% c("median", "maximum", "minimum", "sum"))) {
        temp$risk <- mean(dat$risk, na.rm = TRUE)
        if (all(is.na(dat$uncertainty)) == FALSE) {
          temp$uncertainty <- mean(dat$uncertainty, na.rm = TRUE)
        }
      }
      ecosystem_risk_2 <- rbind(ecosystem_risk_2, temp)

    }
  }
  ecosystem_risk_all <- rbind(ecosystem_risk, ecosystem_risk_type, ecosystem_risk_2)
  ecosystem_risk_all <- ecosystem_risk_all[order(ecosystem_risk_all$type), ]


  # End ecosystem-wide risk

  output_aggr_risk <- list(
    multi_indicator_risk = multi_ind_all,
    multi_pressure_risk = multi_press_all,
    ecosystem_risk = ecosystem_risk_all
  )

  return(output_aggr_risk)

}
