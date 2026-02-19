#' Generate Radar Charts Displaying Pressure-Specific and Overall Risks
#' for Each State Indicator
#'
#' The `plot_radar()`function creates per indicator a ggplot object.
#' The plot shows the risks of all types and effect directions. The associated
#' uncertainty can optionally be displayed. In the middle the plot displays
#' the multi pressure score of a chosen effect type.
#'
#' @param risk_scores output from the \code{\link{risk}} function.
#' @param aggregated_scores output from the \code{\link{aggregate_risk}} function.
#' @param type character string, type used for the multi-pressure score, can be
#'        any type that has been evaluated. The default is \code{combined}.
#' @param pathway character string specifying the multi-pressure score, should
#'        be plotted for each pathway individual \code{individual} or as a
#'        combined score \code{combined}. The default is \code{combined}.
#'
#' @return a list of ggplot2 objects one for each indicator, the order depends
#'  on the order in the risk_score data set. Each plot shows the risks
#'  for one state indicator for each pressure and type of assessment. In the
#'  center of the plot the multi-pressure score (either in blue or in red)
#'  and the associated aggregated uncertainty (in black) is shown.
#'  If one indicator has been assessed with both pathways, one plot is generated
#'  for each of the pathways.
#'  The uncertainty of each individual risk is shown as a ring around the risks in grey.
#'
#'
#' @seealso \code{\link{risk}}, \code{\link{aggregate_risk}} to generate result tables/output
#'          that serve here as input
#'
#' @export
#'
#' @examples
#' ### Demo with output data from the risk() and aggregate_risk() functions
#' #   based on expert scores
#' # The examples can run for a longer time, thus they are in dontrun{}.
#' # Using default settings for the indicator-specific overall risk score (coloured value)
#' # and associated uncertainty score (black value) (i.e., combined across both types)
#' \dontrun{
#' p_radar <- plot_radar(
#'   risk_scores = ex_output_risk_expert,
#'   aggregated_scores = ex_output_aggregate_risk_expert
#' )
#' p_radar[[1]] # display radar chart for first indicator
#'
#' # Show overall risk score based on direct effects only
#' p_radar_direct <- plot_radar(
#'   risk_scores = ex_output_risk_expert,
#'   aggregated_scores = ex_output_aggregate_risk_expert,
#'   type = "direct"
#' )
#' p_radar_direct[[1]]
#'
#'
#'
#' ### Demo with combined expert-based and model-based pathways
#' 
#' combined_risk <- rbind(ex_output_risk_expert, ex_output_risk_model)
#' aggr_risk <- aggregate_risk(risk_results = combined_risk)
#' 
#'
#' # Default settings (combined type and pathway)
#' 
#' p_radar_comb <- plot_radar(
#'   risk_scores = combined_risk,
#'   aggregated_scores = aggr_risk
#' )
#' p_radar_comb[[1]]
#' 
#'
#' # Show overall risk score based on direct/indirect effects only for both
#' # pathways combined
#' 
#' p_radar_comb_dindi <- plot_radar(
#'   risk_scores = ex_output_risk_expert,
#'   aggregated_scores = ex_output_aggregate_risk_expert,
#'   type = "direct_indirect"
#' )
#' p_radar_comb_dindi[[1]]
#' }

plot_radar <- function(risk_scores, aggregated_scores,
  type = "combined", pathway = "combined") {

  ### check if suggested packages are installed: 
  if(!requireNamespace(c("geomtextpath"), quietly = TRUE)){
    stop(
      "Package \"geomtextpath\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if(!requireNamespace(c("ggnewscale"), quietly = TRUE)){
    stop(
      "Package \"ggnewscale\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if(!requireNamespace(c("ggplot2"), quietly = TRUE)){
    stop(
      "Package \"ggplot2\" must be installed to use this function.",
      call. = FALSE
    )
  }

  
  multi_press <- aggregated_scores$multi_pressure_risk
  multi_ind <- aggregated_scores$multi_indicator_risk
  system_risk <- aggregated_scores$ecosystem_risk

  indicators <- unique(risk_scores$indicator)
  pathways <- unique(risk_scores$pathway)
  plot_range <- c(0, 10)

  # test how many state indicators have been assessed with both pathways and
  # add plots to the list accordingly
  unique_pairs <- unique(risk_scores[c("indicator", "pathway")])
  # Count unique combinations for each indicator
  n_pathway_ind  <- as.data.frame(table(unique_pairs$indicator))
  names(n_pathway_ind) <- c("indicator", "n_unique_combinations")
  n_ind_2_pathways <- nrow(n_pathway_ind[n_pathway_ind$n_unique_combinations == 2, ])

  radar_plots <- vector("list", length = length(indicators) + n_ind_2_pathways)

  # Get maximum range for plot to standardize for all indicators
  max_risk <- abs(max(plot_range))

  # Check for different effect directions -> use transparency to
  # highlight them
  risk_scores$direction <- 0
  for (i in 1:nrow(risk_scores)) {
    if (is.na(risk_scores$risk[i]) == TRUE) {
      risk_scores$direction[i] <- 0
    } else if (risk_scores$risk[i] <= 0) {
      risk_scores$direction[i] <- "negative effect"
    } else {
      risk_scores$direction[i] <- "positive effect"
    }
  }
  risk_scores$direction <- as.factor(risk_scores$direction)
  risk_scores$abs_risk <- abs(risk_scores$risk)


  if (all(is.na(risk_scores$uncertainty)) == FALSE) { # Radar plots with uncertainties
    for (i in 1:length(indicators)) {
      for(j in 1:length(pathways)){

        dat <- risk_scores[risk_scores$indicator == indicators[i], ]
        dat <- dat[dat$pathway == pathways[j], ]
        dat <- dat[is.na(dat$risk) == FALSE, ]
        multi_press_helper <- multi_press[multi_press$indicator == indicators[i], ]
        status_helper <- unique(dat$status)
        pathway_helper <- unique(dat$pathway)

        if(nrow(dat) == 0){
          next
        }

        dat$uncertainty_cut <- cut(
          dat$uncertainty,
          breaks = seq(0, max(dat$uncertainty, na.rm = TRUE) + 1, 1),
          labels = seq(1, max(dat$uncertainty, na.rm = TRUE) + 1, by = 1)
        )
        unc_colours <- c(`1` = "gray90", `2` = "gray55", `3` = "black")

        if (length(multi_press_helper$risk[multi_press_helper$type == type &
                                           multi_press_helper$pathway == pathway]) != 0) {
          if (multi_press_helper$risk[multi_press_helper$type == type &
                                      multi_press_helper$pathway == pathway] < 0) {
            col_multi_press <- "red"
          } else {
            col_multi_press <- "blue"
          }
        } else {
          multi_press_helper$risk <- 0
          col_multi_press <- "black"
          warning(paste("No multi-pressure score detected for indicator:",
                        indicators[i],
                        "\nCheck if specified type and pathway are available for this indicator.
          Score is set to 0."))
        }

        # Calculate from i and j the position of the plot in the plot list
        k <- (i - 1) * length(pathways) + j


        # Create a radar plot for each indicator
        radar_plots[[k]] <- ggplot2::ggplot(data = dat) +

          # Uncertainty
          ggplot2::geom_col(ggplot2::aes(x = !!rlang::sym("pressure"), y = max_risk + 4,
                                         fill = !!rlang::sym("uncertainty_cut")), position = "dodge") +
          ggplot2::scale_y_continuous(limits = c(plot_range[1] - 5, plot_range[2] + 4),
                                      breaks = c(0, plot_range[2]/4*1, plot_range[2]/4*2, plot_range[2]/4*3, plot_range[2]/4*4)) +
          ggplot2::scale_fill_manual(breaks = c(1:(max(dat$uncertainty) +1)),
                                     values = unc_colours, na.value = "white") +

          # Add "blank" to distinguish uncertainty and risk columns
          ggplot2::geom_col(ggplot2::aes(x = !!rlang::sym("pressure"), y = max_risk + 3),
                            fill = "white", position = "dodge") +
          ggplot2::labs(fill = "Uncertainty") +

          # Add reference lines
          ggplot2::geom_hline(yintercept = 0,  colour = "grey80") +
          ggplot2::geom_hline(yintercept = plot_range[2]/4*1,  colour = "grey80") +
          ggplot2::geom_hline(yintercept = plot_range[2]/4*2, colour = "grey80") +
          ggplot2::geom_hline(yintercept = plot_range[2]/4*3,  colour = "grey80") +
          ggplot2::geom_hline(yintercept = plot_range[2]/4*4, colour = "grey80") +

          # Risk: set new fill scale
          ggnewscale::new_scale("fill") +
          ggplot2::geom_col(ggplot2::aes(x =!!rlang::sym("pressure"),
                                         y = !!rlang::sym("abs_risk"), fill = !!rlang::sym("direction"),
                                         alpha = !!rlang::sym("type")),
                            position = "dodge") +
          # Colour risk, NA same as background blank
          ggplot2::scale_fill_manual(values = c("red", "blue"), na.value = "white",
                                     guide = ggplot2::guide_legend(order = 1)) +
          {if (length(unique(dat$type)) == 1)
            ggplot2::scale_alpha_manual(values = 1, guide = "none")} +
          {if (length(unique(dat$type)) > 1)
            ggplot2::scale_alpha_manual(values = c(seq(1, 0.1,
                                                       length.out = length(unique(dat$type)))), na.translate = FALSE,
                                        guide = ggplot2::guide_legend(order = 2, override.aes = list(fill = "red")))} +

          # Add multi-pressure risk, colour refers to the status
          ggplot2::annotate("point", x = 0, y = plot_range[1] - 5, size = 25,
                            colour = "white") +
          ggplot2::annotate("text", x = 0, y = plot_range[1] - 2.5,
                            label = round(multi_press_helper[multi_press_helper$type == type &
                                                               multi_press_helper$pathway == pathway, ]$risk, 2),
                            colour = col_multi_press, size = 9) +
          ggplot2::annotate("text", x = 0, y = plot_range[1] - 5,
                            label = round(multi_press_helper[multi_press_helper$type == type &
                                                               multi_press_helper$pathway == pathway, ]$uncertainty, 2),
                            colour = "black", size = 9) +
          ggplot2::annotate("text", x = 0.25,
                            y = c(0, plot_range[2]/4*1, plot_range[2]/4*2, plot_range[2]/4*3, plot_range[2]/4*4),
                            label = c(0, paste0("\u00B1 ", plot_range[2]/4*1), paste0("\u00B1 ", plot_range[2]/4*2),
                                      paste0("\u00B1 ", plot_range[2]/4*3), paste0("\u00B1 ", plot_range[2]/4*4)),
                            color = "grey40") +
          # Convert it to a pie
          # ggplot2::coord_polar() +
          geomtextpath::coord_curvedpolar() +

          # Modify the theme
          ggplot2::labs(
            title = paste0(indicators[i], ", ", status_helper, ", ", pathway_helper),
            x = "", y = "",
            fill = "Direction of effect", alpha = "Type of effect") +
          ggplot2::theme(
            axis.ticks = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_line(linewidth = 0.25,
                                                       colour = "grey80")
          )
      }


    }
  } else {  # Alternatively, radar plots without uncertainties

    for (i in 1:length(indicators)) {
      for(j in 1:length(pathways)){

        dat <- risk_scores[risk_scores$indicator == indicators[i], ]
        dat <- dat[dat$pathway == pathways[j], ]
        dat <- dat[is.na(dat$risk) == FALSE, ]
        multi_press_helper <- multi_press[multi_press$indicator == indicators[i], ]
        status_helper <- unique(dat$status)
        pathway_helper <- unique(dat$pathway)

        if(nrow(dat) == 0){
          next
        }

        if (length(multi_press_helper[multi_press_helper$type == type &
          multi_press_helper$pathway == pathway, ]$risk) != 0) {
          if (multi_press_helper[multi_press_helper$type == type &
                                 multi_press_helper$pathway == pathway, ]$risk < 0) {
            col_multi_press <- "red"
            } else {
            col_multi_press <- "blue"
            }
          } else {
            multi_press_helper$risk <- 0
            col_multi_press <- "black"
            warning(paste("No multi pressure score detected for indicator ",
              indicators[i],
              "\nCheck if specified type and pathway are available for this indicator.
              Score is set to 0."))
          }

        # Calculate from i and j the position of the plot in the plot list
        k <- (i - 1) * length(pathways) + j

        # Create a radar plot for each indicator
        radar_plots[[k]] <- ggplot2::ggplot(data = dat) +
          ggplot2::geom_col(ggplot2::aes(x = !!rlang::sym("pressure"),
            y = max_risk + 4, fill = NA)) +
          ggplot2::scale_y_continuous(
            limits = c(plot_range[1] - 3, plot_range[2] + 4),
            breaks = c(0, plot_range[2]/4*1, plot_range[2]/4*2, plot_range[2]/4*3,
              plot_range[2]/4*4)) +

          # Add "blank" to distinguish uncertainty and risk columns
          ggplot2::geom_col(ggplot2::aes(x = !!rlang::sym("pressure"), y = max_risk + 4),
            fill = "white", position = "dodge") +
          ggplot2::labs(fill = "Uncertainty") +

          # Add reference lines
          ggplot2::geom_hline(yintercept = 0,  colour = "grey80") +
          ggplot2::geom_hline(yintercept = plot_range[2]/4*1,  colour = "grey80") +
          ggplot2::geom_hline(yintercept = plot_range[2]/4*2, colour = "grey80") +
          ggplot2::geom_hline(yintercept = plot_range[2]/4*3,  colour = "grey80") +
          ggplot2::geom_hline(yintercept = plot_range[2]/4*4, colour = "grey80") +

          # Risk
          ggplot2::geom_col(ggplot2::aes(x =!!rlang::sym("pressure"),
            y = !!rlang::sym("abs_risk"), fill = !!rlang::sym("direction"),
            alpha = !!rlang::sym("type")),
            position = "dodge") +
          ggplot2::scale_alpha_manual(values = c(seq(1, 0.1,
            length.out = length(unique(dat$type)))), na.translate = FALSE,
            guide = ggplot2::guide_legend(order = 2, override.aes = list(fill = "red"))) +

          # Colour risk, NA same as background blank
          ggplot2::scale_fill_manual(values = c("red", "blue"), na.value = "white",
            guide = ggplot2::guide_legend(order = 1)) +

          # Add multi-pressure risk, colour refers to the status
          ggplot2::annotate("point", x = 0, y = plot_range[1] - 3, size = 15,
            colour = "white") +
          ggplot2::annotate("text", x = 0, y = plot_range[1] - 3,
            label = round(multi_press_helper[multi_press_helper$type == type &
                multi_press_helper$pathway == pathway, ]$risk, 2),
            colour = col_multi_press, size = 4) +

          # Convert it to a pie
          # ggplot2::coord_polar() +
          geomtextpath::coord_curvedpolar() +

          # Modify the theme
          ggplot2::labs(title = paste0(indicators[i], ", ", status_helper, ", ", pathway_helper),
            x = "", y = "",
            fill = "Direction of effect", alpha = "Type of effect") +
          ggplot2::theme(
            axis.ticks = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_line(linewidth = 0.25,
              colour = "grey80")
          )
      }

    }  # end of for loop
  }  # end of else statement

  # delete empty plots from list
  clean_list <- Filter(Negate(is.null), radar_plots)


  return(clean_list)

}

