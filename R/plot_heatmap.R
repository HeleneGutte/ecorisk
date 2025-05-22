#' Generate a Heatmap Overview of Individual Risk Scores, Aggregated Risk Scores, and
#' Overall Ecosystem Risk
#'
#' The function `plot_heatmap()` creates for each effect type an aggregated plot
#' with a heatmap of the risk scores of each state indicator - pressure combination.
#' The aggregated multi pressure and multi state indicator scores are shown
#' to the left and below the heatmap. In the bottom left corner the ecosystem
#' risk is displayed. Uncertainty can be plotted as frame of the heatmap tiles on a gray scale.
#'
#' @param risk_scores output from the \code{\link{risk}} function.
#' @param aggregated_scores output from the \code{\link{aggregate_risk}} function.
#' @param order_ind character value defining the order of state indicators shown
#'        on the y-axis from top to bottom. If \code{NULL} (default), order is alphabetic.
#' @param order_press character value defining the order of pressures shown on the x-axis
#'        from left to right. If \code{NULL} (default), order is alphabetic.
#' @param pathway a character string specifying the pathway which should be used
#'        for the multi pressure and multi indicator scores. Default is "combined".
#' @param uncertainty logical, determines whether uncertainty should be plotted or not,
#'        if uncertainty scores are provided by the risk scores. Default is \code{TRUE}.
#' @param output_2_pathway_indicators Optionally. An integer value specifying whether
#'        for those state indicators that have been assessed with both pathways
#'        two plots should generated, one for each pathway (\code{2}), or only
#'        one plot, where the risk scores are averaged from both pathways
#'        (\code{1}). The default is \code{NULL}.
#' @param title a string specifying the title of the heatmap. If \code{NULL} (default),
#'        the type of the output data frame from the \code{\link{risk}} function is
#'        displayed.
#' @param risk_scale_steps integer value representing the step size for the risk
#'        scale in the legend. Can only take the value 1 (default), 2 and 5.
#' @param text_size_axis_text integer value specifying text size of axis text.
#'        If \code{NULL} (default), ggplot2 default settings are used.
#' @param text_size_axis_title integer value specifying text size of axis title.
#'        If \code{NULL} (default), ggplot2 default settings are used.
#'
#'
#' @return a list of ggplot objects, one for each type of effect.
#'
#' @seealso \code{\link{risk}}, \code{\link{aggregate_risk}} to generate result tables/output
#'          that serve here as input.
#'
#' @export
#'
#' @examples
#' ### Demo with output data from the risk() and aggregate_risk() functions
#' #   based on expert scores.
#'
#' # Using default settings for the overall risk scores and associated uncertainty
#' # scores (i.e. in this case, combined across both types)
#' p_heat <- plot_heatmap(
#'   risk_scores = ex_output_risk_expert,
#'   aggregated_scores = ex_output_aggregate_risk_expert
#' )
#' # For each type in both input datasets, a heatmap is generated
#' p_heat[[1]] # display direct effects
#' p_heat[[2]] # display direct/indirect effects
#'
#' # Hide uncertainty results and order indicators and pressures manually
#' \donttest{
#'   p_heat_mod <- plot_heatmap(
#'     risk_scores = ex_output_risk_expert,
#'     aggregated_scores = ex_output_aggregate_risk_expert,
#'     order_ind = c("phytoplankton", "herring", "cod", "seabirds"),
#'     order_press = c("temperature", "salinity", "oxygen", "nutrient",
#'       "fishing"),
#'     uncertainty = FALSE
#'   )
#'   p_heat_mod[[1]]
#' }
#'
#'
#' ### Demo with combined expert-based and model-based pathways
#'
#' combined_risk <- rbind(ex_output_risk_expert, ex_output_risk_model)
#' aggr_risk <- aggregate_risk(risk_results = combined_risk)
#'
#' # Default settings (combined type and pathway)
#' p_heat_comb <- plot_heatmap(
#'   risk_scores = combined_risk,
#'   aggregated_scores = aggr_risk
#' )
#' p_heat_comb[[1]]
#'
#'
#' ### Demo with two indicators assessed with both pathways
#' risk_model <- ex_output_risk_model[c(1, 3, 5, 7, 8, 9, 12, 14:16), ]
#' risk_model$pressure <- c(
#'  "nutrient", "temperature", "salinity", "oxygen", "fishing",   # for zooplankton
#'  "nutrient", "temperature", "salinity", "oxygen", "fishing")   # for cod
#'
#' dummy_model <- risk_model |>
#'  dplyr::mutate(indicator = dplyr::case_when(
#'    indicator == "zooplankton_mean_size" ~ "phytoplankton",
#'    .default = "cod"
#'  ))
#'
#' risk_comb <- rbind(ex_output_risk_expert, dummy_model)
#' aggr_risk_comb <- aggregate_risk(risk_results = risk_comb)
#'
#' # show results from both types and pathways individually and order the state
#' # indicators manually
#' p_heat_2_paths <- plot_heatmap(risk_scores = risk_comb,
#'                        aggregated_scores = aggr_risk_comb,
#'                        output_2_pathway_indicators = 2,
#'                        order_ind = c("phytoplankton", "herring", "cod", "seabirds"))
#' p_heat_2_paths
#'
#' # show one plot per type and average across the pathways
#' p_heat_mean_path <- plot_heatmap(risk_scores = risk_comb,
#'                        aggregated_scores = aggr_risk_comb,
#'                        output_2_pathway_indicators = 1,
#'                        order_ind = c("phytoplankton", "herring", "cod", "seabirds"))
#' p_heat_mean_path



plot_heatmap <- function(risk_scores, aggregated_scores,
  order_ind = NULL, order_press = NULL, pathway = "combined",
  uncertainty = TRUE, output_2_pathway_indicators = NULL, title = NULL,
  risk_scale_steps = 1, text_size_axis_text = NULL, text_size_axis_title = NULL) {

  ### Data input validation ----
  if (!is.null(order_ind)) {
    if (!(all(order_ind %in% unique(risk_scores$indicator)) |
        all(unique(risk_scores$indicator) %in% order_ind))) {
      stop("The indicator names in 'order_ind' do not match the indicators in
           the 'risk_scores' data frame.")
    } else {
      order_ind <- rev(order_ind)
    }
  } else {
    order_ind <- rev(sort(unique(risk_scores$indicator)))
  }

  if (!is.null(order_press)) {
    if (!(all(order_press %in% unique(risk_scores$pressure)) |
        all(unique(risk_scores$pressure) %in% order_press))) {
      stop("The pressure names in 'order_press' do not match the pressures in the
           'risk_scores' data frame.")
    }
  } else {
    order_press <- sort(unique(risk_scores$pressure))
  }

  if (!risk_scale_steps %in% c(1, 2, 5)) {
    stop("The 'risk_scale_steps' argument can only take the values 1, 2, or 5.")
  }


  if(is.null(output_2_pathway_indicators) == TRUE){
    two_paths <- 0
  }else if(output_2_pathway_indicators %in% c(1, 2)){
    two_paths <- output_2_pathway_indicators
  }else{
    stop("The 'output_2_pathway_indicators' can only take the values 1, 2, or NULL")
  }

  ### ----

  multi_press <- aggregated_scores$multi_pressure_risk
  multi_press <- multi_press[multi_press$pathway == pathway, ]
  multi_ind <- aggregated_scores$multi_indicator_risk
  multi_ind <- multi_ind[multi_ind$pathway == pathway, ]
  ecosystem_risk <- aggregated_scores$ecosystem_risk
  ecosystem_risk <- ecosystem_risk[ecosystem_risk$pathway == pathway, ]

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

  multi_press$direction <- 0
  multi_ind$direction <- 0
  for (i in 1:nrow(multi_press)) {
    if (multi_press$risk[i] <= 0) {
      multi_press$direction[i] <- "negative effect"
    } else {
      multi_press$direction[i] <- "positive effect"
    }
  }
  for (i in 1:nrow(multi_ind)) {
    if (multi_ind$risk[i] <= 0) {
      multi_ind$direction[i] <- "negative effect"
    } else {
      multi_ind$direction[i] <- "positive effect"
    }
  }

  type_names <- unique(risk_scores$type)
  heat_plots <- vector("list", length = length(type_names))

  # Define colors
  positive_colors <- colorspace::sequential_hcl(n = 13, palette = "Blues 3",
    rev = TRUE)
  negative_colors <- colorspace::sequential_hcl(n = 13, palette = "Reds 3",
    rev = TRUE)

  my_colors <- c(
    `1` = positive_colors[4], `2` = positive_colors[5],
    `3` = positive_colors[6], `4` = positive_colors[7], `5` = positive_colors[8],
    `6` = positive_colors[9], `7` = positive_colors[10], `8` = positive_colors[11],
    `9` = positive_colors[12], `10` = positive_colors[13], `0` = "gray80",
    `-1` = negative_colors[4], `-2` = negative_colors[5], `-3` = negative_colors[6],
    `-4` = negative_colors[7], `-5` = negative_colors[8], `-6` = negative_colors[9],
    `-7` = negative_colors[10], `-8` = negative_colors[11], `-9` = negative_colors[12],
    `-10` = negative_colors[13]
  )
  unc_colours <- c(`1` = "gray90", `2` = "gray55", `3` = "black")

  # Helper function to get legend
  get_legend <- function(myggplot) {
    tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }

  if(two_paths == 1){
    risk_scores <- stats::aggregate(cbind(risk, vulnerability, uncertainty) ~ indicator +
                                   pressure + type, data = risk_scores, FUN = mean)
  }else if(two_paths == 2){
    heat_plots <- vector("list", length = length(type_names)*2)

    # get pathway names
    pathway_names <- unique(risk_scores$pathway)
    for(i in 1:length(pathway_names)){
      for(j in 1:length(type_names)){
        dat <- risk_scores[risk_scores$pathway == pathway_names[i], ]
        dat <- dat[dat$type == type_names[j], ]
        if(nrow(dat) == 0){
          next
        }

        ### Workaround settings to be able to separate tiles in the
        ### heatmaps for uncertainty framing ---

        # Get number of pressures and indicators in the dataset
        n_press <- length(unique(dat$pressure))
        n_ind <- length(unique(dat$indicator))

        # Calculate x and y positions for the tiles
        x_settings <- rep(1, n_press)
        y_settings <- rep(1, n_ind)
        if (n_press > 1) {
          for (k in 2:n_press) {
            x_settings[k] <- x_settings[k - 1] + 2
          }
        }
        if (n_ind > 1) {
          for (l in 2:n_ind) {
            y_settings[l] <- y_settings[l - 1] + 2
          }
        }

        # Attach the x and y positions to the dataset, sort indicators and
        # pressures
        dat$indicator <- factor(dat$indicator, levels = order_ind)
        dat$pressure <- factor(dat$pressure, levels = order_press)

        dat <- dat[order(dat$pressure), ]
        x_labels <- order_press[order_press %in% dat$pressure]
        exclude_press <- order_press[!order_press %in% dat$pressure]
        dat$x <- rep(x_settings, times = table(dat$pressure, exclude = exclude_press))

        dat <- dat[order(dat$indicator), ]
        y_labels <- order_ind[order_ind %in% dat$indicator]
        exclude_inds <- order_ind[!order_ind %in% dat$indicator]
        dat$y <- rep(y_settings, times = table(dat$indicator, exclude = exclude_inds))

        ### End of workaround settings

        # Create heatmap with uncertainty as a frame around each tile
        if (uncertainty == TRUE) {
          # Plotting of heatmap tiles
          heat_risk <-
            ggplot2::ggplot(data = dat, ggplot2::aes(x = !!rlang::sym("x"),
                                                     y = !!rlang::sym("y"))) +
            ggplot2::geom_tile(ggplot2::aes(
              fill = cut(!!rlang::sym("risk"),
                         breaks = c(rev(-0.01:-11), 0:10), labels = -10:10),
              colour = cut(!!rlang::sym("uncertainty"),
                           breaks = seq(0, max(!!rlang::sym("uncertainty")) + 1, 1),
                           labels = seq(1, max(!!rlang::sym("uncertainty")) + 1, by = 1))),
              linewidth = 1, width = 1.5, height = 1.5) +
            ggplot2::scale_fill_manual(values = my_colors,
                                       name = "Risk", na.value = "gray80") +
            ggplot2::scale_color_manual(values = unc_colours, na.value = "white",
                                        name = "Uncertainty") +
            ggplot2::scale_x_continuous(breaks = x_settings) +
            ggplot2::scale_y_continuous(breaks = y_settings) +
            ggplot2::labs(x = "", y = "") +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              axis.text = ggplot2::element_blank(),
              axis.title = ggplot2::element_blank()
            )

        } else { # Alternatively without uncertainty as a frame

          # Plotting of heatmap tiles
          heat_risk <- ggplot2::ggplot(data = dat,
                                       mapping = ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"))) +
            ggplot2::geom_tile(ggplot2::aes(fill = cut(!!rlang::sym("risk"),
                                                       breaks = c(rev(-0.01:-11), 0:10), labels = -10:10)),
                               width = 1.5,height = 1.5) +
            ggplot2::scale_fill_manual(values = my_colors,
                                       name = "risk", na.value = "gray80") +
            ggplot2::scale_x_continuous(breaks = x_settings) +
            ggplot2::scale_y_continuous(breaks = y_settings) +
            ggplot2::labs(x = "", y = "") +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              axis.text = ggplot2::element_blank(),
              axis.title = ggplot2::element_blank()
            )
        }

        # Create multi-press scores
        dat_multi_press <- multi_press[multi_press$type == type_names[i], ]
        dat_multi_press$indicator <- factor(dat_multi_press$indicator,
                                            levels = order_ind)
        dat_multi_press <- dat_multi_press[order(dat_multi_press$indicator), ]
        # attach x and y position to data
        dat_multi_press$y <- y_settings
        dat_multi_press$x <- 1

        # Create multi-press score margin in plot
        heat_multi_press <- ggplot2::ggplot(data = dat_multi_press,
                                            mapping = ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"))) +
          ggplot2::geom_tile(ggplot2::aes(fill = cut(!!rlang::sym("risk"),
                                                     breaks = c(rev(-0.01:-11), 0:10), labels = -10:10)),
                             width = 1.5, height = 1.5) +
          ggplot2::scale_fill_manual(values = my_colors,
                                     name = "risk", na.value = "gray80") +
          ggplot2::scale_y_continuous(breaks = y_settings,
                                      labels = y_labels) +
          ggplot2::labs(y = "Indicator") +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.ticks.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(size = text_size_axis_text),
            axis.title.y = ggplot2::element_text(size = text_size_axis_title),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank(),
            legend.position = "none"
          )

        # Create multi-ind scores
        dat_multi_ind <- multi_ind[multi_ind$type == type_names[i], ]
        dat_multi_ind$pressure <- factor(dat_multi_ind$pressure,
                                         levels = order_press)
        dat_multi_ind <- dat_multi_ind[order(dat_multi_ind$pressure), ]
        dat_multi_ind$y <- 1
        dat_multi_ind$x <- x_settings

        # Create multi-ind score margin in plot
        heat_multi_ind <- ggplot2::ggplot(data = dat_multi_ind,
                                          mapping = ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"))) +
          ggplot2::geom_tile(ggplot2::aes(fill = cut(!!rlang::sym("risk"),
                                                     breaks = c(rev(-0.01:-11), 0:10), labels = -10:10)),
                             width = 1.5, height = 1.5) +
          ggplot2::scale_fill_manual(values = my_colors,
                                     name = "risk", na.value = "gray80") +
          ggplot2::scale_x_continuous(breaks = x_settings,
                                      labels = x_labels) +
          ggplot2::labs(x = "Pressure") +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.ticks.y = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1,
                                                size = text_size_axis_text),
            axis.title.x = ggplot2::element_text(size = text_size_axis_title),
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.grid.minor.x = ggplot2::element_blank(),
            legend.position = "none"
          )


        # Add ecosystem-wide risk score
        dat_sys <- ecosystem_risk[ecosystem_risk$type == type_names[i], ]
        dat_sys$round_risk <- round(dat_sys$risk, 2)
        sys_risk <- ggplot2::ggplot(data = dat_sys,
                                    mapping = ggplot2::aes(x = 0, y = 0)) +
          ggplot2::geom_label(ggplot2::aes(label = !!rlang::sym("round_risk")),
                              colour = "white", size = 5, fontface = "bold", fill = "grey40") +
          ggplot2::theme_void()


        # Create dummy plot for a full legend
        dummy_data <- data.frame(
          x = letters[1:21],
          y = c(-10:0, 1:10),
          unc = rep(1:3, length.out = 21)
        )
        if (uncertainty == TRUE) {
          dummy_plot <- ggplot2::ggplot(data = dummy_data,
                                        mapping = ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"))) +
            ggplot2::geom_tile(ggplot2::aes(
              fill = cut(!!rlang::sym("y"), breaks = c(rev(-0.01:-11), 0:10), labels = -10:10),
              colour = cut(!!rlang::sym("unc"), breaks = seq(0, max(!!rlang::sym("unc")) + 1, 1),
                           labels = seq(1, max(!!rlang::sym("unc")) + 1, by = 1))),
              width = 1.5, linewidth = 1) +
            ggplot2::scale_fill_manual(values = my_colors,
                                       name = "Risk", na.value = "gray80") +
            ggplot2::scale_color_manual(values = unc_colours,
                                        na.value = "white", name = "Uncertainty",
                                        guide = ggplot2::guide_legend(override.aes = list(fill = "white"))) +
            ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))
        } else {
          dummy_plot <- ggplot2::ggplot(data = dummy_data,
                                        mapping = ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"))) +
            ggplot2::geom_tile(ggplot2::aes(
              fill = cut(!!rlang::sym("y"), breaks = c(rev(-0.01:-11), 0:10), labels = -10:10)),
              width = 1.5, linewidth = 1) +
            ggplot2::scale_fill_manual(values = my_colors,
                                       name = "Risk", na.value = "gray80") +
            ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))
        }


        # Add legend to heatmap
        legend <- get_legend(dummy_plot)
        heat_risk <- heat_risk + ggplot2::theme(legend.position = "none")

        # Calculate from i and j the position of the plot in the plot list
        k <- (i - 1) * length(pathway_names) + j

        # Put all pieces together
        heat_plots[[k]] <- gridExtra::arrangeGrob(
          grobs = list(heat_multi_press, heat_risk, legend, heat_multi_ind, sys_risk),
          layout_matrix = rbind(c(1, 2, 3), c(5, 4, NA)),
          widths = c(1, 4, 0.75),
          heights = c(4, 1),
          top = paste0(type_names[j], " effects, ", pathway_names[i] )
        )

        heat_plots[[k]] <- ggpubr::as_ggplot(heat_plots[[k]])



      }
    }

  }

  if(two_paths == 1 | two_paths == 0){
    for (i in 1:length(type_names)) {
      dat <- risk_scores[risk_scores$type == type_names[i], ]

      ### Workaround settings to be able to separate tiles in the
      ### heatmaps for uncertainty framing ---

      # Get number of pressures and indicators in the dataset
      n_press <- length(unique(dat$pressure))
      n_ind <- length(unique(dat$indicator))

      # Calculate x and y positions for the tiles
      x_settings <- rep(1, n_press)
      y_settings <- rep(1, n_ind)
      if (n_press > 1) {
        for (j in 2:n_press) {
          x_settings[j] <- x_settings[j - 1] + 2
        }
      }
      if (n_ind > 1) {
        for (k in 2:n_ind) {
          y_settings[k] <- y_settings[k - 1] + 2
        }
      }

      # Attach the x and y positions to the dataset, sort indicators and
      # pressures
      dat$indicator <- factor(dat$indicator, levels = order_ind)
      dat$pressure <- factor(dat$pressure, levels = order_press)

      dat <- dat[order(dat$pressure), ]
      x_labels <- order_press[order_press %in% dat$pressure]
      exclude_press <- order_press[!order_press %in% dat$pressure]
      dat$x <- rep(x_settings, times = table(dat$pressure, exclude = exclude_press))

      dat <- dat[order(dat$indicator), ]
      y_labels <- order_ind[order_ind %in% dat$indicator]
      exclude_inds <- order_ind[!order_ind %in% dat$indicator]
      dat$y <- rep(y_settings, times = table(dat$indicator, exclude = exclude_inds))

      ### End of workaround settings

      # Create heatmap with uncertainty as a frame around each tile
      if (uncertainty == TRUE) {
        # Plotting of heatmap tiles
        heat_risk <-
          ggplot2::ggplot(data = dat, ggplot2::aes(x = !!rlang::sym("x"),
                                                   y = !!rlang::sym("y"))) +
          ggplot2::geom_tile(ggplot2::aes(
            fill = cut(!!rlang::sym("risk"),
                       breaks = c(rev(-0.01:-11), 0:10), labels = -10:10),
            colour = cut(!!rlang::sym("uncertainty"),
                         breaks = seq(0, max(!!rlang::sym("uncertainty")) + 1, 1),
                         labels = seq(1, max(!!rlang::sym("uncertainty")) + 1, by = 1))),
            linewidth = 1, width = 1.5, height = 1.5) +
          ggplot2::scale_fill_manual(values = my_colors,
                                     name = "Risk", na.value = "gray80") +
          ggplot2::scale_color_manual(values = unc_colours, na.value = "white",
                                      name = "Uncertainty") +
          ggplot2::scale_x_continuous(breaks = x_settings) +
          ggplot2::scale_y_continuous(breaks = y_settings) +
          ggplot2::labs(x = "", y = "") +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank()
          )

      } else { # Alternatively without uncertainty as a frame

        # Plotting of heatmap tiles
        heat_risk <- ggplot2::ggplot(data = dat,
                                     mapping = ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"))) +
          ggplot2::geom_tile(ggplot2::aes(fill = cut(!!rlang::sym("risk"),
                                                     breaks = c(rev(-0.01:-11), 0:10), labels = -10:10)),
                             width = 1.5,height = 1.5) +
          ggplot2::scale_fill_manual(values = my_colors,
                                     name = "risk", na.value = "gray80") +
          ggplot2::scale_x_continuous(breaks = x_settings) +
          ggplot2::scale_y_continuous(breaks = y_settings) +
          ggplot2::labs(x = "", y = "") +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            axis.text = ggplot2::element_blank(),
            axis.title = ggplot2::element_blank()
          )
      }

      # Create multi-press scores
      dat_multi_press <- multi_press[multi_press$type == type_names[i], ]
      dat_multi_press$indicator <- factor(dat_multi_press$indicator,
                                          levels = order_ind)
      dat_multi_press <- dat_multi_press[order(dat_multi_press$indicator), ]
      # attach x and y position to data
      dat_multi_press$y <- y_settings
      dat_multi_press$x <- 1

      # Create multi-press score margin in plot
      heat_multi_press <- ggplot2::ggplot(data = dat_multi_press,
                                          mapping = ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"))) +
        ggplot2::geom_tile(ggplot2::aes(fill = cut(!!rlang::sym("risk"),
                                                   breaks = c(rev(-0.01:-11), 0:10), labels = -10:10)),
                           width = 1.5, height = 1.5) +
        ggplot2::scale_fill_manual(values = my_colors,
                                   name = "risk", na.value = "gray80") +
        ggplot2::scale_y_continuous(breaks = y_settings,
                                    labels = y_labels) +
        ggplot2::labs(y = "Indicator") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(size = text_size_axis_text),
          axis.title.y = ggplot2::element_text(size = text_size_axis_title),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
          panel.grid.minor.y = ggplot2::element_blank(),
          legend.position = "none"
        )

      # Create multi-ind scores
      dat_multi_ind <- multi_ind[multi_ind$type == type_names[i], ]
      dat_multi_ind$pressure <- factor(dat_multi_ind$pressure,
                                       levels = order_press)
      dat_multi_ind <- dat_multi_ind[order(dat_multi_ind$pressure), ]
      dat_multi_ind$y <- 1
      dat_multi_ind$x <- x_settings

      # Create multi-ind score margin in plot
      heat_multi_ind <- ggplot2::ggplot(data = dat_multi_ind,
                                        mapping = ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"))) +
        ggplot2::geom_tile(ggplot2::aes(fill = cut(!!rlang::sym("risk"),
                                                   breaks = c(rev(-0.01:-11), 0:10), labels = -10:10)),
                           width = 1.5, height = 1.5) +
        ggplot2::scale_fill_manual(values = my_colors,
                                   name = "risk", na.value = "gray80") +
        ggplot2::scale_x_continuous(breaks = x_settings,
                                    labels = x_labels) +
        ggplot2::labs(x = "Pressure") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.ticks.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1,
                                              size = text_size_axis_text),
          axis.title.x = ggplot2::element_text(size = text_size_axis_title),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor.y = ggplot2::element_blank(),
          panel.grid.minor.x = ggplot2::element_blank(),
          legend.position = "none"
        )


      # Add ecosystem-wide risk score
      dat_sys <- ecosystem_risk[ecosystem_risk$type == type_names[i], ]
      dat_sys$round_risk <- round(dat_sys$risk, 2)
      sys_risk <- ggplot2::ggplot(data = dat_sys,
                                  mapping = ggplot2::aes(x = 0, y = 0)) +
        ggplot2::geom_label(ggplot2::aes(label = !!rlang::sym("round_risk")),
                            colour = "white", size = 5, fontface = "bold", fill = "grey40") +
        ggplot2::theme_void()


      # Create dummy plot for a full legend
      dummy_data <- data.frame(
        x = letters[1:21],
        y = c(-10:0, 1:10),
        unc = rep(1:3, length.out = 21)
      )
      if (uncertainty == TRUE) {
        dummy_plot <- ggplot2::ggplot(data = dummy_data,
                                      mapping = ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"))) +
          ggplot2::geom_tile(ggplot2::aes(
            fill = cut(!!rlang::sym("y"), breaks = c(rev(-0.01:-11), 0:10), labels = -10:10),
            colour = cut(!!rlang::sym("unc"), breaks = seq(0, max(!!rlang::sym("unc")) + 1, 1),
                         labels = seq(1, max(!!rlang::sym("unc")) + 1, by = 1))),
            width = 1.5, linewidth = 1) +
          ggplot2::scale_fill_manual(values = my_colors,
                                     name = "Risk", na.value = "gray80") +
          ggplot2::scale_color_manual(values = unc_colours,
                                      na.value = "white", name = "Uncertainty",
                                      guide = ggplot2::guide_legend(override.aes = list(fill = "white"))) +
          ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))
      } else {
        dummy_plot <- ggplot2::ggplot(data = dummy_data,
                                      mapping = ggplot2::aes(x = !!rlang::sym("x"), y = !!rlang::sym("y"))) +
          ggplot2::geom_tile(ggplot2::aes(
            fill = cut(!!rlang::sym("y"), breaks = c(rev(-0.01:-11), 0:10), labels = -10:10)),
            width = 1.5, linewidth = 1) +
          ggplot2::scale_fill_manual(values = my_colors,
                                     name = "Risk", na.value = "gray80") +
          ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1))
      }


      # Add legend to heatmap
      legend <- get_legend(dummy_plot)
      heat_risk <- heat_risk + ggplot2::theme(legend.position = "none")

      # Put all pieces together
      heat_plots[[i]] <- gridExtra::arrangeGrob(
        grobs = list(heat_multi_press, heat_risk, legend, heat_multi_ind, sys_risk),
        layout_matrix = rbind(c(1, 2, 3), c(5, 4, NA)),
        widths = c(1, 4, 0.75),
        heights = c(4, 1),
        top = paste0(type_names[i], " effects")
      )
      heat_plots[[i]] <- ggpubr::as_ggplot(heat_plots[[i]])

    }
  }

  # delete empty plots from list
  clean_list <- Filter(Negate(is.null), heat_plots)

  return(clean_list)

}

