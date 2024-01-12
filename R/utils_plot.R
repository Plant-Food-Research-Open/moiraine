#' Plot correlation matrix
#'
#' Plots a correlation matrix using the corrplot package.
#'
#' @param cormat A correlation matrix.
#' @param ... Other arguments to be passed to the [corrplot::corrplot()]
#'   function.
#' @returns A correlation plot
#' @export
plot_correlation_matrix <- function(cormat, ...) {
  if (!requireNamespace("corrplot", quietly = TRUE)) {
    stop(
      "Package \"corrplot\" must be installed to use this function.",
      call. = FALSE
    )
  }
  if (requireNamespace("grDevices", quietly = TRUE)) {
    colours <- grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(9, "RdBu"))
    )(100)
  } else {
    colours <- NULL
  }

  corrplot::corrplot(
    cormat,
    method = "circle",
    addCoef.col = "black",
    col = colours,
    type = "upper",
    diag = FALSE,
    order = "original",
    tl.col = "black",
    tl.srt = 45,
    ...
  )
}


#' Plots a full correlation matrix (corrplot-style)
#'
#' Generates a plot of a correlation matrix in the style of the corrplot
#' package.
#'
#' @param mat Correlation matrix to plot.
#' @param rows_title Character, title for rows. Default value is `NULL`.
#' @param cols_title Character, title for cols. Default value is `NULL`.
#' @param title Character, title of the plot. Default value is `NULL`.
#' @param show_cor Logical, should the correlation values be added to the plot?
#'   Default value is `TRUE`.
#' @param min_show_cor Numeric, the minimum value below which correlation
#'   coefficients values are not added to the plot (i.e. only a circle will
#'   appear for these values but no text). Ignored if `show_cor` is `FALSE`.
#'   Default value is 0.2.
#' @param round_cor Integer, how many decimal places to show for the correlation
#'   coefficients. Ignored if `show_cor` is `FALSE`. Default value is 2.
#' @returns a ggplot.
#' @export
plot_correlation_matrix_full <- function(mat,
                                         rows_title = NULL,
                                         cols_title = NULL,
                                         title = NULL,
                                         show_cor = TRUE,
                                         min_show_cor = 0.2,
                                         round_cor = 2) {
  ## For devtools::check
  from <- to <- from_num <- to_num <- rad <- label <- value <- NULL

  toplot <- mat |>
    tibble::as_tibble(rownames = "from") |>
    tidyr::pivot_longer(
      cols = -from,
      names_to = "to",
      values_to = "value"
    ) |>
    dplyr::mutate(
      rad = sqrt(0.25 * abs(value)),
      from = factor(from, levels = rev(rownames(mat))),
      to = factor(to, levels = colnames(mat)),
      from_num = as.numeric(from),
      to_num = as.numeric(to)
    )

  if (show_cor) {
    toplot <- toplot |>
      dplyr::mutate(
        label = round(value, round_cor),
        label = dplyr::case_when(
          abs(value) >= min_show_cor ~ paste0(label),
          TRUE ~ ""
        )
      )
  }

  p <- toplot |>
    ggplot2::ggplot() +
    ggforce::geom_circle(
      aes(
        x0 = to_num,
        y0 = from_num,
        r = rad,
        fill = value
      ),
      na.rm = TRUE,
      colour = "white",
      linewidth = 0.1
    )

  if (show_cor) {
    p <- p + ggplot2::geom_text(
      aes(
        x = to_num,
        y = from_num,
        label = label
      ),
      hjust = 0.5,
      vjust = 0.5
    )
  }

  p <- p +
    ggplot2::scale_x_continuous(
      breaks = seq_len(ncol(mat)),
      labels = colnames(mat),
      limits = c(0.5, ncol(mat) + 0.5),
      position = "top",
      expand = ggplot2::expansion()
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(nrow(mat)),
      labels = rev(rownames(mat)),
      limits = c(0.5, nrow(mat) + 0.5),
      expand = ggplot2::expansion()
    ) +
    ggplot2::scale_fill_distiller(
      palette = "RdBu",
      limits = c(-1, 1),
      direction = -1,
      guide = ggplot2::guide_colourbar(
        title.vjust = 0.8,
        ticks.colour = "black"
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(hjust = 0.5, size = 12),
      axis.title.y = ggplot2::element_text(hjust = 0.5, size = 12),
      axis.text.x.top = ggplot2::element_text(
        angle = 45,
        hjust = 0,
        vjust = 0,
        size = 8
      ),
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_line(colour = "gray80"),
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = title,
      x = cols_title,
      y = rows_title,
      fill = "Correlation"
    )

  return(p)
}

#' ggpairs plot with custom colours
#'
#' Creates a ggpairs plot (see [GGally::ggpairs()]) where the colours and shapes
#' can differ between the upper triangle, lower triangle and diagonal plots.
#'
#' @param toplot Tibble in wide format, with observations as rows and variables
#'   as columns.
#' @param vars Character vector, names of the columns in `toplot` that
#'   correspond to variables to be used for the plot matrix.
#' @param colour_upper Character, name of column in `toplot` to use for
#'   colouring observations in the upper triangle plots. Default value is
#'   `NULL`.
#' @param colour_diag Character, name of column in `toplot` to use for colouring
#'   observations in the diagonal plots. By default, will follow `colour_upper`.
#' @param colour_lower Character, name of column in `toplot` to use for
#'   colouring observations in the lower triangle plots. By default, will follow
#'   `colour_upper`.
#' @param shape_upper Character, name of column in `toplot` to use for shaping
#'   observations in the upper triangle plots. Default value is `NULL`.
#' @param shape_lower Character, name of column in `toplot` to use for shaping
#'   observations in the lower triangle plots. By default, will follow
#'   `shape_upper`.
#' @param scale_colour_upper ggplot2 colour scale to use for the upper triangle
#'   plots. Default value is `NULL` (if `colour_upper` is not `NULL`, will use
#'   ggplot2 default colour scales).
#' @param scale_colour_diag ggplot2 colour scale to use for the diagonal plots.
#'   If `NULL` (default), the colour scale used for the upper triangle plots
#'   will be used if `colour_diag` is equal to `colour_upper`; or the colour
#'   scale used for the lower triangle plots will be used if `colour_diag` is
#'   equal to `colour_lower`.
#' @param scale_colour_lower ggplot2 colour scale to use for the lower triangle
#'   plots. If `NULL` (default), the colour scale used for the upper triangle
#'   plots will be used.
#' @param scale_shape_upper ggplot2 shape scale to use for the upper triangle
#'   plots. Default value is `NULL` (if `shape_upper` is not `NULL`, will use
#'   ggplot2 default shape scale).
#' @param scale_shape_lower ggplot2 shape scale to use for the lower triangle
#'   plots. If `NULL` (default), the shape scale used for the upper triangle
#'   plots will be used.
#' @param title Character, title of the plot. Default value is `NULL` (no title
#'   added to the plot).
#' @param point_size Numeric, the size of points (in pt) in the plot. Default
#'   value is 1.5.
#' @returns A `ggmatrix` plot.
#' @examples
#' \dontrun{
#' library(palmerpenguins)
#' library(ggplot2)
#' data("penguins")
#'
#' vars <- c(
#'   "bill_length_mm",
#'   "bill_depth_mm",
#'   "flipper_length_mm"
#' )
#'
#' toplot <- penguins |>
#'   dplyr::filter(!is.na(bill_length_mm))
#'
#' # simple scatterplots of the variables
#' ggpairs_custom(toplot, vars)
#'
#' # colouring points by species, using custom colour palette
#' ggpairs_custom(
#'   toplot,
#'   vars,
#'   colour_upper = "species",
#'   scale_colour_upper = scale_colour_brewer(palette = "Set1")
#' )
#'
#' # now adding the sex variable as shape of the points
#' ggpairs_custom(
#'   toplot,
#'   vars,
#'   colour_upper = "species",
#'   scale_colour_upper = scale_colour_brewer(palette = "Set1"),
#'   shape_upper = "sex"
#' )
#'
#' # using the lower plots to show the island as colour
#' ggpairs_custom(
#'   toplot,
#'   vars,
#'   colour_upper = "species",
#'   scale_colour_upper = scale_colour_brewer(palette = "Set1"),
#'   shape_upper = "sex",
#'   colour_lower = "island",
#'   scale_colour_lower = scale_colour_viridis_d()
#' )
#'
#' # showing species and sex in upper plots, body mass and island
#' # in lower plots
#' ggpairs_custom(
#'   toplot,
#'   vars,
#'   colour_upper = "species",
#'   scale_colour_upper = scale_colour_brewer(palette = "Set1"),
#'   shape_upper = "sex",
#'   shape_lower = "island",
#'   colour_lower = "body_mass_g",
#'   scale_colour_lower = scale_colour_viridis_c(option = "plasma")
#' )
#'
#' # same as above, but the diagonal plots show density per year
#' ggpairs_custom(
#'   toplot,
#'   vars,
#'   colour_upper = "species",
#'   scale_colour_upper = scale_colour_brewer(palette = "Set1"),
#'   shape_upper = "sex",
#'   shape_lower = "island",
#'   colour_lower = "body_mass_g",
#'   scale_colour_lower = scale_colour_viridis_c(option = "plasma"),
#'   colour_diag = "year"
#' )
#'
#' # common legend if the diagonal follows the colour of the lower plots
#' ggpairs_custom(
#'   toplot,
#'   vars,
#'   colour_upper = "species",
#'   scale_colour_upper = scale_colour_brewer(palette = "Set1"),
#'   colour_lower = "island",
#'   scale_colour_lower = scale_colour_brewer(palette = "Accent"),
#'   colour_diag = "island"
#' )
#' }
#' @export
ggpairs_custom <- function(toplot,
                           vars,
                           colour_upper = NULL,
                           colour_diag = colour_upper,
                           colour_lower = colour_upper,
                           shape_upper = NULL,
                           shape_lower = shape_upper,
                           scale_colour_upper = NULL,
                           scale_colour_diag = NULL,
                           scale_colour_lower = NULL,
                           scale_shape_upper = NULL,
                           scale_shape_lower = NULL,
                           title = NULL,
                           point_size = 1.5) {
  legend_val <- NULL

  ## Copying scales to avoid changing stuff outside of the function
  scale_colour_upper <- .make_copy_scales(scale_colour_upper)
  scale_colour_diag <- .make_copy_scales(scale_colour_diag)
  scale_colour_lower <- .make_copy_scales(scale_colour_lower)
  scale_shape_upper <- .make_copy_scales(scale_shape_upper)
  scale_shape_lower <- .make_copy_scales(scale_shape_lower)

  ## Checking for continuous variables
  if (!is.null(shape_upper)) {
    if (is.numeric(toplot[[shape_upper]])) {
      shape_upper <- NULL
      message("Cannot use 'shape_upper' for numeric variable.")
    }
  }

  if (!is.null(shape_lower)) {
    if (is.numeric(toplot[[shape_lower]])) {
      shape_lower <- NULL
      message("Cannot use 'shape_lower' for numeric variable.")
    }
  }

  if (!is.null(colour_diag)) {
    if (is.double(toplot[[colour_diag]])) {
      colour_diag <- NULL
      message("Cannot colour diagonal according to continuous variable.")
    }
  }

  any_vars <- c(
    "colour_upper", "colour_diag", "colour_lower",
    "shape_upper", "shape_lower"
  ) |>
    purrr::map_lgl(~ !is.null(get(.x))) |>
    any()

  ## Checking for the presence of custom aesthetics
  if (any_vars) legend_val <- 2

  ## ============== ##
  ## Upper triangle ##
  ## ============== ##

  if (!is.null(colour_upper) || !is.null(shape_upper)) {
    ## Need to set up default scales so they can be duplicated
    ## in the diagonal and lower triangle
    if (!is.null(colour_upper) && is.null(scale_colour_upper)) {
      if (is.numeric(toplot[[colour_upper]])) {
        scale_colour_upper <- ggplot2::scale_colour_continuous()
      } else {
        scale_colour_upper <- ggplot2::scale_colour_discrete()
      }
    }

    if (!is.null(shape_upper) && is.null(scale_shape_upper)) {
      scale_shape_upper <- ggplot2::scale_shape(na.value = 4)
    }

    ggplot_aes_upper <- .make_ggplot_aes(
      colour = colour_upper,
      shape = shape_upper
    )

    ggally_upper <- function(data, mapping, ...) {
      p <- ggplot2::ggplot(
        data,
        mapping
      ) +
        ggplot2::geom_point(
          mapping = eval(ggplot_aes_upper),
          size = point_size
        ) +
        scale_colour_upper +
        scale_shape_upper

      return(p)
    }
  } else {
    ggally_upper <- "points"
  }

  ## ============== ##
  ## Lower triangle ##
  ## ============== ##

  if (!is.null(colour_lower) || !is.null(shape_lower)) {
    same_colour_lower <- all(
      !is.null(colour_lower),
      colour_lower == colour_upper,
      is.null(scale_colour_lower),
      !is.null(scale_colour_upper)
    )

    same_shape_lower <- all(
      !is.null(shape_lower),
      shape_lower == shape_upper,
      is.null(scale_shape_lower),
      !is.null(scale_shape_upper)
    )

    ## If the colour is the same as for upper plots and we are not given a scale
    if (same_colour_lower) {
      scale_colour_lower <- ggplot2::ggproto(NULL, scale_colour_upper)
    }

    ## If we still don't have a scale
    ## (so it's not the same colour as upper plots)
    if (!is.null(colour_lower) && is.null(scale_colour_lower)) {
      if (is.numeric(toplot[[colour_lower]])) {
        scale_colour_lower <- ggplot2::scale_colour_continuous()
      } else {
        scale_colour_lower <- ggplot2::scale_colour_discrete()
      }
    }

    ## If the shape is the same as for lower plots and we are not given a scale
    if (same_shape_lower) {
      scale_shape_lower <- ggplot2::ggproto(NULL, scale_shape_upper)
    }

    ## If we still don't have a scale
    ## (so it's not the same shape as upper plots)
    if (!is.null(shape_lower) && is.null(scale_shape_lower)) {
      scale_shape_lower <- ggplot2::scale_shape(na.value = 4)
    }

    ggplot_aes_lower <- .make_ggplot_aes(
      colour = colour_lower,
      shape = shape_lower
    )

    ggally_lower <- function(data, mapping, ...) {
      p <- ggplot2::ggplot(
        data,
        mapping
      ) +
        ggplot2::geom_point(
          mapping = eval(ggplot_aes_lower),
          size = point_size
        ) +
        scale_colour_lower +
        scale_shape_lower

      return(p)
    }
  } else {
    same_colour_lower <- is.null(colour_upper)
    same_shape_lower <- is.null(shape_upper)
    ggally_lower <- "points"
  }

  ## ============== ##
  ##     Diagonal   ##
  ## ============== ##

  if (!is.null(colour_diag)) {
    same_colour_diag <- all(
      colour_diag == colour_upper,
      is.null(scale_colour_diag),
      !is.null(scale_colour_upper)
    )

    ## If we have the same colour as upper plots and no scale
    if (same_colour_diag) {
      scale_colour_diag <- ggplot2::ggproto(NULL, scale_colour_upper)
    }

    same_colour_ld <- all(
      colour_diag == colour_lower,
      is.null(scale_colour_diag),
      !is.null(scale_colour_lower),
      !is.null(scale_colour_lower)
    )

    ## Otherwise if we have the same colour as lower plots (will not be
    ## fulfilled if we have the same colour as upper plots as we created
    ## a scale)
    if (same_colour_ld) {
      scale_colour_diag <- ggplot2::ggproto(NULL, scale_colour_lower)
    }

    ## If we still don't have a scale
    if (is.null(scale_colour_diag)) {
      if (is.numeric(toplot[[colour_diag]])) {
        scale_colour_diag <- ggplot2::scale_colour_continuous()
      } else {
        scale_colour_diag <- ggplot2::scale_colour_discrete()
      }
    }

    ## Adding the fill legend
    scale_fill_diag <- ggplot2::ggproto(NULL, scale_colour_diag)
    scale_fill_diag$aesthetics <- "fill"

    ggally_diag <- function(data, mapping, ...) {
      p <- ggplot2::ggplot(
        data,
        mapping
      ) +
        ggplot2::geom_density(
          mapping = ggplot2::aes(
            colour = !!sym(colour_diag),
            fill = !!sym(colour_diag),
            group = !!sym(colour_diag)
          ),
          alpha = 0.6
        ) +
        scale_colour_diag +
        scale_fill_diag

      return(p)
    }
  } else {
    same_colour_diag <- is.null(colour_upper)
    same_colour_ld <- is.null(colour_lower)
    ggally_diag <- GGally::wrap("densityDiag", alpha = 0.6, fill = "gray")
  }

  ## ============== ##
  ##    Final plot  ##
  ## ============== ##

  p <- GGally::ggpairs(
    toplot,
    columns = vars,
    upper = list(
      continuous = ggally_upper
    ),
    diag = list(
      continuous = ggally_diag
    ),
    lower = list(
      continuous = ggally_lower
    ),
    progress = FALSE,
    legend = legend_val
  ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    ggplot2::ggtitle(title)


  ## ============== ##
  ## Adding legends ##
  ## ============== ##

  if (!is.null(legend_val)) {
    combined_lower_legend <- all(
      all(!is.null(colour_lower), !is.null(shape_lower)),
      shape_lower == colour_lower,
      !same_shape_lower,
      !same_colour_lower
    )

    combined_legend_ud <- all(
      !is.null(colour_upper),
      !is.null(colour_diag),
      same_colour_diag
    )
    combined_legend_ld <- all(
      !is.null(colour_lower),
      !is.null(colour_diag),
      same_colour_ld
    )

    ## Adding label to legend title
    if (!same_colour_lower) {
      if (!is.null(colour_lower)) {
        if (inherits(scale_colour_lower$name, "waiver")) {
          scale_colour_lower$name <- colour_lower
        }

        suffix1 <- ifelse(
          combined_legend_ld,
          "(bottom-left & diagonal)",
          "(bottom-left)"
        )

        scale_colour_lower$name <- paste0(
          scale_colour_lower$name,
          "\n",
          suffix1
        )
      }

      if (!is.null(colour_upper)) {
        if (inherits(scale_colour_upper$name, "waiver")) {
          scale_colour_upper$name <- colour_upper
        }

        suffix2 <- ifelse(
          combined_legend_ud,
          "(top-right & diagonal)",
          "(top-right)"
        )

        scale_colour_upper$name <- paste0(
          scale_colour_upper$name,
          "\n",
          suffix2
        )
      }
    }

    if (!same_shape_lower) {
      if (!is.null(shape_lower)) {
        if (inherits(scale_shape_lower$name, "waiver")) {
          scale_shape_lower$name <- shape_lower
        }

        suffix1 <- ifelse(
          combined_legend_ld & combined_lower_legend,
          "(bottom-left & diagonal)",
          "(bottom-left)"
        )

        scale_shape_lower$name <- paste0(scale_shape_lower$name, "\n", suffix1)
      }

      if (!is.null(shape_upper)) {
        if (inherits(scale_shape_upper$name, "waiver")) {
          scale_shape_upper$name <- shape_upper
        }

        suffix2 <- ifelse(
          combined_legend_ud & all(shape_upper == colour_upper),
          "(top-right & diagonal)",
          "(top-right)"
        )

        scale_shape_upper$name <- paste0(scale_shape_upper$name, "\n", suffix2)
      }
    }

    ## If we need to add a colour legend for the bottom
    ## (different from shape variable)
    if (!is.null(colour_lower) && !same_colour_lower && !combined_lower_legend) {
      ## If discrete legend increase size point,
      ## otherwise simply set it to appear first
      scale_colour_lower_legend <- ggplot2::ggproto(NULL, scale_colour_lower)
      if (scale_colour_lower_legend$guide == "legend") {
        scale_colour_lower_legend$guide <- ggplot2::guide_legend(
          override.aes = list(size = 3, shape = 15),
          order = 5
        )
      } else {
        scale_colour_lower_legend$guide <- ggplot2::guide_colourbar(
          order = 5
        )
      }

      ggplot_aes_lower_clegend <- .make_ggplot_aes(colour = colour_lower)

      p[1, 2] <- p[1, 2] +
        ggnewscale::new_scale_colour() +
        ggplot2::geom_point(
          mapping = eval(ggplot_aes_lower_clegend),
          size = NA,
          na.rm = TRUE
        ) +
        scale_colour_lower_legend
    }

    if (!is.null(shape_lower) && !same_shape_lower && !combined_lower_legend) {
      ## Increase point size and make it appear after lower colour (if any)
      scale_shape_lower_legend <- ggplot2::ggproto(NULL, scale_shape_lower)
      scale_shape_lower_legend$guide <- ggplot2::guide_legend(
        override.aes = list(size = point_size),
        order = 6
      )

      ggplot_aes_lower_slegend <- .make_ggplot_aes(shape = shape_lower)

      p[1, 2] <- p[1, 2] +
        ggnewscale::new_scale("shape") +
        ggplot2::geom_point(
          mapping = eval(ggplot_aes_lower_slegend),
          size = NA,
          na.rm = TRUE
        ) +
        scale_shape_lower_legend
    }

    if (combined_lower_legend) {
      ## Matching order of scale and colour legends
      scale_colour_lower_legend <- ggplot2::ggproto(NULL, scale_colour_lower)
      if (scale_colour_lower_legend$guide == "legend") {
        scale_colour_lower_legend$guide <- ggplot2::guide_legend(
          order = 5
        )
      } else {
        scale_colour_lower_legend$guide <- ggplot2::guide_colourbar(
          order = 5
        )
      }

      ## Only need to increase point size for one
      ## since the legends will be combined
      scale_shape_lower_legend <- ggplot2::ggproto(NULL, scale_shape_lower)
      scale_shape_lower_legend$guide <- ggplot2::guide_legend(
        override.aes = list(size = point_size),
        order = 5
      )

      ggplot_aes_lower_comblegend <- .make_ggplot_aes(
        colour = colour_lower,
        shape = shape_lower
      )

      p[1, 2] <- p[1, 2] +
        ggnewscale::new_scale_colour() +
        ggnewscale::new_scale("shape") +
        ggplot2::geom_point(
          mapping = eval(ggplot_aes_lower_comblegend),
          size = NA,
          na.rm = TRUE
        ) +
        scale_colour_lower_legend +
        scale_shape_lower_legend
    }

    if (!same_colour_diag && !same_colour_ld && !is.null(colour_diag)) {
      ## if discrete legend increase point size, otherwise set order
      scale_colour_diag_legend <- ggplot2::ggproto(NULL, scale_colour_diag)
      if (scale_colour_diag_legend$guide == "legend") {
        scale_colour_diag_legend$guide <- ggplot2::guide_legend(
          override.aes = list(size = 3, shape = 15, alpha = 0.6),
          order = 8
        )
      } else {
        scale_colour_diag_legend$guide <- ggplot2::guide_colourbar(
          order = 8
        )
      }

      if (inherits(scale_colour_diag_legend$name, "waiver")) {
        scale_colour_diag_legend$name <- colour_diag
      }
      scale_colour_diag_legend$name <- paste0(
        scale_colour_diag_legend$name,
        "\n(diagonal)"
      )

      ggplot_aes_diag_clegend <- .make_ggplot_aes(colour = colour_diag)

      p[1, 2] <- p[1, 2] +
        ggnewscale::new_scale_colour() +
        ggplot2::geom_point(
          mapping = eval(ggplot_aes_diag_clegend),
          size = NA,
          na.rm = TRUE
        ) +
        scale_colour_diag_legend
    }
  }



  return(p)
}

#' Scatter plot function
#'
#' Creates a scatter plot
#'
#' The function adds a loess curve to summarise the trend between the covariate
#' and the samples score. If `colour_by` is used, and the corresponding variable
#' is numeric, the loess curve will not take into account this variable. If
#' instead the `colour_by` variable is a character or factor, a loess curve will
#' be fitted separately for each category.
#'
#' @inheritParams plot_x_wrapper
#' @returns A ggplot.
#' @export
plot_x_continuous <- function(toplot,
                              x,
                              y,
                              colour,
                              shape,
                              point_alpha = 1,
                              add_se = TRUE) {
  ## If the variable used for colour is numeric, will result in continuous
  ## colour scale. In that case, smoothing or boxplots should ignore the colour
  ## Otherwise, smoothing or boxplots will be split by colour (but it needs to
  ## be a factor)
  colour_sum <- NULL

  if (!is.null(colour)) {
    if (!is.numeric(toplot[[colour]])) colour_sum <- colour
  }

  gg_aes <- .make_ggplot_aes(colour = colour, shape = shape, fill = colour)
  gg_aes_sum <- .make_ggplot_aes(colour = colour_sum, fill = colour_sum)

  p <- toplot |>
    ggplot2::ggplot(aes(x = !!sym(x), y = !!sym(y))) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    ggplot2::geom_smooth(
      mapping = eval(gg_aes_sum),
      method = "loess",
      formula = y ~ x,
      alpha = 0.2,
      se = add_se
    ) +
    ggplot2::geom_point(
      mapping = eval(gg_aes),
      alpha = point_alpha
    )

  return(p)
}

#' Violin plot function function
#'
#' Creates a violin plot
#'
#' If `colour_by` is used, and the corresponding variable is numeric, the
#' violins and boxplots will not take into account this variable. If instead the
#' `colour_by` variable is a character or factor, a separate violin and boxplot
#' will be drawn for each category.
#'
#' @inheritParams plot_x_wrapper
#' @returns A ggplot.
#' @export
plot_x_discrete <- function(toplot,
                            x,
                            y,
                            colour,
                            shape,
                            point_alpha = 1,
                            add_boxplot = TRUE) {
  if (!rlang::is_installed("ggbeeswarm")) {
    stop("Please install the 'ggbeeswarm' package to use this function.")
  }

  ## If the variable used for colour is numeric, will result in continuous
  ## colour scale. In that case, smoothing or boxplots should ignore the colour
  ## Otherwise, smoothing or boxplots will be split by colour (but it needs to
  ## be a factor)
  colour_sum <- NULL

  if (!is.null(colour)) {
    if (!is.numeric(toplot[[colour]])) colour_sum <- colour
  }

  gg_aes <- .make_ggplot_aes(colour = colour, shape = shape, group = colour_sum)
  gg_aes_sum <- .make_ggplot_aes(colour = colour_sum, fill = colour_sum)

  p <- toplot |>
    ggplot2::ggplot(aes(x = !!sym(x), y = !!sym(y))) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    )

  dodge_width <- 0.95
  dodge <- ggplot2::position_dodge(width = 0.95)

  if (is.null(colour_sum)) {
    dodge_width <- NULL

    p <- p +
      ggplot2::geom_violin(
        mapping = eval(gg_aes_sum),
        colour = NA,
        fill = "grey80",
        alpha = 0.2,
        position = dodge
      )
  } else {
    p <- p +
      ggplot2::geom_violin(
        mapping = eval(gg_aes_sum),
        colour = NA,
        alpha = 0.2,
        position = dodge
      )
  }

  p <- p +
    ggbeeswarm::geom_quasirandom(
      mapping = eval(gg_aes),
      dodge.width = dodge_width,
      alpha = point_alpha
    )

  if (add_boxplot) {
    p <- p +
      ggplot2::geom_boxplot(
        mapping = eval(gg_aes_sum),
        fill = NA,
        position = dodge,
        outlier.shape = NA
      )
  }

  return(p)
}

#' Wrapper to create plot
#'
#' Wrapper around [plot_x_continuous()] and [plot_x_discrete()], will choose
#' which one to use depending on whether the x-axis variable is continuous or
#' discrete.
#'
#' @param toplot Tibble with data to plot.
#' @param x Character, name of the column in `toplot` to use as x-axis.
#' @param y Character, name of the column in `toplot` to use as y-axis.
#' @param colour Character, name of the column in `toplot` to use as colour.
#' @param shape Character, name of the column in `toplot` to use as shape.
#' @param point_alpha Numeric between 0 and 1, the opacity of the points in the
#'   plot (with 1 = fully opaque, and 0 = fully transparent). Default value is
#'   `1`.
#' @param add_se Logical, should a confidence interval be drawn around the
#'   smoothing curve for scatterplots? Default value is `TRUE`.
#' @param add_boxplot Logical, should a boxplot be drawn on top of the points
#'   for violin plots? Default value is `TRUE`.
#' @param facet_wrap Character, name of the column in `toplot` to use for
#'   faceting (using [ggplot2::facet_wrap()]). Default is `NULL`.
#' @param ncol_wrap Integer, number of columns in the faceted plot if using
#'   `facet_wrap`. Default value is `NULL`.
#' @param facet_grid Character vector of length 2, name of the columns in
#'   `toplot` to use for row (first element) and column (second element)
#'   faceting (using [ggplot2::facet_grid()]). Will be ignored if `facet_wrap`
#'   is not `NULL`. Default is `NULL`.
#' @param scales_facet Character, value to use for the `scales` argument of
#'   [ggplot2::facet_wrap()] or [ggplot2::facet_grid()]. Default value is
#'   `'free_y'`.
#' @returns A ggplot.
#' @export
plot_x_wrapper <- function(toplot,
                           x,
                           y,
                           colour,
                           shape,
                           point_alpha = 1,
                           add_se = TRUE,
                           add_boxplot = TRUE,
                           facet_wrap = NULL,
                           ncol_wrap = NULL,
                           facet_grid = NULL,
                           scales_facet = "free_y") {
  if (is.numeric(toplot[[x]])) {
    p <- plot_x_continuous(
      toplot,
      x = x,
      y = y,
      colour = colour,
      shape = shape,
      point_alpha = point_alpha,
      add_se = add_se
    )
  } else {
    p <- plot_x_discrete(
      toplot,
      x = x,
      y = y,
      colour = colour,
      shape = shape,
      point_alpha = point_alpha,
      add_boxplot = add_boxplot
    )
  }

  if (!is.null(facet_wrap)) {
    p <- p +
      ggplot2::facet_wrap(
        ggplot2::vars(!!sym(facet_wrap)),
        scales = scales_facet,
        ncol = ncol_wrap
      )
  } else if (!is.null(facet_grid)) {
    p <- p +
      ggplot2::facet_grid(
        rows = ggplot2::vars(!!sym(facet_grid[[1]])),
        cols = ggplot2::vars(!!sym(facet_grid[[2]])),
        scales = scales_facet
      )
  }

  return(p)
}


.make_ggplot_aes <- function(...) {
  args <- list(...)

  ggplot_aes <- purrr::imap(
    args,
    ~ dplyr::if_else(
      !is.null(.x),
      paste0(.y, " = ", .x),
      ""
    )
  )

  if (all(ggplot_aes == "")) {
    return(NULL)
  }

  ggplot_aes <- ggplot_aes[ggplot_aes != ""]

  ggplot_aes <- ggplot_aes |>
    paste0(collapse = ", ")

  ggplot_aes <- paste0("ggplot2::aes(", ggplot_aes, ")") |>
    str2expression()

  return(ggplot_aes)
}

.make_copy_scales <- function(scale) {
  if (is.null(scale)) {
    return(NULL)
  }

  return(ggplot2::ggproto(NULL, scale))
}
