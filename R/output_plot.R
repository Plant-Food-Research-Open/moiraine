#' Plots sample scores as scatterplot matrix
#'
#' Plots the samples score from a dimension reduction analysis as a matrix
#' of scatterplots. If there is only one latent dimension, will be plotted as boxplot
#' instead.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @param latent_dimensions Character vector giving the latent dimensions to display.
#' If `NULL` (default value), all latent dimensions will be shown.
#' @param mo_data A `MultiDataSet` object (will be used to extract samples information).
#' @param colour_upper Character, name of column in one of the samples metadata tables
#' from `mo_data` to use for colouring observations in the upper triangle plots.
#' Default value is `NULL`.
#' @param colour_diag Character, name of column in one of the samples metadata tables
#' from `mo_data` to use for colouring observations in the diagonal plots. By default,
#' will follow `colour_upper`.
#' @param colour_lower Character, name of column in one of the samples metadata tables
#' from `mo_data` to use for colouring observations in the lower triangle plots.
#' By default, will follow `colour_upper`.
#' @param shape_upper Character, name of column in one of the samples metadata tables
#' from `mo_data` to use for shaping
#' observations in the upper triangle plots. Default value is `NULL`.
#' @param shape_lower Character, name of column in one of the samples metadata tables
#' from `mo_data` to use for shaping
#' observations in the lower triangle plots. By default, will follow `shape_upper`.
#' @param scale_colour_upper ggplot2 colour scale to use for the upper triangle plots.
#' Default value is `NULL` (if `colour_upper` is not `NULL`, will use ggplot2 default
#' colour scales).
#' @param scale_colour_diag ggplot2 colour scale to use for the diagonal plots.
#' If `NULL` (default), the colour scale used for the upper triangle plots will be used
#' if `colour_diag` is equal to `colour_upper`; or the colour scale used for the
#' lower triangle plots will be used if `colour_diag` is equal to `colour_lower`.
#' @param scale_colour_lower ggplot2 colour scale to use for the lower triangle plots.
#' If `NULL` (default), the colour scale used for the upper triangle plots will be used.
#' @param scale_shape_upper ggplot2 shape scale to use for the upper triangle plots.
#' Default value is `NULL` (if `shape_upper` is not `NULL`, will use ggplot2 default
#' shape scale).
#' @param scale_shape_lower ggplot2 shape scale to use for the lower triangle plots.
#' If `NULL` (default), the shape scale used for the upper triangle plots will be used.
#' @param title Character, title of the plot. If `NULL` (default value), the method
#' name from `method_output` will be used to construct the plot title.
#' @param point_size Numeric, the size of points (in pt) in the plot. Default value
#' is 1.5.
#' @returns A `ggmatrix` plot.
#' @examples
#' \dontrun{
#' ## Let's say we've already prepared a MultiDataSet mo_data, in which the
#' ## datasets have samples metadata with columns treatment (discrete),
#' ## weeks (continuous), tissue_type (discrete), disease_score (continuous).
#'
#' library(ggplot2)
#'
#' pca_res <- run_pca(mo_data, "metabolome")
#' output_pca <- get_output_pca(output_pca)
#' pcs <- paste0("Principal component ", 1:4)
#'
#' # Simple matrix of scatterplot to visualised PCs two by two
#' plot_samples_score(
#'   output_pca,
#'   pcs
#' )
#'
#' # Colouring points according to weeks
#' plot_samples_score(
#'   output_pca,
#'   pcs,
#'   colour_upper = "weeks"
#' )
#' # Adding a custom colour palette
#' plot_samples_score(
#'   output_pca,
#'   pcs,
#'   colour_upper = "weeks",
#'   scale_colour_upper = scale_colour_viridis_c()
#' )
#'
#' # Adding the treatment as shape
#' plot_samples_score(
#'   output_pca,
#'   pcs,
#'   colour_upper = "weeks",
#'   shape_upper = "treatment"
#' )
#'
#' # Using the lower triangle of the plots to display disease score
#' # Again can pass custom colour scale through scale_colour_lower
#' plot_samples_score(
#'   output_pca,
#'   pcs,
#'   colour_upper = "weeks",
#'   shape_upper = "treatment",
#'   colour_lower = "disease_score"
#' )
#'
#' # By default the diagonal plots follow the colour of the upper plots,
#' # but can follow the lower plots instead
#' plot_samples_score(
#'   output_pca,
#'   pcs,
#'   colour_upper = "weeks",
#'   shape_upper = "treatment",
#'   colour_lower = "disease_score",
#'   colour_diag = "disease_score"
#' )
#'
#' # or diagonal can show a different variable
#' plot_samples_score(
#'   output_pca,
#'   pcs,
#'   colour_upper = "weeks",
#'   shape_upper = "treatment",
#'   colour_lower = "tissue_type"
#' )
#'
#' # also the lower plots can have a different shape than the upper plots
#' plot_samples_score(
#'   output_pca,
#'   pcs,
#'   colour_upper = "weeks",
#'   shape_upper = "treatment",
#'   colour_lower = "disease_score",
#'   shape_lower = "tissue_type"
#' )
#' }
#' @export
plot_samples_score <- function(method_output,
                               latent_dimensions = NULL,
                               mo_data = NULL,
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
  ## For devtools::check
  latent_dimension <- score <- id <- NULL

  if (is.null(latent_dimensions)) {
    latent_dimensions <- get_latent_dimensions(method_output)
  } else {
    method_output <- .filter_output_dimensions(method_output, latent_dimensions)
  }

  aes_vars <- c("colour_upper", "colour_diag", "colour_lower", "shape_upper", "shape_lower")

  toplot <- method_output$samples_score

  ## Extracting the additional aesthetics from samples metadata
  if (any(purrr::map_lgl(aes_vars, ~ !is.null(get(.x))))) {
    if (is.null(mo_data)) {
      stop("Need to provide a MultiDataSet object through 'mo_data' argument in order to use one of the aesthetics arguments.")
    }

    mo_data <- check_input_multidataset(mo_data)

    aes_vars |>
      purrr::iwalk(
        ~ paste0(".check_input_var_smetadata_common(", .x, ", mo_data)") |>
          str2expression() |>
          eval()
      )

    aes_cols <- aes_vars |>
      purrr::map(~ get(.x)) |>
      purrr::reduce(union)

    smeta <- get_samples_metadata_combined(mo_data, only_common_cols = FALSE) |>
      tibble::as_tibble() |>
      dplyr::select(sample_id = id, tidyselect::any_of(aes_cols))

    toplot <- toplot |>
      dplyr::left_join(
        smeta,
        by = "sample_id"
      )
  }

  ## Default title
  if (is.null(title)) {
    title <- paste0(
      "Sample scores - ",
      attr(method_output, "method")
    )
  }

  if (length(latent_dimensions) == 1) {

    cause_warning <- !all(
      c(
        is_equal_or_null(colour_diag, colour_upper),
        is_equal_or_null(colour_lower, colour_upper),
        is_equal_or_null(shape_lower, shape_upper)
      )
    )

    if (cause_warning) {
      warning("Only one latent dimension to plot; 'colour_diag', 'colour_lower' and 'shape_lower' argument will be ignored.")
    }

    toplot <- toplot |>
      dplyr::mutate(Vx = "x")

    if (!is.null(scale_colour_upper)) {
      scale_colour_upper <- .make_copy_scales(scale_colour_upper)
      scale_colour_upper$aesthetics <- c("colour", "fill")
    }

    p <- plot_x_discrete(
      toplot,
      x = "Vx",
      y = "score",
      colour = colour_upper,
      shape = shape_upper,
      add_boxplot = FALSE
    ) +
      ggplot2::facet_wrap(~ latent_dimension, scales = "free_y") +
      ggplot2::labs(
        x = NULL,
        title = title
      ) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      scale_colour_upper +
      scale_shape_upper

  } else {

    toplot <- toplot |>
      tidyr::pivot_wider(
        names_from = latent_dimension,
        values_from = score
      )

    p <- ggpairs_custom(
      toplot,
      latent_dimensions,
      colour_upper,
      colour_diag,
      colour_lower,
      shape_upper,
      shape_lower,
      scale_colour_upper,
      scale_colour_diag,
      scale_colour_lower,
      scale_shape_upper,
      scale_shape_lower,
      title,
      point_size
    )
  }

  return(p)
}



#' Plots sample scores as a scatterplot
#'
#' Plots the samples score from a pair of latent dimensions from a dimension
#' reduction analysis as a scatterplot.
#'
#' @param method_output A single Integration method output generated via the
#' `get_output()` function, or a list of two integration method outputs.
#' @param latent_dimensions Character vector of length 2 (if `method_output` is a
#' single output object), or named list of length 2 (if `method_output` is a list
#' of two output objects). In the first case, gives the name of the two latent dimensions
#' that should be represented. In the second case, the names of the list should
#' correspond to the names of the methods, and the values should each be a character
#' giving for the corresponding method the name of the latent dimension to display.
#' @param mo_data A `MultiDataSet` object (will be used to extract samples information).
#' @param colour_by Character, name of column in one of the samples metadata tables
#' from `mo_data` to use for colouring samples in the plot. Default value is `NULL`.
#' @param shape_by Character, name of column in one of the samples metadata tables
#' from `mo_data` to use a shape for the samples in the plot. Default value is `NULL`.
#' @returns A ggplot.
#' @examples
#' \dontrun{
#' ## Let diablo_res be the output from a DIABLO analysis, and mofa_res the
#' ## output from a MOFA analysis. Let mo_set be the corresponding MultiDataSet object.
#' output_diablo <- get_output_diablo(diablo_res)
#' output_mofa <- get_output_mofa2(mofa_res)
#'
#' ## Scatterplot of the first two DIABLO components
#' plot_samples_score_pair(output_diablo, c("Component 1", "Component 2"))
#'
#' ## Adding samples information to the plot - here 'Time' and 'Treatment' should
#' ## two columns in the samples metadata of one of the datasets in mo_set
#' plot_samples_score_pair(
#'   output_diablo,
#'   c("Component 1", "Component 2"),
#'   mo_data <- mo_set,
#'   colour_by = "Time",
#'   shape_by = "Treatment"
#' )
#'
#' ## Comparing the first MOFA factor to the first DIABLO component
#' plot_samples_score_pair(
#'   list(output_diablo, output_mofa),
#'   list("DIABLO" = "Component 1", "MOFA" = "Factor 1"),
#'   mo_data <- mo_set,
#'   colour_by = "Time",
#'   shape_by = "Treatment"
#' )
#'
#' ## Giving custom names to the methods
#' plot_samples_score_pair(
#'   list("DIABLO prefiltered" = output_diablo, "MOFA full" = output_mofa),
#'   list("DIABLO prefiltered" = "Component 1", "MOFA full" = "Factor 1")
#' )
#' }
#' @export
plot_samples_score_pair <- function(method_output,
                                    latent_dimensions,
                                    mo_data = NULL,
                                    colour_by = NULL,
                                    shape_by = NULL) {
  ## For devtools::check
  latent_dimension <- method <- score <- id <- NULL

  if (inherits(method_output, "output_dimension_reduction")) {
    method_output <- .filter_output_dimensions(
      method_output,
      latent_dimensions,
      fixed_length = 2
    )

    toplot <- method_output$samples_score

    ## Default title
    title <- paste0(
      "Sample scores - ",
      attr(method_output, "method")
    )
  } else {
    if (length(method_output) != 2) {
      stop("'method_output' argument: expecting either a single 'output_dimension_reduction' object or a list of two such objects.")
    }

    method_output <- method_output |>
      .check_names_output_list() |>
      .filter_output_dimensions_list(
        latent_dimensions,
        all_present = TRUE,
        fixed_length = 1
      )

    toplot <- method_output[names(latent_dimensions)] |>
      purrr::map_dfr(
        ~ .x$samples_score,
        .id = "method"
      ) |>
      dplyr::mutate(
        latent_dimension = paste0(method, " ", latent_dimension)
      ) |>
      dplyr::select(-method)

    latent_dimensions <- unique(toplot$latent_dimension)

    ## Default title
    title <- paste0(
      "Sample scores - ",
      paste0(names(method_output), collapse = " vs ")
    )
  }

  toplot <- toplot |>
    tidyr::pivot_wider(
      names_from = latent_dimension,
      values_from = score
    )

  aes_vars <- c("colour_by", "shape_by")

  ## Extracting the additional aesthetics from samples metadata
  if (any(purrr::map_lgl(aes_vars, ~ !is.null(get(.x))))) {
    if (is.null(mo_data)) {
      stop("Need to provide a MultiDataSet object through 'mo_data' argument in order to use one of the aesthetics arguments.")
    }

    mo_data <- check_input_multidataset(mo_data)

    aes_vars |>
      purrr::iwalk(
        ~ paste0(".check_input_var_smetadata_common(", .x, ", mo_data)") |>
          str2expression() |>
          eval()
      )

    aes_cols <- aes_vars |>
      purrr::map(~ get(.x)) |>
      purrr::reduce(union)

    smeta <- get_samples_metadata_combined(mo_data, only_common_cols = FALSE) |>
      tibble::as_tibble() |>
      dplyr::select(sample_id = id, tidyselect::any_of(aes_cols))

    toplot <- toplot |>
      dplyr::left_join(
        smeta,
        by = "sample_id"
      )
  }

  gg_aes <- .make_ggplot_aes(colour = colour_by, shape = shape_by)

  p <- ggplot2::ggplot(
    toplot,
    aes(
      x = !!sym(latent_dimensions[[1]]),
      y = !!sym(latent_dimensions[[2]])
    )
  ) +
    ggplot2::geom_point(mapping = eval(gg_aes)) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )

  return(p)
}

#' Plots sample scores against covariate
#'
#' Plots the samples score from the result of an integration method against
#' a covariate from the samples metadata.
#'
#' If the covariate is numeric, the function creates a scatter plot, with a
#' loess curve to summarise the trend between the covariate and the samples score.
#' If `colour_by` is used, and the corresponding variable is numeric, the loess curve
#' will not take into account this variable. If instead the `colour_by` variable is
#' a character or factor, a loess curve will be fitted separately for each category.
#'
#' If the covariate is not numeric, the function creates a violin/boxplot. If `colour_by`
#' is used, and the corresponding variable is numeric, the violins and boxplots
#' will not take into account this variable. If instead the `colour_by` variable is
#' a character or factor, a separate violin and boxplot will be drawn for each category.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @param mo_data A `MultiDataSet` object (will be used to extract samples information).
#' @param covariate Character, name of column in one of the samples metadata tables
#' from `mo_data` to use as x-axis in the plot.
#' @param latent_dimensions Character vector giving the latent dimensions to display.
#' Default value is `NULL`, i.e. all latent dimensions will be shown.
#' @param colour_by Character, name of column in one of the samples metadata tables
#' from `mo_data` to use to colour the samples in the plot. Default value is `NULL`.
#' @param shape_by Character, name of column in one of the samples metadata tables
#' from `mo_data` to use as shape for the samples in the plot. Default value is `NULL`.
#' @param point_alpha Numeric between 0 and 1, the opacity of the points
#' in the plot (with 1 = fully opaque, and 0 = fully transparent). Default value
#' is `1`.
#' @param add_se Logical, should a confidence interval be drawn around the smoothing
#' curves for numerical covariates? Default value is `TRUE`.
#' @param add_boxplot Logical, should a boxplot be drawn on top of the points for
#' categorical covariates? Default value is `TRUE`.
#' @param ncol Integer, number of columns in the faceted plot. Default value is `NULL`.
#' @returns a ggplot.
#' @export
plot_samples_score_covariate <- function(method_output,
                                         mo_data,
                                         covariate,
                                         latent_dimensions = NULL,
                                         colour_by = NULL,
                                         shape_by = NULL,
                                         point_alpha = 1,
                                         add_se = TRUE,
                                         add_boxplot = TRUE,
                                         ncol = NULL) {
  ## for devtools::check
  score <- id <- NULL

  method_output <- .filter_output_dimensions(method_output, latent_dimensions)

  ## Checking inputs
  mo_data <- check_input_multidataset(mo_data)
  .check_input_var_smetadata_common(covariate, mo_data)
  .check_input_var_smetadata_common(colour_by, mo_data)
  .check_input_var_smetadata_common(shape_by, mo_data)

  aes_vars <- c(covariate, colour_by, shape_by) |>
    unique()

  toplot <- method_output$samples_score |>
    dplyr::left_join(
      get_samples_metadata_combined(mo_data, FALSE) |>
        dplyr::select(sample_id = id, tidyselect::all_of(aes_vars)),
      by = "sample_id"
    )

  p <- plot_x_wrapper(
    toplot,
    x = covariate,
    y = "score",
    facet_wrap = "latent_dimension",
    colour = colour_by,
    shape = shape_by,
    point_alpha = point_alpha,
    add_se = add_se,
    add_boxplot = add_boxplot,
    ncol_wrap = ncol
  ) +
    ggplot2::labs(
      y = "Samples score",
      title = paste0("Samples score - ", attr(method_output, "method"))
    )

  return(p)
}

#' Plots features weight distribution
#'
#' Plots the features weight or importance per dataset and latent dimension.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @param latent_dimensions Character vector giving the latent dimensions to display.
#' Default value is `NULL`, i.e. all latent dimensions will be shown.
#' @param datasets Character vector giving the datasets to display.
#' Default value is `NULL`, i.e. all datasets will be shown.
#' @param features_metric Character, which attribute should be plotted: can be
#' `'signed_importance'` (i.e. importance value but with the weight sign),
#' `'importance'` or `'weight'`. Default value is `'signed_importance'`.
#' @param top_n Integer, number of top features (in terms of importance) for which
#' the label should be shown. Default value is `0`.
#' @inheritParams .add_features_labels_toplot
#' @param text_size Numeric, size of the feature labels.
#' @returns A \code{\link[patchwork]{patchwork}} of plots.
#' @export
plot_features_weight_distr <- function(method_output,
                                       latent_dimensions = NULL,
                                       datasets = NULL,
                                       features_metric = c("signed_importance", "weight", "importance"),
                                       top_n = 0,
                                       mo_data = NULL,
                                       label_cols = NULL,
                                       truncate = NULL,
                                       text_size = 2.5) {
  ## for devtools::check
  dataset <- latent_dimension <- label <- feature_id <- importance <- weight <- data <- NULL

  method_output <- method_output |>
    .filter_output_dimensions(latent_dimensions) |>
    .filter_output_datasets(
      datasets
    )
  toplot <- method_output$features_weight

  ## Handling features_metric parameter
  features_metric <- rlang::arg_match(features_metric)
  x_label <- features_metric |>
    stringr::str_replace("_", " ") |>
    stringr::str_to_sentence()

  ## Preparing plotting dataframe
  toplot <- toplot |>
    dplyr::mutate(
      signed_importance = sign(weight) * importance
    ) |>
    dplyr::group_by(latent_dimension, dataset) |>
    dplyr::arrange(!!sym(features_metric)) |>
    dplyr::mutate(
      rank = seq_len(dplyr::n())
    ) |>
    dplyr::ungroup()

  if (top_n > 0) {
    ## Add features label
    toplot <- toplot

    ## Select top features per latent dimension and dataset
    temp <- toplot |>
      dplyr::group_by(dataset, latent_dimension) |>
      dplyr::slice_max(order_by = abs(importance), n = top_n) |>
      dplyr::ungroup() |>
      .add_features_labels_toplot(label_cols, mo_data, truncate) |>
      dplyr::select(feature_id, dataset, latent_dimension, label)

    ## Add column for the top features
    toplot <- toplot |>
      dplyr::left_join(
        temp,
        by = c("feature_id", "dataset", "latent_dimension")
      )
  }

  plot_fct <- function(.x, .y) {
    p <- .x |>
      ggplot2::ggplot(aes(x = !!sym(features_metric), y = rank)) +
      ggplot2::geom_point(size = 1, colour = "gray")

    if (top_n > 0) {
      p <- p +
        ggplot2::geom_point(
          data = dplyr::filter(.x, !is.na(label)),
          size = 1.5,
          colour = "black"
        ) +
        ggrepel::geom_text_repel(
          aes(label = label),
          na.rm = TRUE,
          size = text_size
        )
    }

    if (features_metric == "signed_importance") {
      p <- p +
        ggplot2::scale_x_continuous(limits = c(-1, 1))
    }

    p <- p +
      ggplot2::facet_wrap(~dataset, scales = "free") +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = 0.01)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.subtitle = ggplot2::element_text(hjust = 0.5),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank(),
      ) +
      ggplot2::labs(
        x = x_label,
        y = "Ranked features",
        subtitle = .y
      )

    return(p)
  }

  p <- toplot |>
    dplyr::arrange(latent_dimension) |>
    dplyr::group_by(latent_dimension) |>
    tidyr::nest() |>
    dplyr::mutate(
      plot = purrr::map2(
        data,
        latent_dimension,
        plot_fct
      )
    ) |>
    dplyr::pull(plot) |>
    patchwork::wrap_plots() +
    patchwork::plot_annotation(
      title = paste0("Features weight distribution - ", attr(method_output, "method")),
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    )

  return(p)
}

#' Plots top features importance
#'
#' Plots the top features importance per dataset and latent dimension.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @param latent_dimensions Character vector giving the latent dimensions to display.
#' Default value is `NULL`, i.e. all latent dimensions will be shown.
#' @param group_latent_dims Logical, for integrations methods that construct datasets-
#' specific versions of each latent dimension, should these be grouped? e.g. when
#' DIABLO constructs a snps- and rnaseq version of component 1, should the two be
#' grouped as "Component 1"? Default value is `TRUE`.
#' @param datasets Character vector giving the datasets to display.
#' Default value is `NULL`, i.e. all datasets will be shown.
#' @param n_features Integer, number of top features to display per dataset and
#' latent dimension.
#' @inheritParams .add_features_labels_toplot
#' @param nrow Integer, number of rows over which the dataset panels should be
#' plotted for each latent dimensions.
#' @returns A \code{\link[patchwork]{patchwork}} of plots.
#' @export
plot_top_features <- function(method_output,
                              latent_dimensions = NULL,
                              group_latent_dims = TRUE,
                              datasets = NULL,
                              n_features = 20,
                              mo_data = NULL,
                              label_cols = NULL,
                              truncate = NULL,
                              nrow = 1) {
  ## for devtools::check
  dataset <- latent_dimension <- contribution <- importance <- data <- weight <-  NULL

  if (group_latent_dims) method_output <- .rm_dataset_from_latentdim(method_output)

  method_output <- method_output |>
    .filter_output_dimensions(latent_dimensions) |>
    .filter_output_datasets(datasets)

  method_output$features_weight |>
    dplyr::filter(weight != 0) |>
    .add_features_labels_toplot(label_cols, mo_data, truncate) |>
    dplyr::group_by(dataset, latent_dimension) |>
    dplyr::mutate(
      contribution = dplyr::case_when(
        sign(weight) > 0 ~ "positive",
        sign(weight) < 0 ~ "negative"
      ),
      contribution = factor(contribution, levels = c("negative", "positive"))
    ) |>
    dplyr::slice_max(importance, n = n_features) |>
    dplyr::ungroup() |>
    dplyr::group_by(latent_dimension) |>
    tidyr::nest() |>
    dplyr::mutate(
      data = purrr::map(
        data,
        ~ .x |>
          dplyr::group_by(dataset) |>
          dplyr::arrange(importance) |>
          dplyr::mutate(
            label = make.unique(label),
            label = factor(label, levels = label)
          )
      ),
      plot = purrr::map2(
        data,
        latent_dimension,
        ~ .x |>
          ggplot2::ggplot(ggplot2::aes(
            x = label,
            y = importance,
            colour = contribution,
            fill = contribution
          )) +
          ggplot2::geom_col(width = 0.05) +
          ggplot2::geom_point(size = 3) +
          ggplot2::facet_wrap(~dataset, scales = "free", nrow = nrow) +
          ggplot2::labs(
            x = NULL,
            y = "Importance",
            title = .y,
            colour = "Contribution",
            fill = "Contribution"
          ) +
          ggplot2::coord_flip() +
          ggplot2::scale_colour_manual(
            values = c(
              "negative" = "royalblue",
              "positive" = "firebrick"
            ),
            drop = FALSE,
            na.value = "white"
          ) +
          ggplot2::scale_fill_manual(
            values = c(
              "negative" = "royalblue",
              "positive" = "firebrick"
            ),
            drop = FALSE,
            na.value = "white"
          )
      )
    ) |>
    dplyr::pull(plot) |>
    patchwork::wrap_plots() +
    patchwork::plot_layout(
      guides = "collect",
      ncol = 1
    ) +
    patchwork::plot_annotation(
      title = paste0("Top ", n_features, " features - ", attr(method_output, "method"))
    ) &
    ggplot2::theme_bw() &
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    )
}


#' Plots features weight as a scatterplot
#'
#' Plots the features weight from a pair of latent dimensions from one or two
#' dimension reduction analysis as a scatterplot.
#'
#' @param method_output A single Integration method output generated via the
#' `get_output()` function, or a list of two integration method outputs.
#' @param latent_dimensions Character vector of length 2 (if `method_output` is a
#' single output object), or named list of length 2 (if `method_output` is a list
#' of two output objects). In the first case, gives the name of the two latent dimensions
#' that should be represented. In the second case, the names of the list should
#' correspond to the names of the methods, and the values should each be a character
#' giving for the corresponding method the name of the latent dimension to display.
#' @param datasets Character vector, names of the datasets for which features weight
#' should be plotted. Default value is `NULL`, i.e. all relevant datasets are shown.
#' @param features_metric Character, the features metric that should be plotted on
#' the y-axis. Should be one of `'signed_importance'` (default value), `'weight'` or
#' `'importance'`.
#' @param include_missing_features Logical, whether to show the features that were in the input
#' of one method but not the other, when comparing the results of two different integration
#' methods. Default value is `TRUE`.
#' @param top_n Integer, number of top features (according to the consensus importance metric)
#' to highlight in the plot. Default value is `5`.
#' @inheritParams compute_consensus_importance
#' @inheritParams .add_features_labels_toplot
#' @param ncol Integer, number of columns of datasets in the combined plot.
#' Default value is `NULL`, i.e. will be picked automatically.
#' @param label_size Integer, size of the features label. Default value is `3`.
#' @returns A ggplot.
#' @export
plot_features_weight_pair <- function(method_output,
                                      latent_dimensions,
                                      datasets = NULL,
                                      features_metric = c("signed_importance", "weight", "importance"),
                                      include_missing_features = TRUE,
                                      top_n = 5,
                                      metric = "geometric",
                                      label_cols = NULL,
                                      mo_data = NULL,
                                      truncate = NULL,
                                      ncol = NULL,
                                      label_size = 3) {
  ## For devtools::check
  name <- latent_dimension <- dataset <- x <- y <- importance <- is_top <- NULL
  feature_id <- label <- weight <- NULL

  features_metric <- rlang::arg_match(features_metric)
  y_label <- paste("Features", stringr::str_replace(features_metric, "_", " "))

  if (inherits(method_output, "output_dimension_reduction")) {

    method_output <- .filter_output_dimensions(
      method_output,
      latent_dimensions,
      fixed_length = 2
    ) |>
      .filter_output_datasets(datasets)

    toplot <- method_output$features_weight
    axis_vars <- latent_dimensions
    title <- paste0("Features ", y_label, " - ", attr(method_output, "method"))
  } else {
    if (length(method_output) != 2) {
      stop(
        "'method_output' should be of length 2; this function ",
        "is for comparing the output of two integration methods."
      )
    }

    method_output <- method_output |>
      .check_names_output_list() |>
      .filter_output_dimensions_list(
        latent_dimensions,
        all_present = TRUE,
        fixed_length = 1
      ) |>
      purrr::imap(
        ~.filter_output_datasets(.x, datasets, method_name = .y)
      )

    axis_vars <- latent_dimensions |>
      purrr::imap_chr(
        ~ paste0(.y, " ", .x)
      ) |>
      unname()

    toplot <- method_output |>
      purrr::map_dfr(
        ~ .x$features_weight,
        .id = "name"
      ) |>
      dplyr::mutate(
        latent_dimension = paste0(name, " ", latent_dimension)
      )

    title <- paste0(
      "Features ",
      y_label,
      " - ",
      paste0(names(method_output), collapse = " vs ")
    )
  }

  toplot <- toplot |>
    dplyr::mutate(
      latent_dimension = dplyr::case_when(
        latent_dimension == axis_vars[[1]] ~ "x",
        latent_dimension == axis_vars[[2]] ~ "y"
      ),
      signed_importance = sign(weight) * importance
    ) |>
    dplyr::select(latent_dimension, dataset, feature_id, tidyselect::all_of(features_metric)) |>
    tidyr::pivot_wider(
      names_from = latent_dimension,
      values_from = !!sym(features_metric)
    )

  if (include_missing_features) {
    toplot <- toplot |>
      tidyr::replace_na(list(x = 0, y = 0))
  } else {
    toplot <- toplot |>
      dplyr::filter(!is.na(x), !is.na(y))
  }

  top_features <- compute_consensus_importance(method_output, latent_dimensions, metric, include_missing_features) |>
    dplyr::filter(feature_id %in% toplot$feature_id) |>
    dplyr::group_by(dataset) |>
    dplyr::slice_max(importance, n = top_n)

  p <- toplot |>
    .add_features_labels_toplot(label_cols, mo_data, truncate) |>
    dplyr::mutate(
      is_top = feature_id %in% top_features$feature_id,
      label = dplyr::case_when(
        is_top ~ label,
        TRUE ~ NA_character_
      )
    ) |>
    ggplot2::ggplot(
      aes(
        x = x,
        y = y,
        colour = is_top
      )
    )

  if (features_metric != "weight") {
    p <- p +
      ggplot2::geom_abline(
        intercept = 0,
        slope = 1,
        colour = "gray",
        linetype = 2
      )
  }

  p <- p +
    ggplot2::geom_point(alpha = 0.6) +
    ggrepel::geom_text_repel(
      aes(label = label),
      colour = "black",
      na.rm = TRUE,
      size = label_size
    ) +
    ggplot2::facet_wrap(
      ~dataset,
      scales = dplyr::if_else(features_metric == "weight", "free", "fixed"),
      ncol = ncol
    ) +
    ggplot2::scale_colour_manual(
      values = c("TRUE" = "red", "FALSE" = "black"),
      guide = "none"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      title = title,
      caption = paste0(
        "Top ", top_n, " features per dataset (in red) selected ",
        "according to ", metric, " consensus importance."
      ),
      x = axis_vars[[1]],
      y = axis_vars[[2]]
    )

  return(p)
}

#' Plots features weight against covariate
#'
#' Plots the features weight or importance from the result of an integration
#' method against a covariate from the features metadata.
#'
#' If the covariate is numeric, the function creates a scatter plot, with a
#' loess curve to summarise the trend between the covariate and the features weight.
#' If `colour_by` is used, and the corresponding variable is numeric, the loess curve
#' will not take into account this variable. If instead the `colour_by` variable is
#' a character or factor, a loess curve will be fitted separately for each category.
#'
#' If the covariate is not numeric, the function creates a violin/boxplot. If `colour_by`
#' is used, and the corresponding variable is numeric, the violins and boxplots
#' will not take into account this variable. If instead the `colour_by` variable is
#' a character or factor, a separate violin and boxplot will be drawn for each category.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @param mo_data A `MultiDataSet` object (will be used to extract samples information).
#' @param covariate Character or named list of character, giving for each dataset the name
#' of the column in the corresponding features metadata to use as x-axis in the plot. If one
#' value, will be used for all datasets. If list, the names must correspond to
#' the names of the datasets in `mo_data`. If a dataset is not present in this list,
#' will be excluded from the plot.
#' @param features_metric Character, the features metric that should be plotted on
#' the y-axis. Should be one of `'signed_importance'` (default value), `'weight'` or
#' `'importance'`.
#' @param remove_null_weight Logical, should features with null weight/importance be
#' removed from the plot? Default value is `FALSE`.
#' @param latent_dimensions Character vector giving the latent dimensions to display.
#' Default value is `NULL`, i.e. all latent dimensions will be shown.
#' @param colour_by Character or named list of character, giving for each dataset the
#' name of column in the corresponding feature metadata to use to colour the features in the plot.
#' If one value, will be used for all datasets. If list, the names must correspond to
#' the names of the datasets in `covariate`. Default value is `NULL`.
#' @param shape_by Character or named list of character, giving for each dataset the
#' name of column in the corresponding feature metadata to use as shape for the features in the plot.
#' If one value, will be used for all datasets. If list, the names must correspond to
#' the names of the datasets in `covariate`. Default value is `NULL`.
#' @param point_alpha Numeric between 0 and 1, the opacity of the points
#' in the plot (with 1 = fully opaque, and 0 = fully transparent). Default value
#' is `0.5`.
#' @param add_se Logical, should a confidence interval be drawn around the smoothing
#' curves for numerical covariates? Default value is `TRUE`.
#' @param add_boxplot Logical, should a boxplot be drawn on top of the points for
#' categorical covariates? Default value is `TRUE`.
#' @param scales Character, value to use for the `scales` argument of [ggplot2::facet_grid()].
#' Default value is `'free_x'`.
#' @export
plot_features_weight_covariate <- function(method_output,
                                           mo_data,
                                           covariate,
                                           features_metric = c("signed_importance", "weight", "importance"),
                                           remove_null_weight = FALSE,
                                           latent_dimensions = NULL,
                                           colour_by = NULL,
                                           shape_by = NULL,
                                           point_alpha = 0.5,
                                           add_se = TRUE,
                                           add_boxplot = TRUE,
                                           scales = "free_x") {
  ## for devtools::check
  dataset <- weight <- importance <- NULL

  features_metric <- rlang::arg_match(features_metric)
  y_label <- paste("Features", stringr::str_replace(features_metric, "_", " "))

  method_output <- .filter_output_dimensions(method_output, latent_dimensions)

  ## Checking inputs
  datasets <- levels(method_output$features_weight$dataset)
  mo_data <- check_input_multidataset(mo_data, datasets)


  ## Names to use for plot wrapper
  aes_vars <- list(
    "colour_by" = colour_by,
    "shape_by" = shape_by
  ) |>
    purrr::imap(
      \(.x, .y) {
        if (is.null(.x)) NULL else .y
      }
    )

  ## can't use a loop otherwise error message is uninformative
  covariate <- .make_var_list(covariate, datasets, fixed_length = 1, allow_subset = TRUE)
  ## restricting further arguments to those in covariate list but keeping order of datasets
  datasets <- intersect(datasets, names(covariate))
  colour_by <- .make_var_list(colour_by, datasets, fixed_length = 1, allow_subset = TRUE)
  shape_by <- .make_var_list(shape_by, datasets, fixed_length = 1, allow_subset = TRUE)

  .check_input_var_fmetadata(covariate, mo_data)
  .check_input_var_fmetadata(colour_by, mo_data)
  .check_input_var_fmetadata(shape_by, mo_data)

  vars <- list(
    "covariate" = covariate,
    "colour_by" = colour_by,
    "shape_by" = shape_by
  )

  ## Columns to collect from each dataset's feature metadata
  col_vars <- vars |>
    purrr::transpose() |>
    purrr::map(unlist)

  ## Labels for the plot
  labels_vars <- vars |>
    purrr::imap(
      ~ .x |>
        purrr::discard(is.null) |>
        unlist() |>
        tibble::enframe(name = "dataset", value = "variable") |>
        dplyr::mutate(
          dataset = factor(dataset, levels = datasets),
          is_full = dplyr::n() == length(datasets),
          variable = as.character(variable) ## for the case when tibble is empty
        ) |>
        dplyr::arrange(dataset) |>
        dplyr::mutate(
          variable = factor(variable, levels = unique(variable))
        ) |>
        dplyr::group_by(variable, is_full) |>
        dplyr::summarise(
          dataset = paste(as.character(dataset), collapse = ", "),
          .groups = "drop"
        ) |>
        dplyr::mutate(
          variable = as.character(variable),
          variable = dplyr::case_when(
            (dplyr::n() == 1) & is_full ~ variable,
            TRUE ~ paste0(variable, " (", dataset, ")")
          )
        ) |>
        dplyr::pull(variable) |>
        paste(collapse = dplyr::if_else(.y == "covariate", ", ", ",\n"))
    ) |>
    purrr::map(~ (if (.x == "") NULL else .x))

  ## Collecting features metadata
  fmeta <- get_features_metadata(mo_data) |>
    purrr::imap_dfr(
      ~ .x |>
        dplyr::select(
          feature_id,
          tidyselect::all_of(col_vars[[.y]])
        ) |>
        dplyr::mutate(
          dplyr::across(
            .cols = dplyr::where(\(x){!is.numeric(x)}),
            .fns = as.character
          )
        ),
      .id = "dataset"
    ) |>
    dplyr::mutate(
      dataset = factor(dataset, levels = datasets)
    )

  toplot <- method_output$features_weight |>
    dplyr::filter(dataset %in% names(col_vars)) |>
    dplyr::left_join(fmeta, by = c("feature_id", "dataset")) |>
    dplyr::mutate(
      signed_importance = sign(weight) * importance
    )

  if (remove_null_weight) {
    toplot <- toplot |>
      dplyr::filter(weight != 0)
  }

  p <- plot_x_wrapper(
    toplot,
    x = "covariate",
    y = features_metric,
    facet_grid = c("latent_dimension", "dataset"),
    colour = aes_vars[["colour_by"]],
    shape = aes_vars[["shape_by"]],
    point_alpha = point_alpha,
    add_se = add_se,
    add_boxplot = add_boxplot,
    scales_facet = scales
  ) +
    ggplot2::labs(
      y = y_label,
      title = paste(y_label, "-", attr(method_output, "method")),
      x = labels_vars[["covariate"]],
      colour = labels_vars[["colour_by"]],
      fill = labels_vars[["colour_by"]],
      shape = labels_vars[["shape_by"]]
    )

  return(p)
}

##Would work if for so2pls we write "rnaseq-specific" rather than "rnaseq specific"
.rm_dataset_from_latentdim <- function(method_output) {

  ## for devtools::check
  dataset <- latent_dimension <- new_latent_dimension <- NULL

  new_ld_levels <- method_output$features_weight |>
    dplyr::select(dataset, latent_dimension) |>
    dplyr::distinct() |>
    dplyr::mutate(
      latent_dimension = as.character(latent_dimension),
      new_latent_dimension = purrr::map2_chr(
        dataset, latent_dimension,
        function(.x, .y) {
          if (stringr::str_detect(.y, "specific component")) {
            return(.y)
          }

          return(stringr::str_remove(.y, paste0("^", .x, " ")))
        }
      )
    ) |>
    dplyr::select(latent_dimension, new_latent_dimension) |>
    tibble::deframe()

  method_output$features_weight <- method_output$features_weight |>
    dplyr::mutate(
      latent_dimension = as.character(latent_dimension),
      latent_dimension = new_ld_levels[latent_dimension],
      latent_dimension = factor(latent_dimension, unique(unname(new_ld_levels)))
    )

  return(method_output)
}
