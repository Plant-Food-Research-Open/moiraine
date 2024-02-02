#' Get initials from sentence
#'
#' Extracts initials from a sentence as well as any number.
#'
#' @param x Character vector.
#' @returns A character vector, each element containing the initials of the words
#' in `x` in upper case plus any number, all pasted.
.get_initials <- function(x) {
  res <- stringr::str_extract_all(x, "((?<=^| )[[:alpha:]])|(\\d+)")
  res <- purrr::map_chr(res, ~ stringr::str_c(.x, collapse = ""))
  res <- stringr::str_to_upper(res)

  return(res)
}


#' Get features weight correlation
#'
#' Constructs the correlation matrix between the features weight of the latent
#' dimensions obtained with different integration methods.
#'
#' If `include_missing_features` is `FALSE` (default behaviour), and some
#' features are present in the output of one integration method but not
#' the other (e.g. because a different pre-filtering was applied to the input
#' data of the two methods), these features will be ignored. This does not mean
#' that features that were selected by one method but not the other are discarded;
#' in that case the feature will be assigned a weight of 0 for the method that did
#' not select it. This is the recommended behaviour, should only be changed in
#' specific scenarios (e.g. to check whether using all features in a dataset vs
#' doing a variance-based preselection affect which features are deemed most important).
#' If `include_missing_features` is `TRUE`, missing features will
#' be assigned a weight of 0.
#'
#' @param output_list List of integration methods output, each generated via the
#' `get_output()` function. If named, the names will be added at the beginning of
#' each latent dimension' label. If unnamed, the name of the integration method
#' will be used instead.
#' @param include_missing_features Logical, whether features missing in some of the
#' output should be included in the calculation (see Details). Default value is `FALSE`.
#' @returns A correlation matrix.
#' @export
get_features_weight_correlation <- function(output_list, include_missing_features = FALSE) {
  ## For devtools::check
  latent_dimension <- weight <- NULL

  if (include_missing_features) {
    values_fill <- 0
  } else {
    values_fill <- NULL
  }

  res <- output_list |>
    .check_names_output_list() |>
    purrr::imap(
      ~ .x$features_weight |>
        dplyr::mutate(
          ## Need to make sure that the dimension names are unique
          ## Using "___" as a separator as user can supply custom labels,
          ## so we don't want a separator that they could use in their labels
          latent_dimension = as.character(latent_dimension),
          latent_dimension = stringr::str_c(.y, "___", latent_dimension)
        ) |>
        dplyr::select(latent_dimension, feature_id, weight)
    ) |>
    purrr::reduce(dplyr::bind_rows) |>
    tidyr::pivot_wider(
      names_from = latent_dimension,
      values_from = weight,
      values_fill = values_fill
    ) |>
    tibble::column_to_rownames("feature_id") |>
    stats::cor(use = "pairwise.complete.obs")

  return(res)
}

#' Get samples score correlation
#'
#' Constructs the correlation matrix between the samples score of the latent
#' dimensions obtained with different integration methods.
#'
#' @param output_list List of integration methods output, each generated via the
#' `get_output()` function. If named, the names will be added at the beginning of
#' each latent dimension' label. If unnamed, the name of the integration method
#' will be used instead.
#' @returns A correlation matrix.
#' @export
get_samples_score_correlation <- function(output_list) {
  ## For devtools::check
  latent_dimension <- score <- NULL

  res <- output_list |>
    .check_names_output_list() |>
    purrr::imap(
      ~ .x$samples_score |>
        dplyr::mutate(
          ## Need to make sure that the dimension names are unique
          ## Using "___" as a separator as user can supply custom labels,
          ## so we don't want a separator that they could use in their labels
          latent_dimension = as.character(latent_dimension),
          latent_dimension = stringr::str_c(.y, "___", latent_dimension)
        ) |>
        dplyr::select(latent_dimension, sample_id, score)
    ) |>
    purrr::reduce(dplyr::bind_rows) |>
    tidyr::pivot_wider(
      names_from = latent_dimension,
      values_from = score
    ) |>
    tibble::column_to_rownames("sample_id") |>
    stats::cor(use = "pairwise.complete.obs")

  return(res)
}

#' Heatmap of features weight correlation
#'
#' Constructs a lower triangle heatmap of features weight correlation between
#' the latent dimensions constructed by several integration methods.
#'
#' @param output_list List of integration methods output generated via the
#'   `get_output()` function.
#' @param include_missing_features Logical, see
#'   [get_features_weight_correlation()] for details. Default value is `FALSE`.
#' @param legend_ncol Integer, number of columns in the legend. Default value is
#'   `1`.
#' @returns a `ComplexHeatmap::Heatmap` (lower triangle only).
.heatmap_features_weight_corr <- function(output_list,
                                          include_missing_features = FALSE,
                                          legend_ncol = 1) {
  ## For devtools::check
  latent_dimension <- long_name <- short_name <- NULL

  ## Construct the matrix of features weight correlation across the
  ## latent dimensions returned by the different methods
  cor_features_weight <- get_features_weight_correlation(output_list, include_missing_features)

  ## Attributes of the latent dimensions
  ld_df <- tibble::tibble(
    long_name = rownames(cor_features_weight)
  ) |>
    tidyr::separate(
      long_name,
      into = c("method", "latent_dimension"),
      sep = "___",
      remove = FALSE
    ) |>
    dplyr::mutate(
      short_name = .get_initials(latent_dimension)
    )

  ## Latent dimensions labels to use in the plot
  ld_labels <- ld_df |>
    dplyr::select(
      long_name, short_name
    ) |>
    tibble::deframe()

  ## Row and columns annotations
  methods <- sort(unique(ld_df$method))
  methods_cols <- rep(
    RColorBrewer::brewer.pal(8, "Set2"),
    length.out = length(methods)
  )
  names(methods_cols) <- methods

  row_ha <- ComplexHeatmap::rowAnnotation(
    method = ld_df$method,
    col = list(method = methods_cols),
    gp = grid::gpar(col = "white"),
    annotation_legend_param = list(
      direction = "horizontal",
      ncol = legend_ncol,
      title_position = "topcenter"
    )
  )
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    method = ld_df$method,
    col = list(method = methods_cols),
    show_legend = FALSE,
    gp = grid::gpar(col = "white")
  )

  ## Clustering of rows and columns
  hclust_fw <- hclust_matrix_rows(abs(cor_features_weight))

  hll <- ComplexHeatmap::Heatmap(
    cor_features_weight,
    col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    ## Change latent dimension labels and position
    row_labels = ld_labels,
    column_labels = ld_labels,
    row_names_side = "left",
    ## Hide dendrograms
    cluster_rows = hclust_fw,
    cluster_columns = hclust_fw,
    show_row_dend = TRUE,
    show_column_dend = FALSE,
    row_dend_side = "left",
    ## Row and columns annotations
    left_annotation = row_ha,
    bottom_annotation = column_ha,
    ## Tweak legend
    name = "Correlation",
    heatmap_legend_param = list(
      direction = "horizontal",
      title_position = "topcenter"
    ),
    ## Adding title
    row_title = "Features weight correlation",
    ## Removing upper triangle
    rect_gp = grid::gpar(type = "none"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (as.numeric(x) <= 1 - as.numeric(y) + 1e-6) {
        grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = fill, col = "white"))
      }
    }
  )

  return(hll)
}

#' Heatmap of samples score correlation
#'
#' Constructs a upper triangle heatmap of samples score correlation between
#' the latent dimensions constructed by several integration methods.
#'
#' @param output_list List of integration methods output generated via the
#' `get_output()` function.
#' @param hclust_fw Dendrogram of latent dimensions according to the features
#' weight correlation (obtained with `.heatmap_features_weight_corr`).
#' @returns a `ComplexHeatmap::Heatmap` (upper triangle only).
.heatmap_samples_score_corr <- function(output_list, hclust_fw) {
  ## For devtools::check
  latent_dimension <- long_name_new <- long_name <- short_name <- NULL

  ## Construct the samples score correlation matrix between the latent dimensions
  ## of the different methods
  cor_samples_scores <- get_samples_score_correlation(output_list)

  ## Clustering of rows and columns
  hclust_sc <- hclust_matrix_rows(abs(cor_samples_scores))

  ## Matching row ordering to ordering provided
  row_ordering_indx <- stats::order.dendrogram(hclust_fw)
  sc_row_ordering_indx <- stats::order.dendrogram(hclust_sc)

  new_row_order <- sc_row_ordering_indx[order(row_ordering_indx)]

  temp <- rownames(cor_samples_scores)
  cor_samples_scores <- cor_samples_scores[new_row_order, new_row_order]
  rownames(cor_samples_scores) <- temp
  colnames(cor_samples_scores) <- temp

  ## Need to re-compute the dendrogram to print on top of the columns
  hclust_sc_new <- hclust_matrix_rows(abs(cor_samples_scores))

  lookup_names <- rlang::set_names(
    rownames(cor_samples_scores)[new_row_order],
    rownames(cor_samples_scores)
  )

  ## Latent dimension attributes
  ld_df_sc <- tibble::tibble(
    long_name_new = rownames(cor_samples_scores),
    long_name = lookup_names[long_name_new]
  ) |>
    tidyr::separate(
      long_name,
      into = c("method", "latent_dimension"),
      sep = "___",
      remove = FALSE
    ) |>
    dplyr::mutate(
      short_name = .get_initials(latent_dimension)
    )

  ld_labels_sc <- ld_df_sc |>
    dplyr::select(
      long_name_new, short_name
    ) |>
    tibble::deframe()

  ## Row and columns annotations
  methods <- sort(unique(ld_df_sc$method))
  methods_cols <- rep(
    RColorBrewer::brewer.pal(8, "Set2"),
    length.out = length(methods)
  )
  names(methods_cols) <- methods

  row_ha <- ComplexHeatmap::rowAnnotation(
    method = ld_df_sc$method,
    col = list(method = methods_cols),
    gp = grid::gpar(col = "white")
  )
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    method = ld_df_sc$method,
    col = list(method = methods_cols),
    show_legend = FALSE,
    gp = grid::gpar(col = "white")
  )

  hur <- ComplexHeatmap::Heatmap(
    cor_samples_scores,
    col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    ## Change latent dimension labels and position
    row_labels = ld_labels_sc,
    column_labels = ld_labels_sc,
    row_names_side = "right",
    column_names_side = "top",
    ## Hide dendrograms
    cluster_rows = hclust_fw,
    cluster_columns = hclust_sc_new,
    show_row_dend = FALSE,
    show_column_dend = TRUE,
    column_dend_side = "top",
    ## Row and columns annotations
    right_annotation = row_ha,
    top_annotation = column_ha,
    ## Tweak legend
    show_heatmap_legend = FALSE,
    ## Adding title
    column_title = "Samples score correlation",
    column_title_side = "top",
    ## Removing upper triangle
    rect_gp = grid::gpar(type = "none"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if (as.numeric(x) + 1e-6 >= 1 - as.numeric(y)) {
        grid::grid.rect(x, y, w, h, gp = grid::gpar(fill = fill, col = "white"))
      }
    }
  )

  return(hur)
}


#' Heatmap of correlation between latent dimensions
#'
#' Constructs a heatmap displaying the correlation between the latent dimensions
#' constructed with several integration methods. The lower triangle of the
#' heatmap displays the correlation between features weight, while the upper
#' triangle shows the correlation between samples score. Each triangle of the
#' matrix is reordered separately to show the most highly correlated dimensions
#' next to each other.
#'
#' If `output_list` is unnamed, the different elements in the list will be
#' differentiated by the name of the method used to produce them (e.g. DIABLO,
#' sO2PLS, etc). In order to compare different results from a same integration
#' method (e.g. DIABLO applied to the full vs pre-filtered data), it is possible
#' to assign names to the elements of `output_list` (see examples). These names
#' will be used in place of the method name in the plot to identify where the
#' latent dimensions come from.
#'
#' @param output_list List of integration methods output each generated via the
#'   `get_output()` function. If named, the names will be used to annotate the
#'   plot. See details.
#' @param latent_dimensions Named list, where each element is a character vector
#'   giving the latent dimensions to retain in the corresponding element of
#'   `output_list`. Names must match those of `output_list`. Can be used to
#'   filter latent dimensions only in certain elements from output_list (see
#'   examples). If `NULL` (default value), all latent dimensions will be used.
#' @param include_missing_features Logical, see
#'   [get_features_weight_correlation()] for details. Default value is `FALSE`.
#' @param legend_ncol Integer, number of columns in the legend. Default value is
#'   `1`.
#' @returns a `ComplexHeatmap::Heatmap`.
#' @examplesIf FALSE
#' ## Comparing the output from DIABLO, sO2PLS and MOFA
#'
#' res <- list(
#'   get_output_diablo(diablo_res), ## diablo_res: output from diablo_run()
#'   get_output_so2pls(so2pls_res), ## so2pls_res: output from so2pls_o2m()
#'   get_output_mofa2(mofa_res) ## mofa_res: output from run_mofa
#' )
#'
#' comparison_heatmap_corr(res)
#'
#' ## Selecting only some factors from a MOFA run for the comparison
#' ## (for the other methods, all latent dimensions will be retained)
#' comparison_heatmap_corr(
#'   res,
#'   latent_dimensions = list(
#'     "MOFA" = paste0("Factor ", 1:3)
#'   )
#' )
#'
#' ## Comparing two different results from a same integration method -
#' ## diablo_run_full and diablo_run_prefiltered would both be output
#' ## from the diablo_run() function.
#'
#' res <- list(
#'   "DIABLO full" = get_output_diablo(diablo_run_full),
#'   "DIABLO prefiltered" = get_output_diablo(diablo_run_prefiltered)
#' )
#'
#' comparison_heatmap_corr(res)
#' @export
comparison_heatmap_corr <- function(output_list,
                                    latent_dimensions = NULL,
                                    include_missing_features = FALSE,
                                    legend_ncol = 1) {
  output_list <- output_list |>
    .check_names_output_list() |>
    .filter_output_dimensions_list(latent_dimensions)

  ## Correlation heatmap of features weight
  hll <- .heatmap_features_weight_corr(
    output_list,
    include_missing_features,
    legend_ncol
  )
  dendro_fw <- hll@row_dend_param$obj

  ## Correlation heatmap of samples score
  hur <- .heatmap_samples_score_corr(output_list, dendro_fw)

  ## Adding padding around the heatmap titles
  ComplexHeatmap::ht_opt(TITLE_PADDING = grid::unit(c(15, 0), "points"))

  ComplexHeatmap::draw(
    hll + hur,
    ht_gap = grid::unit(-70, "mm"),
    heatmap_legend_side = "bottom"
  )
}


#' Correlation plot between latent components
#'
#' Plots the correlation between either the samples score or the features
#' weight of the latent components obtained from two different integration methods.
#'
#' If `output_list` is unnamed, the different elements in the list will be differentiated
#' by the name of the method used to produce them (e.g. DIABLO, sO2PLS, etc). In
#' order to compare different results from a same integration method (e.g. DIABLO
#' applied to the full vs pre-filtered data), it is possible to assign names to
#' the elements of `output_list`. These names will be used in place
#' of the method name in the plot to identify where the latent dimensions come from.
#'
#' @param output_list List of length 2 of integration methods output, each generated via the
#' `get_output()` function. If named, the names will be used to annotate the plot.
#' See details.
#' @param by Character, should the correlation be calculated based on samples score (
#' `by = 'samples'`) or on features weight (`by = 'features'`), or both (`by = 'both'`,
#' i.e. two matrices will be plotted). Default value is `'both'`.
#' @param latent_dimensions Named list, where each element is a character vector
#' giving the latent dimensions to retain in the corresponding element of `output_list`.
#' Names must match those of `output_list`. Can be used to filter latent dimensions
#' only in certain elements from output_list (see examples). If `NULL` (default value),
#' all latent dimensions will be used.
#' @param include_missing_features Logical, see \code{\link{get_features_weight_correlation}}
#' for details. Default value is `FALSE`.
#' @inheritParams plot_correlation_matrix_full
#' @returns A ggplot.
#' @export
comparison_plot_correlation <- function(output_list,
                                        by = "both",
                                        latent_dimensions = NULL,
                                        include_missing_features = FALSE,
                                        show_cor = TRUE,
                                        min_show_cor = 0.2,
                                        round_cor = 2) {
  if (by == "both") {
    res <- patchwork::wrap_plots(
      comparison_plot_correlation(
        output_list,
        by = "samples",
        latent_dimensions,
        include_missing_features,
        show_cor,
        min_show_cor,
        round_cor
      ),
      comparison_plot_correlation(
        output_list,
        by = "features",
        latent_dimensions,
        include_missing_features,
        show_cor,
        min_show_cor,
        round_cor
      )
    ) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x.top = ggplot2::element_text(angle = 90)
      )

    return(res)
  }

  ## For devtools::check
  ld <- method <- data <- NULL

  if (!(by %in% c("samples", "features"))) {
    stop("'by' argument should be one of 'samples' or 'features'.")
  }

  if (length(output_list) != 2) {
    stop("'output_list' should be of length 2; this function is for comparing the output of two integration methods.")
  }

  output_list <- output_list |>
    .check_names_output_list() |>
    .filter_output_dimensions_list(latent_dimensions)

  cormat <- switch(
    by,
    samples = get_samples_score_correlation(output_list),
    features = get_features_weight_correlation(output_list, include_missing_features)
  )

  ld_list <- tibble::tibble(
    ld = rownames(cormat)
  ) |>
    tidyr::separate(
      ld,
      into = c("method", "label"),
      sep = "___",
      remove = FALSE
    ) |>
    dplyr::group_by(method) |>
    tidyr::nest() |>
    dplyr::mutate(
      data = purrr::map(data, tibble::deframe)
    ) |>
    tibble::deframe()

  cormat <- cormat[names(ld_list[[1]]), names(ld_list[[2]]), drop = FALSE]
  rownames(cormat) <- unname(ld_list[[1]])
  colnames(cormat) <- unname(ld_list[[2]])

  p <- plot_correlation_matrix_full(
    cormat,
    names(ld_list)[[1]],
    names(ld_list)[[2]],
    dplyr::if_else(
      by == "samples",
      "Samples score correlation",
      "Features weight correlation"
    ),
    show_cor = show_cor,
    min_show_cor = min_show_cor,
    round_cor = round_cor
  )

  return(p)
}

#' Calculate features importance score
#'
#' Given a vector of features importance (between 1 and 0), returns a consensus importance
#' score.
#'
#' @param x Numeric vector of importance values (between 0 and 1).
#' @param metric Character, one of the metrics to use to compute the consensus score.
#' Can be one of `'min'`, `'max'`, `'average'`, `'product'`, `'l2'` (for L2-norm),
#' `'geometric'` (for geometric mean) or `'harmonic'` (for harmonic mean).
#' @param offset Numeric (strictly positive), used to replace zero values to compute the geometric and
#' harmonic mean. If `0` (default value), the zero values will be ignored when
#' calculating the geometric and harmonic mean. Accepting vector of values to facilitate use
#' whithin `dplyr::summarise()`.
#' @returns A numeric value, the importance consensus score.
#' @export
consensus_importance_metric <- function(x, metric, offset = 0) {
  if (metric %in% c("geometric", "harmonic")) {
    x[x == 0] <- min(offset)
  }

  res <- switch(metric,
                min = min(x, na.rm = TRUE),
                max = max(x, na.rm = TRUE),
                average = mean(x, na.rm = TRUE),
                product = prod(x, na.rm = TRUE),
                l2 = sqrt(sum(x^2, na.rm = TRUE)),
                geometric = exp(mean(log(x), na.rm = TRUE)),
                harmonic = 1 / mean(1 / x, na.rm = TRUE),
                stop("'metric' should be one of 'min', 'max', 'average', 'product', 'l2', 'geometric' or 'harmonic'.")
  )

  return(res)
}

#' Computes consensus feature importance
#'
#' Computes the consensus feature importance from features weight obtained with
#' different integration methods (considering features importance for one latent component
#' per integration method), or for different latent dimensions constructed by an
#' integration method.
#'
#' If `include_missing_features` is `FALSE` (default behaviour), and some
#' features are present in the output of one integration method but not
#' the other (e.g. because a different pre-filtering was applied to the input
#' data of the two methods), these features will be ignored. This does not mean
#' that features that were selected by one method but not the other are discarded;
#' in that case the feature will be assigned a weight of 0 for the method that did
#' not select it. This is the recommended behaviour, should only be changed in
#' specific scenarios (e.g. to check whether using all features in a dataset vs
#' doing a variance-based preselection affect which features are deemed most important).
#' If `include_missing_features` is `TRUE`, missing features will
#' be assigned a weight of 0.
#' Note that the geometric and harmonic means only work for strictly positive values.
#' Therefore, all importance scores of 0 are replaced with an offset when computing
#' these metrics. The offset is calculated per dataset, and corresponds to the minimum
#' non-null importance score observed across all features in the dataset (and across
#' all latent dimensions), divided by 2. The calculation of the offset is done before
#' removing missing features (if `include_missing_features = FALSE`) so that results
#' are consistent between the two options for `include_missing_features`.
#'
#' @param output_list List of integration methods output, each generated via the
#' `get_output()` function, or a single integration method output (from `get_output()`).
#' @param latent_dimensions Named list (if `output_list` is a list), where each element is a character
#' giving the latent dimension to retain in the corresponding element of `output_list` (1 value).
#' If `output_list` is a single output object, needs to be instead a character vector giving the
#' latent dimensions to retain.
#' @param metric Character, one of the metrics to use to compute the consensus score.
#' Can be one of `'min'`, `'max'`, `'average'`, `'product'`, `'l2'` (for L2-norm),
#' `'geometric'` (for geometric mean) or `'harmonic'` (for harmonic mean). Default
#' value is `'geometric'`. Names must match those of `output_list`.
#' @param include_missing_features Logical, whether features missing in some of the
#' output should be included in the calculation (see Details). Default value is `FALSE`.
#' @returns A tibble giving the consensus importance of each feature.
#' @export
compute_consensus_importance <- function(output_list, latent_dimensions, metric = "geometric", include_missing_features = FALSE) {

  ## For devtools::check
  importance <- dataset <- feature_id <- latent_dimension <- method <- NULL

  ## If input is only one method
  if (inherits(output_list, "output_dimension_reduction")) {
    output_list <- output_list |>
      .filter_output_dimensions(latent_dimensions)

    df <- output_list$features_weight |>
      dplyr::select(method = latent_dimension, dataset, feature_id, importance)

  } else { ## is input is a list of methods
    output_list <- output_list |>
      .check_names_output_list() |>
      .filter_output_dimensions_list(
        latent_dimensions,
        all_present = TRUE,
        fixed_length = 1
      )

    df <- output_list |>
      purrr::map_dfr(
        ~ .x$features_weight |>
          dplyr::select(dataset, feature_id, importance),
        .id = "method"
      )
  }

  ## Extracting minimum non-null value from each dataset to use as offset
  min_vals <- df |>
    dplyr::filter(importance > 0) |>
    dplyr::group_by(dataset) |>
    dplyr::summarise(
      min = min(importance),
      .groups = "drop"
    ) |>
    dplyr::mutate(min = min / 2) |>
    tibble::deframe()

  ## We want to complete so that if include_missing_features is:
  ## - FALSE, we can get rid of the features with missing scores for some methods
  ## - TRUE, we replace the missing values with 0
  df <- df |>
    tidyr::complete(tidyr::nesting(dataset, feature_id), method)

  if (include_missing_features) {
    df <- df |>
      tidyr::replace_na(list(importance = 0))
  } else {
    ## Getting rid of all features that have some missing importance score
    missing_features <- df |>
      dplyr::filter(is.na(importance)) |>
      dplyr::pull(feature_id) |>
      unique()

    df <- df |>
      dplyr::filter(!(feature_id %in% missing_features))
  }

  res <- df |>
    dplyr::group_by(dataset, feature_id) |>
    dplyr::summarise(
      importance = consensus_importance_metric(importance, metric, offset = min_vals[dataset]),
      .groups = "drop_last"
    ) |>
    dplyr::mutate(
      importance = importance / max(importance, na.rm = TRUE)
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(dataset, dplyr::desc(importance))

  return(res)
}



#' Illustrates importance consensus metrics
#'
#' Plots a heatmap to illustrate the behaviour of different importance consensus metrics.
#'
#' @param metrics Character vector of metrics to show. Should be valid values for
#' the `metric` argument of [consensus_importance_metric()], i.e. "min", "max",
#' "average", "product", "l2", "geometric", "harmonic".
#' @returns A ggplot.
#' @export
show_consensus_metrics <- function(metrics = c("min", "harmonic", "geometric", "product", "average", "l2", "max")) {

  ## for devtools::check
  x <- y <- axis <- importance <- metric <- value <- NULL

  toplot <- tidyr::expand_grid(
    x = seq(0, 1, length.out = 50),
    y = seq(0, 1, length.out = 50)
  ) |>
    dplyr::mutate(
      rowline = seq_len(dplyr::n())
    ) |>
    tidyr::pivot_longer(
      cols = c(x, y),
      names_to = "axis",
      values_to = "importance"
    )

  metrics |>
    purrr::map_dfr(
      ~ toplot |>
        dplyr::group_by(rowline) |>
        dplyr::mutate(
          metric = .x,
          value = consensus_importance_metric(importance, .x)
        ) |>
        dplyr::ungroup()
    ) |>
    tidyr::pivot_wider(
      names_from = axis,
      values_from = importance
    ) |>
    dplyr::mutate(metric = factor(metric, levels = metrics)) |>
    dplyr::group_by(metric) |>
    dplyr::mutate(value = value / max(value, na.rm = TRUE)) |>
    ggplot2::ggplot(aes(x = x, y = y, z = value)) +
    ggplot2::geom_tile(aes(fill = value)) +
    ggplot2::geom_contour(colour = "hotpink") +
    ggplot2::facet_wrap(~metric, nrow = 2) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion()) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion()) +
    ggplot2::scale_fill_viridis_c(na.value = "gray") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      title = "Overview of metrics for consensus feature importance",
      x = "Features importance - method 1",
      y = "Features importance - method 2",
      fill = "Consensus importance"
    )

}
