#' Generate MOFA input data
#'
#' Creates an object that can be used as input for the MOFA analysis implemented in the MOFA2 package. It contains the omics datasets
#' as well as features and samples metadata.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param datasets Character vector, the names of the datasets from \code{mo_data} to include in the analysis.
#' @param groups Character, the column name in the samples metadata data-frames to use as groups
#' (use \code{\link{get_samples_metadata}} to view the samples metadata data-frame for each dataset).
#' WARNING: only use if you are familiar with MOFA and its use of groups. See \url{https://biofam.github.io/MOFA2/faq.html},
#' section "FAQ on the multi-group functionality".
#' @param options_list A named list. Should contain at most 3 elements, named 'data_options', 'model_options' and 'training_options'. Provide respectively the
#' data, model and training options to apply for the MOFA run. See \code{\link[MOFA2]{get_default_data_options}}, \code{\link[MOFA2]{get_default_model_options}}
#' and \code{\link[MOFA2]{get_default_training_options}}.
#' @param only_common_samples Logical, whether only the samples present in all datasets should be returned.
#' Default value is \code{FALSE}.
#' @return A \code{\link[MOFA2]{MOFA}} object.
#' @export
get_input_mofa <- function(mo_data, datasets = names(mo_data), groups = NULL, options_list = NULL, only_common_samples = FALSE) {
  get_input_mofa2(mo_data, datasets, covariates = NULL, groups, options_list, only_common_samples)
}


#' Get feature weights from MOFA object
#'
#' Extracts the feature weights from a trained MOFA or MEFISTO model (from the MOFA2 package).
#'
#' Scaling options:
#' * `scale = 'none'`: no scaling performed;
#' * `scale = 'by_view'`: weights are divided by the maximum absolute weight in the corresponding view/dataset;
#' * `scale = 'by_factor'`: weights are divided by the maximum absolute weight in the corresponding factor;
#' * `scale = 'overall'`: weights are divided by the maximum absolute weight across all views/datasets and
#' factors considered.
#'
#' @param object A trained MOFA object.
#' @param views Character or integer vector, name or index of the views (i.e. datasets) for which the feature weights
#' should be extracted. Default value is "all", i.e. all datasets are considered.
#' @param factors Character or integer vector, name or index of the factors for which the feature weights
#' should be extracted. Default value is "all", i.e. all factors are considered.
#' @param abs Logical, should the absolute value of the weights be returned? Default value is `FALSE`.
#' @param scale Character, the type of scaling that should be performed on the feature weights. Possible
#' values are `'none'`, `'by_view'`, `'by_factor'` or `'overall'` (see Details). Default value is `'none'`.
#' @param as.data.frame Logical, whether the function should return a long data-frame instead of a list
#' of matrices. Default value is `TRUE`.
#' @return By default, returns a tibble with columns `view`, `feature`, `factor`, `value`. Alternatively, if
#' `as.data.frame = FALSE`, returns a list of matrices, one per view, with features as rows and factors as columns.
#' @export
mofa_get_weights <- function(object, views = "all", factors = "all", abs = FALSE, scale = "none", as.data.frame = TRUE) {
  if (!is(object, "MOFA")) {
    stop("'object' has to be an instance of MOFA")
  }

  if (!(scale %in% c("none", "by_view", "by_factor", "overall"))) stop("scale should be one of: 'none', 'by_view', 'by_factor', 'overall'.")

  views <- MOFA2:::.check_and_get_views(object, views)
  factors <- MOFA2:::.check_and_get_factors(object, factors)

  ## Get the raw weights as a long data-frame by default,
  ## and filter the relevant views and factors
  view <- factor <- value <- NULL
  weights <- MOFA2::get_expectations(object, "W", as.data.frame = TRUE) |>
    dplyr::filter(
      view %in% views,
      factor %in% factors
    )

  ## Transform to absolute value if needed
  if (abs) weights$value <- abs(weights$value)

  ## If some scaling must be performed
  if (scale != "none") {
    if (scale == "by_view") weights <- dplyr::group_by(weights, view)
    if (scale == "by_factor") weights <- dplyr::group_by(weights, factor)

    ## When scaling by view or factor, the grouping ensures that the maximum
    ## value is selected within the relevant group (i.e. either view or factor)
    ## otherwise if no grouping is performed (when scale = 'overall'), takes
    ## the max absolute value across the entire dataset.
    weights <- weights |>
      dplyr::mutate(value = value / max(abs(value))) |>
      dplyr::ungroup()
  }

  ## If we want to return a list, transform the long data-frame as a list
  ## of tibbles (can add as.data.frame() at the end if we don't want a tibble)
  if (!as.data.frame) {
    weights <- purrr::map(
      views,
      ~ weights |>
        dplyr::filter(view == .x) |>
        tidyr::pivot_wider(
          names_from = factor,
          values_from = value
        ) |>
        dplyr::select(-view)
    )

    names(weights) <- views
  }

  return(weights)
}


#' Plots the correlation between factors and covariates for MOFA
#'
#' Plots the Pearson correlation between MOFA latent factors and covariates
#' (obtained from the samples metadata). This function provides the `ggplot2`
#' version of the plot created with \code{\link[MOFA2]{correlate_factors_with_covariates}}
#' (with the `plot` parameter set to `'r'`).
#'
#' @param mofa_output The output from \code{\link[MOFA2]{run_mofa}}.
#' @param covariates Character vector, the covariates to use for the plot. If `NULL`,
#' all covariates retrieved via `colnames(MOFA2::samples_metadata(mofa_output))`
#' (except `group`, `id` and `sample`) will be used. Default value is `NULL`.
#' @inheritParams plot_correlation_matrix_full
#' @param factor_as_col Logical, should the factors be represented in columns?
#' If `FALSE`, will be represented in rows instead. Default value is `TRUE`.
#' @return A ggplot.
#' @export
mofa_plot_cor_covariates <- function(mofa_output, covariates = NULL, show_cor = TRUE, min_show_cor = 0.2, round_cor = 2, factor_as_col = TRUE) {
  if (!requireNamespace("ggforce", quietly = TRUE)) {
    stop(
      "Package \"ggforce\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!("MOFA" %in% class(mofa_output))) stop("Expecting MOFA object")

  ## for devtools::check
  Factor <- Covariate <- corcoef <- rad <- label <- NULL

  if (is.null(covariates)) {
    covariates <- setdiff(
      colnames(MOFA2::samples_metadata(mofa_output)),
      c("group", "id", "sample")
    )
  }

  covariates <- sort(covariates)

  cor_mat <- suppressWarnings(MOFA2::correlate_factors_with_covariates(
    mofa_output,
    covariates = covariates,
    plot = "r",
    return_data = TRUE
  ))
  rownames(cor_mat) <- stringr::str_remove(rownames(cor_mat), "Factor")

  if (factor_as_col) {
    cor_mat <- t(cor_mat)
    rtitle <- NULL
    ctitle <- "Factors"
  } else {
    rtitle <- "Factors"
    ctitle <- NULL
  }

  p <- plot_correlation_matrix_full(
    cor_mat,
    rows_title = rtitle,
    cols_title = ctitle,
    title = "Pearson correlation between MOFA factors and covariates"
  ) +
    ggplot2::theme(
      axis.text.x.top = ggplot2::element_text(angle = 0, hjust = 0.5)
    )

  return(p)
}


