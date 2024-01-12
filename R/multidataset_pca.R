#' Run PCA on MultiDataSet
#'
#' Runs a Principal Component Analysis on an omics dataset from a MultiDataSet
#' object. This is a wrapper function around the [get_dataset_matrix()] and
#' [run_pca_matrix()] functions.
#'
#' To facilitate the use of dynamic branching with the `targets` package, the
#' `dataset_name` attribute of the resulting object is set as the value of the
#' `dataset_name` parameter, and can be accessed via `attr(res_pca,
#' "dataset_name")`.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param dataset_name Character, name of the omics dataset on which a PCA
#'   should be run.
#' @inheritParams run_pca_matrix
#' @returns A [pcaMethods::pcaRes-class] object containing the result from the
#'   PCA analysis. The attribute `dataset_name` specifies the name of the
#'   dataset analysed.
#' @export
run_pca <- function(mo_data,
                    dataset_name,
                    n_pcs = 10,
                    scale = "none",
                    center = TRUE,
                    method = NULL) {
  mat <- get_dataset_matrix(mo_data, dataset_name, keep_dataset_name = TRUE)
  res <- run_pca_matrix(mat, n_pcs, scale, center, method)

  return(res)
}

#' Run PCA on matrix
#'
#' Runs a Principal Component Analysis on an omics matrix, using the
#' [pcaMethods::pca()] function.
#'
#' @param mat Matrix of omics measurement, with features as rows and samples as
#'   columns.
#' @param n_pcs numeric, number of Principal Components to compute. Default
#'   value is 10.
#' @param scale character, type of scaling that should be applied to the dataset
#'   before running the PCA. Should be one of `'none'`, `'pareto'`, `'vector'`,
#'   `'uv'` (see [pcaMethods::pca()]). Default value is \code{'none'}.
#' @param center boolean, should the dataset be centred prior to running the
#'   PCA? Default value is \code{TRUE}.
#' @param method character, type of PCA that should be applied to the dataset.
#'   See [pcaMethods::listPcaMethods()]. for a list of available methods.
#'   Default value is `'svd'` for datasets with no missing value, and `'nipals'`
#'   for datasets with missing values.
#' @returns A [pcaMethods::pcaRes-class] object containing the result from the
#'   PCA analysis. The attribute `dataset_name` specifies the name of the
#'   dataset analysed.
#' @export
run_pca_matrix <- function(mat,
                           n_pcs = 10,
                           scale = "none",
                           center = TRUE,
                           method = NULL) {
  if (n_pcs < 1) {
    stop("'n_pcs' argument should be an integer value >= 1.")
  }
  .check_names(
    scale,
    c("none", "pareto", "vector", "uv"),
    "'scale' argument should be one of '_C_'."
  )
  if (!is.logical(center)) {
    stop("'center' argument should be either TRUE or FALSE.")
  }
  if (!is.null(method)) {
    .check_names(
      method,
      pcaMethods::listPcaMethods(),
      "'method' argument should be one of: '_C_'."
    )
  }
  if (n_pcs > ncol(mat)) {
    warning(
      "'n_pcs' argument larger than the number of observations ",
      "- setting n_pcs to ", ncol(mat)
    )
    n_pcs <- ncol(mat)
  }

  ## Choosing appropriate method if not provided
  if (is.null(method)) {
    method <- ifelse(any(is.na(mat)), "nipals", "svd")
  }

  ## Want to keep track of the dataset name without this attribute being
  ## used in the PCA results
  dataset_name <- attr(mat, "dataset_name")
  attr(mat, "dataset_name") <- NULL

  ## Running the PCA
  res <- pcaMethods::pca(
    t(mat),
    method = method,
    nPcs = n_pcs,
    scale = scale,
    center = center
  )

  ## to keep track of which dataset was analysed
  attr(res, "dataset_name") <- dataset_name

  return(res)
}

#' Target factory for PCA run and missing values imputation on each omics
#' dataset
#'
#' Creates a list of targets that perform a PCA run for each omics dataset from
#' a `MultiDataSet` object using dynamic branching, and imputes the missing
#' values in those datasets using the results of the PCA runs.
#'
#' @param mo_data_target Symbol, the name of the target containing the
#'   `MultiDataSet` object.
#' @param dataset_names Character vector, the names of the datasets on which a
#'   PCA should be run. If `NULL`, a PCA will be run on all datasets. Default
#'   value is `NULL`.
#' @param target_name_prefix Character, prefix to add to the name of the targets
#'   created by the factory. Default value is `""`.
#' @param complete_data_name Character, the name of the target containing the
#'   `MultiDataSet` with missing data imputed to be created. If `NULL`, will be
#'   selected automatically. Default value is `NULL`.
#' @param ... Further arguments passed to the [run_pca_matrix()] function.
#' @returns A List of targets. If `target_name_prefix = ""` and
#'   `complete_data_name = NULL`, the following targets are created:
#' * `dataset_names_pca`: target containing a character vector that gives the
#' names of the datasets on which a PCA should be run.
#' * `dataset_mats_pca`: a dynamic branching target that applies the
#' [get_dataset_matrix()] function to each dataset specified in `dataset_names`.
#' The results are saved in a list. Note that because it is using dynamic
#' branching, the names of the list are not meaningful. Rather, use
#' `sapply(pca_pca_runs_listruns_list, attr, "dataset_name")` to assess which
#' element of the list corresponds to which omics dataset.
#' * `pca_runs_list`: a dynamic branching target that applies the
#' [run_pca_matrix()] function to each matrix in `dataset_mats_pca`. The results
#' are saved in a list. Note that because it is using dynamic branching, the
#' names of the list are not meaningful. Rather, use `sapply(pca_runs_list,
#' attr, "dataset_name")` to assess which element of the list corresponds to
#' which omics dataset.
#' * `complete_set`: a target that returns a `MultiDataSet` in which missing
#' values have been imputed.
#' @examples
#' \dontrun{
#' ## in the _targets.R
#' library(moiraine)
#' list(
#'   # ... code for importing datasets etc
#'
#'   ## mo_set is the target containing the MultiDataSet object
#'   ## Example 1: running a PCA on all datasets
#'   run_pca_factory(mo_set),
#'
#'   ## Example 2: running a PCA on 'rnaseq' and 'metabolome' datasets
#'   run_pca_factory(
#'     mo_set,
#'     c("rnaseq", "metabolome"),
#'     complete_data_name = "mo_data_complete"
#'   )
#' )
#' }
#' @export
pca_complete_data_factory <- function(mo_data_target,
                                      dataset_names = NULL,
                                      target_name_prefix = "",
                                      complete_data_name = NULL,
                                      ...) {
  ## can't use is.null(dataset_names) otherwise the expression in dataset_names
  ## gets evaluated
  if (is.null(substitute(dataset_names))) {
    dsn <- substitute(names(mo_data_target))
  } else {
    dsn <- substitute(dataset_names)
  }

  ## Target names
  dsn_name <- paste0(target_name_prefix, "dataset_names_pca")
  pca_mat_name <- paste0(target_name_prefix, "pca_mats_list")
  pca_run_name <- paste0(target_name_prefix, "pca_runs_list")
  if (is.null(complete_data_name)) {
    complete_data_name <- paste0(target_name_prefix, "complete_set")
  }

  ## Target symbols
  dsn_target <- as.symbol(dsn_name)
  pca_mat_target <- as.symbol(pca_mat_name)
  pca_run_target <- as.symbol(pca_run_name)

  targets <- list(
    targets::tar_target_raw(
      dsn_name,
      dsn
    ),
    targets::tar_target_raw(
      pca_mat_name,
      substitute(
        get_dataset_matrix(
          mo_data_target,
          dsn_target,
          keep_dataset_name = TRUE
        )
      ),
      pattern = substitute(map(dsn_target)),
      iteration = "list"
    ),
    targets::tar_target_raw(
      pca_run_name,
      substitute(run_pca_matrix(pca_mat_target, ...)),
      pattern = substitute(map(pca_mat_target)),
      iteration = "list"
    ),
    targets::tar_target_raw(
      complete_data_name,
      substitute(get_complete_data(mo_data_target, pca_run_target))
    )
  )

  return(targets)
}

#' Extract arguments used in PCA run
#'
#' Extracts the list of arguments used for each PCA run from a list of PCA
#' results, and formats them into a tibble.
#'
#' @param pca_result The result of a PCA run on each of the datasets, computed
#'   with the [pcaMethods::pca()] function.
#' @returns A tibble with the following columns: `"Omics dataset"`, `"PCA method
#'   used"`, `"Number of Principal Components computed"`, `"Scaling applied"`
#'   and `"Dataset centered"`.
#' @export
get_pca_arguments <- function(pca_result) {
  if (!is.list(pca_result) || !all(sapply(pca_result, class) == "pcaRes")) {
    stop("Input should be a list of pcaRes objects (from pcaMethods package).")
  }

  ## Make sure that the names of the pca_result list are the name of the
  ## datasets
  names(pca_result) <- purrr::map_chr(pca_result, attr, "dataset_name")

  res <- pca_result |>
    purrr::map(function(x) {
      tibble::tibble(
        "Omics dataset" = attr(x, "dataset_name"),
        "PCA method used" = x@method,
        "Number of Principal Components computed" = x@nPcs,
        "Scaling applied" = x@scaled,
        "Dataset centered" = x@centered
      )
    }) |>
    purrr::reduce(dplyr::bind_rows)

  return(res)
}

#' Screeplots for single-omics PCA
#'
#' Produces a scree plot (percentage of variance explained by each principal
#' component) for the PCA run on each omics dataset, using ggplot2.
#'
#' @param pca_result List of PCA runs result on each of the datasets, each
#'   computed with the [run_pca()] function.
#' @param cumulative Logical, should the cumulative variance be plotted? Default
#'   is `FALSE`.
#' @param datasets Optional, character vector with the names of the datasets to
#'   plot. If `NULL` (default value), all datasets will be plotted.
#' @returns A ggplot2 plot.
#' @export
plot_screeplot_pca <- function(pca_result,
                               cumulative = FALSE,
                               datasets = NULL) {
  if (!is.list(pca_result) || !all(sapply(pca_result, class) == "pcaRes")) {
    stop("Input should be a list of pcaRes objects (from pcaMethods package).")
  }

  ## Make sure that the names of the pca_result list are the name of the
  ## datasets
  names(pca_result) <- purrr::map_chr(pca_result, attr, "dataset_name")

  if (is.null(datasets)) datasets <- names(pca_result)
  .check_names(
    datasets,
    names(pca_result),
    "'datasets' argument: '_W_' not existing datasets. Possible values are: '_C_'."
  )

  pca_result <- pca_result[datasets]

  ## To avoid check() problems
  dataset <- pc <- r2 <- label <- NULL

  lapply(datasets, function(i) {
    x <- pca_result[[i]]

    tibble::tibble(
      dataset = i,
      pc = 1:x@nPcs,
      r2 = dplyr::case_when(
        !cumulative ~ x@R2,
        TRUE ~ x@R2cum
      ),
      label = paste0(round(100 * r2, 1), "%")
    )
  }) |>
    purrr::reduce(dplyr::bind_rows) |>
    dplyr::mutate(
      dataset = factor(dataset, levels = datasets),
      pc = factor(pc, levels = 1:max(pc))
    ) |>
    ggplot2::ggplot(aes(x = pc, y = r2)) +
    ggplot2::geom_col(aes(fill = dataset)) +
    ggplot2::geom_text(aes(label = label), vjust = -0.15, size = 3) +
    ggplot2::facet_wrap(~dataset, ncol = 2, scales = "free_y") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(
      title = ifelse(
        cumulative,
        "Principal Component Analysis - Cumulative scree plots",
        "Principal Component Analysis - Scree plots"
      ),
      x = "Principal Components",
      y = ifelse(
        cumulative,
        "Cumulative percentage of variance explained",
        "Percentage of variance explained"
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )
}


#' Samples score plots for single-omics PCA
#'
#' Produces a pairwise samples score plot for the PCA run on each omics dataset,
#' using [GGally::ggpairs()].
#'
#' @param pca_result List of PCA results on each of the datasets, each computed
#'   with the [run_pca()] function.
#' @param pcs Integer vector or named list of integer vectors, the principal
#'   components to display for each dataset. If integer vector (e.g. `1:5`),
#'   will be used for all datasets. Alternatively, a different set of PCs can be
#'   specified through a named list (e.g. `list('snps' = 1:4, 'rnaseq' = 1:5)`).
#'   The length of the list must match the number of datasets to be displayed,
#'   and the names must match the dataset names. Default value is `NULL`, i.e.
#'   all principal components will be plotted for each dataset.
#' @param datasets Optional, character vector of datasets for which the plots
#'   should be created.
#' @param ... Other arguments passed to [plot_samples_score()].
#' @returns A list of ggmatrix plots (or a single ggmatrix plot if `pcs` has
#'   length 1).
#' @examples
#' \dontrun{
#' ## Default: plotting all PCs for all datasets
#' plot_samples_coordinates_pca(pca_result)
#'
#' ## Plotting only the first 3 PCs for each dataset
#' plot_samples_coordinates_pca(
#'   pca_result,
#'   pcs = 1:3
#' )
#'
#' ## Plotting the first 3 PCs for the genomics dataset, 4 PCs for the
#' ## transcriptomics dataset, 5 PCs for the metabolomics dataset
#' plot_samples_coordinates_pca(
#'   pca_result,
#'   pcs = list(
#'     "snps" = 1:3,
#'     "rnaseq" = 1:4,
#'     "metabolome" = 1:5
#'   )
#' )
#'
#' ## Plotting the first 3 PCs for the genomics and transcriptomics datasets
#' plot_samples_coordinates_pca(
#'   pca_result,
#'   pcs = 1:3,
#'   datasets = c("snps", "rnaseq")
#' )
#'
#' # Plotting the first 3 PCs for the genomics dataset and 4 PCs for the
#' ## transcriptomics dataset (no plot for the metabolomics dataset)
#' plot_samples_coordinates_pca(
#'   pca_result,
#'   pcs = list(
#'     "snps" = 1:3,
#'     "rnaseq" = 1:4
#'   ),
#'   datasets = c("snps", "rnaseq")
#' )
#' }
#' @export
plot_samples_coordinates_pca <- function(pca_result,
                                         pcs = NULL,
                                         datasets = NULL,
                                         ...) {
  purrr::walk(
    pca_result,
    function(.x) {
      if (!inherits(.x, "pcaRes")) {
        stop("Expecting a pcaRes object (from run_pca() or pcaMethods::pca() function).")
      }
      return(NULL)
    }
  )

  names(pca_result) <- sapply(pca_result, attr, "dataset_name")

  if (is.null(datasets)) {
    datasets <- names(pca_result)
  } else {
    .check_names(
      datasets,
      names(pca_result),
      "'datasets' argument: '_W_' not existing datasets. Possible values are: '_C_'."
    )
    pca_result <- pca_result[datasets]
  }

  ## Get number of PCs per dataset
  if (is.null(pcs)) {
    pcs <- pca_result |>
      purrr::map(pcaMethods::nPcs) |>
      purrr::map(seq_len)
  }

  ## Get PCA output in standard format
  pca_result <- pca_result |>
    purrr::map(get_output_pca) |>
    purrr::map(
      function(.x) {
        attr(.x, "method") <- paste0(
          attr(.x, "method"),
          ", ",
          attr(.x, "dataset"),
          " dataset"
        )
        return(.x)
      }
    )

  pcs <- .make_var_list(pcs, datasets) |>
    purrr::map(~ paste0("Principal component ", .x))

  fct <- function(.x, ...) {
    plot_samples_score(
      method_output = pca_result[[.x]],
      latent_dimensions = pcs[[.x]],
      ...
    )
  }

  res <- purrr::map(
    datasets,
    fct,
    ...
  )
  names(res) <- datasets

  if (length(res) == 1) res <- res[[1]]

  return(res)
}


#' Get MultiDataSet object with imputed values
#'
#' Replace missing values with imputed values for each dataset of a MultiDataSet
#' object, based on the results of a Principal Component Analysis applied to the
#' corresponding dataset.
#'
#' Uses the [pcaMethods::completeObs()] function to impute missing values.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param pca_result A list in which each element is the result of a PCA run on
#'   a different dataset, computed with the [run_pca()] function.
#' @return A [MultiDataSet::MultiDataSet-class] object, for which the assay of
#'   each dataset is the imputed dataset.
#'
#' @export
get_complete_data <- function(mo_data, pca_result) {
  check_is_multidataset(mo_data)
  if (!is.list(pca_result) || !all(sapply(pca_result, class) == "pcaRes")) {
    stop("'pca_result' should be a list of pcaRes objects (from pcaMethods package).")
  }

  ## Make sure that the names of the pca_result list are the name of the
  ## datasets
  names(pca_result) <- purrr::map_chr(pca_result, attr, "dataset_name")

  ## Checking which datasets have missing data
  datasets <- get_datasets(mo_data)
  any_nas <- sapply(datasets, function(x) {
    any(is.na(x))
  })

  if (!any(any_nas)) { ## if all datasets are complete
    message("All datasets are complete, no imputation to perform.")
    return(mo_data)
  }

  with_nas <- names(any_nas)[any_nas]

  if (length(setdiff(with_nas, names(pca_result)))) {
    warning(
      "PCA results not available for dataset(s) '",
      paste0(setdiff(with_nas, names(pca_result)), collapse = "', '"),
      "'. NA values will not be imputed."
    )
  }

  with_nas <- intersect(with_nas, names(pca_result))

  res <- mo_data

  for (i in with_nas) {
    new_assay <- t(pcaMethods::completeObs(pca_result[[i]]))

    res <- replace_dataset(res, i, new_assay)
  }

  return(res)
}
