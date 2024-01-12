#' Generate mixomics input data for unsupervised methods
#'
#' Creates an object that can be used as input for the MixOmics package.
#' It contains the omics datasets restricted to common samples.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param datasets Character vector, the names of the datasets from \code{mo_data} to include in the analysis.
#' @return A list, in which each element corresponds to one omics dataset, with samples as rows and features as columns.
#' @export
get_input_mixomics_unsupervised <- function(mo_data, datasets = names(mo_data)) {
  mo_data <- check_input_multidataset(mo_data, datasets)

  common_samples <- MultiDataSet::commonIds(mo_data)

  ## Get the datasets
  ds <- get_datasets(mo_data)

  res <- lapply(datasets, function(i) {
    mat <- ds[[i]][, common_samples, drop = FALSE]
    mat <- t(mat)
    attr(mat, "datatype") <- class(mo_data[[i]])

    mat
  })
  names(res) <- datasets ## should not be necessary

  res
}

#' Generate mixomics input data for supervised methods
#'
#' Creates an object that can be used as input for the MixOmics package.
#' It contains the omics datasets restricted to common samples (with no missing group information) and
#' an outcome group for each sample.
#'
#' Samples with missing values in the `group` column of the sample metadata will be removed from the dataset.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param group Character, the column name in the samples metadata data-frames to use as samples group
#' (use \code{\link{get_samples_metadata}} to view the samples information data-frame for each dataset).
#' This column should be either of type factor, character or integer.
#' @param datasets Character vector, the names of the datasets from \code{mo_data} to include in the analysis.
#' @return A list, in which each element corresponds to one omics dataset, with samples as rows and features as columns.
#' The `Y` element is a factor vector with the outcome group of each sample.
#' @export
get_input_mixomics_supervised <- function(mo_data, group, datasets = names(mo_data)) {
  ## for devtools::check
  id <- NULL

  mo_data <- check_input_multidataset(mo_data, datasets)
  mo_data <- MultiDataSet::commonSamples(mo_data)

  ## Get the samples group
  .check_input_var_smetadata_common(group, mo_data)
  grp <- get_samples_metadata_combined(mo_data, only_common_cols = FALSE) |>
    dplyr::select(id, !!sym(group)) |>
    tibble::deframe()

  ## The group should be factor, character or integer, but not numeric
  if (!(typeof(grp) %in% c("factor", "character", "integer"))) stop(group, " column should be a factor, character or integer column.")

  ## remove samples with missing group information
  if (any(is.na(grp))) {
    warning("Some missing values in samples group. Sample with missing group value will be removed from the dataset.")
    grp <- grp[!is.na(grp)]
  }

  ## Get the datasets
  ds <- get_datasets(mo_data)

  res <- lapply(datasets, function(i) {
    mat <- ds[[i]][, names(grp), drop = FALSE]
    mat <- t(mat)
    attr(mat, "datatype") <- class(mo_data[[i]])

    mat
  })
  names(res) <- datasets ## should not be necessary

  res$Y <- factor(grp)

  res
}


.mixomics_get_perf_plsda <- function(plsda_perf) {
  res <- lapply(names(plsda_perf$error.rate), function(i) {
    df1 <- tibble::as_tibble(plsda_perf$error.rate[[i]], rownames = "comp") |>
      dplyr::mutate(measure = i) |>
      tidyr::pivot_longer(
        cols = tidyselect::ends_with(".dist"),
        names_to = "distance",
        values_to = "classification_error_rate"
      )

    if (!is.null(plsda_perf$error.rate.sd)) {
      df2 <- tibble::as_tibble(plsda_perf$error.rate.sd[[i]], rownames = "comp") |>
        dplyr::mutate(measure = i) |>
        tidyr::pivot_longer(
          cols = tidyselect::ends_with(".dist"),
          names_to = "distance",
          values_to = "classification_error_rate_sd"
        )

      return(dplyr::full_join(df1, df2, by = c("comp", "measure", "distance")))
    } else {
      return(df1)
    }
  }) |>
    purrr::reduce(dplyr::bind_rows)

  return(res)
}

.mixomics_cpus_to_bparam <- function(cpus) {
  if (is.null(cpus)) {
    return(BiocParallel::SerialParam())
  }

  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    stop(
      "Package \"BiocParallel\" must be installed to use the 'cpus' argument.",
      call. = FALSE
    )
  }

  if (.Platform$OS.type == "windows") {
    BPPARAM <- BiocParallel::SnowParam(workers = cpus)
  } else {
    BPPARAM <- BiocParallel::MulticoreParam(workers = cpus)
  }

  return(BPPARAM)
}
