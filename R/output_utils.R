#' Get latent dimensions levels from dimension reduction output
#'
#' Extracts the latent dimension levels from the output of a dimension reduction method.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @returns A character vector giving the labels of the latent dimensions.
#' @export
get_latent_dimensions <- function(method_output) {
  if (!inherits(method_output, "output_dimension_reduction")) stop("Expecting an object of class 'output_dimension_reduction' (from get_output()).")

  sc_levels <- levels(method_output$samples_score$latent_dimension)
  fw_levels <- levels(method_output$features_weight$latent_dimension)

  if (!identical(sc_levels, fw_levels)) {
    stop("Latent dimension levels are not the same in the samples score and features weight tables.")
  }

  return(sc_levels)
}

#' Extract selected features
#'
#' Extracts selected features from the output of an integration method. Only
#' features with a non-null weight for at least one latent dimension will be returned.
#' If a `MultiDataSet` object is supplied, information about the features from
#' the features metadata will be added to the resulting table.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @param latent_dimensions Character vector of latent dimensions name. Default
#' value is `NULL` (top contributing features will be returned for all latent
#' dimensions).
#' @param datasets Character vector of datasets name. Default value is `NULL`
#' (top contributing features will be returned for all datasets).
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns A tibble containing one row per feature and latent dimension, giving
#' the weight and importance score of the feature for the corresponding latent dimension.
#' If `mo_data` is supplied, information about the features from
#' the features metadata will be added to the resulting table.
#' @export
get_selected_features <- function(method_output,
                                  latent_dimensions = NULL,
                                  datasets = NULL,
                                  mo_data = NULL) {

  ## for devtools::check
  weight <- latent_dimension <- dataset <- importance <- NULL

  if (!inherits(method_output, "output_dimension_reduction")) {
    stop("Expecting an object of class 'output_dimension_reduction' (from get_output()).")
  }

  method_output <- method_output |>
    .filter_output_dimensions(latent_dimensions) |>
    .filter_output_datasets(datasets)

  res <- method_output$features_weight |>
    dplyr::filter(weight != 0) |>
    join_features_metadata(mo_data) |>
    dplyr::arrange(latent_dimension, dataset, dplyr::desc(importance))

  return(res)
}

#' Extract top features
#'
#' Extracts the features with highest contribution to the latent dimensions
#' constructed by an integration method. Can retain a specific number of top
#' contributing features for each dataset and latent dimension, or all features
#' above a minimum importance score.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @param n_features Integer, the number of features to extract for each latent
#' dimension and dataset. Ignored if `min_importance` is set. Default value is `10`.
#' Will include all ties.
#' @param min_importance Numeric value between 0 and 1, minimum importance score
#' used to select features. Default value is `NULL`, i.e. the top `n_features`
#' features are selected instead.
#' @param latent_dimensions Character vector of latent dimensions name. Default
#' value is `NULL` (top contributing features will be returned for all latent
#' dimensions).
#' @param datasets Character vector of datasets name. Default value is `NULL`
#' (top contributing features will be returned for all datasets).
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns A tibble containing one row per feature and latent dimension, giving
#' the weight and importance score of the feature for the corresponding latent dimension.
#' If `mo_data` is supplied, information about the features from
#' the features metadata will be added to the resulting table.
#' @export
get_top_features <- function(method_output,
                             n_features = 10,
                             min_importance = NULL,
                             latent_dimensions = NULL,
                             datasets = NULL,
                             mo_data = NULL) {

  ## for devtools::check
  latent_dimension <- dataset <- importance <- NULL

  if (!inherits(method_output, "output_dimension_reduction")) {
    stop("Expecting an object of class 'output_dimension_reduction' (from get_output()).")
  }

  method_output <- method_output |>
    .filter_output_dimensions(latent_dimensions) |>
    .filter_output_datasets(datasets)

  res <- method_output$features_weight |>
    dplyr::group_by(latent_dimension, dataset)

  if (is.null(min_importance)) {
    res <- res |>
      dplyr::slice_max(order_by = importance, n = n_features, with_ties = TRUE)
  } else {
    res <- res |>
      dplyr::filter(importance >= min_importance)
  }

  res <- join_features_metadata(res, mo_data) |>
    dplyr::arrange(latent_dimension, dataset, dplyr::desc(importance))

  return(res)
}

#' Check names of output list
#'
#' Checks that the names of a list of outputs from several integration methods
#' are unique. If not named, the name of the method will be used as name.
#'
#' @param output_list List of integration methods output, each generated via one of the
#' `get_output_*()` function.
#' @returns output_list (but named if it wasn't).
.check_names_output_list <- function(output_list) {
  check_class <- output_list |>
    purrr::map_lgl(~ inherits(.x, "output_dimension_reduction"))

  if (!all(check_class)) stop("'output_list': expecting a list of 'output_dimension_reduction' objects (from get_output()).")


  if (is.null(names(output_list))) {
    m <- purrr::map_chr(output_list, ~ attr(.x, "method"))

    if (any(duplicated(m))) {
      stop("'output_list' contains several results from a same integration method. Add unique names to the list to differentiate them (see function help for an example).",
           call. = FALSE
      )
    }

    names(output_list) <- m
  } else {
    if (any(duplicated(names(output_list)))) {
      stop("'output_list' has duplicated names. Please use unique names.", call. = FALSE)
    }
  }

  return(output_list)
}


#' Filter latent dimensions
#'
#' Filters latent dimensions by name in the output of an integration method.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @param latent_dimensions Character vector giving the latent dimensions to retain
#' in the method's output.
#' @param fixed_length Integer, expected length of `latent_dimensions`.
#' If `NULL` (default value), the length of `latent_dimensions` will not be checked.
#' @param method_name Character, name of the method to use in the error message.
#' @returns Similar to `method_output`, but the samples score table
#' or features weight table have been filtered.
.filter_output_dimensions <- function(method_output, latent_dimensions, fixed_length = NULL, method_name = attr(method_output, "method")) {
  ## for devtools::check
  latent_dimension <- NULL

  if (!inherits(method_output, "output_dimension_reduction")) stop("Expecting an object of class 'output_dimension_reduction' (from get_output()).")

  if (is.null(latent_dimensions)) {
    return(method_output)
  }

  error_suppl <- paste0(" for method ", method_name)


  if (!is.null(fixed_length)) {
    if (length(latent_dimensions) != fixed_length) {
      stop("'latent_dimensions' should be of length ", fixed_length, ".")
    }
  }

  ld_levels <- get_latent_dimensions(method_output)
  .check_names(
    latent_dimensions,
    ld_levels,
    paste0(
      "'_W_' are not valid latent dimension names",
      error_suppl,
      ". Possible values are '_C_'."
    )
  )

  res <- method_output |>
    purrr::map(
      ~ .x |>
        dplyr::filter(latent_dimension %in% latent_dimensions) |>
        dplyr::mutate(latent_dimension = droplevels(latent_dimension))
    )
  attr(res, "class") <- attr(method_output, "class")
  attr(res, "method") <- attr(method_output, "method")

  return(res)
}

#' Filter latent dimensions in list
#'
#' Filters latent dimensions by name in a list of outputs from integration methods.
#'
#' @param output_list List of integration method outputs each generated via one of the
#' `get_output()` function.
#' @param latent_dimensions Named list, where each element is a character vector
#' giving the latent dimensions to retain in the corresponding element of `output_list`.
#' Names must match those of `output_list`.
#' @param all_present Logical, whether there should be one element in `latent_dimensions`
#' for each element of `output_list`. If `TRUE`, an error will be
#' returned if the length and names of `output_list` and `latent_dimensions` do
#' not match. Default value is `FALSE`.
#' @param fixed_length Integer, expected length of each element of `latent_dimensions`.
#' If `NULL` (default value), the length of elements in `latent_dimensions` can vary.
#' @returns A list of output similar to `output_list`, but the samples score table
#' or features weight table have been filtered.
.filter_output_dimensions_list <- function(output_list, latent_dimensions, all_present = FALSE, fixed_length = NULL) {
  if (is.null(latent_dimensions)) {
    return(output_list)
  }

  ## Making sure that the names of the latent_dimensions list match the
  ## names of output_list
  if (is.null(names(latent_dimensions))) {
    stop(
      "'latent_dimensions' argument should be a named list, names should be: '",
      paste0(names(output_list), collapse = "', '"),
      "'."
    )
  }

  .check_names(
    names(latent_dimensions),
    names(output_list),
    "'latent_dimensions' argument: list names '_W_' do not match names from 'output_list'. Possible values are '_C_'."
  )

  if (all_present) {
    .check_names(
      names(output_list),
      names(latent_dimensions),
      "'latent_dimensions' argument should match 'output_list' elements; missing value for '_W_'."
    )
  }

  if (!is.null(fixed_length)) {
    if (any(purrr::map_int(latent_dimensions, length) != fixed_length)) {
      stop("Elements in 'latent_dimensions' list should be of length ", fixed_length, ".")
    }
  }

  for (i in names(latent_dimensions)) {
    output_list[[i]] <- .filter_output_dimensions(
      output_list[[i]],
      latent_dimensions[[i]],
      fixed_length = fixed_length,
      method_name = i
    )
  }

  return(output_list)
}

#' Filter datasets
#'
#' Filters datasets by name in the output of an integration method.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @param datasets Character vector giving the datasets to retain in the features
#' weight table of the method's output.
#' @param fixed_length Integer, expected length of `datasets`.
#' If `NULL` (default value), the length of `datasets` will not be checked.
#' @param method_name Character, name of the method to use in the error message.
#' @returns Similar to `method_output`, but the features weight table has been filtered.
.filter_output_datasets <- function(method_output, datasets, fixed_length = NULL, method_name = attr(method_output, "method")) {
  ## for devtools::check
  dataset <- NULL

  if (!inherits(method_output, "output_dimension_reduction")) stop("Expecting an object of class 'output_dimension_reduction' (from get_output()).")

  if (is.null(datasets)) {
    return(method_output)
  }

  if (is.null(method_name)) {
    error_suppl <- ""
  } else {
    error_suppl <- paste0(" for method ", method_name)
  }


  if (!is.null(fixed_length)) {
    if (length(datasets) != fixed_length) {
      stop("'datasets' should be of length ", fixed_length, ".")
    }
  }

  ds_levels <- levels(method_output$features_weight$dataset)
  .check_names(
    datasets,
    ds_levels,
    paste0(
      "'_W_' are not valid dataset names",
      error_suppl,
      ". Possible values are '_C_'."
    )
  )

  res <- method_output
  res$features_weight <- res$features_weight |>
    dplyr::filter(dataset %in% datasets) |>
    dplyr::mutate(dataset = droplevels(dataset))

  return(res)
}
