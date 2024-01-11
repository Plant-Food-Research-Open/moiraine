#' Check null or equality
#'
#' Tests whether an object is `NULL` or equal to some value.
#'
#' @param x The object to test.
#' @param val Value to compare to
#' @returns `TRUE` or `FALSE`.
#' @export
is_equal_or_null <- function(x, val) {

  if (is.null(x)) {
    return(TRUE)
  }

  return(identical(x, val))
}

#' Hierarchical clustering of matrix rows
#'
#' Performs a hierarchical clustering on the rows of a matrix. Code inspired
#' by ComplexHeatmap package.
#'
#' @param x Matrix.
#' @returns A dendrogram.
#' @export
hclust_matrix_rows <- function(x) {
  ## Initial hierarchical clustering
  dmat <- stats::dist(x, method = "euclidean")
  hcres <- stats::hclust(dmat, method = "complete")

  ## Reordering dendrogram
  order_list <- hcres$order
  rmat <- -rowMeans(x, na.rm = TRUE)
  dendro <- stats::reorder(
    stats::as.dendrogram(hcres),
    rmat[sort(order_list)],
    mean
  )

  return(dendro)
}


#' Returns options list as a tibble
#'
#' Transforms a list of options (or parameters) into a tibble with the name of
#' the options (parameters) in one column, and their value in a second column.
#' Vector values are collapsed to span only one column.
#'
#' @param options_list A named list, where each element corresponds to one
#'   option or parameter and the name of the element corresponds to the name of
#'   the option/parameter.
#' @return a tibble, with the `Parameter` column giving the list of the options
#'   or parameters, and the `Value` column giving the values of the
#'   corresponding option or parameter.
#' @export
options_list_as_tibble <- function(options_list) {
  list <- sapply(options_list, function(x) {
    if (length(x) == 1) {
      return(x)
    }
    if (is.null(names(x))) {
      return(paste0("'", paste0(x, collapse = "', '"), "'"))
    }
    return(paste0("'", names(x), "' = '", x, "'", collapse = ", "))
  })

  tibble::tibble(
    Parameter = names(list),
    Value = unname(list)
  )
}

#' Check a character vector against reference vector
#'
#' Check values in a character vector against values in a reference vector.
#' Generates a useful error message listing the discrepancies if needed.
#'
#' @param names_to_check Character vector to be checked.
#' @param correct_names Character vector used as reference.
#' @param message Character, the error message to display in case of
#'   discrepancy.
#' @param wrong_names_code Character, the code used in `message` to be replaced
#'   with the values from `names_to_check` that are not in `correct_names`.
#' @param correct_names_code Character, the code used in `message` to be
#'   replaced with the values from `correct_names`.
#'
#' @noRd
.check_names <- function(names_to_check,
                         correct_names,
                         message,
                         wrong_names_code = "_W_",
                         correct_names_code = "_C_") {
  wrong_names <- setdiff(names_to_check, correct_names)

  if (length(wrong_names)) {
    msg <- stringr::str_replace(
      message,
      wrong_names_code,
      paste0(wrong_names, collapse = "', '")
    )
    msg <- stringr::str_replace(
      msg,
      correct_names_code,
      paste0(correct_names, collapse = "', '")
    )
    stop(msg, call. = FALSE)
  }

  invisible(NULL)
}


#' Check for missing values in vector
#'
#' Checks for the presence of `NULL`s or `NA`s in a vector.
#'
#' @param x Vector.
#' @returns Logical, whether there are any `NULL` or `NA` values in `x`.
#'
#' @noRd
.check_missing <- function(x) {
  any(is.null(x), is.na(x))
}

#' Turn input argument into list
#'
#' Takes an input argument and turns it into a list with one element per dataset
#'
#' @param x The input parameter. Can either be a vector (will be used for all
#'   datasets) or a list with one value per dataset.
#' @param datasets Character vector, names of datasets.
#' @param fixed_length Integer, expected length of the elements in `x`. If
#'   `NULL` (default value), their length will not be checked.
#' @param allow_subset Logical, is `x` allowed to contain elements for a subset
#'   of the datasets? If `FALSE` and `x` is a list, its length must match that
#'   of `datasets`. Default value is `FALSE`.
#' @returns A named list, where each element gives the column name to use in the
#'   samples metadata of the corresponding dataset, and the name of the element
#'   is the name of the dataset. If `x` is `NULL`, each element will be `NULL`
#'   as well.
#'
#' @noRd
.make_var_list <- function(x,
                           datasets,
                           fixed_length = NULL,
                           allow_subset = FALSE) {
  error_prefix <- paste0("'", deparse(substitute(x)), "' argument")

  ## If the argument is NULL, return a list of NULL
  if (is.null(x)) {
    res <- datasets |>
      rlang::set_names() |>
      purrr::map(~NULL)

    return(res)
  }

  ## If only one value given, use for all datasets
  ## Otherwise must match number of datasets
  if (!is.list(x)) {
    x <- purrr::map(seq_along(datasets), ~x) |>
      rlang::set_names(datasets)
  } else if ((length(x) != length(datasets)) && !allow_subset) {
    stop(
      error_prefix,
      " should either be a vector or a named list of length ",
      length(datasets), ".",
      call. = FALSE
    )
  }

  ## Names must match datasets
  if (is.null(names(x))) {
    stop(
      error_prefix,
      " should be named; names should be: '",
      paste0(datasets, collapse = "', '"), "'.",
      call. = FALSE
    )
  }
  .check_names(
    names(x),
    datasets,
    paste0(
      error_prefix, ": '_W_' are not existing datasets or do not match ",
      "the 'datasets' argument. Possible values are: '_C_'.")
  )

  if (!is.null(fixed_length)) {
    x_len <- purrr::map_int(x, length)
    if (!all(x_len == fixed_length)) {
      stop(error_prefix, " values should have length ", fixed_length, ".")
    }
  }

  return(as.list(x))
}

#' Check that variable names corresponds to columns in features metadata
#'
#' Checks whether a variable name corresponds to a column in the features
#' metadata of the corresponding dataset.
#'
#' @param x Named character list, with one element per dataset, each element
#'   giving the name of the column from the features metadata of the
#'   corresponding dataset. The names should correspond to dataset names in
#'   `mo_data`. Should be checked with `.make_var_list()`.
#' @param mo_data A `MultiDataSet` object containing features information for
#'   the datasets. Should be checked with `.check_input_multidataset()`.
#' @returns Nothing. Will throw an error if need be.
#'
#' @noRd
.check_input_var_fmetadata <- function(x, mo_data) {
  error_prefix <- paste0("'", deparse(substitute(x)), "' argument")

  ## Checking whether values are columns in samples metadata
  smetadata_list <- get_features_metadata(mo_data) |>
    purrr::map(colnames)

  purrr::iwalk(
    x,
    function(.x, .y) {
      if (is.null(.x)) {
        return(invisible(NULL))
      }

      .check_names(
        .x,
        smetadata_list[[.y]],
        paste0(
          error_prefix,
          ": '_W_' is not a column in the features metadata for the ",
          .y,
          " dataset. Possible values are: '_C_'."
        )
      )
    }
  )

  return(invisible(NULL))
}


#' Check that variable names corresponds to columns in samples metadata
#'
#' Checks whether a variable name corresponds to a column in the samples
#' metadata of the corresponding dataset. If one value is provided, will be used
#' for all datasets.
#'
#' @param x Named character list, with one element per dataset, each element
#'   giving the name of the column from the samples metadata of the
#'   corresponding dataset. The names should correspond to dataset names in
#'   `mo_data`. Should be checked with `.make_var_list()`.
#' @param mo_data A `MultiDataSet` object containing samples information for the
#'   datasets. Should be checked with `.check_input_multidataset()`.
#' @returns Nothing. Will throw an error if need be.
.check_input_var_smetadata <- function(x, mo_data) {
  error_prefix <- paste0("'", deparse(substitute(x)), "' argument")

  purrr::iwalk(
    x,
    function(.x, .y) {
      if (is.null(.x)) {
        return(invisible(NULL))
      }

      smetadata_df <- get_samples_metadata(mo_data)[[.y]]

      .check_names(
        .x,
        colnames(smetadata_df),
        paste0(
          error_prefix,
          ": '_W_' is not a column in the samples metadata for the ",
          .y,
          " dataset. Possible values are: '_C_'."
        )
      )
    }
  )

  return(invisible(NULL))
}

#' Check that variable names corresponds to columns in samples metadata
#'
#' Checks whether a variable name corresponds to a column in the samples
#' metadata of the corresponding dataset. If one value is provided, will be used
#' for all datasets.
#'
#' @param x Character, name of the column from the samples metadata.
#' @param mo_data A `MultiDataSet` object containing samples information for the
#'   datasets. Should be checked with `.check_input_multidataset()`.
#' @returns Nothing. Will throw an error if need be.
.check_input_var_smetadata_common <- function(x, mo_data) {
  if (is.null(x)) {
    return(invisible(NULL))
  }

  error_prefix <- paste0("'", deparse(substitute(x)), "' argument")
  smetadata_common <- get_samples_metadata_combined(
    mo_data,
    only_common_cols = FALSE
  )

  .check_names(
    x,
    colnames(smetadata_common),
    paste0(
      error_prefix,
      ": '_W_' is not a column in the samples metadata of any dataset. ",
      "Possible values are: '_C_'."
    )
  )

  return(invisible(NULL))
}

#' Adds features label to data-frame
#'
#' Adds the features label to a data-frame for plotting. Can be extracted from
#' the features metadata of a `MultiDataSet` object; otherwise use the feature
#' IDs as label. If some labels are missing, feature IDs will be used instead.
#'
#' @param toplot The data-frame to which the labels should be added.
#' @param label_cols Character or named list of character, giving for each
#'   dataset the name of the column in the corresponding features metadata to
#'   use as label. If one value, will be used for all datasets. If list, the
#'   names must correspond to the names of the datasets in `mo_data`. If a
#'   dataset is missing from the list or no value is provided, feature IDs will
#'   be used as labels. Alternatively, use `feature_id` to get the feature IDs
#'   as labels.
#' @param mo_data A `MultiDataSet` object. Only used if `label_cols` is not
#'   `NULL`.
#' @param truncate Integer, width to which the labels should be truncated (to
#'   avoid issues with very long labels in plots). If `NULL` (default value), no
#'   truncation will be performed.
#' @returns the `toplot` data-frame with an additional column `label`.
#'
#' @noRd
.add_features_labels_toplot <- function(toplot,
                                        label_cols,
                                        mo_data,
                                        truncate = NULL) {
  ## for devtools::check
  dataset <- feature_id <- NULL

  if (!is.null(label_cols)) {
    if (is.null(mo_data)) {
      stop("Need to provide a MultiDataSet object through 'mo_data' argument ",
           "in order to use 'label_cols' argument.")
    }

    datasets <- levels(toplot$dataset)
    if (is.null(datasets)) datasets <- unique(toplot$dataset)
    mo_data <- .check_input_multidataset(mo_data, datasets)

    res <- get_features_labels(mo_data, label_cols, truncate)

    if (!is.null(levels(toplot$dataset))) {
      res <- res |>
        dplyr::mutate(
          dataset = factor(dataset, levels = levels(toplot$dataset))
        )
    }

    toplot <- toplot |>
      dplyr::left_join(
        res,
        by = c("dataset", "feature_id")
      )
  } else {
    toplot <- toplot |>
      dplyr::mutate(
        label = feature_id
      )
  }

  return(toplot)
}

#' Make Quarto report template from Rmd template
#'
#' Generates the Quarto report templates from corresponding Rmd report
#' templates.
#'
#' @noRd
.make_quarto_template <- function() {
  ## Constructing template for report header and first part
  main_rmd <- here::here("inst/templates/report_template.Rmd")

  main_rmd_con <- file(main_rmd, "r")
  main_rmd_lines <- readLines(main_rmd_con)
  close(main_rmd_con)

  yaml_header_end <- which(main_rmd_lines == "---")[2]
  yaml_header <- main_rmd_lines[seq_len(yaml_header_end)]
  content <- main_rmd_lines[-seq_len(yaml_header_end)]

  yaml_header[yaml_header == "date: '`r format(Sys.Date(), \"%B %d, %Y\")`'"] <- "date: today"
  yaml_header[yaml_header == "output:"] <- "format:"
  yaml_header[yaml_header == "  html_document:"] <- "  html:"
  yaml_header <- yaml_header[!(yaml_header %in% c("     toc_float: true", "     theme: flatly"))]

  tmp <- which(yaml_header == "  html:")
  yaml_header <- c(
    yaml_header[seq_len(tmp)],
    "     toc-location: right",
    "     embed-resources: true",
    yaml_header[-seq_len(tmp)]
  )

  content <- content |>
    stringr::str_remove_all(" \\{\\.tabset \\.tabset-pills\\}") |>
    stringr::str_remove_all(" \\{\\.active\\}")

  main_qmd <- here::here("inst/templates/quarto_report_template.qmd")
  main_qmd_con <- file(main_qmd, "w")
  writeLines(c(yaml_header, content), con = main_qmd_con)
  close(main_qmd_con)

  dir("inst/templates/", pattern = "fragment.+Rmd") |>
    purrr::walk(
      function(.x) {
        file_con <- file(here::here("inst/templates", .x), "r")
        content <- readLines(file_con)
        close(file_con)

        content <- content |>
          stringr::str_remove_all(" \\{\\.tabset \\.tabset-pills\\}") |>
          stringr::str_remove_all(" \\{\\.active\\}")

        new_file <- paste0(
          "quarto_",
          stringr::str_replace(.x, "\\.Rmd$", ".qmd")
        )
        new_file_con <- file(here::here("inst/templates", new_file), "w")
        writeLines(content, con = new_file_con)
        close(new_file_con)
      }
    )
}
