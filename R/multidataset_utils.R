#' Checks whether object is MultiDataSet
#'
#' Checks whether an input object is a [MultiDataSet::MultiDataSet-class]
#' object.
#'
#' @param x An object.
#' @returns nothing.
#' @export
check_is_multidataset <- function(x) {
  if (!inherits(x, "MultiDataSet")) {
    stop("'mo_data' argument: Expecting MultiDataSet object", call. = FALSE)
  }

  return(invisible(NULL))
}

#' Check a MultiDataSet input
#'
#' Checks a MultiDataSet object provided as an input. In particular, checks that
#' 1) the input object is a `MultiDataSet` object, 2) the datasets stored match
#' the datasets named provided (if any). Will restrict the MultiDataSet object
#' to only necessary datasets.
#'
#' @param x A input object that hopefully is  a `MultiDataSet` object.
#' @param datasets Character vector of dataset names that should be in `x`. If
#'   `NULL` (default value), dataset names will not be checked
#' @returns The MultiDataSet object restricted to the datasets required (if
#'   `datasets` is not `NULL`.
#' @export
check_input_multidataset <- function(x, datasets = NULL) {

  check_is_multidataset(x)

  if (!is.null(datasets)) {
    .check_names(
      datasets,
      names(x),
      "'_W_' datasets are not present in mo_data. Possible dataset names are: '_C_'."
    )

    x <- x[, datasets, drop = FALSE]
  }

  return(x)
}

#' Get multi-omics measurement datasets
#'
#' Returns the multi-omics datasets as a list of matrices from a MultiDataSet
#' object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns A named list of matrices, each with features as rows and samples as
#'   columns.
#' @export
get_datasets <- function(mo_data) {
  check_is_multidataset(mo_data)

  res <- MultiDataSet::as.list(mo_data)

  return(res)
}

#' Get multi-omics dataset as matrix
#'
#' Extracts an omics dataset as a matrix of measurements from a MultiDataSet
#' object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param dataset_name Character, name of the omics dataset to extract.
#' @param keep_dataset_name Logical, should the dataset name be stored in the
#'   `'dataset_name'` attribute of the resulting matrix? Default value is
#'   `FALSE`.
#' @returns A matrix of measurements with features as rows and samples as
#'   columns. The name of the dataset is stored in the `'dataset_name'`
#'   attribute if `keep_dataset_name` is `TRUE`.
#' @examples
#' \dontrun{
#' ## mo_data is a MultiDataSet object with a dataset called "rnaseq"
#' mat <- get_dataset_matrix(mo_data, "rnaseq", keep_dataset_name = TRUE)
#' ## with keep_dataset_name = TRUE, can recover dataset name as follows:
#' attr(mat, "dataset_name")
#' }
#' @export
get_dataset_matrix <- function(mo_data,
                               dataset_name,
                               keep_dataset_name = FALSE) {
  if (length(dataset_name) > 1) {
    stop("'dataset_name' argument should be of length 1.")
  }

  mo_data <- check_input_multidataset(mo_data, dataset_name)

  mat <- get_datasets(mo_data)[[dataset_name]]

  if (keep_dataset_name) attr(mat, "dataset_name") <- dataset_name

  return(mat)
}

#' Get samples metadata dataframes from MultiDataSet
#'
#' Extracts the samples metadata data-frame (phenoData field) from each dataset
#' for a MultiDataSet object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns A named list of data-frames, one per dataset in the `mo_data`
#'   object.
#' @export
get_samples_metadata <- function(mo_data) {
  check_is_multidataset(mo_data)
  Biobase::pData(mo_data)
}


#' Get combined samples metadata data-frame from MultiDataSet
#'
#' Extracts the samples metadata data-frame (phenoData field) from each dataset
#' for a MultiDataSet object and combine them into one dataframe.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param only_common_cols Logical, whether to retain only common columns. If
#'   `TRUE` (default value), only retain the columns that are present in the
#'   samples metadata of all datasets. If `FALSE`, retain all columns from each
#'   datasets' sample metadata.
#' @returns A data-frame of samples metadata.
#' @export
get_samples_metadata_combined <- function(mo_data, only_common_cols = TRUE) {
  ## Testing of mo_data class done in get_samples_metadata
  smeta <- get_samples_metadata(mo_data)

  common_cols <- smeta |>
    purrr::map(colnames) |>
    purrr::reduce(intersect)

  if (only_common_cols) {
    res <- smeta |>
      purrr::map(
        \(.x) dplyr::select(.x, tidyselect::all_of(common_cols))
      ) |>
      purrr::reduce(dplyr::full_join, by = common_cols)
  } else {
    res <- smeta |>
      purrr::reduce(
        function(.x, .y) {
          cc <- intersect(colnames(.x), colnames(.y))
          dplyr::full_join(.x, .y, by = cc)
        }
      )
  }

  res <- dplyr::distinct(res)

  duplicated_samples <- duplicated(res$id)
  if (any(duplicated_samples)) {
    stop(
      "Conflicting information in samples metadata for samples ",
      paste0(res$id[duplicated_samples], collapse = ", "),
      "."
    )
  }

  rownames(res) <- res$id

  return(res)
}

#' Get features metadata dataframes from MultiDataSet
#'
#' Extracts the features metadata dataframe (featureData field) from each
#' dataset for a MultiDataSet object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns A named list of data-frames, one per dataset in the `mo_data`
#'   object.
#' @export
get_features_metadata <- function(mo_data) {
  check_is_multidataset(mo_data)
  Biobase::fData(mo_data)
}

#' Get sample IDs from MultiDataSet
#'
#' Extract the list of sample IDs from each dataset in a MultiDataSet object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns A named list, with one element per dataset, and each element is a
#'   character vector of sample IDs.
#' @export
get_samples <- function(mo_data) {
  check_is_multidataset(mo_data)
  res <- MultiDataSet::as.list(mo_data)
  purrr::map(res, colnames)
}

#' Get feature IDs from MultiDataSet
#'
#' Extract the list of feature IDs from each dataset in a MultiDataSet object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns A named list, with one element per dataset, and each element is a
#'   character vector of feature IDs.
#' @export
get_features <- function(mo_data) {
  check_is_multidataset(mo_data)
  res <- MultiDataSet::as.list(mo_data)
  purrr::map(res, rownames)
}

#' Join feature metadata to table
#'
#' Adds features metadata information to a table containing feature IDs.
#'
#' @param df Data-frame of tibble with a column `feature_id` containing feature
#'   IDs.
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns The `df` table with additional columns containing information about
#'   the features from the features metadata table.
#' @export
join_features_metadata <- function(df, mo_data) {
  ## for devtools check
  feature_id <- NULL

  if (is.null(mo_data)) {
    return(df)
  }
  check_is_multidataset(mo_data)

  if (!("feature_id" %in% colnames(df))) {
    stop("`df` must contain a 'feature_id' column with feature IDs.")
  }

  fmeta <- get_features_metadata(mo_data) |>
    purrr::reduce(dplyr::bind_rows) |>
    ## Removing empty columns (for example when there are no features
    ## from a given dataset) - need to do this here so we don't remove
    ## columns from `df` that are empty
    dplyr::filter(feature_id %in% df$feature_id) |>
    dplyr::select(
      feature_id,
      where(~ any(!is.na(.x)))
    )

  if (nrow(fmeta) == 0) {
    stop("No features matching IDs present in `df$feature_id`.")
  }

  res <- df |>
    dplyr::left_join(fmeta, by = "feature_id")

  return(res)
}


#' Join samples metadata to table
#'
#' Adds samples metadata information to a table containing sample IDs.
#'
#' @param df Data-frame or tibble with a column `id` containing sample IDs.
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param datasets Character vector, name(s) of datasets from which the samples
#'   metadata should be extracted. If `NULL` (default value), information from
#'   all datasets will be used.
#' @returns The `df` table with additional columns containing information about
#'   the samples from the samples metadata table.
#' @export
join_samples_metadata <- function(df, mo_data, datasets = NULL) {
  ## For devtools check
  id <- NULL

  if (is.null(mo_data)) {
    return(df)
  }

  if (!("id" %in% colnames(df))) {
    stop("`df` must contain an 'id' column with sample IDs.")
  }

  mo_data <- check_input_multidataset(mo_data, datasets)

  smeta <- get_samples_metadata_combined(mo_data, only_common_cols = FALSE) |>
    ## Removing empty columns (for example when there are no samples
    ## from a given dataset) - need to do this here so we don't remove
    ## columns from `df` that are empty
    dplyr::filter(id %in% df$id) |>
    dplyr::select(
      id,
      where(~ any(!is.na(.x)))
    )

  if (nrow(smeta) == 0) {
    stop("No samples matching IDs present in `df$id`.")
  }

  res <- df |>
    dplyr::left_join(smeta, by = "id")

  return(res)
}

#' Adding data-frame to features metadata
#'
#' Adds information from a data-frame to the features metadata of a MultiDataSet
#' object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param df A tibble or data-frame of features information, with at least
#'   columns `feature_id` (giving the feature IDs) and `dataset` (giving the
#'   name of the dataset to which the features belong), and one row per feature
#'   ID. Can contain info about features from different datasets.
#' @returns A MultiDataSet object, with info from `df` adding to its
#'   corresponding features metadata.
#' @export
add_features_metadata <- function(mo_data, df) {
  ## for devtools::check
  dataset <- NULL

  check_input_multidataset(mo_data)

  if (!inherits(df, "data.frame")) {
    stop("'df' argument should be a tibble or data-frame.")
  }

  if (!all(c("feature_id", "dataset") %in% colnames(df))) {
    stop("'df' argument should have columns 'feature_id' and 'dataset'.")
  }

  datasets <- unique(df$dataset)
  .check_names(
    datasets,
    names(mo_data),
    "'df' argument: in 'dataset' column, '_W_' are not datasets in 'mo_data'. Possible dataset names are: '_C_'."
  )

  fmeta_list <- get_features_metadata(mo_data)
  features_id <- get_features(mo_data)

  for (i in datasets) {
    error_msg <- paste0(i, " dataset: ")
    df_add <- df |>
      dplyr::filter(dataset == i) |>
      dplyr::select(-dataset)

    fids <- unique(df_add$feature_id)

    wrong_ids <- setdiff(fids, features_id[[i]]) |> length()
    missing_ids <- setdiff(features_id[[i]], fids) |> length()

    if (wrong_ids > 0) {
      warning(
        error_msg,
        wrong_ids,
        " feature IDs not in 'mo_data', will be removed from features metadata.",
        call. = FALSE
      )
    }
    if (missing_ids > 0) {
      warning(
        error_msg,
        missing_ids,
        " feature IDs missing from 'df' data-frame.",
        call. = FALSE
      )
    }

    check_unique <- max(table(df_add$feature_id))
    if (check_unique > 1) {
      stop(
        error_msg,
        "should have only one row per feature ID in 'df' data-frame."
      )
    }

    df_ex <- fmeta_list[[i]]

    existing_cols <- intersect(colnames(df_ex), colnames(df_add)) |>
      setdiff("feature_id")
    if (length(existing_cols) > 0) {
      stop(
        error_msg,
        "'",
        paste0(existing_cols, collapse = "', '"),
        "' columns already in features metadata table."
      )
    }

    df_new <- df_ex |>
      dplyr::left_join(df_add, by = "feature_id") |>
      Biobase::AnnotatedDataFrame()
    rownames(df_new) <- df_new$feature_id

    mo_data@featureData[[i]]$main <- df_new
  }

  return(mo_data)
}

#' Adding data-frame to samples metadata
#'
#' Adds information from a data-frame to the samples metadata of a MultiDataSet
#' object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param df A tibble or data-frame of samples information, with at least column
#'   `id` (giving the sample IDs), and one row per sample ID.
#' @param datasets Character vector, name of the datasets to which the samples
#'   information should be added. If `NULL` (default value), the information
#'   will be added to the samples metadata of all datasets.
#' @returns A MultiDataSet object, with info from `df` added to its
#'   corresponding samples metadata.
#' @export
add_samples_metadata <- function(mo_data, df, datasets = NULL) {
  check_is_multidataset(mo_data)

  if (!inherits(df, "data.frame")) {
    stop("'df' argument should be a tibble or data-frame.")
  }

  if (!("id" %in% colnames(df))) {
    stop("'df' argument should have columns 'id'.")
  }

  check_unique <- max(table(df$id))
  if (check_unique > 1) {
    stop("'df' argument: should have only one row per sample ID.")
  }

  if (!is.null(datasets)) {
    .check_names(
      datasets,
      names(mo_data),
      "'datasets' argument: '_W_' are not datasets in 'mo_data'. Possible dataset names are: '_C_'."
    )
  } else {
    datasets <- names(mo_data)
  }

  sids <- unique(df$id)
  smeta_list <- get_samples_metadata(mo_data)
  samples_id <- get_samples(mo_data)

  for (i in datasets) {
    error_msg <- paste0(i, " dataset: ")

    wrong_ids <- setdiff(sids, samples_id[[i]]) |> length()
    missing_ids <- setdiff(samples_id[[i]], sids) |> length()

    if (wrong_ids > 0) {
      warning(
        error_msg,
        wrong_ids,
        " sample IDs not in 'mo_data', will be removed from samples metadata.",
        call. = FALSE
      )
    }
    if (missing_ids > 0) {
      warning(
        error_msg,
        missing_ids,
        " sample IDs missing from 'df' data-frame.",
        call. = FALSE
      )
    }

    df_ex <- smeta_list[[i]]

    existing_cols <- intersect(colnames(df_ex), colnames(df)) |>
      setdiff("id")
    if (length(existing_cols) > 0) {
      stop(
        error_msg,
        "'",
        paste0(existing_cols, collapse = "', '"),
        "' columns already in samples metadata table."
      )
    }

    df_new <- df_ex |>
      dplyr::left_join(df, by = "id") |>
      Biobase::AnnotatedDataFrame()
    rownames(df_new) <- df_new$id

    mo_data@phenoData[[i]]$main <- df_new
  }

  return(mo_data)
}


#' Get feature labels
#'
#' Extracts the feature labels from a MultiDataSet object given the name of the
#' column from the feature metadata of each dataset containing the feature
#' labels for this dataset.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param label_cols Character or named list of character, giving for each
#'   dataset the name of the column in the corresponding features metadata to
#'   use as label. If one value, will be used for all datasets. If list, the
#'   names must correspond to the names of the datasets in `mo_data`. Default
#'   value is `feature_id`, i.e. by default the ID of the features will be used
#'   as label.
#' @param truncate Integer, width to which the labels should be truncated (to
#'   avoid issues with very long labels in plots). If `NULL` (default value), no
#'   truncation will be performed.
#' @returns A tibble with columns `dataset`, `feature_id` and `label`.
#' @examples
#' \dontrun{
#' ## This works if each dataset in mo_data has in their features metadata table
#' ## a column called `name` that contains the feature labels.
#' get_features_labels(mo_data, label_cols = "name")
#'
#' ## If instead we want to use a different column for each dataset:
#' get_features_labels(
#'   mo_data,
#'   label_cols = list(
#'     "snps" = "feature_id",
#'     "rnaseq" = "gene_name",
#'     "metabolome" = "comp_formula"
#'   )
#' )
#'
#' ## If we want to use the feature IDs as labels for the genomics dataset,
#' ## we can simply remove it from the list (this is equivalent to the example
#' ## above):
#' get_features_labels(
#'   mo_data,
#'   label_cols = list(
#'     "rnaseq" = "gene_name",
#'     "metabolome" = "comp_formula"
#'   )
#' )
#' }
#' @export
get_features_labels <- function(mo_data, label_cols = "feature_id", truncate = NULL) {
  ## for devtools::check
  label <- feature_id <- dataset <- NULL

  check_is_multidataset(mo_data)
  datasets <- names(mo_data)

  ## Make a list of column name to use per dataset
  label_cols <- .make_var_list(
    label_cols,
    datasets,
    fixed_length = 1,
    allow_subset = TRUE
  )
  ## adding the missing datasets
  label_cols <- c(
    label_cols,
    setdiff(datasets, names(label_cols)) |>
      rlang::set_names() |>
      purrr::map(~"feature_id")
  )

  ## Checking that label_cols has valid input
  .check_input_var_fmetadata(label_cols, mo_data)

  fmeta <- get_features_metadata(mo_data)

  res <- purrr::map_dfr(
    datasets,
    \(.x) {
      fmeta[[.x]] |>
        tibble::as_tibble() |>
        dplyr::mutate(
          dataset = .x,
          label = !!sym(label_cols[[.x]]),
          label = dplyr::coalesce(label, feature_id)
        ) |>
        dplyr::select(dataset, feature_id, label)
    }
  )

  if (!is.null(truncate)) {
    res <- res |>
      dplyr::mutate(
        label = stringr::str_trunc(label, width = truncate)
      )
  }

  return(res)
}


#' Check for missing values in MultiDataSet
#'
#' Checks if there are missing values in each omics dataset of a MultiDataSet
#' object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns Invisible logical vector indicating whether missing values are
#'   present in each dataset.
#' @export
check_missing_values <- function(mo_data) {
  check_is_multidataset(mo_data)
  datasets <- get_datasets(mo_data) ## get matrices from each dataset

  msg <- names(datasets) |>
    rlang::set_names() |>
    purrr::map_chr(function(i) {
      missing_vals <- is.na(datasets[[i]])

      if (any(missing_vals)) {
        where_nas <- which(missing_vals, arr.ind = TRUE)
        n_nas <- nrow(where_nas)
        perc_nas <- 100 * n_nas / length(datasets[[i]])
        n_features <- length(unique(where_nas[, "row"]))
        n_samples <- length(unique(where_nas[, "col"]))
        msg <- paste0(
          n_nas, " (", round(perc_nas, 2), "%) missing values in ",
          i, " dataset, across ",
          n_features, " feature", ifelse(n_features > 1, "s", ""), " and ",
          n_samples, " sample", ifelse(n_samples > 1, "s", ""), "."
        )
      } else {
        msg <- paste("No missing values in", i, "dataset.")
      }
      msg
    })

  res <- purrr::map_lgl(msg, function(i) {
    !stringr::str_detect(i, "No missing values")
  })

  msg <- paste0(msg, collapse = "\n")

  message(msg)

  return(invisible(res))
}

## Because somehow I can't use the subset method from MultiDataSet
.subset_mds <- getMethod("subset", "MultiDataSet")

#' Number of features in each dataset of MultiDataSet object
#'
#' Gives the number of features in each dataset of a MultiDataSet object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns A named integer vector, where each element is the number of features
#'   in the corresponding dataset.
#' @export
n_features <- function(mo_data) {
  check_is_multidataset(mo_data)
  return(vapply(mo_data@featureData, function(y) nrow(y[[1]]), numeric(1)))
}

#' Number of samples in each dataset of MultiDataSet object
#'
#' Gives the number of samples in each dataset of a MultiDataSet object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns A named integer vector, where each element is the number of sample
#'   in the corresponding dataset.
#' @export
n_samples <- function(mo_data) {
  check_is_multidataset(mo_data)
  return(vapply(mo_data@phenoData, function(y) nrow(y[[1]]), numeric(1)))
}

#' Subset a MultiDataSet object by feature
#'
#' Subsets a MultiDataSet object based on a list of feature IDs provided.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param features_id Character vector, a vector of feature IDs (from across the
#'   datasets) to select. Also accepts lists (e.g. list with a vector of feature
#'   IDs per dataset).
#' @returns A [MultiDataSet::MultiDataSet-class] object with only features
#'   specified.
#' @examples
#' \dontrun{
#' ## works with a vector of feature IDs:
#' subset_features(mo_data, c("featureA", "featureB", "featureC"))
#'
#' ## or with a list of feature IDs (typically one per dataset, but doesn't
#' ## have to be):
#' subset_features(
#'   mo_data,
#'   list(
#'     c("omics1_featureA", "omics1_featureB", "omics1_featureC"),
#'     c("omics2_featureA", "omics2_featureB"),
#'     c("omics3_featureA", "omics3_featureB", "omics3_featureC"),
#'   )
#' )
#' }
#' @export
subset_features <- function(mo_data, features_id) {
  features_id <- unique(unlist(features_id))
  cmd <- paste0(
    ".subset_mds(mo_data, feature_id %in% c('",
    paste0(features_id, collapse = "', '"),
    "'))"
  )
  res <- eval(str2expression(cmd))

  return(res)
}

#' Replace matrix dataset within a MultiDataSet object
#'
#' Replaces the matrix for an omics dataset by a new matrix in a MultiDataSet
#' object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param dataset_name Character, the name of the dataset for which the matrix
#'   data should be changed.
#' @param new_data Matrix, the new data. Should have features as rows and
#'   samples as columns. Rownames should match the corresponding feature IDs,
#'   colnames should match the corresponding sample IDs.
#' @returns A MultiDataSet object.
#' @export
replace_dataset <- function(mo_data, dataset_name, new_data) {
  check_is_multidataset(mo_data)
  if (length(dataset_name) > 1) {
    stop("'dataset_name' argument should be of length 1.")
  }
  .check_names(
    dataset_name,
    names(mo_data),
    "'dataset_name' argument: '_W_' is not an existing dataset. Possible dataset names are: '_C_'."
  )

  ## Compare dimensions, rownames and column names
  old_data <- get_datasets(mo_data[, dataset_name])[[1]]
  if (any(dim(old_data) != dim(new_data))) {
    stop(
      "'new_data' argument has incorrect dimensions. Should have ",
      nrow(old_data), " rows (features) and ",
      ncol(old_data), " columns (samples)."
    )
  }
  if (is.null(rownames(new_data)) || is.null(colnames(new_data))) {
    stop("'new_data' should have feature IDs as rownames and sample IDs as column names.")
  }
  if (any(sort(rownames(new_data)) != sort(rownames(old_data)))) {
    stop("'new_data' rownames do not match feature IDs in existing dataset.")
  }
  if (any(sort(colnames(new_data)) != sort(colnames(old_data)))) {
    stop("'new_data' colnames do not match sample IDs in existing dataset.")
  }
  new_data <- new_data[rownames(old_data), colnames(old_data)]

  ## Create a new omics dataset with the new data matrix
  new_eset <- mo_data[[dataset_name]]
  assay_name <- Biobase::assayDataElementNames(new_eset)[1]

  ## We'll replace the first assay but we need to copy the other assays to
  ## the new object
  other_assays <- Biobase::assayDataElementNames(new_eset)[-1] |>
    purrr::map_chr(
      function(j) {
        paste0(", ", j, " = Biobase::assayDataElement(new_eset, '", j, "')")
      }
    ) |>
    paste0(collapse = "")
  assay_expr <- str2expression(
    paste0(
      "Biobase::assayDataNew(storage.mode = 'lockedEnvironment', ",
      assay_name, " = new_data",
      other_assays, ")")
  )
  Biobase::assayData(new_eset) <- eval(assay_expr)

  ds_name <- stringr::str_extract(dataset_name, "(?<=\\+).+$")
  if (is.na(ds_name)) ds_name <- NULL

  ## We know that add_omics_set will return a warning
  ## as we are overwriting the set
  res <- suppressWarnings(
    add_omics_set(
      mo_data,
      new_eset,
      ds_name,
      overwrite = TRUE
    )
  )

  return(res)
}


#' Round values in omics dataset from MultiDataSet object
#'
#' Rounds the values of a given omics dataset within a MultiDataSet object. Can
#' also limit the range of possible values.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param dataset_name Character, the name of the dataset for which the matrix
#'   data should be changed.
#' @param ndecimals Integer, the number of decimals to keep in the dataset.
#'   Default value is 0.
#' @param min_val Numeric, the minimum value allowed in the dataset. Values
#'   below `min_val` will be set to `min_val`.
#' @param max_val Numeric, the maximum value allowed in the dataset. Values
#'   above `max_val` will be set to `max_val`.
#' @returns A MultiDataSet object.
#' @examples
#' \dontrun{
#' ## Let's imagine that we imputed missing values in the genomics dataset from
#' ## mo_data using NIPALS-PCA. The imputed values are continuous, but the
#' ## dataset contains dosage values for a diploid organism (i.e. values can
#' be 0, 1, 2). We'll round the imputed values and make sure they can't be
#' ## negative or higher than 2.
#' round_dataset(mo_data, "snps", min_val = 0, max_val = 2)
#' }
#' @export
round_dataset <- function(mo_data, dataset_name, ndecimals = 0, min_val = -Inf, max_val = Inf) {
  check_is_multidataset(mo_data)
  .check_names(
    dataset_name,
    names(mo_data),
    "'dataset_name' argument: '_W_' is not an existing dataset. Possible dataset names are: '_C_'."
  )

  mat <- get_datasets(mo_data)[[dataset_name]]
  mat <- round(mat, ndecimals)
  mat[mat < min_val] <- min_val
  mat[mat > max_val] <- max_val

  res <- replace_dataset(mo_data, dataset_name, mat)

  return(res)
}
