#' Create a Biobase set object to store omics data
#'
#' Creates a Biobase object to store an omics dataset and associated samples and
#' features metadata.
#'
#' @param dataset Matrix, the omics dataset in matrix form with features as rows
#'   and samples as columns.
#' @param omics_type Character, which type of omics data is being stored?
#'   Possible values are `'genomics'`, `'transcriptomics'`, `'metabolomics'` and
#'   `'phenomics'`. Use `'phenomics'` for any other omics.
#' @param features_metadata Data.frame, a feature annotation data-frame with
#'   features as rows and information about the features as columns. The number
#'   of rows and row names must match those of the `dataset` matrix.
#' @param samples_metadata Data.frame, a samples information data-frame with
#'   samples as rows and information about the samples as columns. The number of
#'   rows and row names must match the number of columns and column names of the
#'   `dataset` matrix.
#' @returns An object derived from [Biobase::eSet-class]:
#' * if `omics_type = 'genomics'`: a [Biobase::SnpSet-class] object;
#' * if `omics_type = 'transcriptomics'`: a [Biobase::ExpressionSet-class]
#'   object.
#' * if `omics_type = 'metabolomics'`: a [MetabolomeSet-class] object.
#' * if `omics_type = 'phenomics'` a [PhenotypeSet-class] object.
#' @examples
#' \dontrun{
#' data_geno <- import_dataset_csv(
#' "genotype_dataset.csv",
#'   col_id = "Marker",
#'   features_as_rows = TRUE
#' )
#' geno_info_features <- import_fmetadata_csv(
#' "genotype_features_info.csv",
#'   col_id = "Marker"
#' )
#' samples_information <- import_smetadata_csv(
#' "samples_information.csv",
#'   col_id = "Sample"
#' )
#' create_omics_set(
#'   dataset = data_geno,
#'   omics_type = "genomics",
#'   features_metadata = geno_info_features,
#'   samples_metadata = samples_information
#' )
#' }
#' @export
create_omics_set <- function(dataset,
                             omics_type = c("phenomics",
                                            "genomics",
                                            "transcriptomics",
                                            "metabolomics"),
                             features_metadata = NULL,
                             samples_metadata = NULL) {

  omics_type <- rlang::arg_match(omics_type)

  res <- switch(
    omics_type,
    "genomics" = new("SnpSet", call = dataset),
    "transcriptomics" = Biobase::ExpressionSet(assayData = dataset),
    "metabolomics" = new("MetabolomeSet", call = dataset),
    "phenomics" = new("PhenotypeSet", call = dataset)
  )

  if (!is.null(features_metadata)) {

    missing_features <- setdiff(rownames(dataset), rownames(features_metadata))
    if (length(missing_features) > 0) {
      warning(
        length(missing_features),
        " features are not present in feature metadata."
      )
      ## adding the missing features as NAs in the metadata data-frame
      features_metadata[missing_features, ] <- NA
      features_metadata[missing_features, "feature_id"] <- missing_features
    }

    irrelevant_features <- setdiff(
      rownames(features_metadata),
      rownames(dataset)
    )
    if (length(irrelevant_features) > 0) {
      warning(
        length(irrelevant_features), " features",
        " in feature metadata not in dataset, will be removed from metadata."
      )
    }

    features_metadata <- features_metadata[rownames(dataset), ]
    feature_data <- Biobase::AnnotatedDataFrame(features_metadata)

    Biobase::featureData(res) <- feature_data
  } else {
    ## we need to have in the featureData a column with feature ID to be able
    ## to perform subsetting later on
    temp <- data.frame("feature_id" = rownames(dataset))
    rownames(temp) <- rownames(dataset)

    Biobase::featureData(res) <- Biobase::AnnotatedDataFrame(temp)
  }

  if (!is.null(samples_metadata)) {

    missing_samples <- setdiff(colnames(dataset), rownames(samples_metadata))
    if (length(missing_samples)) {
      warning(
        length(missing_samples),
        " samples are not present in samples metadata."
      )
      ## adding the missing samples as NAs in the metadata data-frame
      samples_metadata[missing_samples, ] <- NA
      samples_metadata[missing_samples, "id"] <- missing_samples
    }

    irrelevant_samples <- setdiff(
      rownames(samples_metadata),
      colnames(dataset)
    )
    if (length(irrelevant_samples)) {
      warning(
        length(irrelevant_samples), " samples",
        " in samples metadata not in dataset, will be removed from metadata."
      )
    }

    samples_metadata <- samples_metadata[colnames(dataset), ]
    samples_data <- Biobase::AnnotatedDataFrame(samples_metadata)

    Biobase::phenoData(res) <- samples_data
  }

  return(res)
}

#' Target factory for omics sets creation
#'
#' Creates a list of targets that generate omics sets from targets containing
#' datasets, features and samples metadata.
#'
#' @param datasets Vector of symbols, the names of the targets containing the
#'   omics datasets.
#' @param omics_types Character vector, which type of omics data is being stored
#'   for each dataset? Possible values are `'genomics'`, `'transcriptomics'`,
#'   `'metabolomics'` and `'phenomics'`.  Use `'phenomics'` for any other omics.
#'   Use `'phenomics'` for any other omics.
#' @param features_metadatas Vector of symbols, the names of the targets
#'   containing the features metadata data-frame associated with each omics
#'   dataset. Use `NULL` if no feature metadata exists for a dataset.
#' @param samples_metadatas Vector of symbols, the names of the targets
#'   containing the samples metadata data-frame associated with each omics
#'   dataset. Use `NULL` if no samples metadata exists for a dataset.
#' @param target_name_suffixes Character vector, a suffix to add to the name of
#'   the targets created by this target factory for each dataset. If none
#'   provided, the suffixes will be extracted from the `datasets` argument.
#'   Default value is NULL.
#' @return A list of target objects, with three datasets provided, and
#'   `target_name_suffixes = c("geno", "transcripto", "metabo")`, the following
#'   targets will be returned: `set_geno`, `set_transcripto` and `set_metabo`.
#' @examples
#' \dontrun{
#' ## in the _targets.R
#' library(moiraine)
#' library(targets)
#'
#' list(
#'   ## targets to import the different datasets
#'
#'   ## Example where genomics dataset has no features metadata information
#'   ## Will generate the following targets: set_geno, set_transcripto
#'   create_omics_set_factory(
#'     datasets = c(data_geno, data_transcripto),
#'     omics_types = c("genomics", "transcriptomics"),
#'     features_metadata = c(NULL, fmeta_transcripto),
#'     samples_metadata = c(smeta_geno, smeta_transcripto)
#'   )
#' )
#' }
#' @export
create_omics_set_factory <- function(datasets,
                                     omics_types,
                                     features_metadatas = NULL,
                                     samples_metadatas = NULL,
                                     target_name_suffixes = NULL) {
  n_datasets <- length(.input_to_symbol(rlang::enquo(datasets)))
  n_fmeta <- length(.input_to_symbol(rlang::enquo(features_metadatas)))
  n_smeta <- length(.input_to_symbol(rlang::enquo(samples_metadatas)))

  if (length(omics_types) != n_datasets) {
    stop("'datasets' and 'omics_types' vectors must have the same length.")
  }

  if (length(target_name_suffixes) != n_datasets &&
        !is.null(target_name_suffixes)) {
    stop("'datasets' and 'target_name_suffixes' vectors must have the same length.")
  }

  if (n_fmeta != n_datasets && !is.null(substitute(features_metadatas))) {
    stop("'datasets' and 'features_metadatas' vectors must have the same length.")
  }

  if (n_smeta != n_datasets && !is.null(substitute(samples_metadatas))) {
    stop("'datasets' and 'samples_metadatas' vectors must have the same length.")
  }

  values <- tibble::tibble(
    ds = .input_to_symbol(rlang::enquo(datasets)),
    ot = omics_types,
    fmet = .input_to_symbol(rlang::enquo(features_metadatas)),
    smet = .input_to_symbol(rlang::enquo(samples_metadatas))
  )

  if (is.null(target_name_suffixes)) {
    target_name_suffixes <- stringr::str_remove(
      .symbol_to_char(rlang::enquo(datasets)),
      "data_"
    )
  }

  values$target_name_suffix <- target_name_suffixes

  targets <- tarchetypes::tar_map(
    values = values,
    names = "target_name_suffix",
    targets::tar_target_raw(
      "set",
      substitute(
        create_omics_set(
          dataset = ds,
          omics_type = ot,
          features_metadata = fmet,
          samples_metadata = smet
        )
      )
    )
  )

  return(targets)
}

#' Adds an omics set to a MultiDataSet object
#'
#' Adds a omics set to an existing MultiDataSet object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param omics_set A [Biobase::eSet-class] object, created via
#'   [create_omics_set()]. Currently accepted objects: [Biobase::SnpSet-class],
#'   [Biobase::ExpressionSet-class], [MetabolomeSet-class],
#'   [PhenotypeSet-class].
#' @param ds_name Character, name of the dataset (will be used as suffix for the
#'   name of the dataset in the resulting MultiDataSet object).
#' @param ... Further arguments passed to `[MultiDataSet::add_snps()],
#'   [MultiDataSet::add_rnaseq()], [add_metabo()] or[add_pheno()] (depending on
#'   `omics_set` class).
#' @returns A [MultiDataSet::MultiDataSet-class] object, the `mo_data` with
#'   `omics_set` as an additional dataset.
#' @examples
#' \dontrun{
#' add_omics_set(mo_data, omics_set, "exp1")
#' }
#' @export
add_omics_set <- function(mo_data, omics_set, ds_name, ...) {
  check_is_multidataset(mo_data)

  type <- class(omics_set)[[1]]
  add_function <- switch(
    type,
    "SnpSet" = MultiDataSet::add_snps,
    "ExpressionSet" = MultiDataSet::add_rnaseq,
    "MetabolomeSet" = add_metabo,
    "PhenotypeSet" = add_pheno,
    stop("object ", type, " cannot be added to MultiDataSet object.")
  )
  add_function(mo_data, omics_set, dataset.name = ds_name, ...)
}

#' Create a MultiDataSet object to store multi-omics data
#'
#' Creates a MultiDataSet object from a list of Biobase Set objects to store the
#' different omics sets.
#'
#' @param sets_list List of [Biobase::eSet-class] objects, created via
#'   [create_omics_set()]. Currently accepted objects: [Biobase::SnpSet-class],
#'   [Biobase::ExpressionSet-class], [MetabolomeSet-class],
#'   [PhenotypeSet-class].
#' @param datasets_names Optional, vector of character, name for each Set
#'   object. Will be appended to the data type in the resulting object. If the
#'   `sets_list` list contains several objects of the same data type (e.g.
#'   several SnpSets), their names must be unique. If "" is provided, no name
#'   will be appended to the data type for the corresponding dataset.
#' @param show_warnings Logical, should warnings be displayed when adding a set
#'   to the MultiDataSet object? Default value is `TRUE`.
#' @returns [MultiDataSet::MultiDataSet-class] object.
#' @examples
#' \dontrun{
#' ## set_geno, set_transcripto and set_metabo are all Set objects
#' ## Generating a MultiDataSet object with standard name
#' create_multiomics_set(
#'   list(set_geno, set_transcripto, set_metabo)
#' )
#'
#' ## Adding custom names for genomics and metabolomics datasets
#' ## but not for the transcriptomics dataset
#' create_multiomics_set(
#'   list(set_geno, set_transcripto, set_metabo),
#'   datasets_names = c("genome1", "", "lcms")
#' )
#' }
#' @export
create_multiomics_set <- function(sets_list,
                                  datasets_names = NULL,
                                  show_warnings = TRUE) {

  if (!length(sets_list)) {
    stop("sets_list list is empty.")
  }

  cl <- purrr::map_chr(sets_list, \(x) class(x)[[1]])
  poss_classes <- c("SnpSet", "ExpressionSet", "MetabolomeSet", "PhenotypeSet")

  if (!all(cl %in% poss_classes)) {
    stop(
      "Elements in sets_list must be ",
      paste0(poss_classes, collapse = ", "),
      " objects."
    )
  }

  if (!is.null(datasets_names)) {
    if (length(datasets_names) != length(sets_list)) {
      stop("dataset_names vector must have same length as sets_list list.")
    }

    if (any(duplicated(paste0(cl, datasets_names)))) {
      stop("Dataset names for objects of a same type must be unique.")
    }
  } else {
    ## if several sets_list elements have the same data type,
    ## must provide unique datasets_names
    datasets_names <- .make_unique_ids(cl)
  }

  res <- MultiDataSet::createMultiDataSet()

  for (i in seq_along(sets_list)) {
    x <- sets_list[[i]]
    ds_name <- datasets_names[[i]]
    if (ds_name == "") ds_name <- NULL

    res <- add_omics_set(
      res,
      x,
      ds_name,
      warnings = show_warnings
    )
  }

  ## Check that the samples metadata is consistent across the datasets
  temp <- Biobase::pData(res) |>
    purrr::reduce(
      function(.x, .y) {
        cc <- intersect(colnames(.x), colnames(.y))
        dplyr::full_join(.x, .y, by = cc)
      }
    )

  duplicated_samples <- duplicated(temp$id)
  if (any(duplicated_samples)) {
    stop(
      "Conflicting information in samples metadata for samples ",
      paste0(temp$id[duplicated_samples], collapse = ", "),
      "."
    )
  }

  return(res)
}

#' Convert symbols to character
#'
#' Converts a vector of symbols into a character vector of their values.
#'
#' @param x Vector of symbols
#' @returns A character vector with the same size as `x`.
#'
#' @noRd
.symbol_to_char <- function(x) {
  string <- rlang::quo_text(x, nlines = Inf)
  unlist(strsplit(gsub("(c\\(|\\)|\\s)", "", string), ","))
}

#' Converts a vector into symbols vector
#'
#' Converts an input vector into a vector of symbols, to be used for targets
#' factories. If input value is `NULL`, will return `expression(NULL)`.
#'
#' @param x Vector.
#' @returns a vector of symbols with same length as `x`.
#'
#' @noRd
.input_to_symbol <- function(x) {
  string <- rlang::quo_text(x, nlines = Inf)

  if (string == "NULL") {
    return(str2expression("NULL"))
  }

  string <- unlist(strsplit(gsub("(c\\(|\\)|\\s)", "", string), ","))
  res <- str2expression(string)

  return(res)
}

#' Generate suffixes to make values unique.
#'
#' For a given vector of values, returns a vector of the same size containing
#' suffixes to make each value in the input vector unique.
#'
#' @param x Character vector of IDs.
#' @returns Character vector of the same size as `x`, containing at each
#'   position:
#' * `''` if the corresponding value in `x` is unique;
#' * A number as character (e.g. `'1'`, `'2`, etc) if the corresponding value
#'   in `x` is not unique, such that by adding this value to the values in `x`,
#'   it becomes unique.
#' @examples
#' \dontrun{
#' .make_unique_ids(1:5)
#' #> "" "" "" "" ""
#'
#' .make_unique_ids(c(1:5, 1))
#' #> "1" ""  ""  ""  ""  "2"
#' }
#'
#' @noRd
.make_unique_ids <- function(x) {
  are_duplicates <- unique(x[duplicated(x)])
  res <- character(length(x))
  for (v in are_duplicates) {
    indx <- x == v
    res[indx] <- seq_len(sum(indx))
  }

  return(res)
}
