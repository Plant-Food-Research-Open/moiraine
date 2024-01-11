#' Import a dataset from a csv file
#'
#' Reads a csv file and returns a matrix in which the rows corresponds to
#' features (e.g. markers, genes, phenotypes...) and the columns correspond to
#' samples/observations.
#'
#' @param file Character, path to the dataset csv file.
#' @param col_id Character, the name of the column in the file that contains the
#'   ID of the rows (i.e. feature IDs if `features_as_rows` is `TRUE`, or sample
#'   IDs if `features_as_rows` is `FALSE`).
#' @param features_as_rows Logical, do the rows in the file correspond to
#'   features? Default value is `TRUE`, i.e. the file contains features as rows
#'   and samples as columns.
#' @param ... Further arguments passed to [readr::read_csv()].
#' @returns A matrix with the samples as columns and the features as rows.
#'   Feature IDs are used as row names and sample IDs as column names.
#' @examples
#' \dontrun{
#' data_geno <- import_dataset_csv(
#'   "genotype_dataset.csv",
#'   col_id = "Marker",
#'   features_as_rows = TRUE
#' )
#' data_pheno <- import_dataset_csv(
#'   "phenotype_dataset.csv",
#'   col_id = "Sample",
#'   features_as_rows = FALSE
#' )
#' }
#' @export
import_dataset_csv <- function(file, col_id, features_as_rows = TRUE, ...) {
  ## for devtools::check
  value <- feature <- NULL

  data <- readr::read_csv(file, ...)

  if (!(col_id %in% names(data))) {
    stop("'", col_id, "' is not a column in the csv file. ",
         "Please specify a valid column name for argument col_id.")
  }

  if (features_as_rows) {
    res <- as.matrix(dplyr::select(data, -!!sym(col_id)))
    rownames(res) <- data[[col_id]]
  } else {
    ## Pivot the tibble to have the samples as columns
    temp <- data |>
      tidyr::pivot_longer(
        cols = -!!sym(col_id),
        names_to = "feature",
        values_to = "value"
      ) |>
      tidyr::pivot_wider(
        names_from = !!sym(col_id),
        values_from = value
      )

    ## Transform the result into a matrix with appropriate feature rownames
    res <- as.matrix(dplyr::select(temp, -feature))
    rownames(res) <- temp$feature
  }

  return(res)
}

#' Import feature metadata from a csv file
#'
#' Reads a csv file and returns a dataframe in which the rows correspond to
#' features (e.g. markers, genes, phenotypes...) and columns correspond to
#' information about the features. Non-ASCII characters are replaced with ASCII
#' equivalents (using the stringi and textclean packages).
#'
#' @param file Character, path to the dataset csv file.
#' @param col_id Character, the name of the column in the file that contains the
#'   feature IDs.
#' @param ... Further arguments passed to [readr::read_csv()].
#' @returns A data-frame with the features as rows and features information as
#'   columns. Feature IDs are used as row names.
#' @examples
#' \dontrun{
#' geno_info_features <- import_fmetadata_csv(
#'   "genotype_features_info.csv",
#'   col_id = "Marker"
#' )
#' }
#' @export
import_fmetadata_csv <- function(file, col_id, ...) {
  ## for devtools::check
  feature_id <- rownames <- NULL

  data <- readr::read_csv(file, ...)

  if (!(col_id %in% names(data))) {
    stop("'", col_id, "' is not a column in the csv file. ",
         "Please specify a valid column name for argument col_id.")
  }

  res <- data |>
    dplyr::rename(feature_id = !!sym(col_id)) |>
    dplyr::select(feature_id, tidyselect::everything()) |>
    dplyr::mutate(rownames = feature_id) |>
    tibble::column_to_rownames("rownames") |>
    as.data.frame() |>
    .replace_non_ascii_chars()

  return(res)
}

#' Import features metadata from a GFF/GTF file
#'
#' Reads a GFF or GTF annotation file and returns a dataframe in which the rows
#' correspond to features (e.g. genes or transcripts) and columns correspond to
#' information about the features. Non-ASCII characters are replaced with ASCII
#' equivalents (using the stringi and textclean packages).
#'
#' @param file Character, path to the dataset GFF or GTF file.
#' @param feature_type Character, the type of feature to extract from the
#' annotation file. Currently supports `'genes'` and `'transcripts'`.
#' @param add_fields Character vector, fields in the GFF/GTF file to extract
#' that are not imported by default (only use if you've run the function once
#' and realised that some fields are not extracted by the function).
#' @returns A data-frame with the features as rows and features information as
#' columns. Feature IDs are used as row names.
#' @examples
#' \dontrun{
#' import_fmetadata_gff(
#'   "bos_taurus_gene_model.gff3",
#'   "genes",
#'   add_fields = c("name", "description")
#' )
#' }
#' @export
import_fmetadata_gff <- function(file, feature_type, add_fields = NULL) {
  ## for devtools::check
  feature_id <- NULL

  ## For devtools check
  ID <- seqnames <- NULL

  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop(
      "Package \"GenomicFeatures\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    stop(
      "Package \"rtracklayer\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!(feature_type) %in% c("genes", "transcripts")) {
    stop("feature_type argument should be 'genes' or 'transcripts'.")
  }

  ## Read the GFF/GTF file into a TxDb object (for easier manipulation)
  annotation <- GenomicFeatures::makeTxDbFromGFF(file, format = "auto")

  ref_df <- switch(feature_type,
    "transcripts" = GenomicFeatures::transcripts(annotation, use.names = TRUE),
    "genes" = GenomicFeatures::genes(annotation)
  )

  res <- as.data.frame(ref_df)
  res$feature_id <- names(ref_df)
  rownames(res) <- res$feature_id

  ## IMPROVISED PART: might need to change based on other GFF/GTF files
  if (!is.null(add_fields)) {
    ## To go around the fact that makeTxDbFromGFF doesn't import all fields
    format <- stringr::str_extract(file, "(?<=\\.)\\w+$")
    if (format == "gff") format <- "gff3"
    gr <- rtracklayer::import(file, format = format)
    ## some of the columns of the gr_df dataframe might be lists, we want to
    ## get rid of that (might be more efficient ways to do that)
    gr_df <- as.data.frame(gr) |>
      ## if there's only 1 field we want to add, prevents 1 column dataframe to
      ## become vector
      tibble::as_tibble() |>
      dplyr::mutate(
        dplyr::across(where(is.list), ~ sapply(.x, paste, collapse = " "))
      )

    .check_names(
      add_fields,
      names(gr_df),
      "add_fields argument: Fields '_W_' are not present in the annotation file."
    )

    if (format == "gff3") {
      gr_df <- gr_df |>
        dplyr::mutate(
          ID = stringr::str_remove(ID, "^((gene)|(CDS)|(transcript)|(region)):")
        )
    } else if (format == "gtf") {
      id_col <- ifelse(feature_type == "genes", "gene_id", "transcript_id")
      gr_df <- gr_df[, c(id_col, add_fields)] |>
        dplyr::rename(ID = !!sym(id_col))
    }

    gr_df <- gr_df |>
      dplyr::filter(ID %in% res$feature_id) |>
      dplyr::distinct()

    ## TO IMPROVE: safeguard against selecting a exon property field
    ## and trying to import it with the genes
    if (nrow(gr_df) != nrow(res)) {
      stop(
        "Problem in importing field(s) '",
        paste0(add_fields, collapse = "', '"),
        "': check that it has one unique value per feature."
      )
    }

    gr_df <- gr_df[match(res$feature_id, gr_df$ID), c("ID", add_fields)]
    res <- dplyr::bind_cols(res, gr_df[, -1]) ## we're removing the ID column
  }

  ## rename the seqnames column, might change in the future
  res <- res |>
    dplyr::rename(chromosome = seqnames)

  ## Trying to avoid duplicate columns
  if ((feature_type == "genes") && ("gene_id" %in% colnames(res))) {
    res <- dplyr::select(res, -!!sym("gene_id"))
  }
  if ((feature_type == "transcripts") && ("tx_name" %in% colnames(res))) {
    res <- dplyr::select(res, -!!sym("tx_name"), -!!sym("tx_id"))
  }

  res <- dplyr::select(res, feature_id, dplyr::everything()) |>
    .replace_non_ascii_chars()

  return(res)
}

#' Import samples metadata from a csv file
#'
#' Reads a csv file and returns a dataframe in which the rows
#' correspond to features (e.g. markers, genes, phenotypes...) and columns
#' correspond to information about the features.
#'
#' @param file Character, path to the dataset csv file.
#' @param col_id Character, the name of the column in the file that contains
#' the ID of the rows (i.e. sample IDs).
#' @param ... Further arguments passed to [readr::read_csv()].
#' @returns A data-frame with the samples as rows and the samples properties as
#' columns. Sample IDs are used as rownames.
#' @examples
#' \dontrun{
#' samples_information <- import_smetadata_csv(
#'   "samples_information.csv",
#'   col_id = "Sample"
#' )
#' }
#' @export
import_smetadata_csv <- function(file, col_id, ...) {
  ## For devtools::check
  id <- NULL

  data <- readr::read_csv(file, ...)

  if (!(col_id %in% names(data))) {
    stop(
      "'", col_id, "' is not a column in the csv file. ",
      "Please specify a valid column name for argument col_id."
    )
  }

  res <- data |>
    tibble::column_to_rownames(col_id) |>
    as.data.frame()

  ## Adding an ID column to get rid of the MultiDataSet warnings when
  ## creating a multi-omics set
  res$id <- rownames(res)
  res <- dplyr::select(res, id, dplyr::everything())

  return(res)
}


#' Target factory for csv datasets import
#'
#' Creates a list of targets that track each file and import the dataset from
#' each csv file.
#'
#' @param files Character vector, vector of paths to the dataset csv files.
#' @param col_ids Character vector, the name of the column in each file that
#' contains the ID of the rows (i.e. feature IDs if value in `features_as_rowss`
#' is `TRUE` for the corresponding dataset, or sample IDs if value in
#' `features_as_rowss` is `FALSE`).
#' @param features_as_rowss Logical vector, do the rows in each file correspond
#' to features?
#' @param target_name_suffixes Character vector, a suffix to add to the name of
#' the targets created by this target factory for each dataset.
#' @returns A list of target objects. For example, with two files to import and
#' `target_name_suffixes = c("geno", "transcripto")`, the factory returns the
#' following targets:
#' * `dataset_file_geno` and `dataset_file_transcripto`: targets tracking the
#' genomics dataset file and the transcriptomics dataset file, respectively.
#' * `data_geno` and `data_transcripto`: targets that import the genomics and
#' transcriptomics dataset, respectively.
#' @examples
#' \dontrun{
#' ## in the _targets.R
#' library(MOITestPipeline)
#'
#' list(
#'   import_dataset_csv_factory(
#'     c(
#'       "data/genotype_data.csv",
#'       "data/rnaseq_data.csv"
#'     ),
#'     col_ids = c("Marker", "Sample"),
#'     features_as_rows = c(TRUE, FALSE),
#'     target_name_suffixes = c("geno", "transcripto")
#'   )
#' )
#' }
#' @export
import_dataset_csv_factory <- function(files,
                                       col_ids,
                                       features_as_rowss,
                                       target_name_suffixes) {
  ## for devtools::check
  target_name_suffix <- NULL

  values <- tibble::tibble(
    file = files,
    col_id = col_ids,
    features_as_rows = features_as_rowss,
    target_name_suffix = target_name_suffixes
  ) |>
    dplyr::mutate(
      file_target = rlang::syms(paste0("dataset_file_", target_name_suffix))
    )

  ## based on static branching
  targets <- tarchetypes::tar_map(
    values = values,
    names = "target_name_suffix",
    targets::tar_target_raw(
      "dataset_file",
      as.symbol("file"),
      format = "file"
    ),
    targets::tar_target_raw(
      "data",
      substitute(
        import_dataset_csv(
          file_target,
          col_id = col_id,
          features_as_rows = features_as_rows
        )
      )
    )
  )

  return(targets)
}

#' Target factory for csv features metadata import
#'
#' Creates a list of targets that track each file and import the features
#' metadata from each csv file.
#'
#' @param files Character vector, vector of paths to the features metadata
#' csv files.
#' @param col_ids Character vector, the name of the column in each file that
#' contains the features ID.
#' @param target_name_suffixes Character vector, a suffix to add to the name of
#' the targets created by this target factory for each dataset.
#' @returns A list of target objects. For example, with two files to import and
#' `target_name_suffixes = c("geno", "transcripto")`,
#' the factory returns the following targets:
#' * `fmetadata_file_geno` and `fmetadata_file_transcripto`: targets tracking
#' the genomics and transcriptomics features metadata files, respectively.
#' * `fmetadata_geno` and `fmetadata_transcripto`: targets that import the
#' genomics and transcriptomics features metadata dataset.
#' @examples
#' \dontrun{
#' ## in the _targets.R
#' library(MOITestPipeline)
#'
#' list(
#'   import_fmetadata_csv_factory(
#'     c(
#'       "data/genotype_fmetadata.csv",
#'       "data/rnaseq_fmetadata.csv"
#'     ),
#'     col_ids = c("Marker", "Info"),
#'     target_name_suffixes = c("geno", "transcripto")
#'   )
#' )
#' }
#' @export
import_fmetadata_csv_factory <- function(files,
                                         col_ids,
                                         target_name_suffixes) {
  ## for devtools::check
  target_name_suffix <- NULL

  values <- tibble::tibble(
    file = files,
    col_id = col_ids,
    target_name_suffix = target_name_suffixes
  ) |>
    dplyr::mutate(
      file_target = rlang::syms(paste0("fmetadata_file_", target_name_suffix))
    )

  ## based on static branching
  targets <- tarchetypes::tar_map(
    values = values,
    names = "target_name_suffix",
    targets::tar_target_raw(
      "fmetadata_file",
      as.symbol("file"),
      format = "file"
    ),
    targets::tar_target_raw(
      "fmetadata",
      substitute(
        import_fmetadata_csv(
          file_target,
          col_id = col_id
        )
      )
    )
  )

  return(targets)
}


#' Target factory for GFF/GTF features metadata import
#'
#' Creates a list of targets that track each file and import the features
#' metadata from each GFF/GTF file.
#'
#' @param files Character vector, vector of paths to the samples metadata GFF
#' or GTF files.
#' @param feature_types Character vector, the type of features to extract from
#' each annotation file. Currently supports `'genes'` and `'transcripts'`.
#' @param add_fieldss List, where each element is a character vector of  field
#' names in each GFF/GTF file to extract that are not imported by default.
#' If a character vector is provided, will be used for all files to read in.
#' @param target_name_suffixes Character vector, a suffix to add to the name of
#' the targets created by this target factory for each dataset.
#' @returns A list of target objects. For example, with two files to import and
#' `target_name_suffixes = c("geno", "transcripto")`, the factory returns the
#' following targets:
#' * `fmetadata_file_geno` and `fmetadata_file_transcripto`: targets tracking
#' the genomics and transcriptomics annotation files, respectively.
#' * `fmetadata_geno` and `fmetadata_transcripto`: targets that import the
#' genomics and transcriptomics features metadata datasets, respectively.
#' @examples
#' \dontrun{
#' ## in the _targets.R
#' library(MOITestPipeline)
#'
#' list(
#'   import_fmetadata_gff_factory(
#'     c(
#'       "data/annotation.gff",
#'       "data/annotationv2.gtf"
#'     ),
#'     feature_types = c("genes", "transcripts"),
#'     add_fieldss = list(
#'       c("gene_name", "gene_custom_ID"),
#'       c("transcript_name")
#'     ),
#'     target_name_suffixes = c("geno", "transcripto")
#'   )
#' )
#' }
#' @export
import_fmetadata_gff_factory <- function(files,
                                         feature_types,
                                         add_fieldss,
                                         target_name_suffixes) {
  ## for devtools::check
  target_name_suffix <- NULL

  if (!is.list(add_fieldss)) add_fieldss <- list(add_fieldss)

  values <- tibble::tibble(
    file = files,
    feature_type = feature_types,
    add_fields = add_fieldss,
    target_name_suffix = target_name_suffixes
  ) |>
    dplyr::mutate(
      file_target = rlang::syms(paste0("fmetadata_file_", target_name_suffix))
    )

  ## based on static branching
  targets <- tarchetypes::tar_map(
    values = values,
    names = "target_name_suffix",
    targets::tar_target_raw(
      "fmetadata_file",
      as.symbol("file"),
      format = "file"
    ),
    targets::tar_target_raw(
      "fmetadata",
      substitute(
        import_fmetadata_gff(
          file_target,
          feature_type = feature_type,
          add_fields = add_fields
        )
      )
    )
  )

  return(targets)
}

#' Target factory for csv samples metadata import
#'
#' Creates a list of targets that track each file and import the samples
#' metadata from each csv file.
#'
#' @param files Character vector, vector of paths to the samples metadata csv
#'   files.
#' @param col_ids Character vector, the name of the column in each file that
#'   contains the ID of the rows (i.e. the sample IDs).
#' @param target_name_suffixes Character vector, a suffix to add to the name of
#'   the targets created by this target factory for each dataset.
#' @returns A list of target objects. For example, with two files to import and
#'   `target_name_suffixes = c("geno", "transcripto")`, the factory returns the
#'   following targets:
#' * `smetadata_file_geno` and `smetadata_file_transcripto`: targets tracking
#'   the genomics and transcriptomics samples metadata files, respectively.
#' * `smetadata_geno` and `smetadata_transcripto`: targets that import the
#'   genomics and transcriptomics samples metadata datasets, respectively.
#' @examples
#' \dontrun{
#' ## in the _targets.R
#' library(MOITestPipeline)
#'
#' list(
#'   import_smetadata_csv_factory(
#'     c(
#'       "data/genotype_smetadata.csv",
#'       "data/rnaseq_smetadata.csv"
#'     ),
#'     col_ids = c("Sample", "SampleIDs"),
#'     target_name_suffixes = c("geno", "transcripto")
#'   )
#' )
#' }
#' @export
import_smetadata_csv_factory <- function(files,
                                         col_ids,
                                         target_name_suffixes) {
  ## for devtools::check
  target_name_suffix <- NULL

  values <- tibble::tibble(
    file = files,
    col_id = col_ids,
    target_name_suffix = target_name_suffixes
  ) |>
    dplyr::mutate(
      file_target = rlang::syms(paste0("smetadata_file_", target_name_suffix))
    )

  ## based on static branching
  targets <- tarchetypes::tar_map(
    values = values,
    names = "target_name_suffix",
    targets::tar_target_raw(
      "smetadata_file",
      as.symbol("file"),
      format = "file"
    ),
    targets::tar_target_raw(
      "smetadata",
      substitute(
        import_smetadata_csv(
          file_target,
          col_id = col_id
        )
      )
    )
  )

  return(targets)
}

#' Replace non-ASCII characters in vector using database
#'
#' Uses an internal database (lookup table) to replace non-ASCII characters
#' in a character vector.
#'
#' @param x Character vector.
#' @returns The character vector `x` in which non-ASCII characters have been
#' replaced.
#'
#' @noRd
.replace_non_ascii_database <- function(x) {
  non_ascii_indx <- which(stringr::str_detect(x, "[^[:ascii:]]"))

  if (length(non_ascii_indx) > 0) {
    x_new <- x[non_ascii_indx]
    for (i in seq_len(nrow(unicode_replacement))) {
      x_new <- stringr::str_replace_all(
        x_new,
        unicode_replacement$Unicode[i],
        unicode_replacement$Replacement[i]
      )
    }

    x[non_ascii_indx] <- x_new
    return(x)
  } else {
    return(x)
  }
}

#' Replace non-ASCII characters in data-frame
#'
#' Uses the stringi and textclean packages as well as an internal database
#' (lookup table) to replace non-ASCII characters in a data-frame.
#'
#' @param df A data-frame.
#' @returns The data-frame `df` in which non-ASCII characters have been
#'   replaced.
#'
#' @noRd
.replace_non_ascii_chars <- function(df) {

  if (!requireNamespace("stringi", quietly = TRUE)) {
    stop(
      "Package \"stringi\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("textclean", quietly = TRUE)) {
    stop(
      "Package \"textclean\" must be installed to use this function.",
      call. = FALSE
    )
  }

  res <- lapply(df, function(x) {
    if (is.character(x)) {
      ## Detect if any non-ascii characters
      non_ascii_indx <- which(stringr::str_detect(x, "[^[:ascii:]]"))

      if (length(non_ascii_indx)) {
        ## Step 1: try to replace non-ASCII characters using the stringi package
        x_new <- stringi::stri_trans_general(x, "latin-ascii")

        ## Step 2: try to replace remaining non-ASCII characters using our
        ## custom translation database
        x_new <- .replace_non_ascii_database(x_new)

        ## Step 3: replace all remaining non-replaceable non-ASCII characters
        ## with a "?"
        x_new <- textclean::replace_non_ascii(x_new, replacement = "?")

        return(x_new)
      } else {
        return(x)
      }
    } else {
      return(x)
    }
  })

  res <- as.data.frame(res)
  rownames(res) <- rownames(df)

  return(res)
}

