#' Evaluate feature selection against features label
#'
#' Compares the selection of features with different feature labels (e.g.
#' result of a DE analysis) for each latent dimension.
#'
#' @param method_output Integration method output generated via the
#' `get_output()` function.
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param col_names Named character vector, giving for each dataset the
#' name of the column in the features metadata table that contains the features
#' label. If a dataset is not present in this vector, will be excluded from the
#' resulting table.
#' @param latent_dimensions Character vector, the latent dimensions to include
#' in the resulting table. If `NULL` (default value), all latent dimensions will
#' be represented.
#' @returns A tibble, with for each dataset and latent dimension the number of
#' selected and non-selected features per feature label.
#' @export
evaluate_feature_selection_table <- function(method_output,
                                             mo_data,
                                             col_names,
                                             latent_dimensions = NULL) {

  ## for devtools::check
  dataset <- latent_dimension <- feature_label <- not_selected <- NULL
  method <- selection_status_intern <- n <- selected <- NULL

  if (!inherits(method_output, "output_dimension_reduction")) {
    stop("'method_output': expecting an 'output_dimension_reduction' object (from get_output()).")
  }

  ## Checking input MultiDataSet object
  datasets <- levels(method_output$features_weight$dataset)
  mo_data <- check_input_multidataset(mo_data, datasets)

  ## Checking list for vector labels
  col_names <- .make_var_list(col_names, datasets, fixed_length = 1, allow_subset = TRUE)
  .check_input_var_fmetadata(col_names, mo_data)

  datasets <- intersect(names(col_names), datasets) ## making sure we retain the correct dataset order

  ## Extracting features metadata
  fmetadata <- get_features_metadata(mo_data)[datasets]

  features_label <- fmetadata |>
    purrr::imap_dfr(
      ~ .x |>
        tibble::as_tibble() |>
        dplyr::select(feature_id, feature_label = !!sym(col_names[[.y]])) |>
        dplyr::mutate(feature_label = as.character(feature_label)),
      .id = "dataset"
    )

  method_output |>
    .filter_output_datasets(datasets) |>
    .filter_output_dimensions(latent_dimensions) |>
    purrr::pluck("features_weight") |>
    dplyr::left_join(features_label, by = c("dataset", "feature_id")) |>
    dplyr::mutate(
      dataset = factor(dataset, levels = datasets),
      selection_status_intern = dplyr::case_when(
        weight == 0 ~ "not_selected",
        TRUE ~ "selected"
      )
    ) |>
    dplyr::group_by(dataset, selection_status_intern) |>
    dplyr::count(latent_dimension, feature_label) |>
    tidyr::complete(
      tidyr::nesting(dataset, latent_dimension), feature_label,
      fill = list(n = 0)
    ) |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(
      names_from = selection_status_intern,
      values_from = n
    ) |>
    dplyr::arrange(latent_dimension, dataset) |>
    dplyr::mutate(method = attr(method_output, "method")) |>
    dplyr::select(method, latent_dimension, dataset, feature_label, selected, not_selected)
}

#' Makes list of feature sets from data-frame
#'
#' Creates a list of feature sets from an annotation data-frame.
#'
#' @param annotation_df A data-frame of feature annotation in long format, with
#' at least a column of feature ID and a column giving the set to which the
#' feature belongs. If a feature belongs to more than one set, there should be a
#' row for each of these sets.
#' @param col_id Character, name of the column in `annotation_df` data-frame that
#' contains the features ID.
#' @param col_set Character, name of the column `annotation_df` data-frame that
#' contains the sets ID.
#' @returns a named list, where each element corresponds to a set, and contains
#' a vector of features ID of all features belonging to that set.
#' @export
make_feature_sets_from_df <- function(annotation_df, col_id, col_set) {
  ## for devtools::check
  set_id <- data <- NULL

  ## Checking input
  .check_names(
    col_id,
    colnames(annotation_df),
    "'col_id' argument: '_W_' is not a column in 'annotation_df'. Possible values are: '_C_'."
  )

  .check_names(
    col_set,
    colnames(annotation_df),
    "'col_set' argument: '_W_' is not a column in 'annotation_df'. Possible values are: '_C_'."
  )


  ## turn data-frame into list of sets
  res <- annotation_df |>
    dplyr::select(
      feature_id = !!sym(col_id),
      set_id = !!sym(col_set)
    ) |>
    dplyr::group_by(set_id) |>
    tidyr::nest() |>
    dplyr::mutate(
      data = purrr::map(
        data,
        ~ .x |>
          dplyr::pull(feature_id) |>
          unique()
      )) |>
    tibble::deframe()

  return(res)
}

#' Makes list of feature sets from features metadata
#'
#' Creates a list of feature sets from features metadata.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param col_names Named list of character, with one element per dataset for which
#' feature sets should be generated. The name of each element should correspond
#' to the name of a dataset, and the value should be a column name from
#' the feature metadata of the corresponding dataset to use as set ID.
#' @param combine_omics_sets Logical, can sets contain features from different
#' omics datasets? If `FALSE` (the default), feature sets are created for each
#' omics separately. If two sets from different omics have the same ID, they will
#' be made unique.
#' @returns a named list, where each element corresponds to a set, and contains
#' a vector of features ID of all features belonging to that set.
#' @export
make_feature_sets_from_fm <- function(mo_data, col_names, combine_omics_sets = FALSE) {
  ## for devtools::check
  set_id <- data <- dataset <- n <- set_id_unique <- NULL

  ## Checking input
  if (!("MultiDataSet" %in% class(mo_data))) stop("'mo_data' argument: Expecting MultiDataSet object", call. = FALSE)

  col_names <- .make_var_list(col_names, names(mo_data), fixed_length = 1, allow_subset = TRUE)
  .check_input_var_fmetadata(col_names, mo_data)

  fmeta <- get_features_metadata(mo_data)

  res <- col_names |>
    purrr::imap_dfr(
      ~ fmeta[[.y]] |>
        tibble::as_tibble() |>
        dplyr::mutate(
          set_id = !!sym(col_names[[.y]]),
          set_id = as.character(set_id)
        ) |>
        dplyr::select(feature_id, set_id),
      .id = "dataset"
    )

  ## if the sets are not combined between the omics datasets, need to make
  ## sure that each set name is unique
  if (!combine_omics_sets) {
    new_set_ids <- res |>
      dplyr::select(dataset, set_id) |>
      dplyr::distinct() |>
      dplyr::group_by(set_id) |>
      dplyr::mutate(n = dplyr::n()) |>
      dplyr::ungroup() |>
      dplyr::mutate(set_id_unique = dplyr::case_when(
        n > 1 ~ paste0(set_id, " - ", dataset),
        TRUE ~ set_id
      )) |>
      dplyr::select(-n)

    res <- res |>
      dplyr::left_join(new_set_ids, by = c("dataset", "set_id")) |>
      dplyr::mutate(set_id = set_id_unique) |>
      dplyr::select(-set_id_unique)
  }

  ## Turn dataframe into list of sets
  res <- res |>
    dplyr::select(-dataset) |>
    dplyr::group_by(set_id) |>
    tidyr::nest() |>
    dplyr::mutate(
      data = purrr::map(
        data,
        ~ .x |>
          dplyr::pull(feature_id) |>
          unique()
      )
    ) |>
    tibble::deframe()

  return(res)
}

#' Reduce feature sets to match multi-omics dataset
#'
#' Removes from a list of feature sets all features that are not present in the
#' multi-omics dataset.
#'
#' @param feature_sets Named list, where each element corresponds to a feature set,
#' and contains a vector of features ID of all features belonging to that set.
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param datasets Character vector, the names of the datasets for which features
#' assignment should be checked. By default, all datasets in `mo_data` are considered.
#' @returns A feature sets list, i.e. named list where each element corresponds to
#' a feature set, containing the ID of all features that belong to this set and
#' that are present in the multi-omics dataset.
#' @export
reduce_feature_sets_data <- function(feature_sets, mo_data, datasets = names(mo_data)) {
  ## for devtools::check
  feature_id <- NULL

  mo_data <- check_input_multidataset(mo_data, datasets)

  ## Getting the list of all features in the multi-omics dataset
  features_list <- get_features(mo_data)[datasets] |>
    unlist() |>
    unname()

  ## Getting the list of all features in the feature sets
  annotated_features <- feature_sets |>
    unlist() |>
    unname() |>
    unique()

  missing_features <- setdiff(annotated_features, features_list)
  n_missing_features <- length(missing_features)

  if (n_missing_features > 0) {
    message(
      format(n_missing_features, big.mark = ","),
      " features in sets but not in the multi-omics dataset."
    )

    ## Restricting the feature sets to only those features that are in the
    ## multi-omics dataset
    res <- feature_sets |>
      purrr::map(
        ~ .x[.x %in% features_list]
      )
    ## Removing sets with no remaining features
    res <- res[purrr::map_int(res, length) > 0]

    empty_sets <- setdiff(names(feature_sets), names(res))
    message("Removing ", format(length(empty_sets), big.mark = ","), " empty feature sets.")
  } else {
    message("All features in sets are in the multi-omics dataset.")
    res <- feature_sets
  }

  return(res)
}

#' Checks features assignment to sets
#'
#' Checks proportion of features from a multi-omics dataset that are assigned to
#' feature sets.
#'
#' @param feature_sets Named list, where each element corresponds to a feature set,
#' and contains a vector of features ID of all features belonging to that set.
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param datasets Character vector, the names of the datasets for which features
#' assignment should be checked. By default, all datasets in `mo_data` are considered.
#' @returns A tibble giving for each dataset the number and fraction of features
#' that are assigned to at least one feature set. The `message` column is meant to
#' facilitate reporting.
#' @export
check_feature_sets <- function(feature_sets, mo_data, datasets = names(mo_data)) {
  ## for devtools::check
  dataset <- annotated <- n_annotated <- n <- feature_id <- NULL
  frac_annotated <- perc_annotated <- NULL

  mo_data <- check_input_multidataset(mo_data, datasets)

  annotated_features <- unlist(feature_sets) |>
    unique()

  res <- get_features_metadata(mo_data)[datasets] |>
    purrr::map_dfr(
      ~ .x |>
        tibble::as_tibble() |>
        dplyr::select(feature_id),
      .id = "dataset"
    ) |>
    dplyr::mutate(annotated = feature_id %in% annotated_features) |>
    dplyr::group_by(dataset) |>
    dplyr::summarise(
      n_annotated = sum(annotated),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      frac_annotated = n_annotated / n,
      perc_annotated = round(100 * frac_annotated, 1),
      message = paste0(
        dataset,
        " dataset: ",
        format(n_annotated, big.mark = ",", trim = TRUE),
        " of ",
        format(n, big.mark = ",", trim = TRUE),
        " (",
        perc_annotated,
        "%) features assigned to at least one set."
      ),
      dataset = factor(dataset, levels = datasets)
    ) |>
    dplyr::arrange(dataset) |>
    dplyr::select(-perc_annotated)

  res |>
    dplyr::pull(message) |>
    paste(collapse = "\n") |>
    message()

  return(res)
}

#' Enrichment analysis for integration results
#'
#' Performs an enrichment analysis for each latent dimension in an integration
#' result, based on user-defined feature sets. The enrichment analysis is done
#' with the [gage::gage()] function from the [`gage`](https://bioconductor.org/packages/gage/)
#' package.
#'
#' When `add_missing_features` is `TRUE` (default behaviour) and a MultiDataSet object
#' is passed through the `mo_data` argument, features present in the multi-omics dataset
#' but absent in the integration method's results will be added
#' to the method's result with a weight of 0. This make sure that if, from a set
#' of 30 features, 25 of these features were removed during the feature pre-selection
#' stage, the enrichment considers that these 25 features were not given high
#' weights by the method. Otherwise, if `add_missing_features` is `FALSE`, these
#' 25 features will be ignored, and so the enrichment analysis may find that one
#' latent dimension is enriched for this particular set, even though there only are
#' 5 features out of 30 from the set that contribute to the latent dimension.
#' Also note that multiple-testing correction is applied at the latent dimension
#' level, and there is no correction across the latent dimensions.
#'
#' When setting `use_abs` to `FALSE`, for each latent dimension, their enrichment
#' for the features test is tested twice: once for enrichment in features with
#' positive weight/importance, and once for features with negative weight/importance
#' score. This will be indicated in the `direction` column of the resulting tibble.
#'
#' Note that we built this function using the gage vignette on
#' [RNA-Seq Data Pathway and Gene-set Analysis Workflow](https://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf), section 7.1.
#'
#'
#' @param method_output Integration method output generated via the `get_output()` function.
#' @param feature_sets Named list, where each element corresponds to a feature set,
#' and contains a vector of features ID of all features belonging to that set.
#' @param datasets Character vector, the names of the datasets to consider in the
#' enrichment analysis. If `NULL` (default value), features from all datasets
#' will be included in the analysis.
#' @param latent_dimensions Character vector, the latent dimensions for which an enrichment
#' analysis should be performed. If `NULL` (default value), all latent dimensions will
#' be analysed.
#' @param use_abs Logical, whether to use the absolute value of the features metric to
#' perform the enrichment. If `TRUE` (default value), it allows to higlight feature sets
#' in which the features have high weight/importance score, both positive and negative.
#' If `FALSE`, it will instead highlight feature sets in which the weights or importance
#' scores all have the same sign (coordinated change).
#' @param rank_test Logical, whether a non-parametric Wilcoxon Mann-Whitney test
#' should be used instead of the default two-sample t-test (i.e. based on features
#' rank rather than their metric). Default value is `FALSE`.
#' @param min_set_size Integer, the minimum number of features in a set required
#' in order to compute an enrichment score for the set. Default value is 5.
#' @param add_missing_features Logical, whether features that are in a multi-omics dataset
#' (provided through the `mo_data` argument) but don't have a weight in the integration results
#' (e.g. because they were not selected in the pre-processing step) should be added in the
#' results. If `TRUE` (default value), they will be added with an importance or weight of 0.
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object. If
#' `add_missing_features` is true, all features in the multi-omics dataset with no
#' weight in the integration method result will be added with a weight of 0.
#' @param sets_info_df Data-frame, information about the feature sets that will
#' be added to the enrichment results. If `NULL` (default value), no information
#' will be added to the results.
#' @param col_set Character, name of the column in `sets_info_df` containing the
#' set IDs. Should match the names of the `feature_sets` list.
#' @returns a tibble of enrichment results.
#' @export
evaluate_method_enrichment <- function(method_output,
                                       feature_sets,
                                       datasets = NULL,
                                       latent_dimensions = NULL,
                                       use_abs = TRUE,
                                       rank_test = FALSE,
                                       min_set_size = 5,
                                       add_missing_features = FALSE,
                                       mo_data = NULL,
                                       sets_info_df = NULL,
                                       col_set = NULL) {


  ## for devtools::check
  latent_dimension <- data <- weight <- direction <- stat.mean <- pvalue <- NULL
  p.val <- q.val <- p.geomean <- set.size <- exp1 <- adj_pvalue <- NULL

  if (!rlang::is_installed("gage")) {
    stop(
      "Package \"gage\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!inherits(method_output, "output_dimension_reduction")) {
    stop("'method_output': expecting an 'output_dimension_reduction' object (from get_output()).")
  }

  method_output <- method_output |>
    .filter_output_dimensions(latent_dimensions) |>
    .filter_output_datasets(datasets)

  latent_dimensions <- levels(method_output$features_weight$latent_dimension)
  datasets <- levels(method_output$features_weight$dataset)

  if (!is.null(sets_info_df)) {
    if (is.null(col_set)) {
      stop("'col_set' argument should be the name of the column in 'sets_info_df' containing the feature sets ID.")
    }
    .check_names(
      col_set, colnames(sets_info_df),
      "'col_set' argument: the value(s) '_W_' are not column names in 'sets_info_df'. Possible values are: '_C_'."
    )
  }

  ## Extracting for each latent dimension the vector of features weight
  features_weight <- method_output$features_weight |>
    dplyr::group_by(latent_dimension) |>
    tidyr::nest() |>
    dplyr::mutate(data = purrr::map(
      data,
      ~ dplyr::select(.x, feature_id, weight) |>
        tibble::deframe()
    )) |>
    tibble::deframe()


  ## If necessary, add to the feature weight vectors the features that are present
  ## in the annotation but not in the methods' result (probably they have been
  ## filtered out before running the method). Will be added with a weight of 0.
  if (add_missing_features) {
    if (is.null(mo_data)) {
      stop("'add_missing_features' is TRUE, expecting a MultiDataSet object to be passed through the 'mo_data' argument.")
    }

    mo_data <- check_input_multidataset(mo_data, datasets)

    all_features <- get_features(mo_data) |>
      unlist() |>
      unname()

    missing_features <- features_weight |>
      purrr::map(~ setdiff(all_features, names(.x)))

    features_weight <- purrr::map2(
      features_weight,
      missing_features,
      ~ c(
        .x,
        rep_len(0, length(.y)) |>
          rlang::set_names(.y)
      )
    )
  }

  max_set_size <- purrr::map_int(feature_sets, length) |>
    max()

  ## Running the enrichment with GAGE
  res <- features_weight |>
    purrr::map(
      ~ gage::gage(
        .x,
        gsets = feature_sets,
        same.dir = !use_abs,
        set.size = c(min_set_size, max_set_size + 1),
        rank.test = rank_test
      )
    ) |>
    purrr::map(
      ~ .x[setdiff(names(.x), "stats")] |>
        purrr::map(tibble::as_tibble, rownames = "set_id") |>
        purrr::list_rbind(names_to = "direction")
    ) |>
    purrr::list_rbind(names_to = "latent_dimension") |>
    dplyr::rename(
      stat_mean = `stat.mean`,
      pvalue = `p.val`,
      adj_pvalue = `q.val`,
      set_size = `set.size`
    ) |>
    dplyr::mutate(
      latent_dimension = factor(latent_dimension,levels = latent_dimensions)
    ) |>
    dplyr::arrange(adj_pvalue, pvalue) |>
    dplyr::select(-`p.geomean`, -`exp1`)

  if (use_abs) {
    res <- dplyr::select(res, -direction)
  } else {
    res <- res |>
      dplyr::mutate(
        direction = dplyr::case_when(
          direction == "greater" ~ "up-regulation test",
          direction == "less" ~ "down-regulation test"
        )
      )
  }

  if (!is.null(sets_info_df)) {
    res <- res |>
      dplyr::left_join(
        sets_info_df,
        by = c("set_id" = col_set)
      )
  }

  return(res)
}


#' Plots features weight in/not in a set
#'
#' Plots the distribution of features weight from an integration method,
#' depending on whether the features belong to a feature set of interest.
#'
#' @param method_output Integration method output generated via the `get_output()` function.
#' @param feature_set Character vector, features ID belonging to the features
#' set of interest.
#' @param set_name Character, name of the set. Default value is `'set'`.
#' @param features_metric Character, the features metric that should be plotted on
#' the y-axis. Should be one of `'signed_importance'` (default value), `'weight'` or
#' `'importance'`.
#' @param add_missing_features Logical, whether features that are in a multi-omics dataset
#' (provided through the `mo_data` argument) but don't have a weight in the integration results
#' (e.g. because they were not selected in the pre-processing step) should be added in the
#' results. If `TRUE` (default value), they will be added with a weight and importance of 0.
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object. If
#' `add_missing_features` is true, all features in the multi-omics dataset with no
#' weight in the integration method result will be added with a weight and importance of 0.
#' @param datasets Character vector, name of the datasets for which the features
#' importance should be plotted. If `NULL` (default value), all datasets will be
#' considered.
#' @param latent_dimensions Character vector, the latent dimensions to represent
#' in the plot. If `NULL` (default value), all latent dimensions will be represented.
#' @param point_alpha Numeric between 0 and 1, the opacity of the points
#' in the plot (with 1 = fully opaque, and 0 = fully transparent). Default value
#' is `0.5`.
#' @param add_boxplot Logical, should a boxplot be drawn on top of the points for
#' categorical covariates? Default value is `TRUE`.
#' @param scales Character, value to use for the `scales` argument of [ggplot2::facet_grid()].
#' Default value is `'free_x'`.
#' @returns a ggplot.
#' @export
plot_features_weight_set <- function(method_output,
                                     feature_set,
                                     set_name = "set",
                                     features_metric = c("signed_importance", "weight", "importance"),
                                     add_missing_features = FALSE,
                                     mo_data = NULL,
                                     datasets = NULL,
                                     latent_dimensions = NULL,
                                     point_alpha = 0.5,
                                     add_boxplot = TRUE,
                                     scales = "free_x") {
  ## for devtools::check
  weight <- importance <- feature_id <- NULL

  if (!inherits(method_output, "output_dimension_reduction")) stop("'method_output': expecting an 'output_dimension_reduction' object (from get_output()).")

  features_metric <- rlang::arg_match(features_metric)
  y_label <- paste("Features", stringr::str_replace(features_metric, "_", " "))

  method_output <- method_output |>
    .filter_output_dimensions(latent_dimensions) |>
    .filter_output_datasets(datasets)

  x <- method_output$features_weight

  if (add_missing_features) {
    mo_data <- check_input_multidataset(mo_data, levels(x$dataset))

    missing_features <- get_features(mo_data) |>
      purrr::map_dfr(
        ~ tibble::tibble(feature_id = .x),
        .id = "dataset"
      ) |>
      dplyr::filter(!(feature_id %in% unique(x$feature_id)))

    x <- dplyr::bind_rows(
      x,
      tidyr::expand_grid(
        missing_features,
        latent_dimension = levels(x$latent_dimension)
      ) |>
        dplyr::mutate(
          weight = 0,
          importance = 0
        )
    )
  }

  ## Checking which features have weight
  features_with_weight <- unique(x$feature_id)
  n_missing_features <- setdiff(
    feature_set,
    features_with_weight
  ) |>
    length()

  toplot <- x |>
    dplyr::mutate(
      signed_importance = sign(weight) * importance,
      in_set = feature_id %in% feature_set,
      in_set = dplyr::case_when(
        in_set ~ "In set",
        !in_set ~ "Not in set"
      )
    )

  plot_x_wrapper(
    toplot,
    x = "in_set",
    y = features_metric,
    facet_grid = c("latent_dimension", "dataset"),
    colour = NULL,
    shape = NULL,
    point_alpha = point_alpha,
    add_boxplot = add_boxplot,
    scales_facet = scales
  ) +
    ggplot2::labs(
      title = paste(y_label, "-", attr(method_output, "method")),
      x = paste0(
        set_name,
        " - ",
        format(length(feature_set), big.mark = ","),
        " features in set, ",
        format(length(features_with_weight), big.mark = ","),
        " features total with a score."
      ),
      y = y_label,
      caption = dplyr::if_else(
        n_missing_features > 0,
        paste0(n_missing_features, " features from the set with no weight/importance score."),
        "All features from the set with weight and importance."
      )
    )
}

# plot_enrichment_results <- function(enrichment_results,
#                                     latent_dimensions = NULL,
#                                     sign_thr = 0.05,
#                                     top_n = 10){
#
#   ## Checking latent dimensions input
#   all_latent_dimensions <- enrichment_results |>
#     dplyr::pull(latent_dimension) |>
#     unique() |>
#     as.character()
#
#   if(!is.null(latent_dimensions)){
#     .check_names(latent_dimensions,
#                  all_latent_dimensions,
#                  "'latent_dimensions' argument: the following values do not correspond to latent dimension labels: '_W_'. Possible values are: '_C_'.")
#   } else {
#     latent_dimensions <- all_latent_dimensions
#   }
#
#   n_sets_to_plot <- enrichment_results |>
#     dplyr::filter(latent_dimension %in% latent_dimensions) |>
#     dplyr::group_by(latent_dimension) |>
#     dplyr::summarise(n_sign = sum(q.val <= sign_thr, na.rm = TRUE)) |>
#     dplyr::mutate(n_sign = dplyr::case_when(n_sign > top_n ~ n_sign,
#                                             TRUE ~ as.integer(top_n))) |>
#     tibble::deframe()
#
#
#   enrichment_results |>
#     dplyr::filter(latent_dimension %in% latent_dimensions) |>
#     dplyr::group_by(latent_dimension) |>
#     tidyr::nest() |>
#     dplyr::mutate(plot = purrr::map2(
#       data,
#       latent_dimension,
#       ~ .x |>
#         dplyr::slice_min(q.val, n = n_sets_to_plot[[.y]]) |>
#         dplyr::mutate(score = -log10(q.val)) |>
#         dplyr::arrange(score) |>
#         dplyr::mutate(set_id = )
#         ggplot2::ggplot(aes(x = score, y))))
#
#
# }


#' Computes samples silhouette score from method output
#'
#' Calculates the samples silhouette width from the results of a dimension
#' reduction method, according to some samples grouping from the samples
#' metadata of a `MultiDataSet` object.
#'
#' The samples silhouette width and groups average width are calculated using
#' the [cluster::silhouette()] function.
#'
#' @param method_output method_output Integration method output generated via the
#' `get_output()` function.
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param group_col Character, name of the column in one of the samples metadata
#' table from `mo_data` containing the samples grouping to be used.
#' @param latent_dimensions Character vector, the latent dimensions to use for
#' computing distance between samples. If `NULL` (default value), all latent
#' dimensions will be used.
#' @param distance_metric Character, name of the metric to use when computing
#' the distance between samples from their coordinates in the latent dimensions.
#' This will be passed to the [stats::dist()] function. Options include: `"euclidean"`
#' (default value), `"maximum"`, `"manhattan"`, `"canberra"`, `"binary"` or `"minkowski"`.
#' @returns A list with two elements:
#' * `samples_silhouette`: a tibble giving for each sample (`sample_id` column)
#' the group to which it belongs (`group` column), its closest (other) group in the
#' space spanned by the latent dimensions (`neighbour_group`), and its silhouette
#' width (`silhouette_width` column).
#' * `groups_average_silhouette`: a tibble giving for each samples group (`group`
#' column) its average silhouette width (`group_average_width` column).
#' @export
compute_samples_silhouette <- function(method_output,
                                       mo_data,
                                       group_col,
                                       latent_dimensions = NULL,
                                       distance_metric = c("euclidean",
                                                           "maximum",
                                                           "manhattan",
                                                           "canberra",
                                                           "binary",
                                                           "minkowski")) {

  if (!rlang::is_installed("cluster")) {
    stop(
      "Package \"cluster\" must be installed to use this function.",
      call. = FALSE
    )
  }

  if (!inherits(method_output, "output_dimension_reduction")) {
    stop("'method_output': expecting an 'output_dimension_reduction' object (from get_output()).")
  }

  ## For devtools::check
  latent_dimension <- score <- cluster <- neighbor <- sil_width <- group <- NULL
  sample_id <- NULL

  distance_metric <- rlang::arg_match(distance_metric)
  .check_input_var_smetadata_common(group_col, mo_data)

  ## Computing distance matrix
  dist_mat <- method_output |>
    .filter_output_dimensions(latent_dimensions) |>
    purrr::pluck("samples_score") |>
    tidyr::pivot_wider(
      names_from = latent_dimension,
      values_from = score
    ) |>
    tibble::column_to_rownames("sample_id") |>
    as.matrix() |>
    stats::dist(method = distance_metric)

  ## Extracting samples group
  smeta <- get_samples_metadata_combined(mo_data)[attr(dist_mat, "Labels"), , drop = FALSE]

  cluster_vec <- smeta[[group_col]] |>
    rlang::set_names(rownames(smeta))

  if (is.double(cluster_vec)) {
    warning("'group_col' argument: should indicate a column of samples grouping (character, factor or integer), not a numeric column. Converting column to integer.", call. = FALSE)

    ## as.integer() looses vector names :(
    cluster_vec <- as.integer(cluster_vec) |>
      rlang::set_names(names(cluster_vec))
  }

  ## To keep levels ordering if the column was already a factor
  if (!is.factor(cluster_vec)) {
    already_fac <- TRUE
    cluster_vec <- factor(cluster_vec)
  }

  grp_levels <- levels(cluster_vec)

  ## Calculate silhouette
  sil_res <- cluster::silhouette(
    as.integer(cluster_vec),
    dist_mat
  )

  if (all(is.na(sil_res))) {
    stop(
      "Silhouette width calculation failed. Check that there is more than one sample in at least one group.",
      call. = FALSE
    )
  }

  ## Extract individual silhouette widths
  sil_df <- sil_res |>
    tibble::as_tibble() |>
    dplyr::mutate(
      sample_id = names(cluster_vec),
      cluster = grp_levels[cluster],
      neighbor = grp_levels[neighbor]
    ) |>
    dplyr::select(
      sample_id,
      group = cluster,
      neighbour_group = neighbor,
      silhouette_width = sil_width
    )

  ## Extract group average widths
  sil_sum_df <- summary(sil_res)$clus.avg.widths |>
    rlang::set_names(grp_levels) |>
    tibble::enframe(
      name = "group",
      value = "group_average_width"
    )

  ## Retaining level ordering for factor groups
  if (already_fac) {
    sil_df <- sil_df |>
      dplyr::mutate(
        dplyr::across(
          tidyselect::contains("group"),
          ~ factor(.x, levels = grp_levels)
        )
      )

    sil_sum_df <- sil_sum_df |>
      dplyr::mutate(
        group = factor(group, levels = grp_levels)
      )
  }

  list(
    samples_silhouette = sil_df,
    groups_average_silhouette = sil_sum_df
  )
}
