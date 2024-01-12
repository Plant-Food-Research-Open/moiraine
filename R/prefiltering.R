.check_tokeep_value <- function(to_keep_n = NULL, to_keep_prop = NULL, dataset_size) {
  if (.check_missing(to_keep_n)) {
    if (.check_missing(to_keep_prop)) stop("Need to give a value for either 'to_keep_n' or 'to_keep_prop'.")

    if ((to_keep_prop <= 0) | (to_keep_prop >= 1)) stop("'to_keep_prop' argument should be < 0 and > 1.")

    to_keep <- to_keep_prop * dataset_size
  } else {
    if (to_keep_n >= dataset_size) stop("'to_keep_n' argument: number of features to retain equals or exceeds the number of features in the dataset.")
    to_keep <- to_keep_n
  }

  to_keep <- round(to_keep, 0)

  return(to_keep)
}

#' Select features based on Median Absolute Deviation from MultiDataSet
#'
#' Computes the Median Absolute Deviation (MAD) for each feature in an omics dataset
#' from a MultiDataSet object, and select features with the highest MAD values. This
#' is a wrapper function around the [get_dataset_matrix()] and [select_features_mad_matrix()] functions.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param dataset_name Character, name of the omics dataset on which to apply feature pre-selection.
#' @inheritParams select_features_mad_matrix
#' @return A tibble with columns `feature_id`, `mad` and `selected` (logical, indicates whether the feature is selected
#' based on its MAD value). In addition, the name of the dataset filtered is stored in the return object attribute
#' `dataset_name` (which can be accessed via `attr(res, "dataset_name")`).
#' @export
select_features_mad <- function(mo_data,
                                dataset_name,
                                to_keep_n = NULL,
                                to_keep_prop = NULL,
                                with_ties = TRUE) {

  mat <- get_dataset_matrix(mo_data, dataset_name, keep_dataset_name = TRUE)
  res <- select_features_mad_matrix(mat, to_keep_n, to_keep_prop, with_ties)

  return(res)
}

#' Select features based on Median Absolute Deviation from matrix
#'
#' Computes the Median Absolute Deviation (MAD) for each feature in an omics dataset
#' from a MultiDataSet object, and select features with the highest MAD values.
#'
#' @param mat Matrix of omics measurement, with features as rows and samples as columns.
#' @param to_keep_n Integer, the number of features to retain in the dataset. Should be less than the number of
#' features in the dataset. If `NULL` or `NA`, `to_keep_prop` will be used instead.
#' @param to_keep_prop Numeric, the proportion of features to retain in the dataset. Will be ignored if `to_keep_n`
#' is supplied. Value should be > 0 and < 1.
#' @param with_ties Should ties be kept together? If \code{TRUE}, may return more features than requested. Default
#' value is `TRUE`.
#' @return A tibble with columns `feature_id`, `mad` and `selected` (logical, indicates whether the feature is selected
#' based on its MAD value). In addition, the name of the dataset filtered is stored in the return object attribute
#' `dataset_name` (which can be accessed via `attr(res, "dataset_name")`).
#' @export
select_features_mad_matrix <- function(mat,
                                       to_keep_n = NULL,
                                       to_keep_prop = NULL,
                                       with_ties = TRUE) {

  to_keep <- .check_tokeep_value(to_keep_n, to_keep_prop, nrow(mat))

  mad_values <- apply(mat, 1, stats::mad, na.rm = TRUE)

  ## for devtools:check()
  mad <- NULL
  res <- tibble::tibble(
    feature_id = names(mad_values),
    mad = mad_values
  ) |>
    dplyr::arrange(dplyr::desc(mad))

  ## get the features with highest MAD values
  temp <- dplyr::slice_max(res, order_by = mad, n = to_keep, with_ties = with_ties)

  ## create column to say whether a feature was selected or not
  res$selected <- res$feature_id %in% temp$feature_id

  attr(res, "dataset_name") <- attr(mat, "dataset_name")

  return(res)
}


#' Get filtered MultiDataSet object based on variability measure
#'
#' Selects most highly variable features from omics datasets based on features'
#' variability (e.g. MAD or COV).
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param var_list A list with the result from the MAD or COV calculation for each dataset to be filtered, as returned by
#' the \code{\link{select_features_mad}} or \code{\link{select_features_cov}} function.
#' @return A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @examples
#' \dontrun{
#' # Goal: keep 20% of features in dataset1, and 50% of features in dataset2
#' to_keep_prop <- c("dataset1" = 0.2, "dataset_2" = 0.5)
#'
#' # 1) compute MAD values and select features for dataset1 and dataset2
#' mad_list <- lapply(names(to_keep_prop), function(i) {
#'   select_features_mad(mo_data, i, to_keep_prop[i])
#' })
#'
#' # 2) Get the filtered dataset
#' mo_data_filtered <- get_filtered_dataset_variability(mo_data, mad_list)
#' }
#' @export
get_filtered_dataset_variability <- function(mo_data,
                                             var_list) {
  if (!("MultiDataSet" %in% class(mo_data))) stop("Expecting MultiDataSet object")

  ## Make sure that the names of var_list are the name of the datasets
  names(var_list) <- sapply(var_list, attr, "dataset_name")

  ## Get the name of the datasets that need to be filtered and the ones that don't
  datasets_filtered <- names(var_list)
  datasets_non_filtered <- setdiff(names(mo_data), names(var_list))

  ## for devtools::check()
  selected <- feature_id <- NULL

  retained_features <- c(
    ## Get the list of retained features
    lapply(datasets_filtered, function(i) {
      dplyr::filter(var_list[[i]], selected)$feature_id
    }),
    ## Get the list of all features from non-filtered datasets
    lapply(datasets_non_filtered, function(i) {
      Biobase::featureData(mo_data[[i]])$feature_id
    })
  )

  mo_data <- subset_features(mo_data, retained_features)

  return(mo_data)
}

#' Target factory for feature preselection based on Median Absolute Deviation
#'
#' Creates a list of targets to perform feature preselection on datasets from a `MultiDataSet`
#' object by retaining features with the highest Median Absolute Deviation (MAD).
#'
#' @param mo_data_target Symbol, the name of the target containing the `MultiDataSet` object.
#' @param to_keep_ns Named integer vector, the number of feature to retain in each dataset to be prefiltered
#' (names should correspond to a dataset name). Value should be less than the number of features in the
#' corresponding dataset. Set to `NULL` in order to use `to_keep_props` instead.
#' @param to_keep_props Named numeric vector, the proportion of features to retain in each dataset
#' to be prefiltered (names should correspond to a dataset name). Value should be > 0 and < 1.
#' Will be ignored if `to_keep_ns` is not `NULL`.
#' @param with_ties Should ties be kept together? If `TRUE`, may return more features than requested. Default value is `TRUE`.
#' @param target_name_prefix Character, a prefix to add to the name of the targets created by this target factory.
#' Default value is `""`.
#' @param filtered_set_target_name Character, the name of the final target containing the filtered `MultiDataSet` object.
#' If NULL, a name will automatically be supplied. Default value is `NULL`.
#' @return A list of target objects. With `target_name_prefix = ""` and `filtered_set_target_name = NULL`,
#' the following targets are created:
#' * `mad_spec`: a target that generates a grouped tibble where each row corresponds to one dataset to be filtered,
#'   with the columns specifying each dataset name, and associated values from `to_keep_ns`, `to_keep_props` and `with_ties`.
#' * `mad_mat`: a dynamic branching target that run the [get_dataset_matrix()] function for each dataset.
#' * `individual_mad_values`: a dynamic branching target that runs the [select_features_mad_matrix()] function for each dataset.
#' * `filtered_set_mad`: a target to retain from the original `MultiDataSet` object only features selected based on their MAD values.
#' @examples
#' \dontrun{
#' ## in the _targets.R
#' library(moiraine)
#'
#' list(
#'   ## add code here to load the different datasets
#'
#'   ## the following target creates a MultiDataSet object from previously
#'   ## created omics sets (geno_set, trans_set, etc)
#'   tar_target(
#'     mo_set,
#'     create_multiomics_set(geno_set, trans_set, metabo_set, pheno_set)
#'   ),
#'   feature_preselection_mad_factory(
#'     mo_set,
#'     to_keep_ns = c("rnaseq" = 1000, "metabolome" = 500),
#'     filtered_set_target_name = "mo_set_filtered"
#'   ),
#'
#'   ## Another example using to_keep_props
#'   feature_preselection_mad_factory(
#'     mo_set,
#'     to_keep_ns = NULL,
#'     to_keep_props = c("rnaseq" = 0.3, "metabolome" = 0.5),
#'     filtered_set_target_name = "mo_set_filtered"
#'   )
#' )
#' }
#' @export
feature_preselection_mad_factory <- function(mo_data_target,
                                             to_keep_ns,
                                             to_keep_props = NULL,
                                             with_ties = TRUE,
                                             target_name_prefix = "",
                                             filtered_set_target_name = NULL) {
  mad_spec_name <- paste0(target_name_prefix, "mad_spec")
  mad_mat_name <- paste0(target_name_prefix, "mad_mat")
  mad_run_name <- paste0(target_name_prefix, "individual_mad_values")
  if (is.null(filtered_set_target_name)) filtered_set_target_name <- paste0(target_name_prefix, "filtered_set_mad")

  mad_spec_target <- as.symbol(mad_spec_name)
  mad_mat_target <- as.symbol(mad_mat_name)
  mad_run_target <- as.symbol(mad_run_name)

  if (!is.null(to_keep_ns)) {
    dataset_names <- names(to_keep_ns)
    to_keep_ns <- unname(to_keep_ns)
  } else {
    dataset_names <- names(to_keep_props)
    to_keep_props <- unname(to_keep_props)
  }

  if (is.null(dataset_names)) {
    stop("'to_keep_ns' or 'to_keep_props' argument should be named.")
  }

  list(
    ## store the MAD specifications (arguments) as a tibble (one row per dataset to prefilter)
    ## and group it by dataset name so that following targets will be applied to each row in turn
    targets::tar_target_raw(
      mad_spec_name,
      substitute(tibble::tibble(dsn = dataset_names, tkn = to_keep_ns, tkp = to_keep_props, wt = with_ties) |>
                   dplyr::group_by(dsn) |>
                   tar_group()),
      iteration = "group"
    ),

    ## Get datasets as matrices
    targets::tar_target_raw(
      mad_mat_name,
      substitute(get_dataset_matrix(mo_data_target, mad_spec_target$dsn, keep_dataset_name = TRUE)),
      pattern = substitute(map(mad_spec_target)),
      iteration = "list"
    ),

    ## perform the MAD preselection on each dataset
    targets::tar_target_raw(
      mad_run_name,
      substitute(select_features_mad_matrix(
        mad_mat_target,
        to_keep_n = mad_spec_target$tkn,
        to_keep_prop = mad_spec_target$tkp,
        with_ties = mad_spec_target$wt
      )),
      pattern = substitute(map(mad_mat_target, mad_spec_target)),
      iteration = "list"
    ),

    ## Subset the full MultiDataSet object to only retain selected features
    targets::tar_target_raw(
      filtered_set_target_name,
      substitute(get_filtered_dataset_variability(mo_data_target, mad_run_target))
    )
  )
}



#' Select features based on Coefficient of Variation from MultiDataSet
#'
#' Computes the Coefficient of Variation (COV) for each feature in an omics dataset
#' from a MultiDataSet object, and select features with the highest COV values. This
#' is a wrapper function around the [get_dataset_matrix()] and [select_features_cov_matrix()] functions.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param dataset_name Character, name of the omics dataset on which to apply feature pre-selection.
#' @inheritParams select_features_cov_matrix
#' @return A tibble with columns `feature_id`, `cov` and `selected` (logical, indicates whether the feature is selected
#' based on its COV value). In addition, the name of the dataset filtered is stored in the return object attribute
#' `dataset_name` (which can be accessed via `attr(res, "dataset_name")`).
#' @export
select_features_cov <- function(mo_data,
                                dataset_name,
                                to_keep_n = NULL,
                                to_keep_prop = NULL,
                                with_ties = TRUE) {

  mat <- get_dataset_matrix(mo_data, dataset_name, keep_dataset_name = TRUE)
  res <- select_features_cov_matrix(mat, to_keep_n, to_keep_prop, with_ties)

  return(res)
}

#' Select features based on Coefficient of Variation from matrix
#'
#' Computes the Coefficient of Variation (COV) for each feature in an omics dataset
#' from a MultiDataSet object, and select features with the highest COV values.
#'
#' @param mat Matrix of omics measurement, with features as rows and samples as columns.
#' @param to_keep_n Integer, the number of features to retain in the dataset. Should be less than the number of
#' features in the dataset. If `NULL` or `NA`, `to_keep_prop` will be used instead.
#' @param to_keep_prop Numeric, the proportion of features to retain in the dataset. Will be ignored if `to_keep_n`
#' is supplied. Value should be > 0 and < 1.
#' @param with_ties Should ties be kept together? If \code{TRUE}, may return more features than requested. Default
#' value is `TRUE`.
#' @return A tibble with columns `feature_id`, `cov` and `selected` (logical, indicates whether the feature is selected
#' based on its COV value). In addition, the name of the dataset filtered is stored in the return object attribute
#' `dataset_name` (which can be accessed via `attr(res, "dataset_name")`).
#' @export
select_features_cov_matrix <- function(mat,
                                       to_keep_n = NULL,
                                       to_keep_prop = NULL,
                                       with_ties = TRUE) {

  to_keep <- .check_tokeep_value(to_keep_n, to_keep_prop, nrow(mat))

  cov_values <- apply(mat, 1, function(x) {
    stats::sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
  })

  ## for devtools:check()
  cov <- NULL
  res <- tibble::tibble(
    feature_id = names(cov_values),
    cov = cov_values
  ) |>
    dplyr::arrange(dplyr::desc(cov))

  ## get the features with highest COV values
  temp <- dplyr::slice_max(res, order_by = cov, n = to_keep, with_ties = with_ties)

  ## create column to say whether a feature was selected or not
  res$selected <- res$feature_id %in% temp$feature_id

  attr(res, "dataset_name") <- attr(mat, "dataset_name")

  return(res)
}


#' Target factory for feature preselection based on Coefficient of Variation
#'
#' Creates a list of targets to perform feature preselection on datasets from a `MultiDataSet`
#' object by retaining features with the highest Coefficient of Variation (COV).
#'
#' @param mo_data_target Symbol, the name of the target containing the `MultiDataSet` object.
#' @param to_keep_ns Named integer vector, the number of feature to retain in each dataset to be prefiltered
#' (names should correspond to a dataset name). Value should be less than the number of features in the
#' corresponding dataset. Set to `NULL` in order to use `to_keep_props` instead.
#' @param to_keep_props Named numeric vector, the proportion of features to retain in each dataset
#' to be prefiltered (names should correspond to a dataset name). Value should be > 0 and < 1.
#' Will be ignored if `to_keep_ns` is not `NULL`.
#' @param with_ties Should ties be kept together? If `TRUE`, may return more features than requested. Default value is `TRUE`.
#' @param target_name_prefix Character, a prefix to add to the name of the targets created by this target factory.
#' Default value is `""`.
#' @param filtered_set_target_name Character, the name of the final target containing the filtered `MultiDataSet` object.
#' If NULL, a name will automatically be supplied. Default value is `NULL`.
#' @return A list of target objects. With `target_name_prefix = ""` and `filtered_set_target_name = NULL`,
#' the following targets are created:
#' * `cov_spec`: a target that generates a grouped tibble where each row corresponds to one dataset to be filtered,
#'   with the columns specifying each dataset name, and associated values from `to_keep_ns`, `to_keep_props` and `with_ties`.
#' * `cov_mat`: a dynamic branching target that run the [get_dataset_matrix()] function for each dataset.
#' * `individual_cov_values`: a dynamic branching target that runs the [select_features_cov_matrix()] function on each dataset.
#' * `filtered_set_cov`: a target to retain from the original `MultiDataSet` object only features selected based on their COV values.
#' @examples
#' \dontrun{
#' ## in the _targets.R
#' library(moiraine)
#'
#' list(
#'   ## add code here to load the different datasets
#'
#'   ## the following target creates a MultiDataSet object from previously
#'   ## created omics sets (geno_set, trans_set, etc)
#'   tar_target(
#'     mo_set,
#'     create_multiomics_set(geno_set, trans_set, metabo_set, pheno_set)
#'   ),
#'   feature_preselection_cov_factory(
#'     mo_set,
#'     to_keep_ns = c("rnaseq" = 1000, "metabolome" = 500),
#'     filtered_set_target_name = "mo_set_filtered"
#'   ),
#'
#'   ## Another example using to_keep_props
#'   feature_preselection_cov_factory(
#'     mo_set,
#'     to_keep_ns = NULL,
#'     to_keep_props = c("rnaseq" = 0.3, "metabolome" = 0.5),
#'     filtered_set_target_name = "mo_set_filtered"
#'   )
#' )
#' }
#' @export
feature_preselection_cov_factory <- function(mo_data_target,
                                             to_keep_ns,
                                             to_keep_props = NULL,
                                             with_ties = TRUE,
                                             target_name_prefix = "",
                                             filtered_set_target_name = NULL) {
  cov_spec_name <- paste0(target_name_prefix, "cov_spec")
  cov_mat_name <- paste0(target_name_prefix, "cov_mat")
  cov_run_name <- paste0(target_name_prefix, "individual_cov_values")
  if (is.null(filtered_set_target_name)) filtered_set_target_name <- paste0(target_name_prefix, "filtered_set_cov")

  cov_spec_target <- as.symbol(cov_spec_name)
  cov_mat_target <- as.symbol(cov_mat_name)
  cov_run_target <- as.symbol(cov_run_name)

  if (!is.null(to_keep_ns)) {
    dataset_names <- names(to_keep_ns)
    to_keep_ns <- unname(to_keep_ns)
  } else {
    dataset_names <- names(to_keep_props)
    to_keep_props <- unname(to_keep_props)
  }

  if (is.null(dataset_names)) {
    stop("'to_keep_ns' or 'to_keep_props' argument should be named.")
  }

  list(
    ## store the COV specifications (arguments) as a tibble (one row per dataset to prefilter)
    ## and group it by dataset name so that following targets will be applied to each row in turn
    targets::tar_target_raw(
      cov_spec_name,
      substitute(tibble::tibble(dsn = dataset_names, tkn = to_keep_ns, tkp = to_keep_props, wt = with_ties) |>
                   dplyr::group_by(dsn) |>
                   tar_group()),
      iteration = "group"
    ),

    ## Get datasets as matrices
    targets::tar_target_raw(
      cov_mat_name,
      substitute(get_dataset_matrix(mo_data_target, cov_spec_target$dsn, keep_dataset_name = TRUE)),
      pattern = substitute(map(cov_spec_target)),
      iteration = "list"
    ),

    ## perform the COV preselection on each dataset
    targets::tar_target_raw(
      cov_run_name,
      substitute(select_features_cov_matrix(
        cov_mat_target,
        to_keep_n = cov_spec_target$tkn,
        to_keep_prop = cov_spec_target$tkp,
        with_ties = cov_spec_target$wt
      )),
      pattern = substitute(map(cov_mat_target, cov_spec_target)),
      iteration = "list"
    ),

    ## Subset the full MultiDataSet object to only retain selected features
    targets::tar_target_raw(
      filtered_set_target_name,
      substitute(get_filtered_dataset_variability(mo_data_target, cov_run_target))
    )
  )
}


#' Generate sPLS-DA input data (for mixomics)
#'
#' Creates an object that can be used as input for the (s)PLS-DA functions from the mixOmics package.
#' It contains the omics dataset as well as the samples group membership in a list.
#'
#' `multilevel` argument: enables the multilevel option (see
#' [mixOmics site](http://mixomics.org/methods/multilevel/)) to deal with repeated measurements.
#' [mixOmics::splsda()] enables one- and two-factor decomposition. For one-factor decomposition,
#' `multilevel` argument should be the name of the column in the samples metadata that gives the
#' ID of the observation units (e.g. the ID of the subjects that were measured several times). The resulting
#' design matrix (stored in the `multilevel` argument of the returned object) will be a data-frame
#' with one column which gives the ID (as integer) of the observation units corresponding to each sample
#' in the omics datasets. For two-factor decomposition, `multilevel` should be of length 3. The
#' first value, similarly to the one-factor decomposition, should be the name of the column in the
#' samples metadata that gives the ID of the observation units (e.g. the ID of the subjects that were
#' measured several times). The second and third values should be the name of the columns in the samples
#' metadata that give the two factors considered. The resulting design matrix (stored in the `multilevel`
#' argument of the returned object) will be a data-frame with three columns: the first column gives the
#' ID (as integer) of the observation units corresponding to each sample in the omics datasets; the
#' second and third columns give the levels of the two factors.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param dataset_name Character, name of the dataset from `mo_data` to analyse.
#' @param group Character, the column name in the samples information data-frame to use as samples group
#' (use \code{\link{get_samples_metadata}} to view the samples information data-frame for a omics dataset).
#' @param multilevel Character vector of length 1 or 3 to be used as information about repeated measurements.
#' See Details. Default value is `NULL` (no repeated measurements).
#' @returns A list, in which the first element corresponds to the omics dataset, with samples as rows and features as columns, and the second element (named `'Y'`) is a named factor vector, giving for each sample its group.
#' The name of the dataset to be analysed is stored in the `dataset_name` attribute of the returned object.
#' @export
get_input_splsda <- function(mo_data, dataset_name, group, multilevel = NULL) {

  if (length(dataset_name) != 1) {
    stop("'dataset_name' argument should be of length 1.")
  }

  res <- get_input_mixomics_supervised(mo_data, group, dataset_name)
  attr(res, "dataset_name") <- dataset_name

  if (!is.null(multilevel)) attr(res, "multilevel") <- .get_multilevel(res, mo_data, multilevel)

  return(res)
}


#' Assess optimal number of components for sPLS-DA on omics dataset from MultiDataSet object
#'
#' Performs cross-validation for a PLS-DA run (implemented in the `mixOmics` package) on an omics dataset from a
#' `MultiDataSet` object. This allows to estimate the optimal number of latent components to construct.
#' This is intended for feature preselection in the omics dataset (see examples below).
#'
#' This function uses the \code{\link[mixOmics]{plsda}} and  \code{\link[mixOmics]{perf}}
#' function from the `mixOmics` package.
#'
#' @param splsda_input Input for the sPLS-DA functions from mixOmics, created with [get_input_splsda()].
#' @param ncomp_max Integer, the maximum number of latent components to test when estimating the number of
#' latent components to use. Default value is `5`.
#' @param validation Character, which cross-validation method to use, can be one of `"Mfold"` or `"loo"`
#' (see [mixOmics::perf()]). Default value is `"Mfold"`.
#' @param folds Integer, number of folds to use in the M-fold cross-validation (see [mixOmics::perf()]).
#' Default value is 5.
#' @param nrepeat Integer, number of times the cross-validation is repeated (see [mixOmics::perf()]).
#' @param measure Performance measure used to select the optimal value of `ncomp`, can be one of `"BER"` or `"overall"`
#' (see [mixOmics::perf()]).
#' Default value is `"BER"`.
#' @param distance Distance metric used to select the optimal value of `ncomp`, can be one of `"max.dist"`,
#' `"centroids.dist"` or `"mahalanobis.dist"` (see [mixOmics::perf()]). Default value is `"centroids.dist"`.
#' @param cpus Integer, number of cpus to use.
#' @param progressBar Logical, whether to display a progress bar during the optimisation of `ncomp`. Default
#' value is `TRUE`.
#' @return A list as per the output of the [mixOmics::perf()] function, with the following additional elements:
#' \itemize{
#' \item `dataset_name`: the name of the dataset analysed;
#' \item `group`: column name in the samples information data-frame used as samples group;
#' \item `optim_ncomp`: the optimal number of latent components as per the `measure` and `distance` specified;
#' \item `optim_measure`: the measure used to select the optimal number of latent components;
#' \item `optim_distance`: the distance metric used to select the optimal number of latent components.
#' }
#' In addition, the name of the dataset analysed and the column name in the samples information data-frame
#' used as samples group as stored as attributes `dataset_name` and `group`, respectively.
#' @export
perf_splsda <- function(splsda_input,
                        ncomp_max = 5,
                        validation = "Mfold",
                        folds = 5,
                        nrepeat = 50,
                        measure = "BER",
                        distance = "centroids.dist",
                        cpus = 1,
                        progressBar = TRUE) {

  dataset_name <- setdiff(names(splsda_input), "Y")
  multilevel <- attr(splsda_input, "multilevel")

  ## Run the PLS-DA with several latent components
  plsda_res <- mixOmics::plsda(
    splsda_input[[dataset_name]],
    splsda_input$Y,
    ncomp = ncomp_max,
    multilevel = multilevel
  )

  ## Use cross-validation to assess performance
  plsda_res_perf <- mixOmics::perf(
    plsda_res,
    validation = validation,
    folds = folds,
    nrepeat = nrepeat,
    cpus = cpus,
    progressBar = progressBar
  )


  res <- c(
    plsda_res_perf,
    list(
      optim_ncomp = plsda_res_perf$choice.ncomp[measure, distance],
      optim_measure = measure,
      optim_distance = distance
    )
  )

  class(res) <- class(plsda_res_perf)

  attr(res, "dataset_name") <- dataset_name

  return(res)
}


#' Performs sPLS-DA on omics dataset from MultiDataSet object
#'
#' Performs a sPLS-DA (implemented in the `mixOmics`) package on a omics dataset from a
#' MultiDataSet object. This is intended for feature preselection in the omics dataset
#' (see \code{\link{get_filtered_dataset_splsda}}).
#'
#' This function uses the \code{\link[mixOmics]{plsda}} function from the `mixOmics` package.
#' Note that the sPLS-DA method can select the same feature for several latent components, so the number of
#' features retained for a dataset might be less than the number specified in the `to_keep_n` argument.
#'
#' @param splsda_input Input for the sPLS-DA functions from mixOmics, created with [get_input_splsda()].
#' @param perf_res Result of the \code{\link{perf_splsda}} function. If not supplied, sPLS-DA will be run on
#' dataset specified by argument `dataset_name` with number of latent components specified by argument `comp`.
#' @param to_keep_n Integer, the number of features to retain in the dataset. Should be less than the number of
#' features in the dataset. If `NULL` or `NA`, `to_keep_prop` will be used instead.
#' @param to_keep_prop Numeric, the proportion of features to retain in the dataset. Will be ignored if `to_keep_n`
#' is supplied. Value should be > 0 and < 1.
#' @param ncomp Integer, number of latent components to construct. Ignored if `perf_res` is supplied.
#' Default value is `NULL`.
#' @return A list as per the output of the \code{\link[mixOmics]{splsda}} function.
#' @export
run_splsda <- function(splsda_input, perf_res, to_keep_n = NULL, to_keep_prop = NULL, ncomp = NULL) {

  if (!missing(perf_res)) ncomp <- perf_res$optim_ncomp

  dataset_name <- setdiff(names(splsda_input), "Y")
  multilevel <- attr(splsda_input, "multilevel")

  ## Check the number of features to retain
  to_keep <- .check_tokeep_value(to_keep_n, to_keep_prop, ncol(splsda_input[[dataset_name]]))

  ## Get the keepX vector: to_keep needs to be split by latent component
  keepX <- rep(round(to_keep / ncomp, 0), ncomp)
  keepX[1] <- to_keep - sum(keepX[-1]) ## make sure that the total adds up to to_keep
  names(keepX) <- paste0("comp", 1:ncomp)

  ## Run the sPLS-DA (with feature selection)
  splsda_res <- mixOmics::splsda(
    splsda_input[[dataset_name]],
    splsda_input$Y,
    ncomp = ncomp,
    keepX = keepX,
    multilevel = multilevel
  )

  attr(splsda_res, "dataset_name") <- dataset_name

  return(splsda_res)
}

.get_splsda_selected_features <- function(splsda_res) {
  unlist(lapply(1:splsda_res$ncomp, function(j) {
    unique(mixOmics::selectVar(splsda_res, comp = j)$name)
  }))
}


#' Get filtered MultiDataSet object based on sPLS-DA runs
#'
#' Selects features most associated with the phenotype of interest from omics datasets based on results from sPLS-DA
#' applied to the corresponding omics datasets.
#'
#' Note that the sPLS-DA method can select the same feature for several latent components, so the number of
#' features retained for a dataset might be less than the number specified in the `to_keep` argument.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param splsda_res_list A list with the result from a sPLS-DA run for each dataset to be filtered, as returned by
#' the \code{\link{run_splsda}} function.
#' @return A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @examples
#' \dontrun{
#' # Goal: keep 20% of features in dataset1, and 50% of features in dataset2
#' # outcome_group is the outcome of interest in the samples metadata
#' to_keep_prop <- c("dataset1" = 0.2, "dataset_2" = 0.5)
#'
#' # 1) assess optimal number of latent components for dataset1 and dataset2
#' splsda_perf_runs <- lapply(names(to_keep_prop), function(i) {
#'   perf_splsda(mo_data, i, "outcome_group")
#' })
#'
#' # 2) run sPLS-DA with optimal number of latent components for dataset1 and dataset2
#' splsda_runs <- lapply(splsda_perf_runs, function(x) {
#'   run_splsda(mo_data, x, to_keep_prop = to_keep_prop[attr(x, "dataset_name")])
#' })
#'
#' # 3) Get the filtered dataset
#' mo_data_filtered <- get_filtered_dataset_splsda(mo_data, splsda_runs)
#' }
#' @export
get_filtered_dataset_splsda <- function(mo_data,
                                        splsda_res_list) {
  if (!("MultiDataSet" %in% class(mo_data))) stop("Expecting MultiDataSet object")

  ## Make sure that the names of splsda_res_list are the name of the datasets
  names(splsda_res_list) <- sapply(splsda_res_list, attr, "dataset_name")

  ## Get the name of the datasets that need to be filtered and the ones that don't
  datasets_filtered <- names(splsda_res_list)
  datasets_non_filtered <- setdiff(names(mo_data), names(splsda_res_list))

  retained_features <- c(
    ## Get the list of retained features from sPLS-DA runs
    lapply(datasets_filtered, function(i) {
      .get_splsda_selected_features(splsda_res_list[[i]])
    }),
    ## Get the list of all features from non-filtered datasets
    lapply(datasets_non_filtered, function(i) {
      Biobase::featureData(mo_data[[i]])$feature_id
    })
  )

  mo_data <- subset_features(mo_data, retained_features)

  return(mo_data)
}

#' Target factory for feature preselection based on sPLS-DA
#'
#' Creates a list of targets to perform feature preselection on datasets from a `MultiDataSet`
#' object with sPLS-DA (from the `mixOmics` package).
#'
#' @param mo_data_target Symbol, the name of the target containing the `MultiDataSet` object.
#' @param group Character, the column name in the samples information data-frame to use as samples group.
#' @param to_keep_ns Named integer vector, the number of feature to retain in each dataset to be prefiltered
#' (names should correspond to a dataset name). Value should be less than the number of features in the
#' corresponding dataset. Set to `NULL` in order to use `to_keep_props` instead.
#' @param to_keep_props Named numeric vector, the proportion of features to retain in each dataset
#' to be prefiltered (names should correspond to a dataset name). Value should be > 0 and < 1.
#' Will be ignored if `to_keep_ns` is not `NULL`.
#' @param target_name_prefix Character, a prefix to add to the name of the targets created by this target factory.
#' Default value is `""`.
#' @param filtered_set_target_name Character, the name of the final target containing the filtered `MultiDataSet` object.
#' If NULL, a name will automatically be supplied. Default value is `NULL`.
#' @param multilevel Character vector of length 1 or 3 to be used as information about repeated measurements.
#' See [get_input_splsda()] for details. Default value is `NULL` (no repeated measurements).
#' @param ... Further arguments passed to the \code{\link{perf_splsda}} function.
#' @return A list of target objects. With `target_name_prefix = ""` and `filtered_set_target_name = NULL`,
#' the following targets are created:
#' * `splsda_spec`: generates a grouped tibble where each row corresponds to one dataset to be filtered,
#'   with the columns specifying each dataset name, and associated values from `to_keep_ns` and `to_keep_props`.
#'   * `individual_splsda_input`: a dynamic branching target that runs the [get_input_splsda()] function for each dataset.
#' * `individual_splsda_perf`: a dynamic branching target that runs the [perf_splsda()] function for each dataset.
#' * `individual_splsda_run`: a dynamic branching target that runs the [run_splsda()] function for each dataset,
#'   using the results from `individual_splsda_perf` to guide the number of latent components to construct.
#' * `filtered_set_slpsda`: a target to retain from the original `MultiDataSet` object only features selected in each sPLS-DA run.
#' @examples
#' \dontrun{
#' ## in the _targets.R
#' library(moiraine)
#'
#' list(
#'   ## add code here to load the different datasets
#'
#'   ## the following target creates a MultiDataSet object from previously
#'   ## created omics sets (geno_set, trans_set, etc)
#'   tar_target(
#'     mo_set,
#'     create_multiomics_set(geno_set, trans_set, metabo_set, pheno_set)
#'   ),
#'   feature_preselection_splsda_factory(
#'     mo_set,
#'     group = "outcome_group",
#'     to_keep_ns = c("rnaseq" = 1000, "metabolome" = 500),
#'     filtered_set_target_name = "mo_set_filtered",
#'     folds = 10 ## example of an argument passed to perf_splsda
#'   ),
#'
#'   ## Another example using to_keep_props
#'   feature_preselection_splsda_factory(
#'     mo_set,
#'     group = "outcome_group",
#'     to_keep_ns = NULL,
#'     to_keep_props = c("rnaseq" = 0.3, "metabolome" = 0.5),
#'     filtered_set_target_name = "mo_set_filtered",
#'     folds = 10 ## example of an argument passed to perf_splsda
#'   )
#' )
#' }
#' @export
feature_preselection_splsda_factory <- function(mo_data_target,
                                                group,
                                                to_keep_ns,
                                                to_keep_props = NULL,
                                                target_name_prefix = "",
                                                filtered_set_target_name = NULL,
                                                multilevel = NULL,
                                                ...) {
  splsda_spec_name <- paste0(target_name_prefix, "splsda_spec")
  splsda_input_name <- paste0(target_name_prefix, "individual_splsda_input")
  splsda_perf_name <- paste0(target_name_prefix, "individual_splsda_perf")
  splsda_run_name <- paste0(target_name_prefix, "individual_splsda_run")
  if (is.null(filtered_set_target_name)) filtered_set_target_name <- paste0(target_name_prefix, "filtered_set_slpsda")

  splsda_spec_target <- as.symbol(splsda_spec_name)
  splsda_input_target <- as.symbol(splsda_input_name)
  splsda_perf_target <- as.symbol(splsda_perf_name)
  splsda_run_target <- as.symbol(splsda_run_name)

  if (!is.null(to_keep_ns)) {
    dataset_names <- names(to_keep_ns)
    to_keep_ns <- unname(to_keep_ns)
  } else {
    dataset_names <- names(to_keep_props)
    to_keep_props <- unname(to_keep_props)
  }

  if (is.null(dataset_names)) {
    stop("'to_keep_ns' or 'to_keep_props' argument should be named.")
  }

  list(
    ## store the splsda specifications (arguments) as a tibble (one row per dataset to prefilter)
    ## and group it by dataset name so that following targets will be applied to each row in turn
    targets::tar_target_raw(
      splsda_spec_name,
      substitute(
        tibble::tibble(dsn = dataset_names, tkn = to_keep_ns, tkp = to_keep_props) |>
          dplyr::group_by(dsn) |>
          tar_group()),
      iteration = "group"
    ),

    ## Generate the input objects
    targets::tar_target_raw(
      splsda_input_name,
      substitute(get_input_splsda(mo_data_target, splsda_spec_target$dsn, group, multilevel)),
      pattern = substitute(map(splsda_spec_target)),
      iteration = "list"
    ),

    ## run the perf function for each row of the specfications dataframe
    targets::tar_target_raw(
      splsda_perf_name,
      substitute(perf_splsda(splsda_input_target, ...)),
      pattern = substitute(map(splsda_input_target)),
      iteration = "list"
    ),

    ## run final sPLS-DA on each dataset
    targets::tar_target_raw(
      splsda_run_name,
      substitute(
        run_splsda(
          splsda_input_target,
          perf_res = splsda_perf_target,
          to_keep_n = splsda_spec_target$tkn,
          to_keep_prop = splsda_spec_target$tkp
        )
      ),
      pattern = substitute(map(splsda_input_target, splsda_perf_target, splsda_spec_target)),
      iteration = "list"
    ),

    ## Subset the full MultiDataSet object to only retain selected features
    targets::tar_target_raw(
      filtered_set_target_name,
      substitute(get_filtered_dataset_splsda(mo_data_target, splsda_run_target))
    )
  )
}

#' Diagnostics plots for MAD-based feature preselection
#'
#' Displays the MAD distribution across all features in the original (i.e. non-filtered) datasets,
#' with a vertical red line showing the cut-off used by the preselection function.
#'
#' @param mad_list A list with the result from the MAD calculation for each dataset
#' to be filtered, as returned by the \code{\link{select_features_mad}} function.
#' @return A ggplot.
#' @export
plot_feature_preselection_mad <- function(mad_list) {
  ## Get the MAD distributions into one data-frame
  df <- lapply(mad_list, function(x) {
    x |>
      dplyr::mutate(dataset = attr(x, "dataset_name"))
  }) |>
    purrr::reduce(dplyr::bind_rows)

  ## Get for each dataset the MAD cutoff used in the preselection
  selected <- NULL ## for devtools::check()
  min_lines <- df |>
    dplyr::filter(selected) |>
    dplyr::group_by(dataset) |>
    dplyr::summarise(min = min(mad))


  ## Make the plot
  mad <- dataset <- NULL ## for devtools::check()
  ggplot2::ggplot(df, aes(x = mad, fill = dataset)) +
    ggplot2::geom_histogram(bins = 40, alpha = 0.8) +
    ggplot2::geom_vline(data = min_lines, aes(xintercept = min), colour = "tomato") +
    ggplot2::facet_wrap(~dataset, scales = "free") +
    ggplot2::labs(
      title = "Median Absolute Deviation distribution for prefiltered datasets",
      x = "Features' Median Absolute Deviation",
      y = "Count (features)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )
}


#' Diagnostics plots for COV-based feature preselection
#'
#' Displays the COV distribution across all features in the original (i.e. non-filtered) datasets,
#' with a vertical red line showing the cut-off used by the preselection function.
#'
#' @param cov_list A list with the result from the COV calculation for each dataset
#' to be filtered, as returned by the \code{\link{select_features_cov}} function.
#' @return A ggplot.
#' @export
plot_feature_preselection_cov <- function(cov_list) {
  ## Get the COV distributions into one data-frame
  df <- lapply(cov_list, function(x) {
    x |>
      dplyr::mutate(dataset = attr(x, "dataset_name"))
  }) |>
    purrr::reduce(dplyr::bind_rows)

  ## Get for each dataset the COV cutoff used in the preselection
  selected <- NULL ## for devtools::check()
  min_lines <- df |>
    dplyr::filter(selected) |>
    dplyr::group_by(dataset) |>
    dplyr::summarise(min = min(cov))

  ## Make the plot
  cov <- dataset <- NULL ## for devtools::check()
  ggplot2::ggplot(df, aes(x = cov, fill = dataset)) +
    ggplot2::geom_histogram(bins = 40, alpha = 0.8) +
    ggplot2::geom_vline(data = min_lines, aes(xintercept = min), colour = "tomato") +
    ggplot2::facet_wrap(~dataset, scales = "free") +
    ggplot2::labs(
      title = "Coefficient of Variation distribution for prefiltered datasets",
      x = "Features' Coefficient of Variation",
      y = "Count (features)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "none"
    )
}


#' Diagnostics plots for sPLS-DA-based feature preselection
#'
#' Displays the PLS-DA classification performance across different number of
#' latent components for each prefiltered dataset. The classification error rates are computed with different measures
#' (column facets) and different distance metrics (colours). A vertical grey bar represents for each dataset the number
#' of latent components selected for the feature preselection step. In addition, a circle highlights the measure and
#' distance metric used to select the number of latent component.
#'
#' @param perf_splsda_res A list with the result from the \code{\link{perf_splsda}} for each dataset
#' to be filtered.
#' @param measure Which measure(s) should be displayed? Can be one of `"BER"`
#' or `"overall"`. If NULL, all measures will be displayed. Default value is `NULL`.
#' @param distance Which measure(s) should be displayed? Can be one of `"max.dist"`,
#' `"centroids.dist"` or `"mahalanobis.dist"`. If NULL, all measures will be displayed. Default value is `NULL`.
#' @return A ggplot.
#' @export
plot_feature_preselection_splsda <- function(perf_splsda_res,
                                             measure = NULL,
                                             distance = NULL) {
  ## For devtools::check
  comp <- classification_error_rate <- classification_error_rate_sd <- selected <- NULL
  ## A bit stupid but avoids renaming the argument
  measure_plot <- measure
  distance_plot <- distance

  ## Get the cross-validation performances into one data-frame
  df <- lapply(perf_splsda_res, function(x) {
    .mixomics_get_perf_plsda(x) |>
      dplyr::mutate(
        dataset = attr(x, "dataset_name"),
        selected = (measure == x$optim_measure) &
          (distance == x$optim_distance) &
          (comp == paste0("comp", x$optim_ncomp))
      )
  }) |>
    purrr::reduce(dplyr::bind_rows) |>
    dplyr::mutate(comp = as.numeric(stringr::str_extract(comp, "\\d+")))

  .check_names(measure_plot, unique(df$measure), "'measure' argument: '_W_' not recognised. Possible measures are: '_C_'.")
  .check_names(distance_plot, unique(df$distance), "'distance' argument: '_W_' not recognised. Possible distances are: '_C_'.")

  if (!is.null(measure)) df <- dplyr::filter(df, measure %in% measure_plot)
  if (!is.null(distance)) df <- dplyr::filter(df, distance %in% distance_plot)

  to_add_ncomp <- lapply(perf_splsda_res, function(x) {
    tidyr::expand_grid(
      dataset = attr(x, "dataset_name"),
      comp = as.numeric(stringr::str_extract(x$optim_ncomp, "\\d+")),
      measure = unique(df$measure)
    )
  }) |>
    purrr::reduce(dplyr::bind_rows)

  res_plot <- ggplot2::ggplot(df, aes(x = comp, colour = distance)) +
    ggplot2::geom_vline(aes(xintercept = comp), data = to_add_ncomp, colour = "grey", alpha = 0.3, size = 5)

  if ("classification_error_rate_sd" %in% names(df)) {
    res_plot <- res_plot +
      ggplot2::geom_errorbar(
        aes(
          ymin = classification_error_rate - classification_error_rate_sd,
          ymax = classification_error_rate + classification_error_rate_sd
        ),
        alpha = 0.5,
        width = 0.2
      )
  }

  res_plot <- res_plot +
    ggplot2::geom_line(aes(y = classification_error_rate)) +
    ggplot2::geom_point(aes(y = classification_error_rate)) +
    ggplot2::geom_point(aes(y = classification_error_rate),
                        data = dplyr::filter(df, selected),
                        shape = 1,
                        size = 4,
                        show.legend = FALSE
    ) +
    ggplot2::facet_grid(dataset ~ measure) +
    ggplot2::scale_x_continuous(breaks = 1:max(df$comp)) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::labs(
      title = "PLS-DA classification performance\nfor prefiltered datasets",
      x = "Number of latent components computed by PLS-DA",
      y = "Classification error rate"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    )

  return(res_plot)
}
