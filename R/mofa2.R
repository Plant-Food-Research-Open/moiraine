#' Generate input data for MOFA2 package
#'
#' Creates an object that can be used as input for the MOFA or MEFISTO analysis implemented in the MOFA2 package. It contains the omics datasets
#' as well as features and samples metadata. Should not be called directly; instead use \link{get_input_mofa} or \link{get_input_mefisto}.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param datasets Character vector, the names of the datasets from \code{mo_data} to include in the analysis.
#' @param covariates Character or character vector of length 2, the column name(s) in the samples metadata
#' data-frames to use as continuous covariates. if NULL, creates an input object for MOFA. If not null,
#' creates an input object for MEFISTO.
#' @param groups Character, the column name in the samples metadata data-frames to use as group.
#' @param options_list A named list. Should contain at most 3 or 4 elements (depending on whether the input is for
#' MOFA or MEFISTO), named 'data_options', 'model_options', 'training_options' and 'mefisto_options' (latter only
#' for MEFISTO input).
#' @param only_common_samples Logical, whether only the samples present in all datasets should be returned.
#' Default value is \code{FALSE}.
#' @return A \code{\link[MOFA2]{MOFA}} object.
#' @export
get_input_mofa2 <- function(mo_data, datasets, covariates, groups, options_list, only_common_samples) {
  ## -------- Checking input --------

  ## Check input MultiDataSet and dataset names
  mo_data <- check_input_multidataset(mo_data, datasets)

  ## Check whether the input is mofa or mefisto
  input_mefisto <- !is.null(covariates)

  ## Check input options_list for format only
  opt_list_names <- c("data_options", "model_options", "training_options")
  if (input_mefisto) opt_list_names <- c(opt_list_names, "mefisto_options")

  if (!is.null(options_list)) {
    ## Must be a list
    if (!is.list(options_list)) {
      stop(
        "Argument 'options_list' should be a named list with length <= ",
        length(opt_list_names),
        "; possible names: '",
        paste0(opt_list_names, collapse = "', '"),
        "'."
      )
    }
    ## Must be named
    if (is.null(names(options_list))) {
      stop(
        "Argument 'options_list' should be a named list; possible names: '",
        paste0(opt_list_names, collapse = "', '"),
        "'."
      )
    }
    ## Must have correct names
    .check_names(
      names(options_list),
      opt_list_names,
      "The following names in 'options_list' arguments are not valid: '_W_'. Possible names are: '_C_'."
    )
  }

  ## -------- Getting data --------

  if (only_common_samples) mo_data <- MultiDataSet::commonSamples(mo_data)
  ds_list <- get_datasets(mo_data)

  ## -------- Getting samples metadata --------

  ## for devtools::check()
  id <- NULL

  samples_info <- get_samples_metadata_combined(mo_data, only_common_cols = FALSE) |>
    dplyr::rename(sample = id) |>  ## need to have a column 'sample' for MOFA2
    dplyr::mutate( ## avoid issues when saving models to HDF5 with date columns
      dplyr::across(
        tidyselect::where(.is_poxist),
        as.character
      )
    )

  ## Checking samples metadata for the presence of a group column
  if (any(colnames(samples_info) == "group")) {
    ## MOFA2 doesn't like if there's a column called group that is not used as group
    colnames(samples_info)[colnames(samples_info) == "group"] <- "group_metadata"
    warning("Renaming 'group' column in samples metadata to 'group_metadata'.", call. = FALSE)
  }

  ## -------- MEFISTO - checking covariates --------

  ## check that the covariates exist and that they are continuous
  if (input_mefisto) {
    if (any(covariates == "group")) covariates[covariates == "group"] <- "group_metadata"

    .check_names(
      covariates,
      colnames(samples_info),
      "In 'covariates' argument: '_W_' are not columns in the samples metadata. Possible values are: '_C_'."
    )

    if (!all(sapply(covariates, function(x) {
      is.numeric(samples_info[[x]])
    }))) {
      stop("'covariates' argument should correspond to numerical variables in the samples metadata.")
    }
  }

  ## -------- Filling samples that are not in all datasets --------

  all_samples_id <- unique(unlist(lapply(ds_list, colnames)))
  ds_list <- lapply(ds_list, function(x) {
    if (length(setdiff(all_samples_id, colnames(x)))) {
      missing_samples <- setdiff(all_samples_id, colnames(x))
      x <- cbind(
        x,
        matrix(NA, nrow = nrow(x), ncol = length(missing_samples), dimnames = list(rownames(x), missing_samples))
      )
    }
    return(x[, all_samples_id])
  })


  ## -------- Getting features metadata --------

  ## for devtools::check()
  feature_id <- feature <- view <- NULL

  ## Getting features metadata
  features_info <- lapply(datasets, function(i) {
    Biobase::fData(mo_data[[i]]) |>
      dplyr::rename(feature = feature_id) |>
      dplyr::mutate(view = i) |>
      dplyr::select(feature, view, dplyr::everything())
  }) |>
    purrr::reduce(dplyr::bind_rows) |>
    dplyr::mutate( ## avoid issues when saving models to HDF5 with date columns
      dplyr::across(
        tidyselect::where(.is_poxist),
        as.character
      )
    )


  ## -------- Checking groups --------

  if (!is.null(groups)) {
    if (groups == "group") groups <- "group_metadata"

    if (length(groups) > 1) stop("Parameter 'groups' should have length 1.")
    .check_names(
      groups,
      colnames(samples_info),
      "Parameter 'groups': '_W_' is not a column in the samples information data-frame. Possibles values for the 'groups' parameter are: '_C_'."
    )

    groups_vect <- samples_info[, groups]
    samples_info$group <- samples_info[, groups]
  } else {
    groups_vect <- NULL
  }


  ## -------- Checking likelihood options --------

  if (!is.null(options_list$model_options$likelihoods)) {
    liks <- options_list$model_options$likelihoods

    if (!identical(sort(names(liks)), sort(datasets))) {
      stop(
        "In 'options_list' parameter: options_list$model_options$likelihoods should have length ",
        length(datasets),
        " and names: '",
        paste0(datasets, collapse = "', '"),
        "'."
      )
    }

    if (any(!(liks %in% c("gaussian", "poisson", "bernoulli")))) stop("In 'options_list' parameter: possible values for options_list$model_options$likelihoods are: 'gaussian','bernoulli', 'poisson'.")

    ## If any dataset is to be modelled with a Poisson likelihood, need to make sure it's rounded first
    for (i in names(liks)[liks == "poisson"]) {
      if (!is.integer(ds_list[[i]])) {
        warning("Dataset ", i, " is to be modelled with a poisson likelihood, but is not integer. Transforming to integer.", call. = FALSE)
        ds_list[[i]] <- round(ds_list[[i]], 0)
        mode(ds_list[[i]]) <- "integer"
      }
    }

    ## If any dataset is to be modelled with a Bernoulli likelihood, need to make sure it's binary
    for (i in names(liks)[liks == "bernoulli"]) {
      if (!all(ds_list[[i]] %in% c(0, 1))) stop("Dataset ", i, " is to be modelled with a bernoulli likelihood, but contains values other than 0 and 1. Please transform the dataset or change the likelihood to be used.")
    }
  }


  ## -------- Creating the input object --------

  res <- MOFA2::create_mofa(ds_list, groups = groups_vect)

  if (input_mefisto) res <- MOFA2::set_covariates(res, covariates = t(samples_info[, covariates, drop = FALSE]))

  ## Add metadata to the object
  MOFA2::samples_metadata(res) <- samples_info
  MOFA2::features_metadata(res) <- features_info


  ## -------- Defining and checking options --------
  opts_list <- opt_list_names |>
    rlang::set_names() |>
    purrr::map(
      ~ eval(str2expression(paste0("MOFA2::get_default_", .x, "(res)")))
    )


  ## For each type of options passed through options list,
  ## Check that they have the correct names and pass the values
  ## to the options list
  purrr::iwalk(
    options_list,
    ~ .check_names(
      names(.x),
      names(opts_list[[.y]]),
      paste0(
        "Argument 'options_list': the following names from options_list$",
        .y,
        " are not ",
        stringr::str_replace(.y, "_", " "),
        " parameters in MOFA2: '_W_'.\nPossible names are: '_C_'."
      )
    )
  )

  for (i in names(options_list)) opts_list[[i]][names(options_list[[i]])] <- options_list[[i]]

  ## -------- MEFISTO - adding automatic covariate values to infer --------
  if (input_mefisto & is.null(opts_list$mefisto_options$new_values)) {
    if (!is.null(groups)) {
      cov_ref_df <- samples_info |>
        dplyr::filter(!!sym(groups) == opts_list$mefisto_options$warping_ref)
    } else {
      cov_ref_df <- samples_info
    }

    cov_range <- purrr::map(covariates, ~ range(cov_ref_df[[.x]], na.rm = TRUE))
    cov_step <- purrr::map(covariates, ~ 10^floor(log10(min(diff(sort(unique(cov_ref_df[[.x]]))), na.rm = TRUE))))


    opts_list$mefisto_options$new_values <- t(purrr::reduce(purrr::map2(cov_range, cov_step, ~ seq(.x[[1]], .x[[2]], by = .y)), expand.grid))
    rownames(opts_list$mefisto_options$new_values) <- covariates
  }

  ## -------- Adding options to MOFA2 input --------
  if (input_mefisto) {
    ## MEFISTO input
    res <- MOFA2::prepare_mofa(
      object = res,
      data_options = opts_list$data_options,
      model_options = opts_list$model_options,
      training_options = opts_list$training_options,
      mefisto_options = opts_list$mefisto_options
    )
  } else {
    ## MOFA input
    res <- MOFA2::prepare_mofa(
      object = res,
      data_options = opts_list$data_options,
      model_options = opts_list$model_options,
      training_options = opts_list$training_options
    )
  }

  return(res)
}

.is_poxist <- function(x) {
  return(inherits(x, "POSIXt"))
}
