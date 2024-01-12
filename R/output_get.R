#' Extract output of integration method in standard format
#'
#' Extract samples score and features weight from the result of an integration
#' method. The `get_output()` function provides a wrapper around the methods'
#' specific `get_output_*()` functions.
#'
#' @param method_output The output of an integration method.
#' @param use_average_dimensions Logical, should the (weighted) average of the
#'   samples scores for each latent dimension across the datasets be used? If
#'   `FALSE`, a separate set of sample scores will be returned for each dataset
#'   for each of the latent dimensions. Only applies to sPLS, DIABLO and sO2PLS
#'   results. Default value is `TRUE`.
#' @returns An S3 object of class `output_dimension_reduction`, i.e. a named
#'   list, with the following elements:
#' * `features_weight`: tibble of features weight (loadings) for each latent
#'   dimension, with columns `feature_id`, `dataset`, `latent_dimension`,
#'   `weight` (unscaled feature weight for the corresponding latent dimension),
#'   `importance` (which corresponds to the scaled absolute weight, i.e. 1 =
#'   feature has the maximum absolute weight for the corresponding latent
#'   dimension and dataset, 0 = the feature was not selected for the
#'   corresponding latent dimension)
#' * `samples_score`: tibble of samples score for each latent component, with
#'   columns `sample_id`, `latent_dimension`, `score` (unscaled samples score
#'   for the corresponding latent dimension)
#' * `variance_explained`: tibble of the fraction of variance explained by each
#'   latent component for the relevant datasets.
#' @export
get_output <- function(method_output, use_average_dimensions = TRUE) {
  if (inherits(method_output, "pcaRes")) {
    return(get_output_pca(method_output))
  }

  if (inherits(method_output, "mixo_splsda")) {
    return(get_output_splsda(method_output))
  }

  if (inherits(method_output, "mixo_pls") ||
      inherits(method_output, "mixo_spls")) {
    return(get_output_spls(method_output, use_average_dimensions))
  }

  if (inherits(method_output, "block.plsda") ||
      inherits(method_output, "block.splsda")) {
    return(get_output_diablo(method_output, use_average_dimensions))
  }

  if (inherits(method_output, "MOFA")) {
    return(get_output_mofa2(method_output))
  }

  if (inherits(method_output, "o2m")) {
    return(get_output_so2pls(method_output, use_average_dimensions))
  }

  stop("Input unknown.")

  return(NULL)
}

#' @rdname get_output
#' @export
get_output_pca <- function(method_output) {
  ## for devtools::check
  latent_dimension <- dataset <- weight <- NULL
  importance <- score <- feature_id <- sample_id <- NULL

  if (!inherits(method_output, "pcaRes")) {
    stop("Expecting a pcaRes object (from run_pca() or pcaMethods::pca() function).")
  }

  ld_levels <- paste0(
    "Principal component ",
    seq_len(pcaMethods::nPcs(method_output))
  )

  ## Features weights
  features_weight <- pcaMethods::loadings(method_output) |>
    tibble::as_tibble(rownames = "feature_id") |>
    tidyr::pivot_longer(
      cols = -feature_id,
      names_to = "latent_dimension",
      values_to = "weight"
    ) |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "PC",
        "Principal component "
      ),
      latent_dimension = factor(latent_dimension, levels = ld_levels),
      dataset = attr(method_output, "dataset_name"),
      dataset = factor(dataset)
    ) |>
    dplyr::group_by(dataset, latent_dimension) |>
    dplyr::mutate(importance = abs(weight) / max(abs(weight))) |>
    dplyr::ungroup() |>
    dplyr::arrange(latent_dimension, dataset, feature_id) |>
    dplyr::select(feature_id, dataset, latent_dimension, weight, importance)

  ## Samples scores
  samples_score <- pcaMethods::scores(method_output) |>
    tibble::as_tibble(rownames = "sample_id") |>
    tidyr::pivot_longer(
      cols = -sample_id,
      names_to = "latent_dimension",
      values_to = "score"
    ) |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "PC",
        "Principal component "
      ),
      latent_dimension = factor(latent_dimension, levels = ld_levels)
    ) |>
    dplyr::select(sample_id, latent_dimension, score) |>
    dplyr::arrange(latent_dimension, sample_id)

  ## Variance explained
  temp <- unname(pcaMethods::R2cum(method_output))
  variance_explained <- tibble::tibble(
    latent_dimension = ld_levels,
    dataset = attr(method_output, "dataset_name"),
    prop_var_expl = c(temp[[1]], diff(temp))
  ) |>
    dplyr::mutate(
      latent_dimension = factor(latent_dimension, levels = ld_levels),
      dataset = factor(dataset)
    )

  res <- list(
    features_weight = features_weight,
    samples_score = samples_score,
    variance_explained = variance_explained
  )
  class(res) <- "output_dimension_reduction"
  attr(res, "method") <- "PCA"
  attr(res, "dataset_name") <- attr(method_output, "dataset_name")

  return(res)
}

#' @rdname get_output
#' @export
get_output_splsda <- function(method_output) {
  ## for devtools::check
  latent_dimension <- dataset <- weight <- feature_id <- sample_id <- NULL
  importance <- score <- prop_var_expl <- NULL

  if (!inherits(method_output, "mixo_splsda")) {
    stop("Expecting a mixo_splsda object (from run_splsda() function).")
  }

  ds_level <- attr(method_output, "dataset_name")
  ld_levels <- paste0("Component ", seq_len(method_output$ncomp))

  features_weight <- method_output$loadings$X |>
    tibble::as_tibble(rownames = "feature_id") |>
    tidyr::pivot_longer(
      cols = -feature_id,
      names_to = "latent_dimension",
      values_to = "weight"
    ) |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "comp",
        "Component "
      ),
      latent_dimension = factor(latent_dimension, levels = ld_levels),
      dataset = ds_level,
      dataset = factor(dataset)
    ) |>
    dplyr::group_by(latent_dimension) |>
    dplyr::mutate(importance = abs(weight) / max(abs(weight))) |>
    dplyr::ungroup() |>
    dplyr::arrange(latent_dimension, dataset, feature_id) |>
    dplyr::select(feature_id, dataset, latent_dimension, weight, importance)

  samples_score <- method_output$variates$X |>
    tibble::as_tibble(rownames = "sample_id") |>
    tidyr::pivot_longer(
      cols = -sample_id,
      names_to = "latent_dimension",
      values_to = "score"
    ) |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "comp",
        "Component "
      ),
      latent_dimension = factor(latent_dimension, levels = ld_levels)
    ) |>
    dplyr::select(sample_id, latent_dimension, score) |>
    dplyr::arrange(latent_dimension, sample_id)

  variance_explained <- method_output$prop_expl_var$X |>
    tibble::enframe(
      name = "latent_dimension",
      value = "prop_var_expl"
    ) |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "comp",
        "Component "
      ),
      latent_dimension = factor(latent_dimension, levels = ld_levels),
      dataset = ds_level,
      dataset = factor(dataset)
    ) |>
    dplyr::select(latent_dimension, dataset, prop_var_expl) |>
    dplyr::arrange(latent_dimension, dataset)

  res <- list(
    features_weight = features_weight,
    samples_score = samples_score,
    variance_explained = variance_explained
  )

  class(res) <- "output_dimension_reduction"
  attr(res, "method") <- "sPLS-DA"
  attr(res, "dataset_name") <- ds_level

  return(res)
}

#' @rdname get_output
#' @export
get_output_spls <- function(method_output, use_average_dimensions = TRUE) {
  ## for devtools::check
  dim_int <- dataset <- latent_dimension <- feature_id <- sample_id <- NULL
  weight <- score <- level <- prop_var_expl <- NULL

  if (!(inherits(method_output, "mixo_pls") ||
        inherits(method_output, "mixo_spls"))) {
    stop("Expecting a mixo_pls or mixo_spls object (from spls_run() function).")
  }

  if (inherits(method_output, "mixo_splsda")) {
    stop(
      "Expecting output from sPLS and not sPLS-DA. ",
      "Please use get_output_splsda() instead."
    )
  }

  ds_levels <- attr(method_output, "datasets_name")

  comps <- seq_len(method_output$ncomp)

  ## Replacing dataset names
  names(method_output$loadings) <- ds_levels
  names(method_output$prop_expl_var) <- ds_levels

  ## Get proportion of variance explained as tibble
  variance_explained <- purrr::map_dfr(
    method_output$prop_expl_var,
    ~ tibble::enframe(
      .x,
      name = "latent_dimension",
      value = "prop_var_expl"
    ),
    .id = "dataset"
  )

  if (use_average_dimensions) {
    ld_levels <- paste0("Component ", comps)

    features_weight <- method_output$loadings
    samples_score <- spls_get_wa_coord(method_output)
  } else {
    ld_levels <- tidyr::expand_grid(
      latent_dimension = comps,
      dataset = ds_levels
    ) |>
      dplyr::mutate(
        level = paste0(dataset, " Component ", latent_dimension)
      ) |>
      dplyr::pull(level)

    features_weight <- method_output$loadings |>
      .add_name_to_cols()
    samples_score <- method_output$variates |>
      rlang::set_names(ds_levels) |>
      .add_name_to_cols() |>
      purrr::reduce(cbind)
    variance_explained <- variance_explained |>
      dplyr::mutate(
        latent_dimension = paste0(dataset, " ", latent_dimension)
      )
  }

  ## Features weights
  features_weight <- features_weight |>
    purrr::map_dfr(
      ~ .x |>
        tibble::as_tibble(rownames = "feature_id") |>
        tidyr::pivot_longer(
          cols = -feature_id,
          names_to = "latent_dimension",
          values_to = "weight"
        ),
      .id = "dataset"
    ) |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "comp",
        "Component "
      ),
      dim_int = stringr::str_extract(latent_dimension, "\\d+$")
    ) |>
    dplyr::arrange(dim_int) |>
    dplyr::mutate(
      latent_dimension = factor(latent_dimension, levels = ld_levels),
      dataset = factor(dataset, levels = ds_levels)
    ) |>
    dplyr::select(
      feature_id,
      dataset,
      latent_dimension,
      weight
    ) |>
    dplyr::group_by(dataset, latent_dimension) |>
    dplyr::mutate(importance = abs(weight) / max(abs(weight))) |>
    dplyr::ungroup() |>
    dplyr::arrange(latent_dimension, dataset, feature_id)

  ## Samples scores
  samples_score <- samples_score |>
    tibble::as_tibble(rownames = "sample_id") |>
    tidyr::pivot_longer(
      cols = -sample_id,
      names_to = "latent_dimension",
      values_to = "score"
    ) |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "comp(?=\\d)",
        "Component "
      ),
      latent_dimension = factor(latent_dimension, levels = ld_levels)
    ) |>
    dplyr::select(
      sample_id,
      latent_dimension,
      score
    ) |>
    dplyr::arrange(latent_dimension, sample_id)

  ## Variance explained
  variance_explained <- variance_explained |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "comp(?=\\d)",
        "Component "
      ),
      latent_dimension = factor(latent_dimension, levels = ld_levels),
      dataset = factor(dataset, levels = ds_levels)
    ) |>
    dplyr::select(
      latent_dimension, dataset, prop_var_expl
    ) |>
    dplyr::arrange(latent_dimension, dataset)

  res <- list(
    features_weight = features_weight,
    samples_score = samples_score,
    variance_explained = variance_explained
  )
  class(res) <- "output_dimension_reduction"
  attr(res, "method") <- "sPLS"

  return(res)
}


#' @rdname get_output
#' @export
get_output_diablo <- function(method_output, use_average_dimensions = TRUE) {
  ## for devtools::check
  level <- prop_var_expl <- feature_id <- sample_id <- NULL

  if (!(inherits(method_output, "block.plsda") ||
        inherits(method_output, "block.splsda"))) {
    stop("Expecting a block.plsda or block.splsda object (from diablo_run() function).")
  }

  ## for devtools::check
  dim_int <- dataset <- latent_dimension <- weight <- score <- NULL

  ds_levels <- method_output$names$blocks
  ds_levels <- setdiff(ds_levels, "Y")

  comps <- seq_len(max(method_output$ncomp))

  ## Get proportion of variance explained as tibble
  variance_explained <- purrr::map_dfr(
    method_output$prop_expl_var[ds_levels],
    ~ tibble::enframe(
      .x,
      name = "latent_dimension",
      value = "prop_var_expl"
    ),
    .id = "dataset"
  )

  if (use_average_dimensions) {
    ld_levels <- paste0(
      "Component ",
      comps
    )

    features_weight <- method_output$loadings[ds_levels]
    samples_score <- diablo_get_wa_coord(method_output)
  } else {
    ld_levels <- tidyr::expand_grid(
      latent_dimension = comps,
      dataset = ds_levels
    ) |>
      dplyr::mutate(
        level = paste0(dataset, " Component ", latent_dimension)
      ) |>
      dplyr::pull(level)

    features_weight <- .add_name_to_cols(method_output$loadings[ds_levels])
    samples_score <- method_output$variates[ds_levels] |>
      .add_name_to_cols() |>
      purrr::reduce(cbind)
    variance_explained <- variance_explained |>
      dplyr::mutate(
        latent_dimension = paste0(dataset, " ", latent_dimension)
      )
  }

  ## Features weights
  features_weight <- features_weight |>
    purrr::map_dfr(
      ~ .x |>
        tibble::as_tibble(rownames = "feature_id") |>
        tidyr::pivot_longer(
          cols = -feature_id,
          names_to = "latent_dimension",
          values_to = "weight"
        ),
      .id = "dataset"
    ) |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "comp",
        "Component "
      ),
      dim_int = stringr::str_extract(latent_dimension, "\\d+$")
    ) |>
    dplyr::arrange(dim_int) |>
    dplyr::mutate(
      latent_dimension = factor(latent_dimension, levels = ld_levels),
      dataset = factor(dataset, levels = ds_levels)
    ) |>
    dplyr::select(
      feature_id,
      dataset,
      latent_dimension,
      weight
    ) |>
    dplyr::group_by(dataset, latent_dimension) |>
    dplyr::mutate(importance = abs(weight) / max(abs(weight))) |>
    dplyr::ungroup() |>
    dplyr::arrange(latent_dimension, dataset, feature_id)

  ## Samples scores
  samples_score <- samples_score |>
    tibble::as_tibble(rownames = "sample_id") |>
    tidyr::pivot_longer(
      cols = -sample_id,
      names_to = "latent_dimension",
      values_to = "score"
    ) |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "comp(?=\\d)",
        "Component "
      ),
      latent_dimension = factor(latent_dimension, levels = ld_levels)
    ) |>
    dplyr::select(
      sample_id,
      latent_dimension,
      score
    ) |>
    dplyr::arrange(latent_dimension, sample_id)

  ## Variance explained
  variance_explained <- variance_explained |>
    dplyr::mutate(
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "comp(?=\\d)",
        "Component "
      ),
      latent_dimension = factor(latent_dimension, levels = ld_levels),
      dataset = factor(dataset, levels = ds_levels)
    ) |>
    dplyr::select(
      latent_dimension, dataset, prop_var_expl
    ) |>
    dplyr::arrange(latent_dimension, dataset)

  res <- list(
    features_weight = features_weight,
    samples_score = samples_score,
    variance_explained = variance_explained
  )
  class(res) <- "output_dimension_reduction"
  attr(res, "method") <- "DIABLO"

  return(res)
}

#' @rdname get_output
#' @export
get_output_mofa2 <- function(method_output) {
  ## For devtools::check
  view <- feature <- factor_int <- dataset <- feature_id <- value <- NULL
  latent_dimension <- weight <- group <- prop_var_expl <- sample_id <- NULL

  if (!inherits(method_output, "MOFA")) {
    stop("Expecting a MOFA object (from MOFA2::run_mofa() function).")
  }

  ## Checking whether this is a MOFA or a MEFISTO output
  method_label <- ifelse(
    length(method_output@mefisto_options) > 0,
    "MEFISTO",
    "MOFA"
  )

  ds_levels <- names(method_output@dimensions$D)

  ## Features weights
  features_weight <- MOFA2::get_weights(method_output,
                                        as.data.frame = TRUE,
                                        scale = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      view = as.character(view),
      feature = as.character(feature),
      factor = stringr::str_replace(factor, "Factor", "Factor "),
      factor_int = as.numeric(stringr::str_extract(factor, "\\d+$"))
    ) |>
    dplyr::arrange(factor_int) |>
    dplyr::mutate(factor = factor(factor, levels = unique(factor))) |>
    dplyr::select(
      feature_id = feature,
      dataset = view,
      latent_dimension = factor,
      weight = value
    ) |>
    dplyr::group_by(dataset, latent_dimension) |>
    dplyr::mutate(
      importance = abs(weight) / max(abs(weight)),
      dataset = factor(dataset, levels = ds_levels)
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(latent_dimension, dataset, feature_id)

  ## Samples scores
  samples_score <- MOFA2::get_factors(method_output,
                                      as.data.frame = TRUE,
                                      scale = FALSE
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      factor = as.character(factor),
      factor = stringr::str_replace(factor, "Factor", "Factor "),
      factor = factor(factor, levels = levels(features_weight$latent_dimension))
    ) |>
    dplyr::select(
      sample_id = sample,
      latent_dimension = factor,
      score = value
    ) |>
    dplyr::arrange(latent_dimension, sample_id)

  ## Variance explained
  variance_explained <- MOFA2::get_variance_explained(
    method_output,
    as.data.frame = TRUE
  )$r2_per_factor |>
    tibble::as_tibble()
  if (dplyr::n_distinct(variance_explained$group) == 1) {
    variance_explained <- variance_explained |>
      dplyr::select(-group)
  }

  variance_explained <- variance_explained |>
    dplyr::rename(
      dataset = view,
      latent_dimension = factor,
      prop_var_expl = value
    ) |>
    dplyr::mutate(
      prop_var_expl = prop_var_expl / 100,
      latent_dimension = as.character(latent_dimension),
      latent_dimension = stringr::str_replace(
        latent_dimension,
        "Factor",
        "Factor "
      ),
      latent_dimension = factor(
        latent_dimension,
        levels = levels(features_weight$latent_dimension)
      )
    ) |>
    dplyr::select(
      latent_dimension,
      dataset,
      dplyr::everything(),
      prop_var_expl
    ) |>
    dplyr::arrange(latent_dimension, dataset)

  res <- list(
    features_weight = features_weight,
    samples_score = samples_score,
    variance_explained = variance_explained
  )
  class(res) <- "output_dimension_reduction"
  attr(res, "method") <- method_label

  return(res)
}

#' @rdname get_output
#' @export
get_output_so2pls <- function(method_output, use_average_dimensions = TRUE) {
  if (!inherits(method_output, "o2m")) {
    stop("Expecting a o2m object (from so2pls_o2m() function).")
  }

  ## for devtools::check
  latent_dimension <- dataset <- weight <- feature_id <- sample_id <- NULL
  score <- ld <- type <- component <- NULL

  names_datasets <- attr(method_output, "datasets_name")

  variance_explained <- so2pls_get_variance_explained(method_output)

  if (use_average_dimensions) {
    ## Linking matrices in sO2PLS to the latent dimensions
    ref_df <- tibble::tibble(
      lname = c("Xjoint", "Yjoint", "Xorth", "Yorth"),
      dataset = rep(names_datasets, 2),
      latent_dimension = rep(c("joint", "specific"), each = 2)
    ) |>
      dplyr::mutate(
        latent_dimension = dplyr::case_when(
          latent_dimension == "specific" ~ paste0(dataset, " ", latent_dimension),
          TRUE ~ latent_dimension
        ),
        latent_dimension = paste0(latent_dimension, " component")
      )

    ## Features weight
    features_weight <- seq_len(nrow(ref_df)) |>
      purrr::map_dfr(
        function(x) {
          lds <- OmicsPLS::loadings(method_output, ref_df$lname[[x]])
          colnames(lds) <- paste0(
            ref_df$latent_dimension[[x]],
            " ",
            seq_len(ncol(lds))
          )

          if (all(lds == 0)) {
            return(NULL)
          }

          lds |>
            tibble::as_tibble(rownames = "feature_id") |>
            tidyr::pivot_longer(
              cols = -feature_id,
              names_to = "latent_dimension",
              values_to = "weight"
            ) |>
            dplyr::mutate(dataset = ref_df$dataset[[x]])
        }
      )

    average_coords <- so2pls_get_wa_coord(method_output)
    colnames(average_coords) <- paste0(
      "joint component ",
      seq_len(ncol(average_coords))
    )

    ## Samples scores
    samples_score <- dplyr::bind_rows(
      ## taking the samples score for the joint components
      average_coords |>
        tibble::as_tibble(rownames = "sample_id") |>
        tidyr::pivot_longer(
          cols = -sample_id,
          names_to = "latent_dimension",
          values_to = "score"
        ),

      ## taking the samples score for the specific components
      purrr::map_dfr(
        3:4,
        function(x) {
          lds <- OmicsPLS::scores(method_output, ref_df$lname[[x]])
          colnames(lds) <- paste0(
            ref_df$latent_dimension[[x]],
            " ",
            seq_len(ncol(lds))
          )

          if (all(lds == 0)) {
            return(NULL)
          }

          lds |>
            tibble::as_tibble(rownames = "sample_id") |>
            tidyr::pivot_longer(
              cols = -sample_id,
              names_to = "latent_dimension",
              values_to = "score"
            )
        }
      )
    )

    ld_levels <- unique(features_weight$latent_dimension)

    variance_explained <- variance_explained |>
      dplyr::mutate(
        latent_dimension = stringr::str_remove(latent_dimension, ".+(?=joint)")
      )
  } else {
    ref_df <- tibble::tibble(
      lname = c("Xjoint", "Yjoint", "Xorth", "Yorth"),
      dataset = rep(names_datasets, 2),
      latent_dimension = rep(c("joint", "specific"), each = 2)
    ) |>
      dplyr::mutate(
        latent_dimension = paste0(dataset, " ", latent_dimension, " component")
      )

    ## Features weight
    features_weight <- seq_len(nrow(ref_df)) |>
      purrr::map_dfr(
        function(x) {
          lds <- OmicsPLS::loadings(method_output, ref_df$lname[[x]])
          colnames(lds) <- paste0(
            ref_df$latent_dimension[[x]],
            " ",
            seq_len(ncol(lds))
          )

          if (all(lds == 0)) {
            return(NULL)
          }

          lds |>
            tibble::as_tibble(rownames = "feature_id") |>
            tidyr::pivot_longer(
              cols = -feature_id,
              names_to = "latent_dimension",
              values_to = "weight"
            ) |>
            dplyr::mutate(dataset = ref_df$dataset[[x]])
        }
      )

    ## Samples score
    samples_score <- seq_len(nrow(ref_df)) |>
      purrr::map_dfr(
        function(x) {
          lds <- OmicsPLS::scores(method_output, ref_df$lname[[x]])
          colnames(lds) <- paste0(
            ref_df$latent_dimension[[x]],
            " ",
            seq_len(ncol(lds))
          )

          if (all(lds == 0)) {
            return(NULL)
          }

          lds |>
            tibble::as_tibble(rownames = "sample_id") |>
            tidyr::pivot_longer(
              cols = -sample_id,
              names_to = "latent_dimension",
              values_to = "score"
            )
        }
      )

    ld_levels <- tibble::tibble(
      ld = unique(samples_score$latent_dimension)
    ) |>
      tidyr::separate(ld,
                      c("dataset", "type", "component"),
                      sep = " (component )?",
                      remove = FALSE
      ) |>
      dplyr::mutate(dataset = factor(dataset, levels = names_datasets))

    ld_levels <- dplyr::bind_rows(
      ld_levels |>
        dplyr::filter(type == "joint") |>
        dplyr::arrange(component, dataset),
      ld_levels |>
        dplyr::filter(type == "specific")
    ) |>
      dplyr::pull(ld)
  }

  features_weight <- features_weight |>
    dplyr::mutate(
      latent_dimension = factor(latent_dimension, levels = ld_levels),
      dataset = factor(dataset, levels = names_datasets)
    ) |>
    dplyr::select(
      feature_id,
      dataset,
      latent_dimension,
      weight
    ) |>
    dplyr::group_by(dataset, latent_dimension) |>
    dplyr::mutate(importance = abs(weight) / max(abs(weight))) |>
    dplyr::ungroup() |>
    dplyr::arrange(latent_dimension, dataset, feature_id)

  samples_score <- samples_score |>
    dplyr::mutate(
      latent_dimension = factor(latent_dimension, levels = ld_levels)
    ) |>
    dplyr::select(
      sample_id,
      latent_dimension,
      score
    ) |>
    dplyr::arrange(latent_dimension, sample_id)

  variance_explained <- variance_explained |>
    dplyr::mutate(
      latent_dimension = factor(latent_dimension, levels = ld_levels),
      dataset = factor(dataset, levels = names_datasets)
    ) |>
    dplyr::arrange(latent_dimension, dataset)

  res <- list(
    features_weight = features_weight,
    samples_score = samples_score,
    variance_explained = variance_explained
  )
  class(res) <- "output_dimension_reduction"
  attr(res, "method") <- "sO2PLS"

  return(res)
}

#' Add name of element to column names for list of tibbles
#'
#' Given a list of tibbles or data-frame, append the name of each element
#' to the name of each column in the corresponding tibble.
#'
#' @param l List of tibbles or data-frames
#' @returns A list of tibbles.
#'
#' @noRd
.add_name_to_cols <- function(l) {
  l |>
    purrr::imap(
      function(.x, .y) {
        colnames(.x) <- paste0(.y, " ", colnames(.x))

        return(.x)
      }
    )
}
