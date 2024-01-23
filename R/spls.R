#' Generate sPLS input data (for mixomics)
#'
#' Creates an object that can be used as input for the (s)PLS functions from the
#' mixOmics package. It contains the two omics datasets selected, restricted to
#' common samples.
#'
#' `multilevel` argument: enables the multilevel option (see [mixOmics
#' site](http://mixomics.org/methods/multilevel/)) to deal with repeated
#' measurements. [mixOmics::spls()] enables one- and two-factor decomposition.
#' For one-factor decomposition, `multilevel` argument should be the name of the
#' column in the samples metadata that gives the ID of the observation units
#' (e.g. the ID of the subjects that were measured several times). The resulting
#' design matrix (stored in the `multilevel` argument of the returned object)
#' will be a data-frame with one column which gives the ID (as integer) of the
#' observation units corresponding to each sample in the omics datasets. For
#' two-factor decomposition, `multilevel` should be of length 3. The first
#' value, similarly to the one-factor decomposition, should be the name of the
#' column in the samples metadata that gives the ID of the observation units
#' (e.g. the ID of the subjects that were measured several times). The second
#' and third values should be the name of the columns in the samples metadata
#' that give the two factors considered. The resulting design matrix (stored in
#' the `multilevel` argument of the returned object) will be a data-frame with
#' three columns: the first column gives the ID (as integer) of the observation
#' units corresponding to each sample in the omics datasets; the second and
#' third columns give the levels of the two factors.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param mode Character, the mode of PLS to use of the analysis (see [sPLS
#'   documentation](http://mixomics.org/methods/spls/)). Should be one of
#'   `'regression'`, `'canonical'`, `'invariant'` or `'classic'`.
#' @param datasets Character vector of length 2, the names of the datasets from
#'   `mo_data` to include in the analysis.
#' @param multilevel Character vector of length 1 or 3 to be used as information
#'   about repeated measurements. See Details. Default value is `NULL`, i.e.
#'   the multilevel option will not be used.
#' @returns A list, in which each element corresponds to one omics dataset, with
#'   samples as rows and features as columns. The mode to use for the analysis
#'   is stored in the `mode` attribute of the returned object.
#' @export
get_input_spls <- function(mo_data,
                           mode,
                           datasets = names(mo_data),
                           multilevel = NULL) {

  if (length(datasets) != 2) {
    stop("'datasets' argument: expecting a vector of length 2 (as (s)PLS integrates 2 datasets).")
  }

  .check_names(
    mode,
    c("regression", "canonical", "invariant", "classic"),
    "'mode' argument should be one of: '_C_'."
  )

  res <- get_input_mixomics_unsupervised(mo_data, datasets)
  attr(res, "mode") <- mode

  if (!is.null(multilevel)) {
    attr(res, "multilevel") <- .get_multilevel(res, mo_data, multilevel)
  }

  return(res)
}




#' Run sPLS algorithm
#'
#' Runs the sPLS algorithm ([mixOmics::spls()]) from the mixOmics package.
#'
#' @param spls_input A mixOmics input object created with [get_input_spls()].
#' @param ... Arguments passed to [mixOmics::spls()].
#' @returns An object of class `mixo.spls` (if `keepX` and/or `keepY` arguments
#'   were provided) or `mix.pls` (if they were not), see [mixOmics::spls()] and
#'   [mixOmics::pls()].
#' @export
spls_run <- function(spls_input, ...) {

  args <- list(...)

  if (any(c("keepX", "keepY") %in% names(args))) {
    spls_res <- mixOmics::spls(
      spls_input[[1]],
      spls_input[[2]],
      mode = attr(spls_input, "mode"),
      multilevel = attr(spls_input, "multilevel"),
      ...
    )
  } else {
    spls_res <- mixOmics::pls(
      spls_input[[1]],
      spls_input[[2]],
      mode = attr(spls_input, "mode"),
      multilevel = attr(spls_input, "multilevel"),
      ...
    )
  }

  attr(spls_res, "datasets_name") <- names(spls_input)

  return(spls_res)
}

#' Select the optimal ncomp from sPLS cross-validation results
#'
#' Select the optimal number of components to compute in a sPLS run from the
#' cross-validation results obtained with [mixOmics::perf()] on a PLS or sPLS
#' result, using the the mean `Q2.total` values. Note that this function is
#' experimental, the corresponding diagnostic plots should be considered when
#' selecting the optimal number of components to use.
#'
#' The selection is made as follows:
#' * If all `Q2` values are below the threshold specified with `thr`, the
#' number of components yielding the highest `Q2` value is selected.
#' * If all `Q2` values are above the threshold, the number of components
#' yielding the lowest `Q2` value is selected.
#' * If the `Q2` values are increasing, the number of components `n` is
#' selected such that `n+1` is the smallest number of components with a `Q2`
#' value above the threshold.
#' * If the `Q2` values are decreasing, the number of components `n` is
#' selected such that `n+1` is the smallest number of components with a `Q2`
#' value below the threshold.
#'
#' @param spls_perf List, the result of [mixOmics::perf()].
#' @param thr Numeric, the threshold to be used for the `Q2` values. Default
#'   value is 0.0975.
#' @param min_ncomp Integer, the minimum `ncomp` value to be returned. Default
#'   value is 1, i.e. this argument does not play a role in selecting the `comp`
#'   value. Can be useful if we want at least 2 latent components for final
#'   plots.
#' @returns An integer, the optimal number of components to use for the sPLS
#'   run.
#' @export
spls_get_optim_ncomp <- function(spls_perf, thr = 0.0975, min_ncomp = 1) {
  df <- spls_perf$measures$Q2.total$summary
  comps <- df$comp
  q2_vals <- df$mean

  ## If all below the threshold, select the one with the highest value
  if (all(q2_vals < thr)) {
    res <- comps[which.max(q2_vals)]
  } else if (all(q2_vals > thr)) {
    ## If all above the threshold, select the one with the lowest value
    res <- comps[which.min(q2_vals)]
  } else {
    ## Otherwise, 2 scenarios (that are mentioned in the mixOmics tutorial)

    ## 1. if the Q2 values decrease, we want to retain the highest comp value
    ##    with a Q2 value above the threshold
    res1 <- which(q2_vals < thr)[1] - 1

    ## 2. if the Q2 values increase, we want to retain the highest comp value
    ##   with a Q2 value below the threshold
    res2 <- which(q2_vals > thr)[1] - 1

    res <- comps[max(res1, res2)]
  }

  res <- max(res, min_ncomp)
  return(res)
}

#' Performs cross-validation for mixomics sPLS to select optimal keepX and keepY
#'
#' Peforms cross-validation to assess the optimal number of features to retain
#' in each dataset for a sPLS run (implemented in the `mixOmics` package).
#'
#' If no value is provided for `keepX` or `keepY`, the sequence `seq(5, 30, 5)`
#' will be used, truncated to only retain values that are inferior or equal to
#' the number of features in the X or Y dataset.
#'
#' @param spls_input A `mixOmics` input object created with
#'   [get_input_mixomics_unsupervised()].
#' @param keepX Numeric vector, values for the number of features to retain from
#'   dataset X to test. Default value is `NULL` (default sequence of values will
#'   be used, see details).
#' @param keepY Numeric vector, values for the number of features to retain from
#'   dataset Y to test. Default value is `NULL` (default sequence of values will
#'   be used, see details).
#' @param cpus Integer, the number of CPUs to use when running the code in
#'   parallel. For advanced users, see the `BPPARAM` argument of
#'   [mixOmics::tune.spls()].
#' @param ... Further arguments passed to [mixOmics::tune.spls()].
#' @returns A list (see [mixOmics::tune.spls()]).
#' @export
spls_tune <- function(spls_input,
                      keepX = NULL,
                      keepY = NULL,
                      cpus = NULL,
                      ...) {


  ## Getting the list of keepX values to test
  if (is.null(keepX)) keepX <- seq(5, 30, 5)
  keepX <- keepX[keepX <= ncol(spls_input[[1]])]

  ## Getting the list of keepY values to test
  if (is.null(keepY)) keepY <- seq(5, 30, 5)
  keepY <- keepY[keepY <= ncol(spls_input[[2]])]

  BPPARAM <- .mixomics_cpus_to_bparam(cpus)

  spls_tune_res <- mixOmics::tune.spls(
    spls_input[[1]],
    spls_input[[2]],
    test.keepX = keepX,
    test.keepY = keepY,
    mode = attr(spls_input, "mode"),
    multilevel = attr(spls_input, "multilevel"),
    BPPARAM = BPPARAM,
    ...
  )

  attr(spls_tune_res, "datasets_name") <- names(spls_input)

  return(spls_tune_res)
}

#' Displays results of sPLS tuning
#'
#' Displays the results of cross-validation to tune the number of components to
#' retain in each dataset for a sPLS run. Similar to
#' [mixOmics::plot.tune.spls()].
#'
#' The plot displays the correlation or RSS between the latent components
#' obtained with the corresponding values of `keepX` (x-axis) and `keepY`
#' (y-axis) and the latent components of the full model (i.e. that retains all
#' features). The colour of the points shows the mean correlation/RSS across the
#' cross-validation folds, and the size of the points' shadow (in gray)
#' represents the coefficient of variation (COV) of the correlation/RSS, i.e.
#' standard error divided by mean. The point corresponding with the optimal
#' value for `keepX` and `keepY` is indicating with a red border.
#'
#' @param spls_tune_res The result of [spls_tune()] or [mixOmics::tune.spls()].
#' @returns A ggplot.
#' @export
spls_plot_tune <- function(spls_tune_res) {
  ## for devtools::check
  measure <- comp <- V <- sd <- mean <- NULL

  datasets_name <- attr(spls_tune_res, "datasets_name") |>
    paste0(" dataset")

  df <- tibble::as_tibble(spls_tune_res$measure.pred) |>
    dplyr::filter(measure == spls_tune_res$call$measure) |>
    dplyr::select(-tidyselect::starts_with("value")) |>
    dplyr::mutate(
      comp = paste0("Component ", comp),
      cov = sd / mean,
      V = dplyr::case_when(
        V == "u" ~ datasets_name[[1]],
        V == "t" ~ datasets_name[[2]]
      )
    )

  keepX <- keepY <- optimum.keepA <- cov <- NULL
  ggplot2::ggplot(df, aes(x = keepX, y = keepY)) +
    ggplot2::geom_point(
      aes(size = cov, colour = optimum.keepA),
      shape = 21,
      fill = "gray70",
      colour = "gray60"
    ) +
    ggplot2::geom_point(
      aes(fill = mean, colour = optimum.keepA),
      size = 3,
      shape = 21
    ) +
    ggplot2::scale_colour_manual(
      values = c("FALSE" = "black", "TRUE" = "red"),
      guide = "none"
    ) +
    ggplot2::scale_fill_viridis_c(option = "plasma", direction = 1) +
    ggplot2::scale_size(range = c(3, 6)) +
    ggplot2::facet_grid(comp ~ V) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "bottom"
    ) +
    ggplot2::labs(
      title = paste0("Measure: '", spls_tune_res$call$measure, "'"),
      x = paste0("Nb of features retained in ", datasets_name[[1]]),
      y = paste0("Nb of features retained in ", datasets_name[[2]]),
      fill = "Mean",
      size = "COV"
    )
}


#' Get parameters from sPLS run
#'
#' Extracts the `ncomp`, `keepX` and `keepY` parameters from a sPLS run and
#' format them into a table.
#'
#' @param spls_res The output from \code{\link[mixOmics]{spls}} or
#'   \code{\link{spls_run}}.
#' @return A tibble.
#' @export
spls_get_params <- function(spls_res) {
  ds_names <- attr(spls_res, "datasets_name")
  if (is.null(ds_names)) ds_names <- c("X", "Y")

  keep_descr <- "Number of features retained in dataset X for each latent component"

  tibble::tibble(
    Parameter = c(
      "ncomp",
      "keepX",
      "keepY"
    ),
    Description = c(
      "Number of latent component",
      stringr::str_replace(keep_descr, "X", ds_names[[1]]),
      stringr::str_replace(keep_descr, "X", ds_names[[2]])
    ),
    Value = c(
      paste0(spls_res$ncomp),
      paste0(spls_res$keepX, collapse = ", "),
      paste0(spls_res$keepY, collapse = ", ")
    )
  )
}


#' Computes average sample coordinates for sPLS components
#'
#' Computes the average sample coordinates for sPLS components across the two
#' datasets.
#'
#' @param spls_res The output from the [spls_run()] or [mixOmics::spls()]
#'   function.
#' @returns A matrix of samples coordinates, with samples in rows and joint
#'   components in columns.
#' @export
spls_get_wa_coord <- function(spls_res) {

  res <- (spls_res$variates$X + spls_res$variates$Y) / 2
  return(res)
}

#' Get design table for mixOmics multilevel option
#'
#' Constructs the design matrix for the multilevel option for the mixOmics
#' package.
#'
#' @param spls_input A mixOmics input object created with [get_input_spls()].
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param multilevel Character vector of length 1 or 3 to be used as information
#'   about repeated measurements.
#' @returns A data-frame, the design matrix to use.
#' @noRd
.get_multilevel <- function(spls_input, mo_data, multilevel) {
  if (is.null(mo_data)) {
    stop("Need to provide 'mo_data' argument to use the 'multilevel' option.")
  }

  if (!(length(multilevel) %in% c(1, 3))) {
    stop(
      "'multilevel' argument should be of length 1 ",
      "(for one-factor decomposition) or 3 (for two-factor decomposition)."
    )
  }
  mo_data <- check_input_multidataset(mo_data, names(spls_input))

  ## extract samples metadata
  .check_input_var_smetadata_common(multilevel, mo_data)
  samples_info <- get_samples_metadata_combined(
    mo_data,
    only_common_cols = FALSE
  )

  res <- samples_info[rownames(spls_input[[1]]), multilevel, drop = FALSE] |>
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = as.factor
      )
    )

  ## Turning first column into integer IDs since that is what we see on the
  ## mixOmics examples (see ?mixOmics::spls)
  res[, 1] <- res[, 1] |>
    as.integer()

  return(res)
}

