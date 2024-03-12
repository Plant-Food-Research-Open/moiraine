#' Pairwise PLS datasets comparison
#'
#' Runs a Projection to Latent Structure (PLS) analysis on a pair of omics datasets, as per the mixOmics package.
#'
#' Note only one latent component is computed as only the first latent component will be used to assess the
#' correlation between datasets.
#'
#' @param mixomics_data A \code{mixOmics} input object created with \code{\link{get_input_mixomics_supervised}}.
#' @param datasets_name Character vector of length 2, the names of the two omics datasets to analyse.
#' @param ... Additional parameters to be passed to the \code{\link[mixOmics]{pls}} function.
#' @param verbose Logical, should details be printed during execution? Default value is `TRUE`.
#' @return A named list; each element is an object of class \code{\link[mixOmics]{pls}},
#' which provides the result of the PLS run. The name of the datasets analysed is stored as a character vector
#' in the `datasets_name` attribute.
#' @examples
#' \dontrun{
#' run_pairwise_pls(mo_set, c("rnaseq", "metabolome"))
#' }
#' @export
run_pairwise_pls <- function(mixomics_data, datasets_name, ..., verbose = TRUE) {
  if (length(datasets_name) != 2) stop("'datasets_name' argument: expecting a character vector of length 2.")
  .check_names(datasets_name, setdiff(names(mixomics_data), "Y"), "'datasets_name' argument: '_W_' is/are not valid dataset name. Possible dataset names are: '_C_'.")

  if (verbose) message("Running PLS on ", datasets_name[1], " and ", datasets_name[2], " datasets...")

  ## we do it this way so that the "Call" of the resulting object is more informative for the user
  cmd <- str2expression(paste0("mixOmics::pls(mixomics_data[['", datasets_name[1], "']], mixomics_data[['", datasets_name[2], "']], ncomp = 1, ...)"))
  res <- eval(cmd)

  if (verbose) message("Done.")

  attr(res, "datasets_name") <- datasets_name

  return(res)
}

#' Get pairwise correlations from PLS run
#'
#' Computes the correlation matrix between datasets based on the first component of a PLS run.
#'
#' The correlation coefficient between two datasets is computed as the correlation coefficient between the first component of each
#' dataset obtained from a Projection to Latent Structure (PLS) run, from the \code{mixOmics} package.
#'
#' @param pairwise_pls_result List containing the results from the pairwise PLS runs, computed by the \code{\link{run_pairwise_pls}} function.
#' @return A matrix of correlation coefficients between datasets.
#' @examples
#' \dontrun{
#' ## get mixomics input from MultiDataSet object
#' mixomics_data <- get_input_mixomics_supervised(mo_set, "outcome_group")
#'
#' ## Get the list of dataset names
#' datasets <- setdiff(names(mixomics_data), "Y")
#'
#' ## Get all possible pairwise combinations of dataset names
#' ds_pairs <- utils::combn(datasets, 2)
#'
#' ## run PLS for each pair of datasets
#' pls_res_list <- lapply(1:ncol(ds_pairs), function(i) {
#'   run_pairwise_pls(mo_set, ds_pairs[, i])
#' })
#'
#' ## extract the pairwise correlation matrix
#' diablo_get_pairwise_pls_corr(pls_res_list)
#' }
#' @export
diablo_get_pairwise_pls_corr <- function(pairwise_pls_result) {
  ds_list <- lapply(pairwise_pls_result, attr, "datasets_name")
  ds_names <- unique(unlist(ds_list))

  res <- diag(1, length(ds_names))
  rownames(res) <- ds_names
  colnames(res) <- ds_names

  for (i in 1:length(ds_list)) {
    ds <- ds_list[[i]]
    x <- pairwise_pls_result[[i]]
    res[ds[1], ds[2]] <- res[ds[2], ds[1]] <- stats::cor(x$variates$X[, "comp1"], x$variates$Y[, "comp1"])
  }

  return(res)
}
#' Generate DIABLO design matrix
#'
#' Generates a design matrix for the DIABLO algorithm, based on the correlation between datasets inferred from pairwise PLS runs.
#'
#' Use a threshold to detect pairs of datasets that are highly correlated. For those pairs of datasets, the corresponding cell in
#' the design matrix is set to \code{high_val}. For pairs of datasets with a correlation below the threshold, the corresponding cell in the
#' design matrix is set to \code{low_val}.
#'
#' @param cormat The correlation matrix between datasets, obtained from \code{diablo_get_pairwise_pls_corr}.
#' @param threshold Numeric, correlation value above which datasets are considered as highly correlated (see Details).
#' Default value is 0.8.
#' @param low_val Numeric, value in the design matrix for datasets that are not highly correlated.
#' Default value is 0.1.
#' @param high_val Numeric, value in the design matrix for datasets that are highly correlated.
#' Default value is 1.
#' @param y_val Numeric, value in the design matrix between datasets and the outcome (Y). Default value is 1.
#' @return A numeric matrix, to be used as design matrix when running \code{\link[mixOmics]{block.plsda}},
#' with one row per dataset and one column per dataset.
#' @export
diablo_generate_design_matrix <- function(cormat, threshold = 0.8, low_val = 0.1, high_val = 1, y_val = 1) {
  res <- matrix(0, nrow = nrow(cormat), ncol = ncol(cormat))

  ## Any correlation above 0.8 is considered strong, this is the recommended value by the mixOmics authors
  res[cormat >= threshold] <- high_val
  res[cormat < threshold] <- low_val

  rownames(res) <- rownames(cormat)
  colnames(res) <- colnames(cormat)

  res <- cbind(
    res,
    matrix(y_val, nrow = nrow(res), ncol = 1, dimnames = list(rownames(res), "Y"))
  )

  res <- rbind(
    res,
    matrix(y_val, nrow = 1, ncol = ncol(res), dimnames = list("Y", colnames(res)))
  )

  diag(res) <- 0

  return(res)
}

#' Target factory for pairwise PLS and design matrix estimation for DIABLO run
#'
#' Creates a list of targets to perform a PLS run on each pair of datasets, and uses the results to
#' assess the correlation between datasets and create a design matrix for the DIABLO algorithm.
#'
#' @param mixomics_data A \code{mixOmics} input object created with \code{\link{get_input_mixomics_supervised}}.
#' @param ... Additional parameters to be passed to the \code{\link{run_pairwise_pls}} function.
#' @param threshold Numeric, correlation value above which datasets are considered as highly correlated (see Details).
#' Default value is 0.8.
#' @param low_val Numeric, value in the design matrix for datasets that are not highly correlated.
#' Default value is 0.1.
#' @param high_val Numeric, value in the design matrix for datasets that are highly correlated.
#' Default value is 1.
#' @param y_val Numeric, value in the design matrix between datasets and the outcome (Y). Default value is 1.
#' @param target_name_prefix Character, prefix to add to the name of the targets created by the
#' factory. Default value is `""`.
#' @return A list of targets. For example, with `target_name_prefix = ""`, the following targets are created:
#' * `diablo_pairs_datasets`: a target that generates a list of all possible pairs of dataset names.
#' * `diablo_pls_runs_list`: a dynamic branching target that runs the PLS algorithm on each possible pair of datasets.
#'   The target returns a list with the PLS results for each pair of datasets.
#' * `diablo_pls_correlation_matrix`: a target that computes from the PLS results list a correlation matrix between the datasets.
#' * `diablo_design_matrix`: a target that constructs from the datasets correlation matrix a design matrix to use for the
#'   DIABLO algorithm.
#' @examples
#' \dontrun{
#' ## in the _targets.R file
#' library(moiraine)
#'
#' list(
#'   ## code to import the datasets, etc
#'   ## mo_set is the target containing the MultiDataSet object
#'   tar_target(
#'     mixomics_input,
#'     get_input_mixomics_supervised(mo_set, "outcome_group")
#'   ),
#'   diablo_pairwise_pls_factory(mixomics_input)
#' )
#' }
#' @export
diablo_pairwise_pls_factory <- function(mixomics_data, ..., threshold = 0.8, low_val = 0.1, high_val = 1, y_val = 1, target_name_prefix = "") {
  ## Target names
  pair_ds_name <- paste0(target_name_prefix, "diablo_pairs_datasets")
  pls_runs_name <- paste0(target_name_prefix, "diablo_pls_runs_list")
  cormat_name <- paste0(target_name_prefix, "diablo_pls_correlation_matrix")
  designmat_name <- paste0(target_name_prefix, "diablo_design_matrix")

  ## Target symbols
  pair_ds_target <- as.symbol(pair_ds_name)
  pls_runs_target <- as.symbol(pls_runs_name)
  cormat_target <- as.symbol(cormat_name)

  targets <- list(
    targets::tar_target_raw(
      pair_ds_name,
      substitute(utils::combn(setdiff(names(mixomics_data), "Y"), 2, simplify = FALSE))
    ),
    targets::tar_target_raw(pls_runs_name,
                            substitute(run_pairwise_pls(mixomics_data, pair_ds_target[[1]], ...)),
                            pattern = substitute(map(pair_ds_target)),
                            iteration = "list"
    ),
    targets::tar_target_raw(
      cormat_name,
      substitute(diablo_get_pairwise_pls_corr(pls_runs_target))
    ),
    targets::tar_target_raw(
      designmat_name,
      substitute(diablo_generate_design_matrix(cormat_target, threshold = threshold, low_val = low_val, high_val = high_val, y_val = y_val))
    )
  )

  return(targets)
}

#' Generate a design matrix for DIABLO
#'
#' Generates a design matrix for DIABLO, following a predesigned pattern as
#' recommended by the mixOmics authors.
#'
#' @param datasets_name Character vector, the names of the datasets to
#'   integrate. Should include the value `"Y"` to represent the samples outcome
#'   groups.
#' @param design_matrix Character, the type of design matrix to generate. Should
#'   be one of `"null"`, `"weighted_full"` or `"full"`.
#' @returns A matrix with as many rows and columns as the length of
#'   `datasets_name`, and filled with either 0 (`design_matrix = "null"`), 0.1
#'   (`design_matrix = "weighted_full"`) or 1 (`design_matrix = "full"`). Values
#'   in the diagonal are set to 0, and values in the `"Y"` row and columns are
#'   set to 1.
#' @export
diablo_predefined_design_matrix <- function(datasets_name,
                                          design_matrix = c("null",
                                                            "weighted_full",
                                                            "full")) {
  design_matrix <- rlang::arg_match(design_matrix)

  datasets_name <- c(setdiff(datasets_name, "Y"), "Y")

  temp <- c("null" = 0, "weighted_full" = 0.1, "full" = 1)
  res <- matrix(
    temp[design_matrix],
    nrow = length(datasets_name),
    ncol = length(datasets_name),
    dimnames = list(datasets_name, datasets_name)
  )

  res[, "Y"] <- 1
  res["Y", ] <- 1

  diag(res) <- 0

  res
}

#' Runs DIABLO algorithm
#'
#' Runs the DIABLO algorithm (\code{\link[mixOmics]{block.splsda}}) from the \code{mixOmics} package.
#'
#' The \code{design_matrix} argument can either be a custom design matrix (for example as constructed via the
#' \code{\link{diablo_generate_design_matrix}} function); or a character indicating the type of design matrix
#' to generate. Possible values include:
#' \itemize{
#'     \item \code{'null'}: Off-diagonal elements of the design matrix are set to 0;
#'     \item \code{'weighted_full'}: Off-diagonal elements of the design matrix are set to 0.1;
#'     \item \code{'full'}: Off-diagonal elements of the design matrix are set to 1.
#' }
#'
#' @param mixomics_data A \code{mixOmics} input object created with \code{\link{get_input_mixomics_supervised}}.
#' @param design_matrix Either numeric matrix created through \code{\link{diablo_generate_design_matrix}}, or character
#' (accepted values are \code{'null'}, \code{'weighted_full'}, \code{'full'}). See Details.
#' @param ... Arguments to be passed to the \code{\link[mixOmics]{block.splsda}} function.
#' @returns An object of class `block.splsda` (if `keepX` argument was provided) or `block.splsda`
#' (if it was not), see [mixOmics::block.splsda()] and [mixOmics::block.plsda()].
#' @export
diablo_run <- function(mixomics_data, design_matrix, ...) {
  datasets_name <- setdiff(names(mixomics_data), "Y")

  if (is.character(design_matrix)) {
    design_matrix <- diablo_predefined_design_matrix(
      datasets_name,
      design_matrix
    )
  }

  args <- list(...)
  if ("keepX" %in% names(args)) {
    res <- mixOmics::block.splsda(
      X = mixomics_data[datasets_name],
      Y = mixomics_data[["Y"]],
      design = design_matrix,
      ...
    )
  } else {
    res <- mixOmics::block.plsda(
      X = mixomics_data[datasets_name],
      Y = mixomics_data[["Y"]],
      design = design_matrix,
      ...
    )
  }

  return(res)
}


#' Plots DIABLO perf results
#'
#' Displays the error rate of a DIABLO run cross-validation to estimate the optimal number of components (\code{ncomp})
#'
#' @param perf_res The cross-validation results, computed with \code{\link[mixOmics]{perf}}.
#' @return A \code{ggplot2} object.
#' @export
diablo_plot_perf <- function(perf_res) {
  x <- perf_res[["WeightedVote.error.rate"]]
  x_sd <- perf_res[["WeightedVote.error.rate.sd"]]

  measure_labels <- c("Overall.ER" = "Overall error rate", "Overall.BER" = "Balanced error rate")

  ## for devtools::check()
  measure <- comp <- error_rate <- distance <- measure <- error_rate_min <- error_rate_max <- sd <- NULL

  toplot <- lapply(names(x), function(i) {
    df <- tibble::as_tibble(x[[i]][c("Overall.ER", "Overall.BER"), ], rownames = "measure") |>
      tidyr::pivot_longer(
        cols = -measure,
        names_to = "comp",
        values_to = "error_rate"
      )

    df_sd <- tibble::as_tibble(x_sd[[i]][c("Overall.ER", "Overall.BER"), ], rownames = "measure") |>
      tidyr::pivot_longer(
        cols = -measure,
        names_to = "comp",
        values_to = "sd"
      )

    dplyr::full_join(df, df_sd, by = c("measure", "comp")) |>
      dplyr::mutate(distance = i)
  }) |>
    purrr::reduce(dplyr::bind_rows) |>
    dplyr::mutate(
      comp = as.numeric(stringr::str_extract(comp, "\\d+")),
      error_rate_min = error_rate - sd,
      error_rate_max = error_rate + sd,
      distance = stringr::str_remove(distance, "\\.dist"),
      measure = measure_labels[measure]
    )

  ggplot2::ggplot(toplot, aes(x = comp, colour = distance)) +
    ggplot2::geom_errorbar(aes(ymin = error_rate_min, ymax = error_rate_max), width = 0.1) +
    ggplot2::geom_line(aes(y = error_rate)) +
    ggplot2::geom_point(aes(y = error_rate)) +
    ggplot2::facet_grid(measure ~ .) +
    ggplot2::scale_colour_brewer(palette = "Set1") +
    ggplot2::labs(
      x = "Components",
      y = "Weighted vote classification error rate"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

#' Get the optimal ncomp value
#'
#' Selects the optimal \code{comp} value (number of components to compute) from a DIABLO cross-validation run,
#' for a given error measure and distance metric.
#'
#' @param perf_res The cross-validation results, computed with \code{\link[mixOmics]{perf}}.
#' @param measure The error measure for which to obtain the optimal value; possible values are \code{'Overall.ER'}
#' and \code{'Overall.BER'}. Default value is \code{'Overall.BER'}.
#' @param distance The distance metric for which to obtain the optimal value; possible values are \code{'max.dist'},
#' \code{'centroids.dist'} and \code{'mahalanobis.dist'}. Default value is \code{'centroids.dist'}.
#' @param min_ncomp Integer, the minimum `ncomp` value to be returned. Default value is 1, i.e. this argument does
#' not play a role in selecting the `comp` value. Can be useful if we want at least 2 latent components for final plots.
#' @return An integer, the optimal value of \code{ncomp} to use for the DIABLO run.
#' @export
diablo_get_optim_ncomp <- function(perf_res, measure = "Overall.BER", distance = "centroids.dist", min_ncomp = 1) {
  max(perf_res$choice.ncomp$WeightedVote[measure, distance], min_ncomp)
}

#' Tunes keepX arg for DIABLO
#'
#' Performs cross-validation to estimate the optimal number of features to retain from each dataset for a DIABLO run.
#'
#' The \code{design_matrix} argument can either be a custom design matrix (for example as constructed via the
#' \code{\link{diablo_generate_design_matrix}} function); or a character indicating the type of design matrix
#' to generate. Possible values include:
#' \itemize{
#'     \item \code{'null'}: Off-diagonal elements of the design matrix are set to 0;
#'     \item \code{'weighted_full'}: Off-diagonal elements of the design matrix are set to 0.1;
#'     \item \code{'full'}: Off-diagonal elements of the design matrix are set to 1.
#' }
#'
#' @param mixomics_data A \code{mixOmics} input object created with \code{\link{get_input_mixomics_supervised}}.
#' @param design_matrix Either numeric matrix created through \code{\link{diablo_generate_design_matrix}}, or character
#' (accepted values are \code{'null'}, \code{'weighted_full'}, \code{'full'}). See Details.
#' @param keepX_list Named list, gives for each omics dataset in the mixOmics input (i.e. excluding the response Y) a vector of values
#' to test (i.e. number of features to return from this dataset). If \code{NULL} (default), a standard grid will
#' be applied for each dataset and latent component, testing values: \code{seq(5, 30, 5)}.
#' @param cpus Integer, the number of CPUs to use when running the code in parallel. For advanced users,
#' see the \code{BPPARAM} argument of \code{\link[mixOmics]{tune.block.splsda}}.
#' @param ... Arguments to be passed to the \code{\link[mixOmics]{tune.block.splsda}} function.
#' @return A list, see \code{\link[mixOmics]{tune.block.splsda}}.
#' @export
diablo_tune <- function(mixomics_data, design_matrix, keepX_list = NULL, cpus = NULL, ...) {
  ## Take care of the design matrix
  datasets_name <- setdiff(names(mixomics_data), "Y")

  if (is.character(design_matrix)) {
    design_matrix <- diablo_predefined_design_matrix(
      datasets_name,
      design_matrix
    )
  }

  ## Create the grid of values to be tested for keepX
  if (is.null(keepX_list)) {
    keepX_list <- lapply(datasets_name, function(i) {
      res <- seq(5, 30, 5)
      return(res[res <= ncol(mixomics_data[[i]])])
    })
    names(keepX_list) <- datasets_name
  } else {
    if (length(keepX_list) != length(datasets_name)) {
      stop("keepX_list argument should have as many elements as the number of omics datasets in mixomics_data.")
    }

    if (length(setdiff(names(keepX_list), datasets_name))) {
      stop("keepX_list names should match omics datasets names in mixomics_data.")
    }

    if (is.null(names(keepX_list))) names(keepX_list) <- datasets_name
  }

  BPPARAM <- .mixomics_cpus_to_bparam(cpus)

  mixOmics::tune.block.splsda(
    X = mixomics_data[datasets_name],
    Y = mixomics_data[["Y"]],
    design = design_matrix,
    test.keepX = keepX_list,
    BPPARAM = BPPARAM,
    ...
  )
}


#' Plots DIABLO tune results
#'
#' Displays the error rate of a DIABLO run cross-validation to estimate the
#' optimal number of features to retain from each dataset (\code{keepX}).
#'
#' @param tune_res The cross-validation results, computed with [diablo_tune()].
#' @return A `ggplot2` object.
#' @export
diablo_plot_tune <- function(tune_res) {
  ## For devtools::check()
  id <- error <- sd <- comp <- dataset <- nb_retained_features <- NULL

  error_label <- c(
    "overall" = "Classification error rate",
    "BER" = "Balanced error rate"
  )
  error_label <- error_label[tune_res$measure]
  datasets_name <- names(tune_res$choice.keepX)
  comps <- colnames(tune_res$error.rate)
  comps_colours <- RColorBrewer::brewer.pal(max(length(comps), 3), "Set2")
  comps_colours <- comps_colours[1:length(comps)]
  names(comps_colours) <- comps

  ## Reading error sd
  if (is.null(rownames(tune_res$error.rate.sd))) {
    rownames(tune_res$error.rate.sd) <- rownames(tune_res$error.rate)
  }
  if (is.null(colnames(tune_res$error.rate.sd))) {
    colnames(tune_res$error.rate.sd) <- colnames(tune_res$error.rate)
  }

  df_sd <- tibble::as_tibble(tune_res$error.rate.sd, rownames = "id") |>
    tidyr::pivot_longer(
      cols = tidyselect::starts_with("comp"),
      names_to = "comp",
      values_to = "sd"
    )

  ## Plot dataframe
  toplot <- tibble::as_tibble(tune_res$error.rate, rownames = "id") |>
    tidyr::pivot_longer(
      cols = tidyselect::starts_with("comp"),
      names_to = "comp",
      values_to = "error"
    ) |>
    dplyr::left_join(df_sd, by = c("id", "comp")) |>
    tidyr::separate(
      id,
      into = datasets_name,
      sep = "_",
      remove = FALSE,
      convert = TRUE
    ) |>
    dplyr::mutate(
      error_min = error - sd,
      error_max = error + sd
    ) |>
    dplyr::select(-sd)

  plots_list <- list()
  for (i in comps) {
    df <- toplot |>
      dplyr::filter(comp == i) |>
      dplyr::select(-comp) |>
      dplyr::arrange(dplyr::desc(error)) |>
      dplyr::mutate(id = factor(id, levels = id))

    ## Only plot the legend if we're plotting the first component
    plot1 <- df |>
      dplyr::select(-tidyselect::starts_with("error")) |>
      tidyr::pivot_longer(
        cols = -id,
        names_to = "dataset",
        values_to = "nb_retained_features"
      ) |>
      dplyr::mutate(dataset = factor(dataset, levels = datasets_name)) |>
      ggplot2::ggplot(aes(x = dataset, y = id, fill = nb_retained_features)) +
      ggplot2::geom_tile() +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_fill_viridis_c(
        option = "plasma",
        guide = ifelse(i == comps[1], "colourbar", "none")
      ) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid = ggplot2::element_blank(),
        legend.position = "bottom",
        plot.margin = grid::unit(c(5.5, 0, 5.5, 5.5), "pt")
      ) +
      ggplot2::labs(
        x = "Dataset",
        y = NULL,
        fill = "Number of retained features"
      )

    plot2 <- df |>
      ggplot2::ggplot(aes(x = id)) +
      ggplot2::geom_col(aes(y = error), fill = comps_colours[i], width = 1) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(0, 0.05))
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5),
        plot.margin = grid::unit(c(5.5, 5.5, 5.5, 0), "pt"),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor.y = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        x = NULL,
        y = error_label,
        title = i
      )

    plots_list[[i]] <- ggpubr::ggarrange(
      plot1,
      plot2,
      common.legend = TRUE,
      legend = "none",
      align = "h",
      widths = c(1, 5)
    )

    ## grab the legend from plot1
    if (i == comps[1]) common_legend <- ggpubr::get_legend(plot1)
  }

  ggpubr::ggarrange(
    plotlist = plots_list,
    legend.grob = common_legend,
    legend = "bottom"
  )
}

#' Formatted table with optimal keepX values
#'
#' Produces a nicely formatted table with the optimal number of features to select from each dataset for a
#' DIABLO run according to the results of a cross-validation analysis. Used mostly for producing a nice report.
#'
#' @param tune_res The cross-validation results, computed with \code{\link{diablo_tune}}.
#' @return A tibble with a Dataset column, one column for each latent component and a Total column giving the number of
#' features to retain from the corresponding dataset across all latent components.
#' @export
diablo_table_optim_keepX <- function(tune_res) {
  ## for the devtools::check()
  comp <- nb_features <- NULL

  res <- lapply(names(tune_res$choice.keepX), function(i) {
    tibble::tibble(
      Dataset = i,
      comp = paste0("Component ", 1:length(tune_res$choice.keepX[[i]])),
      nb_features = tune_res$choice.keepX[[i]]
    ) |>
      tidyr::pivot_wider(
        names_from = comp,
        values_from = nb_features
      )
  }) |>
    purrr::reduce(dplyr::bind_rows) |>
    dplyr::mutate(Total = rowSums(dplyr::across(where(is.numeric))))

  return(res)
}

## From the .add_average_blocks function in mixOmics

#' Get weighted average coordinates
#'
#' Computes the samples coordinates in the weighted average latent components space from a DIABLO result object.
#'
#' @param diablo_res The output from \code{\link[mixOmics]{block.splsda}} or \code{\link{diablo_run}}.
#' @return A matrix with one row per sample and one column per latent component.
#' @export
diablo_get_wa_coord <- function(diablo_res) {
  datasets_name <- names(diablo_res$X)
  arrays <- diablo_res$variates[datasets_name]
  weights <- diablo_res$weights

  weighted_arrays <- lapply(datasets_name, function(i) {
    variates <- arrays[[i]]
    w <- diag(weights[which(datasets_name == i), ])
    weighted_variates <- variates %*% w
    dimnames(weighted_variates) <- dimnames(variates)

    weighted_variates
  })

  wtd_sum <- Reduce(f = "+", weighted_arrays)

  ## weighted mean = weighted sum / sum(weight)
  res <- base::sweep(wtd_sum, MARGIN = 2, colSums(weights), FUN = "/")

  return(res)
}

#' Plots DIABLO output
#'
#' Displays the samples coordinates for a given latent component across the datasets.
#' This is a copy of the \code{\link[mixOmics]{plotDiablo}} function, the only difference is that the
#' margins have been increased to accomodate for a title.
#'
#' @param diablo_res The output from \code{\link[mixOmics]{block.splsda}} or \code{\link{diablo_run}}.
#' @param ncomp Integer, the latent component to plot.
#' @param legend Logical, should the legend be added to the plot? Default value is `TRUE`.
#' @param legend.ncol Integer, the number of columns for the legend. If none specified, will be calculated as
#' `min(5, nlevels(diablo_res$Y))`.
#' @param col.per.group Named character vector, provides the colours to use for each phenotypic group. Names
#' must match the levels of `diablo_res$Y`.
#' @return None.
#' @export
diablo_plot <- function(diablo_res,
                        ncomp = 1,
                        legend = TRUE,
                        legend.ncol,
                        col.per.group = NULL) {
  object <- diablo_res
  Y <- object$Y
  if (is.null(col.per.group)) col.per.group <- mixOmics::color.mixo(1:nlevels(Y))
  col.per.group <- mixOmics:::.get.cols.and.group(
    col.per.group = col.per.group,
    group = Y,
    object = object,
    n_ind = length(object$names$sample)
  )
  col.per.group <- col.per.group$col.per.group
  opar <- graphics::par()[!names(graphics::par()) %in% c("cin", "cra", "csi", "cxy", "din", "page")]
  indY <- object$indY
  object$variates <- c(object$variates[-indY], object$variates[indY])
  object$loadings <- c(object$loadings[-indY], object$loadings[indY])
  VarX <- do.call(cbind, lapply(object$variates, function(i) i[, ncomp]))
  datNames <- colnames(VarX)
  if (ncol(VarX) <= 2) {
    stop("This function is only available when there are more than 3 blocks")
  }
  if (length(ncomp) != 1 | ncomp > min(object$ncomp)) {
    stop(paste0("'ncomp' must be a numeric value lower than ", min(object$ncomp), ", which is min(object$ncomp)"))
  }
  if (missing(legend.ncol)) legend.ncol <- min(5, nlevels(Y))
  numberOfCols <- ncol(VarX) - 1
  numberOfRows <- numberOfCols
  mat <- matrix(0, nrow = numberOfRows, ncol = numberOfRows)
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) mat[i, j] <- paste(i, j, sep = "_")
  }
  plotType <- list(
    cor = mat[lower.tri(mat)],
    scatter = mat[upper.tri(mat)],
    lab = diag(mat)
  )
  graphics::par(
    mfrow = c(numberOfRows + 1, numberOfCols),
    mar = rep.int(1 / 2, 4),
    oma = c(2, 2, 6, 2)
  )
  graphics::layout(
    matrix(c(1:(numberOfCols)^2, rep((numberOfCols)^2 + 1, numberOfCols)),
           numberOfRows + 1,
           numberOfCols,
           byrow = TRUE
    ),
    heights = c(rep(1, numberOfRows), 0.25 * floor(nlevels(Y) / legend.ncol))
  )
  for (i in 1:numberOfRows) {
    for (j in 1:numberOfCols) {
      ptype <- unlist(lapply(plotType, function(x) {
        intersect(paste(i, j, sep = "_"), x)
      }))
      mixOmics:::splotMatPlot(
        x = VarX[, j],
        y = VarX[, i],
        datNames,
        Y,
        ptype,
        col.per.group
      )
      if (i == 1 & j %in% seq(2, numberOfRows, 1)) graphics::Axis(side = 3, x = VarX[, i])
      if (j == numberOfRows & i %in% seq(1, numberOfRows - 1, 1)) graphics::Axis(side = 4, x = VarX[, i])
    }
  }
  plot(1:3, 1:3, type = "n", axes = FALSE, xlab = "", ylab = "")
  if (legend) {
    graphics::legend("center",
                     legend = levels(Y), col = col.per.group,
                     pch = 19, ncol = legend.ncol, cex = 1.5
    )
  }
  graphics::par(opar)
}

#' Plots DIABLO features correlation circle
#'
#' Displays the DIABLO correlation circle plot, but uses available feature metadata
#' to display feature names.
#'
#' @param diablo_res The output from \code{\link[mixOmics]{block.splsda}} or \code{\link{diablo_run}}.
#' @inheritParams get_features_labels
#' @param ... Additional arguments passed to \code{\link[mixOmics]{plotVar}}.
#' @returns A plot (see \code{\link[mixOmics]{plotVar}}).
#' @examples
#' \dontrun{
#' # Use the default features ID for the plot
#' diablo_plot_var(
#'   diablo_final_run,
#'   mo_data,
#'   "feature_id",
#'   overlap = FALSE,
#'   cex = rep(2, 3),
#'   comp = 1:2
#' )
#'
#' # Using a different column from the feature metadata of each omics dataset
#' diablo_plot_var(
#'   diablo_final_run,
#'   mo_presel_supervised,
#'   c(
#'     "snps" = "feature_id",
#'     "rnaseq" = "gene_name",
#'     "metabolome" = "comp_name"
#'   ),
#'   overlap = FALSE,
#'   cex = rep(2, 3),
#'   comp = 1:2
#' )
#' }
#' @export
diablo_plot_var <- function(diablo_res,
                            mo_data,
                            label_cols = "feature_id",
                            truncate = NULL,
                            ...) {
  ## for devtools::check
  dataset <- data <- feature_id <- NULL

  ## Get the names of the datasets to plot
  ## (if not provided via the blocks argument, plot all datasets)
  datasets <- list(...)$blocks
  if (is.null(datasets)) datasets <- names(diablo_res$X)

  mo_data <- check_input_multidataset(mo_data, datasets)

  ## Need to get feature names only for features in input data
  features_list <- datasets |>
    purrr::map(
      ~rownames(diablo_res$loadings[[.x]])
    ) |>
    unlist()

  ## Extract labels to use for features
  feature_names <- get_features_labels(mo_data, label_cols, truncate) |>
    dplyr::filter(feature_id %in% features_list) |>
    dplyr::group_by(dataset) |>
    tidyr::nest() |>
    dplyr::mutate(
      data = purrr::map(data, tibble::deframe)
    ) |>
    tibble::deframe()


  mixOmics::plotVar(diablo_res, var.names = feature_names, ...)
}

#' Plots DIABLO circos plot
#'
#' Displays the DIABLO circos plot, but uses available feature metadata
#' to display feature names.
#'
#' @param diablo_res The output from \code{\link[mixOmics]{block.splsda}} or \code{\link{diablo_run}}.
#' @inheritParams get_features_labels
#' @param ... Additional arguments passed to \code{\link[mixOmics]{circosPlot}}.
#' @export
diablo_plot_circos <- function(diablo_res, mo_data, label_cols, truncate = NULL, ...) {
  ## for devtools::check
  dataset <- data <- feature_id <- NULL

  ## Get the names of the datasets to plot
  ## (if not provided via the blocks argument, plot all datasets)
  datasets <- list(...)$blocks
  if (is.null(datasets)) datasets <- names(diablo_res$X)

  mo_data <- check_input_multidataset(mo_data, datasets)

  ## Need to get feature names only for features in input data
  features_list <- datasets |>
    purrr::map(
      ~rownames(diablo_res$loadings[[.x]])
    ) |>
    unlist()

  ## Extract labels to use for features
  feature_names <- get_features_labels(mo_data, label_cols, truncate) |>
    dplyr::filter(feature_id %in% features_list) |>
    dplyr::group_by(dataset) |>
    tidyr::nest() |>
    dplyr::mutate(
      data = purrr::map(data, tibble::deframe)
    ) |>
    tibble::deframe()

  mixOmics::circosPlot(diablo_res, var.names = feature_names, ...)
}

#' Get parameters from DIABLO run
#'
#' Extracts the `ncomp` and `keepX` parameters from a DIABLO run and format them into a table.
#'
#' @param diablo_res The output from \code{\link[mixOmics]{block.splsda}} or \code{\link{diablo_run}}.
#' @return A tibble.
#' @export
diablo_get_params <- function(diablo_res) {
  tibble::tibble(
    Parameter = c(
      "ncomp",
      "keepX"
    ),
    Description = c(
      "Number of latent component",
      "Number of features retained in each X for each latent component"
    ),
    Value = c(
      paste0(max(diablo_res$ncomp)),
      paste0(sapply(names(diablo_res$keepX), function(i) {
        paste0(i, ": ", paste0(diablo_res$keepX[[i]], collapse = ", "))
      }), collapse = "<br>")
    )
  )
}
