#' Get functions used for each integration method
#'
#' @returns A named list of character vectors, each element corresponding to the
#'   name of the functions that are used for data integration through a
#'   particular method (e.g. DIABLO, MOFA, etc).
#' @export
get_method_functions <- function() {
  list(
    sPLS = c(
      "get_input_spls",
      "spls_run",
      "mixOmics::perf",
      "spls_get_optim_ncomp",
      "spls_tune"
    ),
    sO2PLS = c(
      "get_input_omicspls",
      "so2pls_crossval_o2m_adjR2",
      "so2pls_get_optim_ncomp_adj",
      "so2pls_crossval_o2m",
      "so2pls_get_optim_ncomp",
      "so2pls_crossval_sparsity",
      "so2pls_get_optim_keep",
      "so2pls_o2m"
    ),
    MOFA = c(
      "get_input_mofa",
      "run_mofa"
    ),
    DIABLO = c(
      "get_input_mixomics_supervised",
      "diablo_run",
      "mixOmics::perf",
      "diablo_get_optim_ncomp",
      "diablo_tune"
    )
  )
}

#' Aggregate regexp patterns for a search.
#'
#' @param x Character vector, the regex patterns.
#' @returns Character, aggregated regex patterns.
#' @export
aggr_patterns_fct <- function(x) {
  paste0("(", x, ")", collapse = "|")
}

#' Extract/plot running time of functions used for different integration
#' methods.
#'
#'
#' @param target_patterns Named character vector, regex patterns used to extract
#'   targets by their names.
#' @param patterns_to_methods Character vector of same length as
#'   `target_patterns`, gives the integration method to which each pattern from
#'   `target_patterns` corresponds.
#' @param target_exclude_patterns Character vector, regex used to exclude some
#'   targets from the analysis.
#' @param method_functions Named list of character vector, giving for each
#'   integration method (i.e. values in `patterns_to_methods`) the name of the
#'   functions used in the pipeline to be considered when computing the running
#'   time. Names should match values in `patterns_to_methods`.
#' @returns Either a tibble of running time for each function used as part of
#'   each integration method (`get_targets_running_time`), or a ggplot
#'   (`plot_running_time`).
#' @examples
#' \dontrun{
#' ## By default uses all targets whose name starts with "spls_" as part of the
#' ## "sPLS" method, etc
#' get_targets_running_time()
#' plot_running_time()
#'
#' ## If we ran two versions of MOFA, e.g. one with supervised input and the
#' ## other with unsupervised input
#' get_targets_running_time(target_exclude_patterns = "^mofa_unsupervised")
#' plot_running_time(target_exclude_patterns = "^mofa_unsupervised")
#'
#' ## If instead we want to compare the two MOFA run times
#' get_targets_running_time(
#'   target_patterns = c(
#'     "MOFA supervised" = "^mofa_(?!=unsupervised)",
#'     "MOFA unsupervised" = "^mofa_unsupervised"
#'   ),
#'   patterns_to_methods = c("MOFA", "MOFA")
#' )
#' plot_running_time(
#'   target_patterns = c(
#'     "MOFA supervised" = "^mofa_(?!=unsupervised)",
#'     "MOFA unsupervised" = "^mofa_unsupervised"
#'   ),
#'   patterns_to_methods = c("MOFA", "MOFA")
#' )
#' }
#' @export
get_targets_running_time <- function(target_patterns = c("sPLS" = "^spls_",
                                                         "sO2PLS" = "^so2pls_",
                                                         "MOFA" = "^mofa_",
                                                         "DIABLO" = "^diablo_"),
                                     patterns_to_methods = c("sPLS", "sO2PLS",
                                                             "MOFA", "DIABLO"),
                                     target_exclude_patterns = NULL,
                                     method_functions = get_method_functions()) {

  ## For devtools::check()
  command <- name <- fct <- pattern <- method <- bytes <- seconds <- NULL

  if (length(patterns_to_methods) != length(target_patterns)) {
    stop("Length of `target_patterns` and `patterns_to_methods` should match.")
  }

  patterns_labels <- names(target_patterns)
  if (is.null(patterns_labels)) patterns_labels <- target_patterns

  if (is.null(names(method_functions))) {
    stop("`method_functions` should be a named list of character vectors.")
  }
  missing_methods <- setdiff(unique(patterns_to_methods), names(method_functions))
  if (length(missing_methods) > 0) {
    warning(
      "The following methods don't have an entry in the `method_functions` list: ",
      paste0(missing_methods, collapse = ", "), ". ",
      "All targets matching the corresponding pattern(s) will be returned."
    )
    method_functions <- c(
      method_functions,
      missing_methods |>
        rlang::set_names() |>
        purrr::map(\(x) {"^[^\\(]+"})
    )
  }

  df_manifest <- targets::tar_manifest()

  df <- targets::tar_meta() |>
    dplyr::select(-command) |>
    dplyr::left_join(df_manifest, by = c("name"))

  if (!is.null(target_exclude_patterns)) {
    df <- df |>
      dplyr::filter(!stringr::str_detect(name, aggr_patterns_fct(target_exclude_patterns)))
  }

  purrr::map(seq_along(target_patterns), \(i){
    mth <- patterns_to_methods[[i]]
    mth_fcts <- paste0(method_functions[[mth]], "\\(")

    df |>
      dplyr::filter(stringr::str_detect(name, target_patterns[[i]])) |>
      dplyr::mutate(
        pattern = patterns_labels[[i]],
        method = mth,
        fct = stringr::str_extract(command, aggr_patterns_fct(mth_fcts)),
        fct = stringr::str_remove(fct, "\\($")
      ) |>
      dplyr::filter(!is.na(fct))
  }) |>
    purrr::list_rbind() |>
    dplyr::select(pattern, method, name, fct, bytes, seconds)
}

#' @rdname get_targets_running_time
#' @export
plot_running_time <- function(target_patterns = c("sPLS" = "^spls_",
                                                  "sO2PLS" = "^so2pls_",
                                                  "MOFA" = "^mofa_",
                                                  "DIABLO" = "^diablo_"),
                              patterns_to_methods = c("sPLS", "sO2PLS",
                                                      "MOFA", "DIABLO"),
                              target_exclude_patterns = NULL,
                              method_functions = get_method_functions()) {
  ## For devtools::check()
  pattern <- seconds <- fct <- fct_pat <- time <- time_accr <- total <- n <- NULL

  toplot <- get_targets_running_time(
    target_patterns,
    patterns_to_methods,
    target_exclude_patterns,
    method_functions
  )  |>
    dplyr::arrange(pattern, dplyr::desc(seconds)) |>
    dplyr::mutate(
      fct_pat = paste0(pattern, "_", fct),
      fct_pat = factor(fct_pat, levels = rev(unique(fct_pat)))
    ) |>
    dplyr::group_by(pattern) |>
    dplyr::mutate(
      time = seconds / 60,
      total = sum(time),
      n = (1:dplyr::n()) %% 2,
    ) |>
    dplyr::ungroup()

  toadd <- toplot |>
    dplyr::group_by(pattern) |>
    dplyr::arrange(dplyr::desc(fct_pat)) |>
    dplyr::mutate(
      time_accr = cumsum(time),
      rank = rank(dplyr::desc(time))
    ) |>
    dplyr::slice_max(time, n = 3, with_ties = TRUE) |>
    dplyr::group_by() |>
    dplyr::filter(!(time < 1 & rank > 1)) |>
    dplyr::mutate(time = time_accr - (time / 2))


  breaks_fct <- function(x) {
    if (max(x) > 60) {
      seq(0, max(x), by = 30)
    } else {
      scales::extended_breaks()(x)
    }
  }

  toplot |>
    ggplot2::ggplot(ggplot2::aes(x = stats::reorder(pattern, total), y = time)) +
    ggplot2::geom_col(
      ggplot2::aes(colour = fct_pat, fill = factor(n)),
      width = 0.5,
      position = "stack",
      linewidth = 0,
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = fct),
      data = toadd,
      nudge_x = 0.5
    ) +
    ggplot2::scale_y_continuous(
      breaks = breaks_fct,
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::scale_fill_manual(values = c("#8350ba", "#aa87d0")) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Running time (min)"
    ) +
    ggplot2::theme_bw()
}
