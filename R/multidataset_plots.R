#' Upset plot of samples
#'
#' Generates an upset plot to compare the samples present in each omics dataset.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @returns An upset plot.
#' @export
plot_samples_upset <- function(mo_data) {
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop(
      "Package \"UpSetR\" must be installed to use this function.",
      call. = FALSE
    )
  }

  sets <- get_samples(mo_data)

  if (length(MultiDataSet::commonIds(mo_data))) {
    query <- list(list(
      query = UpSetR::intersects,
      params = names(sets),
      color = "firebrick",
      active = TRUE
    ))
  } else {
    query <- NULL
  }

  UpSetR::upset(
    UpSetR::fromList(sets),
    nsets = length(sets),
    mainbar.y.label = "Common samples between datasets",
    sets.x.label = "Number of samples per dataset",
    order.by = c("freq", "degree"),
    decreasing = c(TRUE, TRUE),
    point.size = 4,
    line.size = 2,
    text.scale = c(1.6, 1.4, 1.6, 1.4, 1.8, 1.6),
    queries = query
  )
}

#' Per-dataset density plot for MultiDataSet object
#'
#' Displays a density plot of the values in each dataset of a MultiDataSet
#' object.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param datasets Character vector, names of the datasets to include in the
#'   plot. By default, all datasets are included.
#' @param combined Logical, should the different datasets be represented in the
#'   same plot? If `FALSE` (default value), each dataset will be represented in
#'   its own subplot. Default value is `TRUE`.
#' @param scales Character, how should the axes be plotted if `combined =
#'   FALSE`. Can be either `'fixed'`, i.e. the same limits will be applied to
#'   the axes of each subplot; or `'free'`, i.e. the axis limits will be adapted
#'   to each subplot. Ignored if `combined = TRUE`. Default value is `'fixed`.
#' @returns A ggplot.
#' @export
plot_density_data <- function(mo_data,
                              datasets = names(mo_data),
                              combined = TRUE,
                              scales = "fixed") {
  mo_data <- check_input_multidataset(mo_data, datasets)
  ds_list <- get_datasets(mo_data)

  ## for devtools::check()
  feature_id <- dataset <- value <- NULL

  df <- lapply(names(ds_list), function(i) {
    tibble::as_tibble(ds_list[[i]], rownames = "feature_id") |>
      tidyr::pivot_longer(
        cols = -feature_id,
        names_to = "sample_id",
        values_to = "value"
      ) |>
      dplyr::mutate(dataset = i)
  }) |>
    purrr::reduce(dplyr::bind_rows) |>
    dplyr::mutate(dataset = factor(dataset, levels = datasets)) |>
    dplyr::filter(!is.na(value))


  the_plot <- ggplot2::ggplot(
    df,
    aes(x = value, colour = dataset, fill = dataset)
  ) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      title = "Density plots - per dataset",
      x = "Values",
      y = "Density",
      colour = "Dataset",
      fill = "Dataset"
    )

  if (!combined) {
    the_plot <- the_plot + ggplot2::facet_wrap(~dataset, scales = scales)
  }

  return(the_plot)
}

#' Per-dataset mean-sd trend plots for MultiDataSet object
#'
#' Displays for each dataset in a MultiDataSet object the trend between features
#' mean and standard deviation across all samples.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param datasets Character vector, names of the datasets to include in the
#'   plot. By default, all datasets are included.
#' @param by_rank Logical, should the x-axis display the rank of the features
#'   (ordered by mean) rather than the features mean? Default value is `FALSE`,
#'   i.e. the x axis represents the mean of the features.
#' @param colour_log10 Should the colour legend be on the log10 scale? Default
#'   value is `TRUE`.
#' @returns A ggplot.
#' @export
plot_meansd_data <- function(mo_data,
                             datasets = names(mo_data),
                             by_rank = FALSE,
                             colour_log10 = TRUE) {
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    stop(
      "Package \"hexbin\" must be installed to use this function.",
      call. = FALSE
    )
  }

  mo_data <- check_input_multidataset(mo_data, datasets)
  ds_list <- get_datasets(mo_data)

  ## for devtools::check()
  feature_id <- dataset <- mean <- sd <- use_hexbin <- use_gam <- value <- NULL

  df <- names(ds_list) |>
    purrr::map(
      function(i) {
        tibble::as_tibble(ds_list[[i]], rownames = "feature_id") |>
          tidyr::pivot_longer(
            cols = -feature_id,
            names_to = "sample_id",
            values_to = "value"
          ) |>
          dplyr::group_by(feature_id) |>
          dplyr::summarise(
            mean = mean(value, na.rm = TRUE),
            sd = stats::sd(value, na.rm = TRUE)
          ) |>
          dplyr::mutate(dataset = i) |>
          dplyr::arrange(mean) |>
          dplyr::mutate(
            rank = seq_len(dplyr::n()),
            use_hexbin = nrow(ds_list[[i]]) > 30,
            use_gam = nrow(ds_list[[i]]) > 10
          )
      }) |>
    purrr::reduce(dplyr::bind_rows) |>
    dplyr::mutate(dataset = factor(dataset, levels = datasets))

  x_var <- ifelse(by_rank, "rank", "mean")

  ggplot2::ggplot(df, ggplot2::aes(x = !!sym(x_var), y = sd)) +
    ggplot2::geom_hex(data = dplyr::filter(df, use_hexbin)) +
    ggplot2::geom_point(data = dplyr::filter(df, !use_hexbin), alpha = 0.7) +
    ggplot2::scale_fill_viridis_c(
      trans = ifelse(colour_log10, "log10", "identity"),
      option = "viridis"
    ) +
    ggplot2::facet_wrap(~dataset, scales = "free") +
    ggplot2::geom_smooth(
      data = dplyr::filter(df, use_gam),
      method = "gam",
      formula = y ~ s(x, bs = "cs"),
      colour = "deeppink"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      title = "Features mean-sd trend per dataset",
      x = ifelse(by_rank, "Features rank (sorted by mean)", "Features mean"),
      y = "Features standard deviation",
      fill = "Feature count"
    )
}


#' Plots omics data as heatmap
#'
#' For a given set of features, plots their value across the samples as a
#' heatmap.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param features Character vector, the ID of features to show in the plot.
#' @param center Logical, whether the data should be centered (feature-wise).
#'   Default value is `FALSE`.
#' @param scale Logical, whether the data should be scaled (feature-wise).
#'   Default value is `FALSE`.
#' @param samples Character vector, the ID of samples to include in the plot. If
#'   `NULL` (default), all samples will be used.
#' @param only_common_samples Logical, whether only samples that are present in
#'   all datasets should be plotted. Default value is `FALSE`.
#' @param samples_info Character vector, column names from the samples metadata
#'   tables of the datasets to be represented in the plot as samples annotation.
#' @param features_info Named list of character vectors, where each element
#'   corresponds to a dataset, and gives the column names from the features
#'   metadata of the dataset to be represented in the plot as features
#'   annotation. The names of the list must correspond to dataset names in the
#'   `mo_data` object.
#' @param colours_list Named list, where each element gives the colour palette
#'   to use in the samples or features annotation. Names must match values in
#'   `samples_info` vector and elements of `features_info` list. For continuous
#'   palettes, must use [circlize::colorRamp2()] function (see the
#'   [ComplexHeatmap reference
#'   book](https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-annotations.html#simple-annotation)).
#' @inheritParams .add_features_labels_toplot
#' @param legend_title_size Integer, size in points of legend title.
#' @param legend_text_size Integer, size in points of legend elements text.
#' @param ... Additional arguments passed to the [ComplexHeatmap::Heatmap()]
#'   function.
#' @returns A [ComplexHeatmap::Heatmap-class] object.
#' @examples
#' \dontrun{
#' ## Selecting at random 3 features from each dataset
#' random_features <- get_features(mo_set) |>
#'   map(sample, size = 3, replace = FALSE) |>
#'   unlist() |>
#'   unname()
#'
#' plot_data_heatmap(
#'   mo_set,
#'   random_features,
#'   center = TRUE,
#'   scale = TRUE,
#'   show_column_names = FALSE,
#'   only_common_samples = TRUE,
#'   samples_info = c("status", "day_on_feed"),
#'   features_info = c("chromosome"),
#'   colours_list = list(
#'     "status" = c("Control" = "gold", "BRD" = "navyblue"),
#'     "day_on_feed" = colorRamp2(c(5, 65), c("white", "pink3"))
#'   ),
#'   label_cols = list(
#'      "rnaseq" = "Name",
#'     "metabolome" = "name"
#'   )
#'  )
#' }
#' @export
plot_data_heatmap <- function(mo_data,
                              features,
                              center = FALSE,
                              scale = FALSE,
                              samples = NULL,
                              only_common_samples = FALSE,
                              samples_info = NULL,
                              features_info = NULL,
                              colours_list = NULL,
                              label_cols = NULL,
                              truncate = NULL,
                              legend_title_size = 10,
                              legend_text_size = 10,
                              ...) {
  ## For devtools::check
  dataset <- values <- label <- sample_id <- feature_id <- NULL

  mo_data <- check_input_multidataset(mo_data)
  mo_data <- subset_features(mo_data, features)

  ## removing datasets with no features
  ds <- n_features(mo_data)
  ds <- ds[ds > 0]
  if (length(ds) == 0) {
    stop("No feature selected.")
  }
  mo_data <- mo_data[, names(ds)]

  if (only_common_samples) mo_data <- MultiDataSet::commonSamples(mo_data)

  ## Get values from the datasets
  toplot <- get_datasets(mo_data) |>
    purrr::map_dfr(
      ~ .x |>
        tibble::as_tibble(rownames = "feature_id") |>
        tidyr::pivot_longer(
          cols = -feature_id,
          names_to = "sample_id",
          values_to = "values"
        ),
      .id = "dataset"
    )

  ## Make sure there are no missing features (would mean wrong input)
  .check_names(
    features,
    unique(toplot$feature_id),
    "The following features are not present in any dataset: '_W_'."
  )

  ## Filter for specific samples (happens after selecting common samples)
  if (!is.null(samples)) {
    .check_names(
      samples,
      unique(toplot$sample_id),
      "The following samples are not present in any dataset: '_W_'."
    )
    toplot <- dplyr::filter(toplot, sample_id %in% samples)
  }

  ## Checking validity of colours_list input
  if (!is.null(colours_list)) {
    .check_names(
      names(colours_list),
      c("dataset", samples_info, unlist(features_info)),
      "'colours_list' argument: '_W_' are not names of features or samples metadata to be plotted. Possible values are: '_C_'."
    )
  } else {
    colours_list <- list()
  }

  ##  -------------------- Main data matrix -------------------------- ##

  ## Main data matrix - do the centring and scaling if need be
  mat <- toplot |>
    dplyr::select(-dataset) |>
    tidyr::pivot_wider(
      names_from = sample_id,
      values_from = values
    ) |>
    tibble::column_to_rownames("feature_id") |>
    as.matrix() |>
    ## the scale function works on columns, but the features are in rows
    t() |>
    scale(center = center, scale = scale) |>
    t()

  ## Getting features labels
  if (!is.null(label_cols)) {
    row_labels <- get_features_labels(
      mo_data,
      label_cols,
      truncate
    ) |>
      dplyr::filter(feature_id %in% rownames(mat)) |>
      dplyr::select(feature_id, label) |>
      tibble::deframe()

    row_labels <- row_labels[rownames(mat)]
  } else {
    row_labels <- rownames(mat)
  }

  ##  -------------------- Features annotation -------------------------- ##

  features_annot <- toplot |>
    dplyr::select(feature_id, dataset) |>
    dplyr::distinct() |>
    tibble::column_to_rownames("feature_id") |>
    as.data.frame()

  ## adding additional info from metadata
  if (!is.null(features_info)) {
    ## extracting info
    fmeta <- get_features_metadata(mo_data) |>
      purrr::map_dfr(
        ~ dplyr::filter(.x, feature_id %in% features)
      )

    ## checking input validity
    .check_names(
      features_info,
      colnames(fmeta),
      "'features_info' argument: '_W_' is not a column in the features metadata of any dataset. Possible values are: '_C_'."
    )

    ## adding info to row annotation data-frame
    features_annot <- cbind(
      features_annot,
      fmeta |>
        dplyr::select(tidyselect::all_of(features_info))
    )
  }

  ## Default dataset colours
  if (!("dataset" %in% names(colours_list))) {
    datasets <- unique(toplot$dataset)
    dataset_colours <- rep(
      RColorBrewer::brewer.pal(12, "Paired"),
      length.out = length(datasets)
    )
    names(dataset_colours) <- datasets

    colours_list[["dataset"]] <- dataset_colours
  }

  ## Creating heatmap annotation
  row_ha <- ComplexHeatmap::rowAnnotation(
    df = features_annot,
    col = colours_list[names(colours_list) %in% colnames(features_annot)],
    gp = grid::gpar(col = "white"),
    annotation_legend_param = list(
      direction = "horizontal",
      ncol = 3,
      title_position = "topcenter",
      labels_gp = grid::gpar(fontsize = legend_text_size),
      title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold")
    )
  )

  ##  -------------------- Samples annotation -------------------------- ##

  if (!is.null(samples_info)) {
    ## checking input validity
    .check_input_var_smetadata_common(samples_info, mo_data)

    ## extracting info
    samples_annot <- get_samples_metadata_combined(
      mo_data,
      only_common_cols = FALSE
    ) |>
      dplyr::select(tidyselect::all_of(samples_info))

    ## creating heatmap annotation
    column_ha <- ComplexHeatmap::HeatmapAnnotation(
      df = samples_annot,
      col = colours_list[names(colours_list) %in% colnames(samples_annot)],
      gp = grid::gpar(col = "white"),
      annotation_legend_param = list(
        direction = "horizontal",
        ncol = 3,
        title_position = "topcenter",
        labels_gp = grid::gpar(fontsize = legend_text_size),
        title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold")
      )
    )
  } else {
    column_ha <- NULL
  }


  ##  -------------------- Final heatmap -------------------------- ##

  ## Checking whether we need to disable features and samples clustering
  args <- list(...) |>
    .check_clustering_heatmap(mat)

  args <- c(
    args,
    list(
      matrix = mat,
      row_labels = row_labels,
      ## Row and columns annotations
      right_annotation = row_ha,
      top_annotation = column_ha,
      ## Legend
      name = "Values"
    )
  )
  args$heatmap_legend_param <- c(
    args$heatmap_legend_param,
    list(
      direction = "horizontal",
      title_position = "topcenter",
      labels_gp = grid::gpar(fontsize = legend_text_size),
      title_gp = grid::gpar(fontsize = legend_title_size, fontface = "bold")
    )
  )

  hm <- do.call(ComplexHeatmap::Heatmap, args)

  ComplexHeatmap::draw(
    hm,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom"
  )
}

#' Plots omics data vs sample covariate
#'
#' For a given set of features, plots their value against a sample covariate
#' from the samples metadata. Depending on whether the covariate is continuous
#' or discrete, will generate either a scatterplot or a violin plot.
#'
#' @param mo_data A [MultiDataSet::MultiDataSet-class] object.
#' @param covariate Character, name of column in one of the samples metadata
#'   tables from `mo_data` to use as x-axis in the plot.
#' @param features Character vector, the ID of features to show in the plot.
#' @param samples Character vector, the ID of samples to include in the plot. If
#'   `NULL` (default), all samples in the corresponding dataset will be used.
#' @param only_common_samples Logical, whether only samples that are present in
#'   all datasets should be plotted. Default value is `FALSE`.
#' @param colour_by Character, name of column in one of the samples metadata
#'   tables from `mo_data` to use to colour the observations in the plot.
#'   Default value is `NULL`.
#' @param shape_by Character, name of column in one of the samples metadata
#'   tables from `mo_data` to use as shape for the observations in the plot.
#' @param point_alpha Numeric between 0 and 1, the opacity of the points in the
#'   plot (with 1 = fully opaque, and 0 = fully transparent). Default value is
#'   `1`.
#' @param add_se Logical, should a confidence interval be drawn around the
#'   smoothing curves for numerical covariates? Default value is `TRUE`.
#' @param add_boxplot Logical, should a boxplot be drawn on top of the points
#'   for categorical covariates? Default value is `TRUE`.
#' @param ncol Integer, number of columns in the faceted plot. Default value is
#'   `NULL`.
#' @inheritParams .add_features_labels_toplot
#' @returns a ggplot.
#' @examples
#' \dontrun{
#' ## Selecting at random 3 features from each dataset
#' random_features <- get_features(mo_set) |>
#'   map(sample, size = 3, replace = FALSE) |>
#'   unlist() |>
#'   unname()
#'
#' ## Plotting features value against a discrete samples covariate
#' plot_data_covariate(
#'   mo_set,
#'   "feedlot",
#'   random_features,
#'   only_common_samples = TRUE,
#'   colour_by = "status",
#'   shape_by = "geno_comp_cluster"
#' )
#'
#' ## Plotting features value against a continuous samples covariate
#' plot_data_covariate(
#'   mo_set,
#'   "day_on_feed",
#'   random_features,
#'   only_common_samples = TRUE,
#'   colour_by = "status",
#'   shape_by = "geno_comp_cluster"
#' )
#' }
#' @export
plot_data_covariate <- function(mo_data,
                                covariate,
                                features,
                                samples = NULL,
                                only_common_samples = FALSE,
                                colour_by = NULL,
                                shape_by = NULL,
                                point_alpha = 1,
                                add_se = TRUE,
                                add_boxplot = TRUE,
                                ncol = NULL,
                                label_cols = NULL,
                                truncate = NULL) {

  ## for devtools::check
  values <- id <- dataset <- label <- sample_id <- feature_id <- NULL

  mo_data <- check_input_multidataset(mo_data)
  mo_data <- subset_features(mo_data, features)

  ## removing datasets with no features
  ds <- n_features(mo_data)
  ds <- ds[ds > 0]
  if (length(ds) == 0) {
    stop("No feature selected.")
  }
  mo_data <- mo_data[, names(ds)]

  if (only_common_samples) mo_data <- MultiDataSet::commonSamples(mo_data)

  .check_input_var_smetadata_common(covariate, mo_data)
  .check_input_var_smetadata_common(colour_by, mo_data)
  .check_input_var_smetadata_common(shape_by, mo_data)

  toplot <- get_datasets(mo_data) |>
    purrr::map_dfr(
      ~ .x |>
        tibble::as_tibble(rownames = "feature_id") |>
        dplyr::filter(feature_id %in% features) |>
        tidyr::pivot_longer(
          cols = -feature_id,
          names_to = "sample_id",
          values_to = "values"
        ),
      .id = "dataset"
    ) |>
    dplyr::filter(!is.na(values))

  .check_names(
    features,
    unique(toplot$feature_id),
    "The following features are not present in any dataset: '_W_'."
  )

  if (!is.null(samples)) {

    .check_names(
      samples,
      unique(toplot$sample_id),
      "The following samples are not present in any dataset: '_W_'."
    )

    toplot <- dplyr::filter(toplot, sample_id %in% samples)
  }

  aes_vars <- c(covariate, colour_by, shape_by) |>
    unique()

  ## Dealing with features label - doing it separately because we need to
  ## make sure that the labels are unique
  df_features <- toplot |>
    dplyr::select(dataset, feature_id) |>
    dplyr::distinct() |>
    .add_features_labels_toplot(label_cols, mo_data, truncate) |>
    dplyr::mutate(
      label = make.unique(label)
    )

  toplot <- toplot |>
    dplyr::left_join(df_features, by = c("dataset", "feature_id")) |>
    ## Making sure features are ordered according to 'features' input vector
    dplyr::mutate(
      feature_id = factor(feature_id, levels = features)
    ) |>
    dplyr::arrange(feature_id) |>
    dplyr::mutate(
      label = factor(label, levels = unique(label))
    ) |>
    ## Adding samples information
    dplyr::left_join(
      get_samples_metadata_combined(mo_data, FALSE) |>
        dplyr::select(sample_id = id, tidyselect::all_of(aes_vars)),
      by = "sample_id"
    )

  p <- plot_x_wrapper(
    toplot,
    x = covariate,
    y = "values",
    facet_wrap = "label",
    colour = colour_by,
    shape = shape_by,
    point_alpha = point_alpha,
    add_se = add_se,
    add_boxplot = add_boxplot,
    ncol_wrap = ncol
  )

  return(p)
}

#' Get value from list of arguments
#'
#' Given an input named list of arguments name/value pairs, extracts the value
#' of a specific argument, or return a default value if the argument is not
#' present in the list.
#'
#' @param args Named list of arguments, where the names corresponds to the names
#'   of the arguments, and the elements to the value for the corresponding
#'   argument.
#' @param param Character, name of the argument to extract from `args`.
#' @param default Value to return if `param` is not in `args`.
#' @returns Value of `param` according to `args` or `default`.
#'
#' @noRd
.getval <- function(args, param, default) {
  if (is.null(args[[param]])) {
    return(default)
  } else {
    return(args[[param]])
  }
}

#' Check whether clustering is possible for heatmap
#'
#' Checks whether it is possible to perform row or column clustering for
#' plotting a matrix as a heatmap.
#'
#' @param args Named list of arguments to be passed on to the plotting function.
#' @param mat Matrix to be plotted.
#' @returns The named list `args` where if needed row and/or column clustering
#'   has been disabled.
#'
#' @noRd
.check_clustering_heatmap <- function(args, mat) {
  cluster_columns <- .getval(args, "cluster_columns", TRUE)
  cluster_rows <- .getval(args, "cluster_rows", TRUE)
  disabling_columns <- FALSE
  disabling_rows <- FALSE

  ## Checking whether columns clustering is possible
  if (cluster_columns) {
    clustering_distance_columns <- .getval(
      args,
      "clustering_distance_columns",
      "euclidean"
    )
    temp <- ComplexHeatmap:::get_dist(t(mat), clustering_distance_columns)

    if (any(is.na(temp))) {
      args$cluster_columns <- FALSE
      disabling_columns <- TRUE
    }
  }

  ## Checking whether rows clustering is possible
  if (cluster_rows) {
    clustering_distance_rows <- .getval(
      args,
      "clustering_distance_rows",
      "euclidean"
    )
    temp <- ComplexHeatmap:::get_dist(mat, clustering_distance_rows)

    if (any(is.na(temp))) {
      args$cluster_rows <- FALSE
      disabling_rows <- TRUE
    }
  }

  ## Throwing appropriate warning
  if (disabling_columns && disabling_rows) {
    warning(
      "Not enough data to calculate distance between features ",
      "and between samples, disabling clustering of rows and columns.",
      call. = FALSE
    )
  } else if (disabling_columns && !disabling_rows) {
    warning(
      "Not enough data to calculate distance between samples, ",
      "disabling clustering of columns.",
      call. = FALSE
    )
  } else if (!disabling_columns && disabling_rows) {
    warning(
      "Not enough data to calculate distance between features, ",
      "disabling clustering of rows.",
      call. = FALSE
    )
  }

  return(args)
}
