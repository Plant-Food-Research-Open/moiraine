test_that(".get_initials works", {
  expect_equal(
    .get_initials("principal component 10"),
    "PC10"
  )

  expect_equal(
    .get_initials(c("principal component 10", "Factor 2")),
    c("PC10", "F2")
  )
})



test_that("get_features_weight_correlation works", {
  res_list <- test_get_output_list()
  res <- get_features_weight_correlation(res_list)

  n_tot <- res_list |>
    purrr::map(
      ~ length(levels(.x$features_weight$latent_dimension))
    ) |>
    purrr::reduce(`+`)


  ld_list <- res_list |>
    purrr::map(
      ~ levels(.x$features_weight$latent_dimension)
    )

  ## Testing default return
  expect_true(is.matrix(res))
  expect_equal(dim(res), rep(n_tot, 2))
  expect_equal(
    rownames(res),
    seq_along(res_list) |>
      purrr::map(
        ~ paste0(attr(res_list[[.x]], "method"), "___", ld_list[[.x]])
      ) |>
      purrr::reduce(c)
  )

  ## Testing calculations
  expect_equal(
    res["MOFA___Factor 1", "DIABLO___Component 2"],
    cor(
      res_list[[1]]$features_weight |>
        dplyr::filter(latent_dimension == "Component 2") |>
        tidyr::replace_na(list(weight = 0)) |>
        dplyr::pull(weight),
      res_list[[2]]$features_weight |>
        dplyr::filter(latent_dimension == "Factor 1") |>
        tidyr::replace_na(list(weight = 0)) |>
        dplyr::pull(weight)
    )
  )

  ## Custom results names
  res_list_named <- res_list
  names(res_list_named) <- c("A", "B", "C")

  res2 <- get_features_weight_correlation(res_list_named)

  expect_equal(
    rownames(res2),
    seq_along(res_list) |>
      purrr::map(
        ~ paste0(names(res_list_named)[[.x]], "___", ld_list[[.x]])
      ) |>
      purrr::reduce(c)
  )
})

test_that("get_samples_score_correlation works", {
  res_list <- test_get_output_list()

  res <- get_samples_score_correlation(res_list)

  n_tot <- res_list |>
    purrr::map(
      ~ length(levels(.x$samples_score$latent_dimension))
    ) |>
    purrr::reduce(`+`)


  ld_list <- res_list |>
    purrr::map(
      ~ levels(.x$samples_score$latent_dimension)
    )

  ## Testing default return
  expect_true(is.matrix(res))
  expect_equal(dim(res), rep(n_tot, 2))
  expect_equal(
    rownames(res),
    seq_along(res_list) |>
      purrr::map(
        ~ paste0(attr(res_list[[.x]], "method"), "___", ld_list[[.x]])
      ) |>
      purrr::reduce(c)
  )

  ## Testing calculations
  expect_equal(
    res["sO2PLS___joint component 1", "DIABLO___Component 2"],
    cor(
      dplyr::full_join(
        res_list[[1]]$samples_score |>
          dplyr::filter(latent_dimension == "Component 2"),
        res_list[[3]]$samples_score |>
          dplyr::filter(latent_dimension == "joint component 1"),
        by = "sample_id"
      ) |>
        dplyr::select(sample_id, tidyselect::starts_with("score")) |>
        tibble::column_to_rownames("sample_id"),
      use = "complete.obs"
    )[1, 2]
  )

  ## Custom results names
  res_list_named <- res_list
  names(res_list_named) <- c("A", "B", "C")

  res2 <- get_samples_score_correlation(res_list_named)

  expect_equal(
    rownames(res2),
    seq_along(res_list) |>
      purrr::map(
        ~ paste0(names(res_list_named)[[.x]], "___", ld_list[[.x]])
      ) |>
      purrr::reduce(c)
  )
})

test_that(".heatmap_features_weight_corr works", {
  res_list <- test_get_output_list()

  expect_error(.heatmap_features_weight_corr(res_list), NA)
  expect_s4_class(.heatmap_features_weight_corr(res_list), "Heatmap")
})

test_that(".heatmap_samples_score_corr works", {
  res_list <- test_get_output_list()

  res1 <- .heatmap_features_weight_corr(res_list)

  expect_error(.heatmap_samples_score_corr(res_list, res1@row_dend_param$obj), NA)

  res <- .heatmap_samples_score_corr(res_list, res1@row_dend_param$obj)
  expect_s4_class(res, "Heatmap")

  ## checking that the ordering of the rows match what we would get by applying the clustering
  comp_mat <- get_samples_score_correlation(res_list)
  comp_names <- rownames(comp_mat) |>
    stringr::str_remove("^\\w+___") |>
    .get_initials()
  comp_order <- hclust_matrix_rows(abs(comp_mat)) |>
    stats::order.dendrogram()

  expect_equal(
    unname(res@row_names_param$labels)[stats::order.dendrogram(res@row_dend_param$obj)],
    comp_names[comp_order]
  )
})


test_that("comparison_heatmap_corr works", {
  res_list <- test_get_output_list()

  grDevices::pdf(NULL)

  expect_error(comparison_heatmap_corr(res_list), NA)
  res <- comparison_heatmap_corr(res_list)
  expect_equal(
    nrow(res@ht_list[[1]]@matrix),
    7
  )
  expect_equal(
    nrow(res@ht_list[[2]]@matrix),
    7
  )

  res <- comparison_heatmap_corr(res_list, list("sO2PLS" = paste0("joint component ", 1:2), "DIABLO" = "Component 2"))
  expect_equal(
    nrow(res@ht_list[[1]]@matrix),
    4
  )
  expect_equal(
    nrow(res@ht_list[[2]]@matrix),
    4
  )

  grDevices::dev.off()
})


test_that("comparison_plot_correlation works", {
  res_list <- test_get_output_list()

  expect_error(
    comparison_plot_correlation(res_list),
    "'output_list' should be of length 2; this function is for comparing the output of two integration methods."
  )
  expect_error(
    comparison_plot_correlation(res_list[1:2]),
    NA
  )

  expect_error(
    comparison_plot_correlation(res_list[1:2], by = "TEST"),
    "'by' argument should be one of 'samples' or 'features'."
  )
  expect_error(
    comparison_plot_correlation(res_list[1:2], by = "samples"),
    NA
  )
  expect_error(
    comparison_plot_correlation(res_list[1:2], by = "features"),
    NA
  )
})


test_that("consensus_importance_metric works", {
  expect_equal(consensus_importance_metric(1:5, "min"), 1)
  expect_equal(consensus_importance_metric(1:5, "max"), 5)
  expect_equal(consensus_importance_metric(1:5, "average"), 3)
  expect_equal(consensus_importance_metric(1:5, "product"), 120)
  expect_equal(consensus_importance_metric(1:5, "l2"), sqrt(1^2 + 2^2 + 3^2 + 4^2 + 5^2))
  expect_equal(consensus_importance_metric(1:5, "geometric"), prod(1:5)^(1 / 5))
  expect_equal(consensus_importance_metric(1:5, "harmonic"), 5 / (1 + 1 / 2 + 1 / 3 + 1 / 4 + 1 / 5))
})

test_that("compute_consensus_importance works", {
  res_list <- test_get_output_list()

  ld_list <- list("DIABLO" = "Component 1", "sO2PLS" = "joint component 1", "MOFA" = "Factor 1")

  ## function to extract the individual importance scores for a feature
  get_raw_importance <- function(id) {
    res_list |>
      purrr::map_dfr(
        ~ .x$features_weight |>
          dplyr::filter(
            latent_dimension %in% unname(ld_list),
            feature_id == id
          )
      ) |>
      tidyr::replace_na(list(importance = 0)) |>
      dplyr::pull(importance)
  }

  ## function to extract the smallest non-negative importance for a dataset
  get_min_importance <- function(ds) {
    res_list |>
      purrr::map_dfr(
        ~ .x$features_weight |>
          dplyr::filter(
            latent_dimension %in% unname(ld_list),
            dataset == ds
          )
      ) |>
      dplyr::filter(!is.na(importance)) |>
      dplyr::filter(importance > 0) |>
      dplyr::pull(importance) |>
      min()
  }

  res <- compute_consensus_importance(res_list, ld_list, "geometric", include_missing_features = FALSE)

  vals_C1 <- get_raw_importance("featureC_1")
  vals_C2 <- get_raw_importance("featureC_2")
  vals_C10 <- get_raw_importance("featureC_10")
  vals_C10_offset <- vals_C10
  vals_C10_offset[vals_C10_offset == 0] <- get_min_importance("metabolome") / 2
  vals_C34 <- get_raw_importance("featureC_34")

  # expect_equal(
  #   res |>
  #     dplyr::filter(feature_id == "featureC_2") |>
  #     dplyr::pull(importance),
  #   1
  # )
  # expect_equal(
  #   res |>
  #     dplyr::filter(feature_id == "featureC_1") |>
  #     dplyr::pull(importance),
  #   exp(mean(log(vals_C1))) / exp(mean(log(vals_C2)))
  # )
  # expect_equal(
  #   res |>
  #     dplyr::filter(feature_id == "featureC_10") |>
  #     dplyr::pull(importance),
  #   exp(mean(log(vals_C10_offset))) / exp(mean(log(vals_C2)))
  # )
  #
  # res <- compute_consensus_importance(res_list, ld_list, "average")
  # expect_equal(
  #   res |>
  #     dplyr::filter(feature_id == "featureC_1") |>
  #     dplyr::pull(importance),
  #   mean(vals_C1) / mean(vals_C34)
  # )
  # expect_equal(
  #   res |>
  #     dplyr::filter(feature_id == "featureC_10") |>
  #     dplyr::pull(importance),
  #   mean(vals_C10) / mean(vals_C34)
  # )


  res_list[[1]]$features_weight <- res_list[[1]]$features_weight |>
    dplyr::filter(
      !(feature_id %in% paste0("featureB_", c(34, 14, 18)))
    )

  res <- compute_consensus_importance(res_list, ld_list, "geometric", include_missing_features = FALSE)
  res1 <- compute_consensus_importance(res_list, ld_list, "geometric", include_missing_features = TRUE)
  temp <- dplyr::left_join(res, res1, by = c("dataset", "feature_id"))

  expect_equal(temp$importance.x, temp$importance.y)
})

