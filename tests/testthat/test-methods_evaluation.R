test_that("make_feature_sets_from_df works", {
  multiomics_set <- test_get_multidataset()

  ## Make a fake annotation data-frame
  features_list <- get_features_metadata(multiomics_set) |>
    purrr::map(~ dplyr::pull(.x, feature_id)) |>
    unlist() |>
    unname()

  set.seed(612)
  fake_sets <- paste0("set", 1:5) |>
    rlang::set_names() |>
    purrr::map(~ sample(features_list, 10, replace = FALSE))

  fake_sets_df <- fake_sets |>
    purrr::map_dfr(~ tibble::tibble(x = .x),
                   .id = "my_sets"
    )

  ## Checking input
  expect_error(
    make_feature_sets_from_df(fake_sets_df, col_id = "TEST", col_set = "my_sets"),
    "'col_id' argument: 'TEST' is not a column in 'annotation_df'. Possible values are: 'my_sets', 'x'."
  )
  expect_error(
    make_feature_sets_from_df(fake_sets_df, col_id = "x", col_set = "TEST"),
    "'col_set' argument: 'TEST' is not a column in 'annotation_df'. Possible values are: 'my_sets', 'x'."
  )
  expect_error(
    make_feature_sets_from_df(fake_sets_df, col_id = "x", col_set = "my_sets"),
    NA
  )

  res <- make_feature_sets_from_df(fake_sets_df, col_id = "x", col_set = "my_sets")
  expect_equal(res, fake_sets)
})


test_that("make_feature_sets_from_fm works", {
  multiomics_set <- test_get_multidataset()

  fmeta <- get_features_metadata(multiomics_set)

  ## Checking input
  expect_error(
    make_feature_sets_from_fm("TEST", "something"),
    "'mo_data' argument: Expecting MultiDataSet object"
  )
  expect_error(
    make_feature_sets_from_fm(multiomics_set, "test"),
    "'col_names' argument: 'test' is not a column in the features metadata for the snps\\+A dataset. Possible values are:.+"
  )
  expect_error(
    make_feature_sets_from_fm(multiomics_set, list("test")),
    "'col_names' argument should be named; names should be:"
  )
  expect_error(
    make_feature_sets_from_fm(multiomics_set, list("test" = "feature_id")),
    "'col_names' argument: 'test' are not existing datasets or do not match the 'datasets' argument. Possible values are:.+"
  )
  expect_error(
    make_feature_sets_from_fm(multiomics_set, list("snps+A" = "TEST")),
    "'col_names' argument: 'TEST' is not a column in the features metadata for the snps\\+A dataset. Possible values are:.+"
  )
  expect_error(
    make_feature_sets_from_fm(multiomics_set, list("snps+A" = "chromosome")),
    NA
  )
  expect_error(
    make_feature_sets_from_fm(multiomics_set, list("snps+A" = "chromosome", "rnaseq" = "chromosome", "metabolome" = "retention_time", "phenotypes" = "unit")),
    NA
  )

  ## Taking sets for one omics dataset only
  res <- make_feature_sets_from_fm(multiomics_set, list("snps+A" = "chromosome"))
  expect_equal(names(res), unique(fmeta$`snps+A`$chromosome))
  expect_equal(
    res[sort(names(res))] |> purrr::map(sort),
    unique(fmeta$`snps+A`$chromosome) |>
      sort() |>
      rlang::set_names() |>
      purrr::map(~ fmeta$`snps+A` |>
                   dplyr::filter(chromosome == .x) |>
                   dplyr::pull(feature_id) |>
                   sort())
  )

  ## Taking sets for two omics with same set names - no combining
  res <- make_feature_sets_from_fm(multiomics_set, list("snps+A" = "chromosome", "rnaseq" = "chromosome"), combine_omics_sets = FALSE)
  ## set names should have been made unique
  expect_equal(
    sort(names(res)),
    fmeta[c("snps+A", "rnaseq")] |>
      purrr::imap(~ paste0(unique(.x$chromosome), " - ", .y)) |>
      unlist() |>
      unname() |>
      sort()
  )

  true_ans <- fmeta[c("snps+A", "rnaseq")] |>
    purrr::map_dfr(
      ~ dplyr::select(.x, feature_id, chromosome) |>
        tibble::as_tibble(),
      .id = "dataset"
    ) |>
    dplyr::arrange(dataset, feature_id)

  expect_equal(
    res |>
      purrr::map_dfr(~ tibble::tibble(feature_id = .x),
                     .id = "dataset"
      ) |>
      dplyr::mutate(
        chromosome = stringr::str_extract(dataset, "^.+(?= - )"),
        dataset = stringr::str_remove(dataset, "^.+( - ){1}")
      ) |>
      dplyr::arrange(dataset, feature_id),
    true_ans
  )

  ## Taking sets for two omics with same set names - combining
  res <- make_feature_sets_from_fm(multiomics_set, list("snps+A" = "chromosome", "rnaseq" = "chromosome"), combine_omics_sets = TRUE)
  ## set names should have been made unique
  expect_equal(
    sort(names(res)),
    fmeta[c("snps+A", "rnaseq")] |>
      purrr::imap(~ unique(.x$chromosome)) |>
      unlist() |>
      unname() |>
      unique() |>
      sort()
  )

  true_ans <- fmeta[c("snps+A", "rnaseq")] |>
    purrr::map_dfr(~ dplyr::select(.x, feature_id, chromosome) |>
                     tibble::as_tibble()) |>
    dplyr::arrange(feature_id, chromosome)

  expect_equal(
    res |>
      purrr::imap_dfr(~ tibble::tibble(
        feature_id = .x,
        chromosome = .y
      )) |>
      dplyr::arrange(feature_id, chromosome),
    true_ans
  )
})


test_that("reduce_feature_sets_data works", {
  multiomics_set <- test_get_multidataset()

  fake_set <- list(
    "set1" = paste0("featureB_", 1:10),
    "set2" = c(
      paste0("featureB_", 11:15),
      paste0("featureC_", 1:5)
    )
  )

  ## Checking input
  expect_error(
    reduce_feature_sets_data(fake_set, "TEST"),
    "'mo_data' argument: Expecting MultiDataSet object"
  )
  expect_error(
    reduce_feature_sets_data(fake_set, multiomics_set),
    NA
  )

  expect_error(
    reduce_feature_sets_data(fake_set, multiomics_set, "TEST"),
    "'TEST' datasets are not present in mo_data. Possible dataset names are:.+"
  )

  ## When all features from set are present in multi-omics dataset
  expect_message(
    reduce_feature_sets_data(fake_set, multiomics_set),
    "All features in sets are in the multi-omics dataset."
  )
  expect_equal(
    reduce_feature_sets_data(fake_set, multiomics_set),
    fake_set
  )

  ## When not all features from the set are present in multi-omics dataset
  expect_message(
    reduce_feature_sets_data(fake_set, multiomics_set, datasets = "rnaseq"),
    "5 features in sets but not in the multi-omics dataset."
  )
  expect_message(
    reduce_feature_sets_data(fake_set, multiomics_set, datasets = "rnaseq"),
    "Removing 0 empty feature sets."
  )
  expect_equal(
    reduce_feature_sets_data(fake_set, multiomics_set, datasets = "rnaseq"),
    list(
      "set1" = paste0("featureB_", 1:10),
      "set2" = c(
        paste0("featureB_", 11:15)
      )
    )
  )

  ## When one entire set is not in the multi-omics dataset
  fake_set1 <- list(
    "set1" = paste0("featureB_", 1:10),
    "set2" = paste0("featureC_", 1:5)
  )
  expect_message(
    reduce_feature_sets_data(fake_set1, multiomics_set, datasets = "rnaseq"),
    "5 features in sets but not in the multi-omics dataset."
  )
  expect_message(
    reduce_feature_sets_data(fake_set1, multiomics_set, datasets = "rnaseq"),
    "Removing 1 empty feature sets."
  )
  expect_equal(
    reduce_feature_sets_data(fake_set1, multiomics_set, datasets = "rnaseq"),
    list(
      "set1" = paste0("featureB_", 1:10)
    )
  )
})

test_that("check_feature_sets works", {
  multiomics_set <- test_get_multidataset()
  fmeta <- get_features(multiomics_set)

  fake_sets <- list(
    "set1" = fmeta[['snps+A']][1:5],
    "set2" = fmeta[['rnaseq']][1:10]
  )

  nft <- n_features(multiomics_set) |> unname()
  nfu <- c(5, 10, 0, 0)

  ## Checking input
  expect_error(
    check_feature_sets(fake_sets, "TEST"),
    "'mo_data' argument: Expecting MultiDataSet object"
  )
  expect_no_error(check_feature_sets(fake_sets, multiomics_set))
  expect_error(
    check_feature_sets(fake_sets, multiomics_set, "TEST"),
    "'TEST' datasets are not present in mo_data. Possible dataset names are:.+"
  )
  expect_no_error(check_feature_sets(fake_sets, multiomics_set, c("rnaseq", "metabolome")))

  true_ans <- tibble::tibble(
    dataset = names(multiomics_set),
    n_annotated = nfu,
    n = nft,
    frac_annotated = n_annotated / n
  ) |>
    dplyr::mutate(
      dataset = factor(dataset, levels = names(multiomics_set)),
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
      )
    ) |>
    dplyr::select(-perc_annotated)

  expect_equal(
    check_feature_sets(fake_sets, multiomics_set),
    true_ans
  )
  expect_equal(
    check_feature_sets(fake_sets, multiomics_set, "metabolome"),
    true_ans |>
      dplyr::filter(dataset == "metabolome") |>
      dplyr::mutate(dataset = droplevels(dataset))
  )
  expect_equal(
    check_feature_sets(fake_sets, multiomics_set, c("rnaseq", "metabolome")),
    true_ans |>
      dplyr::filter(dataset %in% c("rnaseq", "metabolome")) |>
      dplyr::mutate(dataset = droplevels(dataset))
  )
})


test_that("evaluate_method_enrichment works", {
  multiomics_set <- test_get_multidataset()

  output_method <- test_get_diablo_run() |>
    get_output_diablo()
  output_method$features_weight <- output_method$features_weight |>
    dplyr::filter(!(feature_id %in% paste0("featureB_", 16:20)))
  fake_set <- list(
    "set1" = paste0("featureB_", 1:10),
    "set2" = paste0("featureB_", 11:20)
  )

  fake_set_info <- tibble::tibble(
    set_id = c("set1", "set2"),
    set_description = c("something", "something else")
  )

  ## checking input
  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "TEST"),
    "'TEST' are not valid dataset names for method DIABLO. Possible values are .+"
  )
  expect_error(
    evaluate_method_enrichment(output_method, fake_set),
    NA
  )
  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq"),
    NA
  )

  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", "TEST"),
    "'TEST' are not valid latent dimension names for method DIABLO. Possible values are 'Component 1', 'Component 2'."
  )
  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", "Component 1"),
    NA
  )

  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", sets_info_df = fake_set_info),
    "'col_set' argument should be the name of the column in 'sets_info_df' containing the feature sets ID."
  )
  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", sets_info_df = fake_set_info, col_set = "TEST"),
    "'col_set' argument: the value\\(s\\) 'TEST' are not column names in 'sets_info_df'. Possible values are: 'set_id', 'set_description'."
  )
  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", sets_info_df = fake_set_info, col_set = "set_id"),
    NA
  )

  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", add_missing_features = TRUE),
    "'add_missing_features' is TRUE, expecting a MultiDataSet object to be passed through the 'mo_data' argument."
  )
  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", add_missing_features = TRUE, mo_data = "TEST"),
    "'mo_data' argument: Expecting MultiDataSet object"
  )
  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", add_missing_features = TRUE, mo_data = multiomics_set[, c("snps+A", "metabolome")]),
    "'rnaseq' datasets are not present in mo_data. Possible dataset names are:.+"
  )
  expect_error(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", add_missing_features = TRUE, mo_data = multiomics_set),
    NA
  )

  ## Testing addition of missing features
  expect_equal(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", add_missing_features = FALSE)$set_size,
    rep(c(10, 5), 2)
  )
  expect_equal(
    evaluate_method_enrichment(output_method, fake_set, "rnaseq", add_missing_features = TRUE, mo_data = multiomics_set)$`set_size`,
    rep(10, 4)
  )

  ## Testing effet of min set size
  expect_equal(
    is.na(evaluate_method_enrichment(output_method, fake_set, "rnaseq", add_missing_features = FALSE, min_set_size = 5)$pvalue),
    rep(FALSE, 4)
  )
  expect_equal(
    is.na(evaluate_method_enrichment(output_method, fake_set, "rnaseq", add_missing_features = FALSE, min_set_size = 6)$pvalue),
    rep(c(FALSE, TRUE), each = 2)
  )
})


test_that("compute_samples_silhouette works", {

  if (!requireNamespace("cluster", quietly = TRUE)) {
    expect_error(transform_vsn(), "Package \"cluster\" must be installed to use this function.")
    skip("Package \"cluster\" must be installed to use this function.")
  }

  multiomics_set <- test_get_multidataset()

  output_method <- test_get_diablo_run() |>
    get_output_diablo()

  ## testing input
  expect_error(
    compute_samples_silhouette("test"),
    "'method_output': expecting an 'output_dimension_reduction' object (from get_output()).",
    fixed = TRUE
  )
  expect_error(
    compute_samples_silhouette(output_method, multiomics_set, "TEST"),
    "'group_col' argument: 'TEST' is not a column in the samples metadata of any dataset. Possible values are: 'id', 'name', 'pheno_group', 'time'."
  )
  expect_no_error(
    compute_samples_silhouette(output_method, multiomics_set, "pheno_group")
  )

  expect_error(
    suppressWarnings(compute_samples_silhouette(output_method, multiomics_set, "time")),
    "Silhouette width calculation failed. Check that there is more than one sample in at least one group."
  )

  w <- compute_samples_silhouette(output_method, multiomics_set, "time") |>
    capture_warning()

  expect_equal(
    w$message,
    "'group_col' argument: should indicate a column of samples grouping (character, factor or integer), not a numeric column. Converting column to integer."
  )

  res <- compute_samples_silhouette(output_method, multiomics_set, "pheno_group")

  expect_equal(names(res), c("samples_silhouette", "groups_average_silhouette"))
  expect_s3_class(res$samples_silhouette, "tbl_df")
  expect_s3_class(res$groups_average_silhouette, "tbl_df")

  expect_equal(
    names(res$samples_silhouette),
    c("sample_id", "group", "neighbour_group", "silhouette_width")
  )
  expect_equal(
    names(res$groups_average_silhouette),
    c("group", "group_average_width")
  )

  expect_equal(nrow(res$samples_silhouette), 5)
  expect_equal(nrow(res$groups_average_silhouette), 2)
})
