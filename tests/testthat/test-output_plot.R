test_that("plot_samples_score works", {
  output_diablo <- test_get_diablo_run() |>
    get_output_diablo()
  multiomics_set <- test_get_multidataset()

  expect_error(
    plot_samples_score(output_diablo, shape_upper = "test"),
    "Need to provide a MultiDataSet object through 'mo_data' argument in order to use one of the aesthetics arguments."
  )
  expect_error(
    plot_samples_score(output_diablo, mo_data = multiomics_set, shape_upper = "test"),
    "'shape_upper' argument: 'test' is not a column in the samples metadata of any dataset. Possible values are:.+"
  )
  expect_no_error(
    plot_samples_score(output_diablo, mo_data = multiomics_set, shape_upper = "pheno_group")
  )
  expect_no_warning(
    plot_samples_score(output_diablo, latent_dimensions = "Component 1", mo_data = multiomics_set, shape_upper = "pheno_group")
  )
  expect_warning(
    plot_samples_score(output_diablo, latent_dimensions = "Component 1", mo_data = multiomics_set, shape_lower = "pheno_group"),
    "Only one latent dimension to plot; 'colour_diag', 'colour_lower' and 'shape_lower' argument will be ignored."
  )
})

test_that("plot_samples_score_pair works", {
  multiomics_set <- test_get_multidataset()
  output_diablo <- test_get_diablo_run() |>
    get_output_diablo()
  output_mofa <- test_get_mofa_run() |>
    get_output_mofa2()

  output_list <- list(output_diablo, output_mofa)

  expect_error(plot_samples_score_pair("TEST"), "'method_output' argument: expecting either a single 'output_dimension_reduction' object or a list of two such objects.")
  expect_error(plot_samples_score_pair(list(output_diablo, output_diablo, output_diablo)), "'method_output' argument: expecting either a single 'output_dimension_reduction' object or a list of two such objects.")

  ## Testing single output
  expect_error(plot_samples_score_pair(output_diablo, "TEST"), "'latent_dimensions' should be of length 2.")
  expect_error(plot_samples_score_pair(output_diablo, paste0("TEST ", 1:2)), "'TEST 1', 'TEST 2' are not valid latent dimension names for method DIABLO. Possible values are 'Component 1', 'Component 2'.")
  expect_no_error(plot_samples_score_pair(output_diablo, paste0("Component ", 1:2)))

  ## Testing two output
  expect_error(plot_samples_score_pair(output_list, "TEST"), "'latent_dimensions' argument should be a named list, names should be: 'DIABLO', 'MOFA'.")
  expect_error(plot_samples_score_pair(output_list, list("DIABLO" = "TEST", "MOFA" = "TEST")), "'TEST' are not valid latent dimension names for method DIABLO. Possible values are 'Component 1', 'Component 2'.")
  expect_no_error(plot_samples_score_pair(output_list, list("DIABLO" = "Component 1", "MOFA" = "Factor 1")))

  ## Testing with samples metadata
  expect_error(
    plot_samples_score_pair(output_diablo, paste0("Component ", 1:2), colour_by = "TEST"),
    "Need to provide a MultiDataSet object through 'mo_data' argument in order to use one of the aesthetics arguments."
  )
  expect_error(
    plot_samples_score_pair(output_diablo, paste0("Component ", 1:2), mo_data = multiomics_set, colour_by = "TEST"),
    "'colour_by' argument: 'TEST' is not a column in the samples metadata of any dataset. Possible values are:.+"
  )
  expect_no_error(
    plot_samples_score_pair(output_diablo, paste0("Component ", 1:2), mo_data = multiomics_set, colour_by = "time")
  )
  expect_no_error(
    plot_samples_score_pair(output_diablo, paste0("Component ", 1:2), mo_data = multiomics_set, colour_by = "pheno_group")
  )
})

test_that("plot_samples_score_covariate works", {
  multiomics_set <- test_get_multidataset()
  output_diablo <- test_get_diablo_run() |>
    get_output_diablo()

  expect_error(
    plot_samples_score_covariate("TEST"),
    "Expecting an object of class 'output_dimension_reduction' (from get_output()).",
    fixed = TRUE
  )
  expect_error(
    plot_samples_score_covariate(output_diablo, "TEST"),
    "'mo_data' argument: Expecting MultiDataSet object"
  )
  expect_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "TEST"),
    "'covariate' argument: 'TEST' is not a column in the samples metadata of any dataset. Possible values are:.+"
  )
  expect_no_error( ## continuous covariate
    plot_samples_score_covariate(output_diablo, multiomics_set, "time")
  )
  expect_no_error( ## discrete covariate
    plot_samples_score_covariate(output_diablo, multiomics_set, "pheno_group")
  )

  ## Testing continuous colour_by
  expect_no_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "time", colour_by = "time")
  )
  expect_no_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "pheno_group", colour_by = "time")
  )

  ## Testing discrete colour_by
  expect_no_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "time", colour_by = "pheno_group")
  )
  expect_no_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "pheno_group", colour_by = "pheno_group")
  )

  ## Testing shape_by
  expect_no_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "time", shape_by = "pheno_group")
  )
  expect_no_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "pheno_group", shape_by = "pheno_group")
  )

  ## Testing both
  expect_no_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "time", colour_by = "time", shape_by = "pheno_group")
  )
  expect_no_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "pheno_group", colour_by = "time", shape_by = "pheno_group")
  )
  expect_no_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "time", colour_by = "pheno_group", shape_by = "pheno_group")
  )
  expect_no_error(
    plot_samples_score_covariate(output_diablo, multiomics_set, "pheno_group", colour_by = "pheno_group", shape_by = "pheno_group")
  )
})

test_that("plot_features_weight_pair works - with one method", {
  diablo_output <- test_get_diablo_run() |>
    get_output()

  expect_no_error(plot_features_weight_pair(diablo_output, paste0("Component ", 1:2)))
  expect_no_error(plot_features_weight_pair(diablo_output, paste0("Component ", 1:2), datasets = "rnaseq"))
})

test_that("plot_features_weight_pair works - with two different methods", {
  res_list <- test_get_output_list()
  multiomics_set <- test_get_multidataset()

  ld_list <- list("DIABLO" = "Component 1", "MOFA" = "Factor 1")

  expect_error(
    plot_features_weight_pair(res_list, ld_list),
    "'method_output' should be of length 2; this function is for comparing the output of two integration methods."
  )
  expect_no_error(plot_features_weight_pair(res_list[1:2], ld_list))
  expect_no_error(plot_features_weight_pair(res_list[1:2], ld_list, datasets = "rnaseq"))
  expect_no_error(
    plot_features_weight_pair(
      res_list[1:2],
      ld_list,
      label_cols = list("snps+A" = "chromosome", "rnaseq" = "name", "metabolome" = "name"),
      mo_data = multiomics_set
    )
  )

  expect_true(
    all(stringr::str_detect(
      plot_features_weight_pair(res_list[1:2], ld_list)$data |>
        dplyr::filter(is_top, dataset == "rnaseq") |>
        dplyr::pull(label),
      "^featureB"
    ))
  )
  expect_true(
    !all(stringr::str_detect(
      plot_features_weight_pair(res_list[1:2],
                                ld_list,
                                label_cols = list("snps+A" = "chromosome", "rnaseq" = "name", "metabolome" = "name"),
                                mo_data = multiomics_set
      )$data |>
        dplyr::filter(is_top, dataset == "rnaseq") |>
        dplyr::pull(label),
      "^featureB"
    ))
  )
})

test_that("plot_features_weight_covariate works", {
  multiomics_set <- test_get_multidataset()
  diablo_output <- test_get_diablo_run() |> get_output()
  ds <- levels(diablo_output$features_weight$dataset)

  for (i in seq(length(ds))) {
    j <- ds[i]
    n <- n_features(multiomics_set)[[j]]

    df <- data.frame(
      feature_id = get_features(multiomics_set)[[j]],
      dataset = j,
      cov_char_common = sample(LETTERS[1:3], n, replace = TRUE),
      cov_num_common = runif(n),
      cov_char = sample(LETTERS[(5 * i):(5 * i + 3)], n, replace = TRUE),
      cov_num = runif(n)
    )
    names(df)[5:6] <- paste0(names(df)[5:6], i)

    multiomics_set <- add_features_metadata(multiomics_set, df)
  }

  cov_list <- paste0("cov_char", seq(length(ds))) |>
    rlang::set_names(ds) |>
    as.list()

  ## Checking input (other than features metadata columns)
  expect_error(
    plot_features_weight_covariate("TEST", multiomics_set, cov_list),
    "Expecting an object of class 'output_dimension_reduction' (from get_output()).",
    fixed = TRUE
  )
  expect_error(
    plot_features_weight_covariate(diablo_output, "TEST", cov_list),
    "'mo_data' argument: Expecting MultiDataSet object"
  )
  expect_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, cov_list, features_metric = "TEST"),
    "`features_metric` must be one of \"signed_importance\", \"weight\", or \"importance\", not \"TEST\"."
  )
  expect_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, cov_list, latent_dimensions = "TEST"),
    "'TEST' are not valid latent dimension names for method DIABLO. Possible values are 'Component 1', 'Component 2'."
  )

  expect_no_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, cov_list)
  )

  ## Testing covariate errors
  expect_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, list("TEST" = "cov_char_common")),
    "'covariate' argument: 'TEST' are not existing datasets or do not match the 'datasets' argument. Possible values are: 'snps\\+A', 'rnaseq', 'metabolome'."
  )

  expect_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, list("snps+A" = "TEST")),
    "'covariate' argument: 'TEST' is not a column in the features metadata for the snps\\+A dataset. Possible values are: 'feature_id', 'chromosome', 'position', 'cov_char_common', 'cov_num_common', 'cov_char1', 'cov_num1'."
  )

  expect_no_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, list("snps+A" = "cov_char_common"))
  )
  res <- plot_features_weight_covariate(
    diablo_output,
    multiomics_set,
    list("snps+A" = "cov_char_common", "rnaseq" = "cov_char_common")
  )
  expect_equal(
    res$labels$x,
    "cov_char_common"
  )
  res <- plot_features_weight_covariate(
    diablo_output,
    multiomics_set,
    list("snps+A" = "cov_char1", "rnaseq" = "cov_char2")
  )
  expect_equal(
    res$labels$x,
    "cov_char1 (snps+A), cov_char2 (rnaseq)"
  )

  ## Testing colour_by
  expect_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, cov_list, colour_by = list("TEST" = "cov_char_common")),
    "'colour_by' argument: 'TEST' are not existing datasets or do not match the 'datasets' argument. Possible values are: 'snps\\+A', 'rnaseq', 'metabolome'."
  )
  expect_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, cov_list, colour_by = list("snps+A" = "TEST")),
    "'colour_by' argument: 'TEST' is not a column in the features metadata for the snps\\+A dataset. Possible values are: 'feature_id', 'chromosome', 'position', 'cov_char_common', 'cov_num_common', 'cov_char1', 'cov_num1'."
  )
  expect_error(
    plot_features_weight_covariate(
      diablo_output, multiomics_set, list("snps+A" = "cov_char1"), colour_by = list("rnaseq" = "cov_char_common")
    ),
    "'colour_by' argument: 'rnaseq' are not existing datasets or do not match the 'datasets' argument. Possible values are: 'snps\\+A'."
  )
  expect_no_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, cov_list, colour_by = list("snps+A" = "cov_char_common"))
  )

  ## Testing shape_by
  expect_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, cov_list, shape_by = list("TEST" = "cov_char_common")),
    "'shape_by' argument: 'TEST' are not existing datasets or do not match the 'datasets' argument. Possible values are: 'snps\\+A', 'rnaseq', 'metabolome'."
  )
  expect_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, cov_list, shape_by = list("snps+A" = "TEST")),
    "'shape_by' argument: 'TEST' is not a column in the features metadata for the snps\\+A dataset. Possible values are: 'feature_id', 'chromosome', 'position', 'cov_char_common', 'cov_num_common', 'cov_char1', 'cov_num1'."
  )
  expect_error(
    plot_features_weight_covariate(
      diablo_output, multiomics_set, list("snps+A" = "cov_char1"), shape_by = list("rnaseq" = "cov_char_common")
    ),
    "'shape_by' argument: 'rnaseq' are not existing datasets or do not match the 'datasets' argument. Possible values are: 'snps\\+A'."
  )
  expect_no_error(
    plot_features_weight_covariate(diablo_output, multiomics_set, cov_list, shape_by = list("snps+A" = "cov_char_common"))
  )

  ## Checking the merger of variable names
  res <- plot_features_weight_covariate(
    diablo_output,
    multiomics_set,
    cov_list,
    shape_by = list("snps+A" = "cov_char_common")
  )
  expect_equal(res$labels$x, "cov_char1 (snps+A), cov_char2 (rnaseq), cov_char3 (metabolome)")
  expect_equal(res$labels$shape, "cov_char_common (snps+A)")

  res <- plot_features_weight_covariate(
    diablo_output,
    multiomics_set,
    cov_list,
    shape_by = list("snps+A" = "cov_char_common", "rnaseq" = "cov_char_common", "metabolome" = "cov_char_common")
  )
  expect_equal(res$labels$x, "cov_char1 (snps+A), cov_char2 (rnaseq), cov_char3 (metabolome)")
  expect_equal(res$labels$shape, "cov_char_common")
})
