test_that("get_latent_dimensions works", {
  so2pls_res <- test_get_so2pls_run()
  output_ave <- get_output_so2pls(so2pls_res, use_average_dimensions = TRUE)
  ouput_noav <- get_output_so2pls(so2pls_res, use_average_dimensions = FALSE)

  expect_error(
    get_latent_dimensions("TEST"),
    "Expecting an object of class 'output_dimension_reduction' (from get_output()).",
    fixed = TRUE
  )

  expect_no_error(
    get_latent_dimensions(output_ave)
  )

  expect_no_error(
    get_latent_dimensions(ouput_noav)
  )

  expect_equal(
    get_latent_dimensions(output_ave),
    c("joint component 1", "joint component 2", "rnaseq specific component 1", "metabolome specific component 1")
  )

  expect_equal(
    get_latent_dimensions(ouput_noav),
    c("rnaseq joint component 1", "metabolome joint component 1", "rnaseq joint component 2", "metabolome joint component 2", "rnaseq specific component 1", "metabolome specific component 1")
  )

})

test_that("get_selected_features works", {
  multiomics_set <- test_get_multidataset()
  so2pls_output <- test_get_so2pls_run() |>
    get_output()

  expect_error(
    get_selected_features("TEST"),
    "Expecting an object of class 'output_dimension_reduction' (from get_output()).",
    fixed = TRUE
  )
  expect_no_error(get_selected_features(so2pls_output))
  expect_no_error(get_selected_features(so2pls_output, mo_data = multiomics_set))

  cols_list <- c("feature_id", "dataset", "latent_dimension", "weight", "importance")
  fmeta_cols <- get_features_metadata(multiomics_set)[c("rnaseq", "metabolome")] |>
    purrr::map(colnames) |>
    unlist() |>
    unname() |>
    unique()

  res <- get_selected_features(so2pls_output)
  expect_s3_class(res, "tbl_df")
  expect_equal(colnames(res), cols_list)
  expect_equal(
    nrow(res),
    3*2 + 5*2 + 35 + 40 ## 2 joint components, 1 specific comp per dataset
  )

  res2 <- get_selected_features(so2pls_output, mo_data = multiomics_set)
  expect_identical(
    res, dplyr::select(res2, tidyselect::all_of(cols_list))
  )
  expect_equal(colnames(res2), union(cols_list, fmeta_cols))

  ## Trying with sPLS
  spls_res <- test_get_spls_run()
  spls_output <- get_output(spls_res)

  res <- get_selected_features(spls_output)

  exp <- c(1, 2) |>
    purrr::map(
      ~ mixOmics::selectVar(spls_res, comp = .x)
    ) |>
    purrr::imap_dfr(
      ~ .x |>
        purrr::map_dfr(purrr::pluck, "value") |>
        tibble::as_tibble(rownames = "feature_id"),
      .id = "latent_dimension"
    ) |>
    dplyr::rename(weight = value.var) |>
    dplyr::mutate(
      latent_dimension = paste0("Component ", latent_dimension),
      latent_dimension = factor(latent_dimension),
      dataset = stringr::str_extract(feature_id, "\\w(?=_)")
    ) |>
    dplyr::arrange(latent_dimension, dataset, dplyr::desc(abs(weight))) |>
    dplyr::select(-dataset)

  expect_identical(
    dplyr::select(res, latent_dimension, feature_id, weight),
    exp
  )
})

test_that("get_top_features works", {
  multiomics_set <- test_get_multidataset()
  diablo_output <- test_get_diablo_run() |>
    get_output()

  expect_error(
    get_top_features("TEST"),
    "Expecting an object of class 'output_dimension_reduction' (from get_output()).",
    fixed = TRUE
  )

  expect_no_error(get_top_features(diablo_output, n_features = 10))
  expect_no_error(get_top_features(diablo_output, min_importance = 0.8))
  expect_no_error(get_top_features(diablo_output, n_features = 10, mo_data = multiomics_set))

  ## expected column names
  cols_list <- c("feature_id", "dataset", "latent_dimension", "weight", "importance")
  fmeta_cols <- get_features_metadata(multiomics_set)[c("snps+A", "rnaseq", "metabolome")] |>
    purrr::map(colnames) |>
    unlist() |>
    unname() |>
    unique()

  ## Testing n_features
  res_n <- get_top_features(diablo_output, n_features = 10)
  expect_s3_class(res_n, "tbl_df")
  expect_equal(colnames(res_n), cols_list)
  expect_equal(nrow(res_n), 10*3*2)

  ## Testing min_importance
  res_m <- get_top_features(diablo_output, min_importance = 0.8)
  expect_equal(colnames(res_m), cols_list)
  expect_true(min(res_m$importance) >= 0.8)

  ## Testing filtering
  res_f <- get_top_features(
    diablo_output,
    n_features = 10,
    latent_dimensions = "Component 2",
    datasets = c("rnaseq", "metabolome")
  )
  expect_equal(unique(as.character(res_f$latent_dimension)), "Component 2")
  expect_equal(unique(as.character(res_f$dataset)), c("rnaseq", "metabolome"))

  ## Testing features metadata
  res_fm <- get_top_features(diablo_output, n_features = 10, mo_data = multiomics_set)
  expect_equal(
    dplyr::select(res_fm, tidyselect::all_of(cols_list)),
    res_n
  )
  expect_equal(colnames(res_fm), union(cols_list, fmeta_cols))
})

test_that(".check_names_output_list works", {
  res_list <- test_get_output_list()

  expect_error(.check_names_output_list(res_list), NA)
  expect_equal(
    names(.check_names_output_list(res_list)),
    c("DIABLO", "MOFA", "sO2PLS")
  )

  ## More than one output from the same integration method
  expect_error(
    .check_names_output_list(list(res_list[[1]], res_list[[1]])),
    "'output_list' contains several results from a same integration method.+"
  )
  expect_error(
    .check_names_output_list(list("R1" = res_list[[1]], "R2" = res_list[[1]])),
    NA
  )

  ## Custom names
  expect_error(
    .check_names_output_list(list("R1" = res_list[[1]], "R1" = res_list[[2]])),
    "'output_list' has duplicated names. Please use unique names."
  )
  expect_error(
    .check_names_output_list(list("R1" = res_list[[1]], "R2" = res_list[[2]])),
    NA
  )
})

test_that(".filter_output_dimensions works", {
  ## Input
  res_list <- test_get_output_list()
  good_ex <- paste0("joint component ", 1:2)


  expect_error(
    .filter_output_dimensions("TEST"),
    "Expecting an object of class 'output_dimension_reduction' (from get_output()).",
    fixed = TRUE
  )

  expect_no_error(.filter_output_dimensions(res_list[[3]], good_ex))


  ## Checking values of latent_dimensions list
  expect_error(
    .filter_output_dimensions(res_list[[1]], c("TEST", "TEST2")),
    "'TEST', 'TEST2' are not valid latent dimension names for method DIABLO. Possible values are 'Component 1', 'Component 2'."
  )
  expect_error(
    .filter_output_dimensions(res_list[[1]], c("TEST", "TEST2"), method_name = "THING"),
    "'TEST', 'TEST2' are not valid latent dimension names for method THING. Possible values are 'Component 1', 'Component 2'."
  )
  expect_error(
    .filter_output_dimensions(res_list[[1]], "Component 1"),
    NA
  )

  ## Checking results
  res <- .filter_output_dimensions(res_list[[1]], "Component 1")
  expect_equal(
    levels(res$features_weight$latent_dimension),
    "Component 1"
  )
  expect_equal(
    levels(res$samples_score$latent_dimension),
    "Component 1"
  )

  ## Optional tests wok
  expect_error(
    .filter_output_dimensions(res_list[[3]], paste0("joint component ", 1:2), fixed_length = 1),
    "'latent_dimensions' should be of length 1."
  )
  expect_error(
    .filter_output_dimensions(res_list[[3]], paste0("joint component ", 1:2), fixed_length = 2),
    NA
  )
})


test_that(".filter_output_dimensions_list works", {
  ## Input
  res_list <- test_get_output_list()
  res_list_named <- .check_names_output_list(res_list)
  good_ex <- list("sO2PLS" = paste0("joint component ", 1:2))

  expect_no_error(.filter_output_dimensions_list(res_list_named, good_ex))

  ## Checking names of latent_dimensions list
  expect_error(
    .filter_output_dimensions_list(res_list_named, list(paste0("joint component ", 1:2))),
    "'latent_dimensions' argument should be a named list, names should be: 'DIABLO', 'MOFA', 'sO2PLS'."
  )
  expect_error(
    .filter_output_dimensions_list(res_list_named, list("TEST" = paste0("joint component ", 1:2))),
    "'latent_dimensions' argument: list names 'TEST' do not match names from 'output_list'. Possible values are 'DIABLO', 'MOFA', 'sO2PLS'."
  )
  expect_no_error(.filter_output_dimensions_list(res_list_named, list("sO2PLS" = paste0("joint component ", 1:2))))

  res_list_named2 <- res_list_named
  names(res_list_named2) <- paste0("Method ", 1:3)
  expect_error(
    .filter_output_dimensions_list(res_list_named2, list("DIABLO" = paste0("joint component ", 1:2))),
    "'latent_dimensions' argument: list names 'DIABLO' do not match names from 'output_list'. Possible values are 'Method 1', 'Method 2', 'Method 3'."
  )

  ## Checking values of latent_dimensions list
  expect_error(
    .filter_output_dimensions_list(res_list_named, list("DIABLO" = c("TEST", "TEST2"))),
    "'TEST', 'TEST2' are not valid latent dimension names for method DIABLO. Possible values are 'Component 1', 'Component 2'."
  )
  expect_no_error(.filter_output_dimensions_list(res_list_named, list("DIABLO" = c("Component 1"))))
  expect_no_error(
    .filter_output_dimensions_list(
      res_list_named,
      list("DIABLO" = c("Component 1"), "sO2PLS" = paste0("joint component ", 1:2)))
  )

  ## Checking results
  res <- .filter_output_dimensions_list(res_list_named, list("DIABLO" = c("Component 1")))
  expect_equal(
    levels(res$DIABLO$features_weight$latent_dimension),
    "Component 1"
  )
  expect_equal(
    levels(res$DIABLO$samples_score$latent_dimension),
    "Component 1"
  )

  ## works with several filtering
  res <- .filter_output_dimensions_list(res_list_named, list("DIABLO" = c("Component 1"), "sO2PLS" = paste0("joint component ", 1:2)))
  expect_equal(
    levels(res$DIABLO$features_weight$latent_dimension),
    "Component 1"
  )
  expect_equal(
    levels(res$sO2PLS$features_weight$latent_dimension),
    paste0("joint component ", 1:2)
  )


  ## Optional tests wok
  expect_error(
    .filter_output_dimensions_list(res_list_named, list("DIABLO" = c("Component 1"), "sO2PLS" = paste0("joint component ", 1:2)), all_present = TRUE),
    "'latent_dimensions' argument should match 'output_list' elements; missing value for 'MOFA'."
  )
  expect_error(
    .filter_output_dimensions_list(res_list_named, list("DIABLO" = c("Component 1"), "sO2PLS" = paste0("joint component ", 1:2), "MOFA" = paste0("Factor 1")), all_present = TRUE),
    NA
  )

  expect_error(
    .filter_output_dimensions_list(res_list_named, list("DIABLO" = c("Component 1"), "sO2PLS" = paste0("joint component ", 1:2)), fixed_length = 1),
    "Elements in 'latent_dimensions' list should be of length 1."
  )
  expect_error(
    .filter_output_dimensions_list(res_list_named, list("DIABLO" = c("Component 1"), "sO2PLS" = c("joint component 1")), fixed_length = 1),
    NA
  )
})


test_that(".filter_output_datasets works", {
  ## Input
  res_list <- test_get_output_list()
  good_ex <- c("snps+A", "rnaseq")

  expect_error(
    .filter_output_datasets("TEST"),
    "Expecting an object of class 'output_dimension_reduction' (from get_output()).",
    fixed = TRUE
  )

  ## Checking values of datasets list
  expect_error(
    .filter_output_datasets(res_list[[1]], c("TEST", "TEST2")),
    "'TEST', 'TEST2' are not valid dataset names for method DIABLO. Possible values are 'snps\\+A', 'rnaseq', 'metabolome'."
  )
  expect_error(
    .filter_output_datasets(res_list[[1]], c("TEST", "TEST2"), method_name = "THING"),
    "'TEST', 'TEST2' are not valid dataset names for method THING. Possible values are 'snps\\+A', 'rnaseq', 'metabolome'."
  )
  expect_error(
    .filter_output_datasets(res_list[[1]], "snps+A"),
    NA
  )

  ## Checking results
  expect_equal(
    .filter_output_datasets(res_list[[1]], NULL),
    res_list[[1]]
  )

  res <- .filter_output_datasets(res_list[[1]], "snps+A")
  expect_equal(
    levels(res$features_weight$dataset),
    "snps+A"
  )

  res <- .filter_output_datasets(res_list[[1]], "rnaseq")
  expect_equal(
    levels(res$features_weight$dataset),
    "rnaseq"
  )

  res <- .filter_output_datasets(res_list[[1]], c("snps+A", "rnaseq"))
  expect_equal(
    levels(res$features_weight$dataset),
    c("snps+A", "rnaseq")
  )

  ## Optional tests work
  expect_error(
    .filter_output_datasets(res_list[[1]], c("snps+A", "rnaseq"), fixed_length = 1),
    "'datasets' should be of length 1."
  )
  expect_error(
    .filter_output_datasets(res_list[[1]], c("snps+A", "rnaseq"), fixed_length = 2),
    NA
  )
})
