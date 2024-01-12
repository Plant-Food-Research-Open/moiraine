test_that("get_input_mofa2 works - errors on input", {
  multiomics_set <- test_get_multidataset()
  names_ds <- c("snps+A", "rnaseq", "metabolome")

  ## Testing that wrong inputs return an error
  expect_error(
    get_input_mofa2("test", names_ds, NULL, NULL, NULL, FALSE),
    "Expecting MultiDataSet object"
  )

  expect_error(
    get_input_mofa2(multiomics_set, "test", NULL, NULL, NULL, FALSE),
    "'test' datasets are not present in mo_data. Possible dataset names are:.+"
  )

  expect_error(
    get_input_mofa2(multiomics_set, names_ds, "test", NULL, NULL, FALSE),
    "In 'covariates' argument: 'test' are not columns in the samples metadata\\. Possible values are: .+"
  )

  expect_error(
    get_input_mofa2(multiomics_set, names_ds, NULL, "test", NULL, FALSE),
    "Parameter 'groups': 'test' is not a column in the samples information data-frame\\. Possibles values for the 'groups' parameter are: .+"
  )

  expect_error(
    get_input_mofa2(multiomics_set, names_ds, NULL, NULL, "test", FALSE),
    "Argument 'options_list' should be a named list with length <= 3; possible names: 'data_options', 'model_options', 'training_options'\\."
  )
  expect_error(
    get_input_mofa2(multiomics_set, names_ds, "pheno_group", NULL, "test"),
    "Argument 'options_list' should be a named list with length <= 4; possible names: 'data_options', 'model_options', 'training_options', 'mefisto_options'\\."
  )

  expect_error(
    get_input_mofa2(multiomics_set, names_ds, NULL, NULL, list("test"), FALSE),
    "Argument 'options_list' should be a named list; possible names: 'data_options', 'model_options', 'training_options'\\."
  )
  expect_error(
    get_input_mofa2(multiomics_set, names_ds, "pheno_group", NULL, list("test"), FALSE),
    "Argument 'options_list' should be a named list; possible names: 'data_options', 'model_options', 'training_options', 'mefisto_options'\\."
  )

  expect_error(
    get_input_mofa2(multiomics_set, names_ds, NULL, NULL, list(test = "test"), FALSE),
    "The following names in 'options_list' arguments are not valid: 'test'. Possible names are: .+"
  )

  expect_error(
    get_input_mofa2(multiomics_set, names_ds, "pheno_group", NULL, NULL, FALSE),
    "'covariates' argument should correspond to numerical variables in the samples metadata."
  )

  expect_error(
    suppressMessages(get_input_mofa2(multiomics_set, names_ds, NULL, NULL, list("model_options" = list("test" = 1)), FALSE)),
    "Argument 'options_list': the following names from options_list\\$model_options are not model options parameters in MOFA2: 'test'.\nPossible names are: .+"
  )

  ## This should be correct:
  expect_error(
    suppressMessages(get_input_mofa2(multiomics_set, names_ds, NULL, NULL, NULL, FALSE)),
    NA
  )
})

test_that("get_input_mofa2 works - samples metadata", {
  multiomics_set <- test_get_multidataset()
  names_ds <- c("snps+A", "rnaseq", "metabolome")

  multiomics_set2 <- multiomics_set
  multiomics_set2@phenoData$rnaseq$main$name <- NA

  ## Datasets with conflicting samples metadata information
  expect_error(
    get_input_mofa2(multiomics_set2, names(multiomics_set2)[-1], NULL, NULL, NULL, FALSE),
    "Conflicting information in samples metadata for samples .+"
  )

  ## Getting input samples metadata for the checks
  res1 <- suppressMessages(get_input_mofa2(multiomics_set, names_ds, NULL, NULL, NULL, FALSE))
  smeta <- get_samples_metadata(multiomics_set)[names_ds]

  ## Checking dimensions and dim names
  expect_equal(setdiff(rownames(res1@samples_metadata), "sample"), setdiff(unique(unlist(lapply(smeta, rownames))), "id"))
  expect_equal(colnames(res1@samples_metadata), c("sample", setdiff(unique(unlist(lapply(smeta, colnames))), "id"), "group"))

  ## check that the information has been conserved (actually not super informative since everything is the same)
  expect_equal(res1@samples_metadata[rownames(smeta[[1]]), setdiff(colnames(smeta[[1]]), "id")], smeta[[1]][, -which(colnames(smeta[[1]]) == "id")])
  expect_equal(res1@samples_metadata$sample, rownames(res1@samples_metadata))

  ## Check with dates
  multiomics_set3 <- multiomics_set
  multiomics_set3@phenoData$rnaseq$main$date <- rep(
    as.POSIXct("2023-01-01"),
    length(multiomics_set3@phenoData$rnaseq$main$time)
  )

  res3 <- suppressMessages(get_input_mofa2(multiomics_set3, names_ds, NULL, NULL, NULL, FALSE))
  expect_true(
    is.character(res3@samples_metadata$date)
  )
})

test_that("get_input_mofa2 works - features metadata", {
  multiomics_set <- test_get_multidataset()
  names_ds <- c("snps+A", "rnaseq", "metabolome")

  ## Getting input samples metadata for the checks
  res1 <- suppressMessages(get_input_mofa2(multiomics_set, names_ds, NULL, NULL, NULL, FALSE))
  fmeta <- get_features_metadata(multiomics_set)[names_ds]

  ## Checking dimensions and dim names
  expect_equal(rownames(res1@features_metadata), unname(unlist(lapply(fmeta, rownames))))
  expect_equal(sort(colnames(res1@features_metadata)), sort(setdiff(c(unique(unlist(lapply(fmeta, colnames))), "feature", "view"), "feature_id")))

  ## check that the information has been conserved (actually not super informative since everything is the same)
  expect_equal(res1@features_metadata[rownames(fmeta[[1]]), setdiff(colnames(fmeta[[1]]), "feature_id")], fmeta[[1]][, -1])
  expect_equal(res1@features_metadata$feature, rownames(res1@features_metadata))
  expect_equal(res1@features_metadata$view, rep(names_ds, times = sapply(fmeta, nrow)))

  ## Check with dates
  multiomics_set3 <- multiomics_set
  multiomics_set3@featureData$rnaseq$main$date <- rep(
    as.POSIXct("2023-01-01"),
    length(multiomics_set3@featureData$rnaseq$main$name)
  )

  res3 <- suppressMessages(get_input_mofa2(multiomics_set3, names_ds, NULL, NULL, NULL, FALSE))
  expect_true(
    is.character(res3@features_metadata$date)
  )
})

test_that("get_input_mofa2 works - likelihood checks", {
  multiomics_set <- test_get_multidataset()
  names_ds <- c("snps+A", "rnaseq", "metabolome")

  expect_error(
    get_input_mofa2(multiomics_set, names_ds, NULL, NULL, list(model_options = list(likelihoods = "test")), FALSE),
    "In 'options_list' parameter: options_list\\$model_options\\$likelihoods should have length 3 and names: .+"
  )

  expect_error(
    get_input_mofa2(multiomics_set, names_ds, NULL, NULL, list(model_options = list(likelihoods = c("test" = "test", "rnaseq" = "test", "metabolome" = "test"))), FALSE),
    "In 'options_list' parameter: options_list\\$model_options\\$likelihoods should have length 3 and names: .+"
  )

  expect_error(
    get_input_mofa2(multiomics_set, names_ds, NULL, NULL, list(model_options = list(likelihoods = c("snps+A" = "test", "rnaseq" = "test", "metabolome" = "test"))), FALSE),
    "In 'options_list' parameter: possible values for options_list\\$model_options\\$likelihoods are: 'gaussian','bernoulli', 'poisson'."
  )

  expect_warning(
    suppressMessages(get_input_mofa2(multiomics_set, names_ds, NULL, NULL, list(model_options = list(likelihoods = c("snps+A" = "gaussian", "rnaseq" = "gaussian", "metabolome" = "poisson"))), FALSE)),
    "Dataset metabolome is to be modelled with a poisson likelihood, but is not integer. Transforming to integer."
  )

  expect_error(
    get_input_mofa2(multiomics_set, names_ds, NULL, NULL, list(model_options = list(likelihoods = c("snps+A" = "gaussian", "rnaseq" = "gaussian", "metabolome" = "bernoulli"))), FALSE),
    "Dataset metabolome is to be modelled with a bernoulli likelihood, but contains values other than 0 and 1. Please transform the dataset or change the likelihood to be used."
  )
})


test_that("get_input_mofa2 works - group column", {
  multiomics_set <- test_get_multidataset()
  names_ds <- c("snps+A", "rnaseq", "metabolome")

  res1 <- suppressMessages(get_input_mofa2(multiomics_set, names_ds, NULL, NULL, NULL, FALSE))

  ## renaming the 'pheno_group' column in samples metadata as 'group'
  for (i in names_ds) {
    Biobase::varLabels(multiomics_set@phenoData[[i]]$main) <- stringr::str_replace(
      Biobase::varLabels(multiomics_set@phenoData[[i]]$main),
      "pheno_group",
      "group"
    )
  }

  expect_warning(
    suppressMessages(get_input_mofa2(multiomics_set, names_ds, NULL, NULL, NULL, FALSE)),
    "Renaming 'group' column in samples metadata to 'group_metadata'."
  )

  res2 <- suppressWarnings(suppressMessages(get_input_mofa2(multiomics_set, names_ds, NULL, NULL, NULL, FALSE)))
  expect_equal(res2@samples_metadata$group_metadata, res1@samples_metadata$pheno_group)

  res3 <- suppressWarnings(suppressMessages(get_input_mofa2(multiomics_set, names_ds, NULL, "group", NULL, FALSE)))
  expect_equal(res3@samples_metadata$group, res3@samples_metadata$group_metadata)
})


test_that("get_input_mofa2 works - covariates", {
  multiomics_set <- test_get_multidataset()
  names_ds <- c("snps+A", "rnaseq", "metabolome")

  res1 <- suppressMessages(get_input_mofa2(multiomics_set, names_ds, "time", NULL, NULL, FALSE))
  time_vals <- get_samples_metadata(multiomics_set)[names_ds] |>
    purrr::reduce(rbind) |>
    dplyr::distinct() |>
    dplyr::select(time) |>
    t()

  expect_equal(res1@covariates$group1, time_vals)

  ## Covariates and groups should work together
  expect_error(get_input_mofa2(multiomics_set, names_ds, "time", "pheno_group", NULL, FALSE), NA)
})

test_that("get_input_mofa2 works - common samples", {
  multiomics_set <- test_get_multidataset()
  names_ds <- c("snps+A", "rnaseq", "metabolome")

  res1 <- suppressMessages(get_input_mofa2(multiomics_set, names_ds, "time", NULL, NULL, only_common_samples = FALSE))

  expect_equal(res1@dimensions$N[[1]], dplyr::n_distinct(unlist(purrr::map(names_ds, ~ get_samples_metadata(multiomics_set)[[.x]]$id))))

  res2 <- suppressWarnings(suppressMessages(get_input_mofa2(multiomics_set, names_ds, "time", NULL, NULL, only_common_samples = TRUE)))

  expect_equal(res2@dimensions$N[[1]], length(MultiDataSet::commonIds(multiomics_set[, names_ds])))
})
