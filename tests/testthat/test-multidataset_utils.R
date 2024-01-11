test_that("get_datasets works", {
  multiomics_set <- test_get_multidataset()

  expect_error(get_datasets("TEST"), "Expecting MultiDataSet object")

  ds_list <- get_datasets(multiomics_set)

  expect_equal(names(ds_list), c("snps+A", "rnaseq", "metabolome", "phenotypes"))
  for (x in ds_list) expect_true(is.matrix(x))
  expect_equal(sapply(ds_list, dim), matrix(c(30, 15, 35, 15, 40, 15, 40, 15),
                                            nrow = 2,
                                            dimnames = list(c(), c("snps+A", "rnaseq", "metabolome", "phenotypes"))
  ))

  ## Making sure it is keeping the features and samples name
  expect_equal(rownames(ds_list$snps), paste0("featureA_", 1:30))
  expect_equal(colnames(ds_list$snps), paste0("sample_", 1:15))
})

test_that("get_samples_metadata works", {
  multiomics_set <- test_get_multidataset()

  expect_error(get_samples_metadata("TEST"), "Expecting MultiDataSet object")

  smeta_list <- get_samples_metadata(multiomics_set)

  expect_equal(names(smeta_list), c("snps+A", "rnaseq", "metabolome", "phenotypes"))
  for (x in smeta_list) expect_s3_class(x, "data.frame")
  expect_equal(sapply(smeta_list, dim), matrix(rep(c(15, 4), 4),
                                               nrow = 2,
                                               dimnames = list(c(), c("snps+A", "rnaseq", "metabolome", "phenotypes"))
  ))

  expect_equal(rownames(smeta_list$snps), paste0("sample_", 1:15))
  expect_equal(colnames(smeta_list$snps), c("id", "name", "pheno_group", "time"))
})

test_that("get_samples_metadata_combined works", {
  multiomics_set <- test_get_multidataset()

  expect_error(
    get_samples_metadata_combined(multiomics_set),
    NA
  )

  res <- get_samples_metadata_combined(multiomics_set, only_common_cols = TRUE)
  expect_equal(
    rownames(res),
    get_samples_metadata(multiomics_set) |>
      purrr::map(rownames) |>
      purrr::reduce(union)
  )
  expect_equal(
    colnames(res),
    c("id", "name", "pheno_group", "time")
  )

  ## Testing with different columns
  smeta_list <- test_get_smeta_list()

  smeta_list[[1]] <- smeta_list[[1]] |>
    dplyr::select(-time)
  smeta_list[[2]] <- smeta_list[[2]] |>
    dplyr::mutate(some_col = 2)

  mo_set <- test_get_multidataset(smeta_list = smeta_list)

  expect_equal(
    get_samples_metadata_combined(mo_set, only_common_cols = TRUE) |>
      colnames(),
    c("id", "name", "pheno_group")
  )
  expect_equal(
    get_samples_metadata_combined(mo_set, only_common_cols = FALSE) |>
      colnames(),
    c("id", "name", "pheno_group", "time", "some_col")
  )
})

test_that("get_samples works", {
  multiomics_set <- test_get_multidataset()

  expect_error(get_samples("TEST"), "Expecting MultiDataSet object")
  expect_no_error(get_samples(multiomics_set))

  expect_equal(
    get_samples(multiomics_set),
    list(
      "snps+A" = paste0("sample_", 1:15),
      "rnaseq" = paste0("sample_", 6:20),
      "metabolome" = paste0("sample_", 11:25),
      "phenotypes" = paste0("sample_", 16:30)
    )
  )
})

test_that("get_features works", {
  multiomics_set <- test_get_multidataset()

  expect_error(get_features("TEST"), "Expecting MultiDataSet object")
  expect_no_error(get_features(multiomics_set))

  expect_equal(
    get_features(multiomics_set),
    list(
      "snps+A" = paste0("featureA_", 1:30),
      "rnaseq" = paste0("featureB_", 1:35),
      "metabolome" = paste0("featureC_", 1:40),
      "phenotypes" = paste0("featureD_", 1:40)
    )
  )
})

test_that("get_features_metadata works", {
  multiomics_set <- test_get_multidataset()

  expect_error(get_features_metadata("TEST"), "Expecting MultiDataSet object")

  fmeta_list <- get_features_metadata(multiomics_set)

  expect_equal(names(fmeta_list), c("snps+A", "rnaseq", "metabolome", "phenotypes"))
  for (x in fmeta_list) expect_s3_class(x, "data.frame")
  expect_equal(sapply(fmeta_list, dim), matrix(c(30, 3, 35, 5, 40, 3, 40, 2),
                                               nrow = 2,
                                               dimnames = list(c(), c("snps+A", "rnaseq", "metabolome", "phenotypes"))
  ))

  expect_equal(rownames(fmeta_list$snps), paste0("featureA_", 1:30))
  expect_equal(colnames(fmeta_list$snps), c("feature_id", "chromosome", "position"))
})

test_that("join_features_metadata works", {
  multiomics_set <- test_get_multidataset()

  df <- tidyr::expand_grid(
    dataset = LETTERS[1:3],
    number = 1:3
  ) |>
    dplyr::mutate(
      feature_id = paste0("feature", dataset, "_", number),
      value = runif(dplyr::n())
    )

  expect_error(
    join_features_metadata(df, "TEST"),
    "Expecting MultiDataSet object"
  )

  expect_no_error(join_features_metadata(df, multiomics_set))

  expect_equal(join_features_metadata(df, NULL), df)

  res <- join_features_metadata(df, multiomics_set)
  fmeta_cols <- get_features_metadata(multiomics_set)[1:3] |>
    purrr::map(colnames) |>
    unlist() |>
    unname() |>
    unique()

  expect_equal(nrow(res), nrow(df))
  expect_equal(colnames(res), union(colnames(df), fmeta_cols))

  df2 <- df |>
    dplyr::mutate(feature_id = paste0(feature_id, "*"))

  expect_error(
    join_features_metadata(df2, multiomics_set),
    "No features matching IDs present in `df$feature_id`.",
    fixed = TRUE
  )
})

test_that("join_samples_metadata works", {
  multiomics_set <- test_get_multidataset()

  df <- tibble::tibble(
    id = paste0("sample_", 1:5),
    value = 1:5
  )

  expect_error(
    join_samples_metadata(df, "TEST"),
    "Expecting MultiDataSet object"
  )

  expect_error(
    join_samples_metadata(df, multiomics_set, "TEST"),
    "'TEST' datasets are not present in mo_data. Possible dataset names are: 'snps\\+A', 'rnaseq', 'metabolome', 'phenotypes'."
  )

  expect_no_error(join_samples_metadata(df, multiomics_set))
  expect_no_error(join_samples_metadata(df, multiomics_set, "snps+A"))

  expect_equal(join_samples_metadata(df, NULL), df)

  res <- join_samples_metadata(df, multiomics_set)
  smeta_cols <- get_samples_metadata_combined(multiomics_set, only_common_cols = FALSE) |>
    colnames()

  expect_equal(nrow(res), nrow(df))
  expect_equal(colnames(res), union(colnames(df), smeta_cols))

  expect_error(
    join_samples_metadata(df, multiomics_set, datasets = "rnaseq"),
    "No samples matching IDs present in `df$id`.",
    fixed = TRUE
  )
})

test_that("add_features_metadata works", {
  multiomics_set <- test_get_multidataset()

  df_ok <- get_features(multiomics_set) |>
    purrr::map_dfr(
      ~ tibble::tibble(feature_id = .x),
      .id = "dataset"
    ) |>
    dplyr::mutate(
      pvalue = runif(dplyr::n()),
      de_outcome = sample(LETTERS[1:3], dplyr::n(), replace = TRUE)
    )

  ## Errors in input class
  expect_error(
    add_features_metadata("TEST", df_ok),
    "'mo_data' argument: Expecting MultiDataSet object"
  )
  expect_error(
    add_features_metadata(multiomics_set, "TEST"),
    "'df' argument should be a tibble or data-frame."
  )

  ## Missing columns from df
  expect_error(
    add_features_metadata(multiomics_set, dplyr::rename(df_ok, TEST = feature_id)),
    "'df' argument should have columns 'feature_id' and 'dataset'."
  )
  expect_error(
    add_features_metadata(multiomics_set, dplyr::rename(df_ok, TEST = dataset)),
    "'df' argument should have columns 'feature_id' and 'dataset'."
  )

  ## Wrong dataset name in df
  df_wrong <- df_ok |>
    dplyr::mutate(
      dataset = dplyr::case_when(
        dataset == "snps+A" ~ "TEST",
        TRUE ~ dataset
      )
    )
  expect_error(
    add_features_metadata(multiomics_set, df_wrong),
    "'df' argument: in 'dataset' column, 'TEST' are not datasets in 'mo_data'. Possible dataset names are:.+"
  )

  ## Wrong and missing features ID
  df_wrong <- df_ok |>
    dplyr::mutate(
      feature_id = dplyr::case_when(
        feature_id == "featureA_1" ~ "TEST1",
        feature_id == "featureA_2" ~ "TEST2",
        TRUE ~ feature_id
      )
    )
  w <- capture_warnings(add_features_metadata(multiomics_set, df_wrong))
  expect_match(
    w,
    "snps\\+A dataset: 2 feature IDs not in 'mo_data', will be removed from features metadata.",
    all = FALSE
  )
  expect_match(
    w,
    "snps\\+A dataset: 2 feature IDs missing from 'df' data-frame.",
    all = FALSE
  )

  ## Duplicated rows
  df_wrong <- df_ok |>
    tibble::add_row(
      feature_id = "featureA_1",
      dataset = "snps+A",
      pvalue = 0.42,
      de_outcome = "A"
    )
  expect_error(
    add_features_metadata(multiomics_set, df_wrong),
    "snps\\+A dataset: should have only one row per feature ID in 'df' data-frame."
  )

  ## Columns already exist
  df_wrong <- df_ok |>
    dplyr::rename(
      chromosome = de_outcome
    )
  expect_error(
    add_features_metadata(multiomics_set, df_wrong),
    "snps\\+A dataset: 'chromosome' columns already in features metadata table."
  )

  expect_no_error(add_features_metadata(multiomics_set, df_ok))

  ## Each dataset's features metadata modified
  res <- add_features_metadata(multiomics_set, df_ok)
  ref_fmeta <- get_features_metadata(multiomics_set)
  res_fmeta <- get_features_metadata(res)

  expect_s4_class(res, "MultiDataSet")
  expect_equal(names(res), names(multiomics_set))
  expect_equal(get_datasets(res), get_datasets(multiomics_set))
  expect_equal(get_samples_metadata(res), get_samples_metadata(multiomics_set))
  expect_equal(names(res_fmeta), names(multiomics_set))
  purrr::walk2(
    res_fmeta,
    ref_fmeta,
    function(.x, .y) {
      expect_equal(nrow(.x), nrow(.y))
      expect_equal(colnames(.x), c(colnames(.y), c("pvalue", "de_outcome")))
      expect_equal(.x[colnames(.y)], .y)
    }
  )

  ## What happens to missing features
  df_wrong <- df_ok |>
    dplyr::filter(!(feature_id %in% c("featureA_2", "featureA_3")))
  res <- suppressWarnings(add_features_metadata(multiomics_set, df_wrong))
  ref_fmeta <- get_features_metadata(multiomics_set)[["snps+A"]]
  res_fmeta <- get_features_metadata(res)[["snps+A"]]

  expect_equal(nrow(res_fmeta), nrow(ref_fmeta))
  expect_equal(colnames(res_fmeta), c(colnames(ref_fmeta), c("pvalue", "de_outcome")))
  expect_equal(res_fmeta[colnames(ref_fmeta)], ref_fmeta)
  expect_equal(
    res_fmeta[c("featureA_2", "featureA_3"), "pvalue"],
    rep(NA_real_, 2)
  )
  expect_equal(
    res_fmeta[c("featureA_2", "featureA_3"), "de_outcome"],
    rep(NA_character_, 2)
  )
})

test_that("add_samples_metadata works", {
  multiomics_set <- test_get_multidataset()

  df_ok <- tibble::tibble(
    id = get_samples(multiomics_set) |>
      unlist() |>
      unname() |>
      unique()
  ) |>
    dplyr::mutate(
      new_info = sample(letters[1:4], dplyr::n(), replace = TRUE)
    )

  ## Errors in input class
  expect_error(
    add_samples_metadata("TEST", df_ok),
    "'mo_data' argument: Expecting MultiDataSet object"
  )
  expect_error(
    add_samples_metadata(multiomics_set, "TEST"),
    "'df' argument should be a tibble or data-frame."
  )

  ## Missing columns from df
  expect_error(
    add_samples_metadata(multiomics_set, dplyr::rename(df_ok, TEST = id)),
    "'df' argument should have columns 'id'."
  )

  ## Wrong dataset name
  expect_error(
    add_samples_metadata(multiomics_set, df_ok, datasets = "TEST"),
    "'datasets' argument: 'TEST' are not datasets in 'mo_data'. Possible dataset names are:.+"
  )

  ## Wrong and missing samples ID
  df_wrong <- df_ok |>
    dplyr::filter(id %in% paste0("sample_", 1:15)) |>
    dplyr::mutate(
      id = dplyr::case_when(
        id == "sample_1" ~ "TEST1",
        id == "sample_2" ~ "TEST2",
        TRUE ~ id
      )
    )
  w <- capture_warnings(add_samples_metadata(multiomics_set, df_wrong, datasets = "snps+A"))
  expect_match(
    w,
    "snps\\+A dataset: 2 sample IDs not in 'mo_data', will be removed from samples metadata.",
    all = FALSE
  )
  expect_match(
    w,
    "snps\\+A dataset: 2 sample IDs missing from 'df' data-frame.",
    all = FALSE
  )

  ## Duplicated rows
  df_wrong <- df_ok |>
    tibble::add_row(
      id = "sample_1",
      new_info = "a"
    )
  expect_error(
    add_samples_metadata(multiomics_set, df_wrong),
    "'df' argument: should have only one row per sample ID."
  )

  ## Columns already exist
  df_wrong <- df_ok |>
    dplyr::rename(
      pheno_group = new_info
    )
  expect_error(
    suppressWarnings(add_samples_metadata(multiomics_set, df_wrong)),
    "snps\\+A dataset: 'pheno_group' columns already in samples metadata table."
  )

  expect_no_error(suppressWarnings(add_samples_metadata(multiomics_set, df_ok)))
  expect_no_error(suppressWarnings(add_samples_metadata(multiomics_set, df_ok, datasets = "rnaseq")))

  ## Each dataset's samples metadata modified
  res <- suppressWarnings(add_samples_metadata(multiomics_set, df_ok))
  ref_smeta <- get_samples_metadata(multiomics_set)
  res_smeta <- get_samples_metadata(res)

  expect_s4_class(res, "MultiDataSet")
  expect_equal(names(res), names(multiomics_set))
  expect_equal(get_datasets(res), get_datasets(multiomics_set))
  expect_equal(get_features_metadata(res), get_features_metadata(multiomics_set))
  expect_equal(names(res_smeta), names(multiomics_set))
  purrr::walk2(
    res_smeta,
    ref_smeta,
    function(.x, .y) {
      expect_equal(nrow(.x), nrow(.y))
      expect_equal(colnames(.x), c(colnames(.y), c("new_info")))
      expect_equal(.x[colnames(.y)], .y)
    }
  )

  ## Specifying dataset
  res <- suppressWarnings(add_samples_metadata(multiomics_set, df_ok, datasets = "rnaseq"))
  res_smeta <- get_samples_metadata(res)
  other_ds <- setdiff(names(multiomics_set), "rnaseq")

  expect_s4_class(res, "MultiDataSet")
  expect_equal(names(res), names(multiomics_set))
  expect_equal(get_datasets(res), get_datasets(multiomics_set))
  expect_equal(get_features_metadata(res), get_features_metadata(multiomics_set))
  expect_equal(res_smeta[other_ds], ref_smeta[other_ds])
  expect_equal(nrow(res_smeta[["rnaseq"]]), nrow(ref_smeta[["rnaseq"]]))
  expect_equal(colnames(res_smeta[["rnaseq"]]), c(colnames(ref_smeta[["rnaseq"]]), c("new_info")))
  expect_equal(res_smeta[["rnaseq"]][colnames(ref_smeta[["rnaseq"]])], ref_smeta[["rnaseq"]])

  ## What happens to missing features
  df_wrong <- df_ok |>
    dplyr::filter(!(id %in% c("sample_2", "sample_3")))
  res <- suppressWarnings(add_samples_metadata(multiomics_set, df_wrong))
  ref_smeta <- get_samples_metadata(multiomics_set)[["snps+A"]]
  res_smeta <- get_samples_metadata(res)[["snps+A"]]

  expect_equal(nrow(res_smeta), nrow(ref_smeta))
  expect_equal(colnames(res_smeta), c(colnames(ref_smeta), c("new_info")))
  expect_equal(res_smeta[colnames(ref_smeta)], ref_smeta)
  expect_equal(
    res_smeta[c("sample_2", "sample_3"), "new_info"],
    rep(NA_character_, 2)
  )
})

test_that("check_missing_values works", {
  multiomics_set <- test_get_multidataset()

  expect_error(check_missing_values("TEST"), "Expecting MultiDataSet object")
  expect_message(check_missing_values(multiomics_set),
                 paste0(
                   c(
                     "50 (11.11%) missing values in snps+A dataset, across 27 features and 14 samples.",
                     "No missing values in rnaseq dataset.",
                     "No missing values in metabolome dataset.",
                     "No missing values in phenotypes dataset."
                   ),
                   collapse = "\n"
                 ),
                 fixed = TRUE
  )
  expect_equal(suppressMessages(check_missing_values(multiomics_set)), c("snps+A" = T, "rnaseq" = F, "metabolome" = F, "phenotypes" = F))
})

test_that("subset_features", {
  multiomics_set <- test_get_multidataset()

  t1 <- subset_features(multiomics_set, paste0("feature", rep(LETTERS[1:2], each = 5), "_", rep(1:5, 2)))

  expect_s4_class(t1, "MultiDataSet")
  expect_equal(n_features(t1), c("snps+A" = 5, "rnaseq" = 5, "metabolome" = 0, "phenotypes" = 0))
  expect_equal(n_samples(t1), c("snps+A" = 15, "rnaseq" = 15, "metabolome" = 15, "phenotypes" = 15))

  ## Testing that it works when one of the features name is not in the dataset
  t2 <- subset_features(multiomics_set, c(paste0("feature", rep(LETTERS[1:2], each = 5), "_", rep(1:5, 2)), "TEST"))
  expect_equal(n_features(t2), n_features(t1)) ## can't compare t1 and t2 directly because they save the datasets as different environments
  expect_equal(n_samples(t2), n_samples(t1))

  ## Testing that subsetting works when a list of features is supplied
  t3 <- subset_features(multiomics_set, lapply(LETTERS[1:2], function(i) {
    paste0("feature", i, "_", 1:5)
  }))
  expect_equal(n_features(t3), n_features(t1)) ## can't compare t1 and t2 directly because they save the datasets as different environments
  expect_equal(n_samples(t3), n_samples(t1))
})


test_that("replace_dataset works", {
  multiomics_set <- test_get_multidataset()

  ## Creating a new data matrix
  new_dat <- matrix(
    42, nrow = 30, ncol = 15,
    dimnames = list(
      paste0("featureA_", 1:30),
      paste0("sample_", 1:15)
    )
  )

  ## Testing input
  expect_error(
    replace_dataset("TEST", "snps+A", new_dat),
    "Expecting MultiDataSet object"
  )
  expect_error(
    replace_dataset(multiomics_set, c("snps+A", "rnaseq"), new_dat),
    "'dataset_name' argument should be of length 1."
  )
  expect_error(
    replace_dataset(multiomics_set, "TEST", new_dat),
    "'dataset_name' argument: 'TEST' is not an existing dataset. Possible dataset names are: 'snps+A', 'rnaseq', 'metabolome', 'phenotypes'.",
    fixed = TRUE
  )

  new_dat2 <- matrix(
    42, nrow = 20, ncol = 15,
    dimnames = list(
      paste0("featureA_", 1:20),
      paste0("sample_", 1:15)
    )
  )
  expect_error(
    replace_dataset(multiomics_set, "snps+A", new_dat2),
    "'new_data' argument has incorrect dimensions. Should have 30 rows (features) and 15 columns (samples).",
    fixed = TRUE
  )

  new_dat2 <- matrix(
    42, nrow = 30, ncol = 15,
    dimnames = list(
      NULL,
      paste0("sample_", 1:15)
    )
  )
  expect_error(
    replace_dataset(multiomics_set, "snps+A", new_dat2),
    "'new_data' should have feature IDs as rownames and sample IDs as column names."
  )

  new_dat2 <- matrix(
    42, nrow = 30, ncol = 15,
    dimnames = list(
      paste0("feature_", 1:30),
      paste0("sample_", 1:15)
    )
  )
  expect_error(
    replace_dataset(multiomics_set, "snps+A", new_dat2),
    "'new_data' rownames do not match feature IDs in existing dataset."
  )

  new_dat2 <- matrix(
    42, nrow = 30, ncol = 15,
    dimnames = list(
      paste0("featureA_", 1:30),
      paste0("sampleA_", 1:15)
    )
  )
  expect_error(
    replace_dataset(multiomics_set, "snps+A", new_dat2),
    "'new_data' colnames do not match sample IDs in existing dataset."
  )

  ## Testing that the different assays are correct
  t1 <- replace_dataset(multiomics_set, "snps+A", new_dat)

  expect_equal(names(t1), names(multiomics_set))
  expect_equal(get_datasets(t1)[["snps+A"]], new_dat)
  expect_equal(
    t1@assayData[["snps+A"]]$callProbability,
    multiomics_set@assayData$snps$callProbability
  )

  ## Testing replacement for a dataset with no additional name
  new_dat <- matrix(
    42, nrow = 35, ncol = 15,
    dimnames = list(
      paste0("featureB_", 1:35),
      paste0("sample_", 6:20)
    )
  )
  t2 <- replace_dataset(multiomics_set, "rnaseq", new_dat)

  expect_equal(names(t1), names(multiomics_set))
  expect_equal(get_datasets(t2)[["rnaseq"]], new_dat)
})


test_that("round_dataset works", {
  multiomics_set <- test_get_multidataset()

  ## Getting the original dataset
  ds <- get_datasets(multiomics_set)[["snps+A"]]

  ## Testing input
  expect_error(round_dataset("TEST", "snps+A"), "Expecting MultiDataSet object")
  expect_error(round_dataset(multiomics_set, "TEST"),
               "'dataset_name' argument: 'TEST' is not an existing dataset. Possible dataset names are: 'snps+A', 'rnaseq', 'metabolome', 'phenotypes'.",
               fixed = TRUE
  )

  t1 <- round_dataset(multiomics_set, "snps+A")
  expect_equal(get_datasets(t1)[["snps+A"]], round(ds, 0))

  t2 <- round_dataset(multiomics_set, "snps+A", ndecimals = 2)
  expect_equal(get_datasets(t2)[["snps+A"]], round(ds, 2))

  t3 <- round_dataset(multiomics_set, "snps+A", min_val = 20, max_val = 80)
  expect_equal(range(get_datasets(t3)[["snps+A"]], na.rm = TRUE), c(20, 80))
})


test_that("get_features_labels works", {
  multiomics_set <- test_get_multidataset()

  expect_no_error(
    get_features_labels(multiomics_set)
  )
  expect_error(
    get_features_labels(multiomics_set, "TEST"),
    "'label_cols' argument: 'TEST' is not a column in the features metadata for the snps\\+A dataset.+"
  )
  expect_error(
    get_features_labels(multiomics_set, list("metabolome" = "TEST")),
    "'label_cols' argument: 'TEST' is not a column in the features metadata for the metabolome dataset.+"
  )
  expect_error(
    get_features_labels(multiomics_set, list("TEST" = "name")),
    "'label_cols' argument: 'TEST' are not existing datasets or do not match the 'datasets' argument. Possible values are:.+"
  )
  expect_no_error(
    get_features_labels(multiomics_set, list("metabolome" = "name"))
  )

  ## By default labels are features ID
  res <- get_features_labels(multiomics_set)
  expect_equal(
    nrow(res),
    145
  )
  expect_equal(
    res$feature_id,
    res$label
  )

  res <- get_features_labels(multiomics_set, list("snps+A" = "chromosome", "metabolome" = "name"))
  ## Datasets with no label_col have feature ID as label
  expect_equal(
    res$feature_id[!(res$dataset %in% c("snps+A", "metabolome"))],
    res$label[!(res$dataset %in% c("snps+A", "metabolome"))]
  )
  ## The other one have specific label
  expect_true(
    all(stringr::str_detect(
      res$label[res$dataset == "snps+A"],
      "^ch"
    ))
  )
  expect_true(
    all(res$label[res$feature_id %in% paste0("featureC_", 1:30)] != paste0("featureC_", 1:30))
  )
  ## NA labels are replaced with feature ID
  expect_true(
    all(res$label[res$feature_id %in% paste0("featureC_", 31:40)] == paste0("featureC_", 31:40))
  )

  ## Testing truncation option
  expect_equal(
    get_features_labels(multiomics_set, truncate = 9)$label,
    rep("featur...", 145)
  )
})
