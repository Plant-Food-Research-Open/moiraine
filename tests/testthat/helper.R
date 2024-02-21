test_get_data_list <- function() {
  ## File names
  n_omics <- 4
  data_files <- paste0(
    test_path("fixtures", "data_omics"),
    LETTERS[1:n_omics], ".csv"
  )
  names(data_files) <- LETTERS[1:n_omics]

  ## Load data
  data_list <- lapply(LETTERS[1:n_omics], function(i) {
    import_dataset_csv(
      data_files[i],
      col_id = ifelse(i == "B", "Sample", "Feature"),
      features_as_rows = (i != "B"),
      show_col_types = FALSE
    )
  })
  names(data_list) <- LETTERS[1:n_omics]

  return(data_list)
}

test_get_fmeta_list <- function() {
  ## File names
  n_omics <- 4
  fmeta_files <- paste0(
    test_path("fixtures", "fmeta_omics"),
    LETTERS[1:n_omics], ".csv"
  )
  names(fmeta_files) <- LETTERS[1:n_omics]

  ## Load features metadata
  fmeta_list <- lapply(LETTERS[1:n_omics], function(i) {
    import_fmetadata_csv(
      fmeta_files[i],
      col_id = "Feature",
      show_col_types = FALSE
    )
  })
  names(fmeta_list) <- LETTERS[1:n_omics]

  return(fmeta_list)
}


test_get_smeta_list <- function() {
  ## File names
  n_omics <- 4
  smeta_files <- paste0(
    test_path("fixtures", "smeta_omics"),
    LETTERS[1:n_omics], ".csv"
  )
  names(smeta_files) <- LETTERS[1:n_omics]

  ## Load samples metadata
  smeta_list <- lapply(LETTERS[1:n_omics], function(i) {
    import_smetadata_csv(
      smeta_files[i],
      col_id = "Sample",
      show_col_types = FALSE
    )
  })
  names(smeta_list) <- LETTERS[1:n_omics]

  return(smeta_list)
}

test_get_omics_list <- function(data_list = NULL,
                                fmeta_list = NULL,
                                smeta_list = NULL) {
  if (is.null(data_list)) data_list <- test_get_data_list()
  if (is.null(fmeta_list)) fmeta_list <- test_get_fmeta_list()
  if (is.null(smeta_list)) smeta_list <- test_get_smeta_list()

  types <- c(
    "A" = "genomics",
    "B" = "transcriptomics",
    "C" = "metabolomics",
    "D" = "phenomics"
  )

  omics_sets <- purrr::map(
    names(data_list),
    ~ suppressWarnings(
      create_omics_set(
        data_list[[.x]],
        omics_type = types[[.x]],
        features_metadata = fmeta_list[[.x]],
        samples_metadata = smeta_list[[.x]]
      )
    )
  )

  return(omics_sets)
}

test_get_multidataset <- function(data_list = NULL,
                                  fmeta_list = NULL,
                                  smeta_list = NULL) {
  if (all(is.null(data_list), is.null(fmeta_list), is.null(smeta_list))) {
    multiomics_set <- readRDS(test_path("fixtures", "multiomics_set.rds"))
  } else {
    omics_sets <- test_get_omics_list(data_list, fmeta_list, smeta_list)
    multiomics_set <- create_multiomics_set(
      omics_sets,
      datasets_names = c("A", "", "", "")
    )
  }

  return(multiomics_set)
}

test_get_pca_run <- function() {
  readRDS(test_path("fixtures", "pca_res.rds"))
}

test_get_splsda_run <- function() {
  readRDS(test_path("fixtures", "splsda_res.rds"))
}

test_get_spls_run <- function() {
  readRDS(test_path("fixtures", "spls_res.rds"))
}
test_get_diablo_run <- function() {
  readRDS(test_path("fixtures", "diablo_res.rds"))
}

test_get_so2pls_run <- function() {
  readRDS(test_path("fixtures", "so2pls_res.rds"))
}

test_get_mofa_run <- function() {
  readRDS(test_path("fixtures", "mofa_res.rds"))
}

test_get_output_list <- function() {
  diablo_res <- test_get_diablo_run()
  mofa_res <- test_get_mofa_run()
  so2pls_res <- test_get_so2pls_run()

  res_list <- list(
    get_output_diablo(diablo_res),
    get_output_mofa2(mofa_res),
    get_output_so2pls(so2pls_res)
  )

  return(res_list)
}


test_helper_output <- function(x, method) {
  expect_s3_class(x, "output_dimension_reduction")

  expect_equal(
    names(x),
    c("features_weight", "samples_score", "variance_explained")
  )
  expect_equal(
    names(x$features_weight),
    c("feature_id", "dataset", "latent_dimension", "weight", "importance")
  )
  expect_equal(
    names(x$samples_score),
    c("sample_id", "latent_dimension", "score")
  )
  expect_equal(
    setdiff(names(x$variance_explained), "group"),
    c("latent_dimension", "dataset", "prop_var_expl")
  )

  expect_s3_class(x$features_weight$latent_dimension, "factor")
  expect_s3_class(x$samples_score$latent_dimension, "factor")
  expect_s3_class(x$variance_explained$latent_dimension, "factor")

  expect_s3_class(x$features_weight$dataset, "factor")
  expect_s3_class(x$variance_explained$dataset, "factor")

  expect_equal(attr(x, "method"), method)

  return(invisible(NULL))
}

test_extract_features_weight <- function(x, ld, ds) {
  x$features_weight |>
    dplyr::filter(latent_dimension == ld, dataset == ds) |>
    dplyr::mutate(temp = readr::parse_number(feature_id)) |>
    dplyr::arrange(temp) |>
    dplyr::select(feature_id, weight) |>
    tibble::deframe()
}

test_extract_samples_score <- function(x, ld) {
  x$samples_score |>
    dplyr::filter(latent_dimension == ld) |>
    dplyr::mutate(temp = readr::parse_number(sample_id)) |>
    dplyr::arrange(temp) |>
    dplyr::select(sample_id, score) |>
    tibble::deframe()
}

test_clean_expression <- function(x) {
  as.character(x) |>
    stringr::str_remove_all("\n") |>
    stringr::str_replace_all("\\s+", " ") |>
    stringr::str_remove_all("(?<=function) (?=\\()")
}
