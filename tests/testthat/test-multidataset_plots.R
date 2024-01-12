test_that("plot_samples_upset works", {
  multiomics_set <- test_get_multidataset()

  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    expect_error(plot_samples_upset(multiomics_set), "Package \"UpSetR\" must be installed to use this function.")
    skip("Package \"UpSetR\" must be installed to use this function.")
  }

  grDevices::pdf(NULL)

  ## Just testing that there are no errors -need to find a way to disable plotting
  expect_error(plot_samples_upset(multiomics_set), regexp = NA)

  grDevices::dev.off()
})

test_that("plot_density_data works", {
  multiomics_set <- test_get_multidataset()

  ## Testing input
  expect_error(plot_density_data("TEST"), "Expecting MultiDataSet object")
  expect_error(
    plot_density_data(multiomics_set, c("rnaseq", "TEST")),
    "'TEST' datasets are not present in mo_data. Possible dataset names are: 'snps+A', 'rnaseq', 'metabolome', 'phenotypes'.",
    fixed = TRUE
  )
  expect_error(plot_density_data(multiomics_set, combined = FALSE), NA)
  expect_error(plot_density_data(multiomics_set, combined = FALSE, scales = "free"), NA)

  ## Testing that it returns a ggplot
  plotres <- plot_density_data(multiomics_set)
  expect_s3_class(plotres, "ggplot")
})

test_that("plot_meansd_data works", {
  multiomics_set <- test_get_multidataset()

  if (!requireNamespace("hexbin", quietly = TRUE)) {
    expect_error(plot_meansd_data(multiomics_set), "Package \"hexbin\" must be installed to use this function.")
    skip("Package \"hexbin\" must be installed to use this function.")
  }

  ## Testing input
  expect_error(plot_meansd_data("TEST"), "Expecting MultiDataSet object")
  expect_error(plot_meansd_data(multiomics_set, c("rnaseq", "TEST")),
               "'TEST' datasets are not present in mo_data. Possible dataset names are: 'snps+A', 'rnaseq', 'metabolome', 'phenotypes'.",
               fixed = TRUE
  )

  ## Testing that it returns a ggplot
  plotres <- plot_meansd_data(multiomics_set)
  expect_s3_class(plotres, "ggplot")
})

test_that("plot_data_heatmap works", {

  multiomics_set <- test_get_multidataset()
  f_list <- c(
    paste0("featureA_", 1:5),
    paste0("featureB_", 1:5)
  )
  s_list <- paste0("sample_", 1:10)

  expect_error(
    plot_data_heatmap("TEST"),
    "'mo_data' argument: Expecting MultiDataSet object"
  )
  expect_error(
    plot_data_heatmap(multiomics_set, "TEST"),
    "No feature selected."
  )
  expect_error(
    plot_data_heatmap(multiomics_set, c(f_list[1], "TEST")),
    "The following features are not present in any dataset: 'TEST'."
  )
  expect_error(
    plot_data_heatmap(multiomics_set, f_list, samples = "TEST"),
    "The following samples are not present in any dataset: 'TEST'."
  )

  expect_no_error(suppressWarnings(plot_data_heatmap(multiomics_set, f_list)))
  expect_warning(
    plot_data_heatmap(multiomics_set, f_list),
    "Not enough data to calculate distance between samples, disabling clustering of columns."
  )
  expect_no_warning(plot_data_heatmap(multiomics_set, f_list, cluster_columns = FALSE))

  expect_no_error(plot_data_heatmap(multiomics_set, f_list[1:5]))
  expect_no_warning(plot_data_heatmap(multiomics_set, f_list[1:5]))

  expect_no_error(plot_data_heatmap(multiomics_set, f_list, samples = s_list))

  ## Selecting common samples should be restricted to selected features
  expect_no_error(plot_data_heatmap(multiomics_set, f_list, only_common_samples = TRUE))

  expect_no_error(
    plot_data_heatmap(
      multiomics_set,
      f_list,
      only_common_samples = TRUE,
      label_cols = list(
        "rnaseq" = "name"
      )
    )
  )

  ## Adding rows and column annotations
  expect_error(
    plot_data_heatmap(multiomics_set, f_list, samples_info = "TEST", cluster_columns = FALSE),
    "'samples_info' argument: 'TEST' is not a column in the samples metadata of any dataset. Possible values are: 'id', 'name', 'pheno_group', 'time'."
  )
  expect_no_error(plot_data_heatmap(multiomics_set, f_list, samples_info = "pheno_group", cluster_columns = FALSE))

  expect_error(
    plot_data_heatmap(multiomics_set, f_list, features_info = "TEST", cluster_columns = FALSE),
    "argument: 'TEST' is not a column in the features metadata of any dataset. Possible values are: 'feature_id', 'chromosome', 'position', 'start', 'end', 'name'."
  )
  expect_no_error(plot_data_heatmap(multiomics_set, f_list, features_info = "chromosome", cluster_columns = FALSE)) ## shared by all features
  expect_no_error(plot_data_heatmap(multiomics_set, f_list, features_info = "position", cluster_columns = FALSE)) ## only for one dataset

  expect_no_error(plot_data_heatmap(multiomics_set, f_list, samples_info = "pheno_group", features_info = c("chromosome"), cluster_columns = FALSE))

  ## Testing custom colours
  expect_error(
    plot_data_heatmap(
      multiomics_set,
      f_list,
      samples_info = "pheno_group",
      features_info = "chromosome",
      cluster_columns = FALSE,
      colours_list = list("TEST" = c("rnaseq" = "blue", "snps+A" = "green"))
    ),
    "'colours_list' argument: 'TEST' are not names of features or samples metadata to be plotted. Possible values are: 'dataset', 'pheno_group', 'chromosome'."
  )
  expect_no_error(
    plot_data_heatmap(
      multiomics_set,
      f_list,
      samples_info = "pheno_group",
      features_info = "chromosome",
      cluster_columns = FALSE,
      colours_list = list("dataset" = c("rnaseq" = "blue", "snps+A" = "green"))
    )
  )
  expect_no_error(
    plot_data_heatmap(
      multiomics_set,
      f_list,
      samples_info = "pheno_group",
      features_info = "chromosome",
      cluster_columns = FALSE,
      colours_list = list(
        "dataset" = c("rnaseq" = "blue", "snps+A" = "green"),
        "pheno_group" = c("group1" = "purple", "group2" = "turquoise")
      )
    )
  )
})

test_that("plot_data_covariate works", {

  multiomics_set <- test_get_multidataset()
  f_list <- paste0("feature", LETTERS[1:3], "_1")
  s_list <- paste0("sample_", 1:20)

  expect_error(
    plot_data_covariate("TEST", "pheno_group", f_list),
    "'mo_data' argument: Expecting MultiDataSet object"
  )

  expect_error(
    plot_data_covariate(multiomics_set, "TEST", f_list),
    "'covariate' argument: 'TEST' is not a column in the samples metadata of any dataset. Possible values are:.+"
  )

  expect_no_error(plot_data_covariate(multiomics_set, "pheno_group", f_list))

  expect_error(
    plot_data_covariate(multiomics_set, "pheno_group", c(f_list, "TEST")),
    "The following features are not present in any dataset: 'TEST'."
  )

  expect_error(
    plot_data_covariate(multiomics_set, "pheno_group", f_list, c(s_list, "sample_100")),
    "The following samples are not present in any dataset: 'sample_100'."
  )

  expect_no_error(plot_data_covariate(multiomics_set, "pheno_group", f_list, s_list))

  expect_error(
    plot_data_covariate(multiomics_set, "pheno_group", f_list, colour_by = "TEST"),
    "'colour_by' argument: 'TEST' is not a column in the samples metadata of any dataset. Possible values are:.+"
  )
  expect_error(
    plot_data_covariate(multiomics_set, "pheno_group", f_list, shape_by = "TEST"),
    "'shape_by' argument: 'TEST' is not a column in the samples metadata of any dataset. Possible values are:.+"
  )
  expect_no_error(plot_data_covariate(multiomics_set, "pheno_group", f_list, colour_by = "time", shape_by = "pheno_group"))
})
