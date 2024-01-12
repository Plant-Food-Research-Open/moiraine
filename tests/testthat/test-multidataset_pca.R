test_that("run_pca works", {
  multiomics_set <- test_get_multidataset()

  ## Testing input
  expect_error(run_pca("TEST", "rnaseq"), "Expecting MultiDataSet object")
  expect_error(run_pca(multiomics_set, c("rnaseq", "snps+A")), "'dataset_name' argument should be of length 1.")
  expect_error(run_pca(multiomics_set, "TEST"),
               "'TEST' datasets are not present in mo_data. Possible dataset names are: 'snps+A', 'rnaseq', 'metabolome', 'phenotypes'.",
               fixed = TRUE
  )
  expect_error(run_pca(multiomics_set, "snps+A", n_pcs = 0), "'n_pcs' argument should be an integer value >= 1.")
  expect_error(run_pca(multiomics_set, "snps+A", scale = "TEST"), "'scale' argument should be one of 'none', 'pareto', 'vector', 'uv'.")
  expect_error(run_pca(multiomics_set, "snps+A", center = "TEST"), "'center' argument should be either TRUE or FALSE.")
  expect_error(run_pca(multiomics_set, "snps+A", method = "TEST"), "'method' argument should be one of:.+")

  ## Testing output
  pca_res <- run_pca(multiomics_set, "metabolome", n_pcs = 5, scale = "pareto", center = FALSE, method = "svd")

  expect_s4_class(pca_res, "pcaRes")
  expect_equal(pca_res@completeObs, t(get_datasets(multiomics_set)[["metabolome"]]))
  expect_equal(pca_res@nPcs, 5)
  expect_equal(pca_res@scaled, "pareto")
  expect_equal(pca_res@centered, FALSE)
  expect_equal(pca_res@method, "svd")
  expect_equal(attr(pca_res, "dataset_name"), "metabolome")

  ## Testing automatic method selection
  expect_equal(run_pca(multiomics_set, "snps+A")@method, "nipals") ## with missing values
  expect_equal(run_pca(multiomics_set, "phenotypes")@method, "svd") ## with no missing values
})

test_that("pca_complete_data_factory works", {
  tar_res <- pca_complete_data_factory(mo_data)

  ## Testing that we get a list of target stems
  expect_type(tar_res, "list")
  expect_s3_class(tar_res[[1]], "tar_stem")

  ## Testing targets name
  expect_equal(
    sapply(tar_res, function(x) {
      x$settings$name
    }),
    c("dataset_names_pca", "pca_mats_list", "pca_runs_list", "complete_set")
  )

  ## Testing targets command
  expect_equal(
    sapply(tar_res, function(x) {
      x$command$expr
    }),
    c(
      expression(names(mo_data)),
      expression(get_dataset_matrix(mo_data, dataset_names_pca, keep_dataset_name = TRUE)),
      expression(run_pca_matrix(pca_mats_list)),
      expression(get_complete_data(mo_data, pca_runs_list))
    )
  )

  ## Testing targets pattern
  expect_equal(
    sapply(tar_res, function(x) {
      x$settings$pattern
    }),
    list(NULL, expression(map(dataset_names_pca)), expression(map(pca_mats_list)), NULL)
  )

  ## Testing iteration method (for dynamic branching it should be a list)
  expect_equal(
    sapply(tar_res, function(x) {
      x$settings$iteration
    }),
    c("vector", "list", "list", "vector")
  )

  ## Testing targets format
  expect_equal(
    sapply(tar_res, function(x) {
      x$settings$format
    }),
    c("rds", "rds", "rds", "rds")
  )

  ## Testing custom dataset names
  tar_res <- pca_complete_data_factory(mo_data, setdiff(names(mo_data), "phenotypes"))
  expect_equal(
    sapply(tar_res, function(x) {
      x$command$expr
    }),
    c(
      expression(setdiff(names(mo_data), "phenotypes")),
      expression(get_dataset_matrix(mo_data, dataset_names_pca, keep_dataset_name = TRUE)),
      expression(run_pca_matrix(pca_mats_list)),
      expression(get_complete_data(mo_data, pca_runs_list))
    )
  )

  ## Testing custom prefix for target names
  tar_res <- pca_complete_data_factory(mo_data, target_name_prefix = "test_")
  expect_equal(
    sapply(tar_res, function(x) {
      x$settings$name
    }),
    paste0("test_", c("dataset_names_pca", "pca_mats_list", "pca_runs_list", "complete_set"))
  )

  ## Testing custom complete target name + custom prefix for target names
  tar_res <- pca_complete_data_factory(mo_data, target_name_prefix = "test_", complete_data_name = "final_tar")
  expect_equal(
    sapply(tar_res, function(x) {
      x$settings$name
    }),
    c("test_dataset_names_pca", "test_pca_mats_list", "test_pca_runs_list", "final_tar")
  )

  ## Testing additional parameters
  tar_res <- pca_complete_data_factory(mo_data, n_pcs = 5, scale = "none")
  expect_equal(
    sapply(tar_res, function(x) {
      x$command$expr
    }),
    c(
      expression(names(mo_data)),
      expression(get_dataset_matrix(mo_data, dataset_names_pca, keep_dataset_name = TRUE)),
      expression(run_pca_matrix(pca_mats_list, n_pcs = 5, scale = "none")),
      expression(get_complete_data(mo_data, pca_runs_list))
    )
  )
})


test_that("get_pca_arguments works", {
  multiomics_set <- test_get_multidataset()

  ## Testing input
  expect_error(get_pca_arguments("TEST"), "Input should be a list of pcaRes objects (from pcaMethods package).", fixed = TRUE)
  expect_error(get_pca_arguments(list("TEST", "TEST")), "Input should be a list of pcaRes objects (from pcaMethods package).", fixed = TRUE)

  pca_list <- sapply(names(multiomics_set), function(i) {
    run_pca(multiomics_set, i, n_pcs = 5)
  })

  t1 <- get_pca_arguments(pca_list)

  expect_s3_class(t1, "tbl_df")
  expect_equal(dim(t1), c(4, 5))
  expect_equal(t1$`Omics dataset`, names(multiomics_set))
})


test_that("plot_screeplot_pca works", {
  multiomics_set <- test_get_multidataset()

  pca_list <- sapply(names(multiomics_set), function(i) {
    run_pca(multiomics_set, i, n_pcs = 5)
  })

  expect_error(plot_screeplot_pca("TEST"), "Input should be a list of pcaRes objects (from pcaMethods package).", fixed = TRUE)
  expect_error(plot_screeplot_pca(pca_list, datasets = "TEST"), "'datasets' argument: 'TEST' not existing datasets. Possible values are:.+")


  plotres <- plot_screeplot_pca(pca_list)
  expect_s3_class(plotres, "ggplot")
  expect_no_error(plot_screeplot_pca(pca_list, cumulative = TRUE))
  expect_no_error(plot_screeplot_pca(pca_list, datasets = names(multiomics_set)[1]))
})


test_that("plot_samples_coordinates_pca works", {
  grDevices::pdf(NULL)

  fct <- function(...) {
    print(plot_samples_coordinates_pca(...))
  }

  multiomics_set <- test_get_multidataset()[, c("snps+A", "rnaseq", "metabolome")]
  pca_list <- names(multiomics_set) |>
    purrr::map(
      ~ run_pca(multiomics_set, .x, )
    )

  expect_error(
    fct("TEST"),
    "Expecting a pcaRes object (from run_pca() or pcaMethods::pca() function).",
    fixed = TRUE
  )
  expect_no_error(
    fct(pca_list)
  )

  expect_error(
    fct(pca_list, datasets = "TEST"),
    "'datasets' argument: 'TEST' not existing datasets. Possible values are:.+"
  )
  expect_no_error(
    fct(pca_list, datasets = "rnaseq")
  )

  ## Testing pcs argument
  expect_error(
    fct(pca_list, pcs = list(1:5)),
    "'pcs' argument should either be a vector or a named list of length 3."
  )
  expect_no_error(
    fct(pca_list, pcs = 1:5)
  )
  expect_no_error(
    fct(pca_list, pcs = list("snps+A" = 1:3, "rnaseq" = 1:4, metabolome = 1:2))
  )

  ## Testing basic arguments passing to the plotting function
  expect_error(
    fct(pca_list, colour_upper = "group"),
    "Need to provide a MultiDataSet object through 'mo_data' argument in order to use one of the aesthetics arguments."
  )
  expect_error(
    fct(pca_list, mo_data = multiomics_set, colour_upper = "group"),
    "'colour_upper' argument: 'group' is not a column in the samples metadata of any dataset. Possible values are:+."
  )
  expect_no_error(
    fct(pca_list, mo_data = multiomics_set, colour_upper = "pheno_group")
  )

  grDevices::dev.off()
})



test_that("get_complete_data works", {
  multiomics_set <- test_get_multidataset()

  pca_list <- lapply(names(multiomics_set), function(i) {
    run_pca(multiomics_set, i, n_pcs = 5)
  })
  names(pca_list) <- names(multiomics_set)

  ## Testing input
  expect_error(get_complete_data("TEST"), "Expecting MultiDataSet object")
  expect_error(get_complete_data(multiomics_set, "TEST"), "'pca_result' should be a list of pcaRes objects (from pcaMethods package).", fixed = TRUE)

  ## Testing correct behaviour
  t1 <- get_complete_data(multiomics_set, pca_list)
  expect_s4_class(t1, "MultiDataSet")
  expect_equal(names(t1), names(multiomics_set))
  expect_equal(get_datasets(t1)[c("snps+A", "rnaseq")], lapply(pca_list[c("snps+A", "rnaseq")], function(x) {
    t(x@completeObs)
  }))
  expect_equal(get_datasets(t1)[c("metabolome", "phenotypes")], get_datasets(multiomics_set)[c("metabolome", "phenotypes")])

  ## Testing all complete datasets
  expect_message(get_complete_data(multiomics_set[, c("metabolome", "phenotypes")], pca_list), "All datasets are complete, no imputation to perform.")
  t1 <- suppressMessages(get_complete_data(multiomics_set[, c("metabolome", "phenotypes")], pca_list))
  expect_equal(t1, multiomics_set[, c("metabolome", "phenotypes")])

  ## Testing missing PCA results
  expect_warning(get_complete_data(multiomics_set, pca_list[c("metabolome", "phenotypes")]),
                 "PCA results not available for dataset(s) 'snps+A'. NA values will not be imputed.",
                 fixed = TRUE
  )
  t1 <- suppressWarnings(get_complete_data(multiomics_set, pca_list[c("metabolome", "phenotypes")]))
  expect_equal(t1, multiomics_set)
})
