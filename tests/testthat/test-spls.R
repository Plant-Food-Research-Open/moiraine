test_that("get_input_spls works", {
  multiomics_set <- test_get_multidataset()

  expect_error(
    get_input_spls(multiomics_set, mode = "regression", datasets = "TEST"),
    "'datasets' argument: expecting a vector of length 2 (as (s)PLS integrates 2 datasets).",
    fixed = TRUE
  )

  expect_error(
    get_input_spls(multiomics_set, mode = "regression", datasets = c("rnaseq", "TEST")),
    "'TEST' datasets are not present in mo_data. Possible dataset names are:.+"
  )

  expect_error(
    get_input_spls(multiomics_set, mode = "TEST", datasets = c("rnaseq", "metabolome")),
    "'mode' argument should be one of: 'regression', 'canonical', 'invariant', 'classic'."
  )

  expect_no_error(
    get_input_spls(multiomics_set, mode = "regression", datasets = c("rnaseq", "metabolome"))
  )

  res <- get_input_spls(multiomics_set, mode = "regression", datasets = c("rnaseq", "metabolome"))

  expect_equal(names(res), c("rnaseq", "metabolome"))
  expect_equal(attr(res, "mode"), "regression")
  expect_equal(dim(res[[1]]), c(10, 35))
  expect_equal(dim(res[[2]]), c(10, 40))
  expect_null(attr(res, "multilevel"))
})

test_that("get_input_spls works - with multilevel option", {
  multiomics_set <- test_get_multidataset()

  expect_error(
    get_input_spls(multiomics_set, mode = "regression", datasets = c("rnaseq", "metabolome"), multilevel = "TEST"),
    "'multilevel' argument: 'TEST' is not a column in the samples metadata of any dataset. Possible values are:.+"
  )
  expect_error(
    get_input_spls(multiomics_set, mode = "regression", datasets = c("rnaseq", "metabolome"), multilevel = rep("TEST", 2)),
    "'multilevel' argument should be of length 1 (for one-factor decomposition) or 3 (for two-factor decomposition).",
    fixed = TRUE
  )

  expect_no_error(
    get_input_spls(multiomics_set, mode = "regression", datasets = c("rnaseq", "metabolome"), multilevel = "pheno_group")
  )
  expect_no_error(
    get_input_spls(multiomics_set, mode = "regression", datasets = c("rnaseq", "metabolome"), multilevel = c("pheno_group", "time", "name"))
  )

  ## One-factor decomposition
  res <- get_input_spls(multiomics_set, mode = "regression", datasets = c("rnaseq", "metabolome"), multilevel = "pheno_group")
  res_multilevel <- attr(res, "multilevel")
  expect_s3_class(res_multilevel, "data.frame")
  expect_equal(dim(res_multilevel), c(nrow(res[[1]]), 1))
  expect_equal(colnames(res_multilevel), "pheno_group")
  expect_equal(
    res_multilevel[, 1],
    get_samples_metadata_combined(multiomics_set)[rownames(res[[1]]), "pheno_group"] |>
      as.factor() |>
      as.integer()
  )

  ## Two-factor decomposition
  res <- get_input_spls(multiomics_set, mode = "regression", datasets = c("rnaseq", "metabolome"), multilevel = c("pheno_group", "time", "name"))
  res_multilevel <- attr(res, "multilevel")
  expect_s3_class(res_multilevel, "data.frame")
  expect_equal(dim(res_multilevel), c(nrow(res[[1]]), 3))
  expect_equal(colnames(res_multilevel), c("pheno_group", "time", "name"))
  expect_equal(
    res_multilevel[, 1],
    get_samples_metadata_combined(multiomics_set)[rownames(res[[1]]), "pheno_group"] |>
      as.factor() |>
      as.integer()
  )
  expect_equal(
    res_multilevel[, 2],
    get_samples_metadata_combined(multiomics_set)[rownames(res[[1]]), "time"] |>
      as.factor()
  )
  expect_equal(
    res_multilevel[, 3],
    get_samples_metadata_combined(multiomics_set)[rownames(res[[1]]), "name"] |>
      as.factor()
  )
})

test_that("spls_run works", {
  multiomics_set <- test_get_multidataset()
  spls_input <- get_input_spls(multiomics_set, datasets = c("rnaseq", "metabolome"), mode = "regression")

  ## For a PLS run
  expect_no_error(
    spls_run(spls_input, ncomp = 2)
  )

  res <- spls_run(spls_input, ncomp = 2)
  res_pls <- mixOmics::pls(
    spls_input[[1]],
    spls_input[[2]],
    ncomp = 2,
    mode = "regression"
  )

  expect_equal(
    attr(res, "datasets_name"),
    c("rnaseq", "metabolome")
  )

  res$call <- "removed from test"
  attr(res, "datasets_name") <- NULL
  res_pls$call <- "removed from test"
  expect_identical(res, res_pls)

  ## For a sPLS run
  expect_no_error(
    spls_run(spls_input, ncomp = 2, keepX = c(10, 10), keepY = c(5, 5))
  )

  res <- spls_run(spls_input, ncomp = 2, keepX = c(10, 10), keepY = c(5, 5))
  res_spls <- mixOmics::spls(
    spls_input[[1]],
    spls_input[[2]],
    ncomp = 2,
    mode = "regression",
    keepX = c(10, 10),
    keepY = c(5, 5)
  )

  expect_equal(
    attr(res, "datasets_name"),
    c("rnaseq", "metabolome")
  )

  res$call <- "removed from test"
  attr(res, "datasets_name") <- NULL
  res_spls$call <- "removed from test"
  expect_identical(res, res_spls)

  ## With multilevel option
  expect_no_error(
    spls_run(
      get_input_spls(multiomics_set, datasets = c("rnaseq", "metabolome"), mode = "regression", multilevel = "pheno_group")
    )
  )
  expect_no_error(
    spls_run(
      get_input_spls(multiomics_set, datasets = c("rnaseq", "metabolome"), mode = "regression", multilevel = c("pheno_group", "time", "name"))
    )
  )
})


# test_that("mixOmics pls and spls functions equivalent when no features selection", {
#   multiomics_set <- test_get_multidataset()
#   spls_input <- get_input_spls(multiomics_set, datasets = c("rnaseq", "metabolome"), mode = "regression")
#
#
#   ## Checking that with no features selection we get the same as if we had run mixOmics::pls
#   res_pls <- mixOmics::pls(
#     spls_input[[1]],
#     spls_input[[2]],
#     ncomp = 2,
#     mode = "regression"
#   )
#   res_spls <- mixOmics::spls(
#     spls_input[[1]],
#     spls_input[[2]],
#     ncomp = 2,
#     mode = "regression"
#   )
#
#   common_vals <- setdiff(
#     intersect(names(res_pls), names(res_spls)),
#     "call"
#   )
#   x_pls <- res_pls[common_vals]
#   class(x_pls) <- NULL
#
#   x_spls <- res_spls[common_vals]
#   class(x_spls) <- NULL
#
#   expect_identical(x_pls, x_spls)
# })


test_that("spls_get_optim_ncomp works", {

  make_input <- function(x){
    res <- list()
    res$measures$Q2.total$summary <- data.frame(
      comp = seq_along(x),
      mean = x
    )
    return(res)
  }

  ## If all `Q2` values are below the threshold specified with `thr`,
  ## the number of components yielding the highest `Q2` value is selected.
  expect_equal(
    spls_get_optim_ncomp(make_input(c(0.05, 0.06, 0.03))),
    2
  )

  ## If all `Q2` values are above the threshold, the number of components
  ## yielding the lowest `Q2` value is selected.
  expect_equal(
    spls_get_optim_ncomp(make_input(c(0.5, 0.6, 0.3))),
    3
  )

  ## If the `Q2` values are increasing, the number of components `n` is
  ## selected such that `n+1` is the smallest number of components with a `Q2`
  ## value above the threshold.
  expect_equal(
    spls_get_optim_ncomp(make_input(c(0.02, 0.07, 0.1, 0.4))),
    2
  )

  ## If the `Q2` values are decreasing, the number of components `n` is
  ## selected such that `n+1` is the smallest number of components with a `Q2`
  ## value below the threshold.
  expect_equal(
    spls_get_optim_ncomp(make_input(c(0.4, 0.1, 0.07, 0.05))),
    2
  )

  ## mountain shape
  expect_equal(
    spls_get_optim_ncomp(make_input(c(0.02, 0.05, 0.4, 0.1, 0.07))),
    2
  )

  ## valley shape
  expect_equal(
    spls_get_optim_ncomp(make_input(c(0.4, 0.1, 0.07, 0.05, 0.06, 0.099, 0.2))),
    2
  )
})

test_that("spls_tune works", {
  multiomics_set <- test_get_multidataset()
  spls_input <- get_input_spls(multiomics_set, datasets = c("rnaseq", "metabolome"), mode = "regression")

  ## Suppressing warnings about failure to converge
  expect_no_error(
    suppressWarnings(
      spls_tune(spls_input, ncomp = 2, keepX = c(5, 8, 10), keepY = c(5, 8, 10), validation = "Mfold", folds = 5, nrepeat = 1)
    )
  )

  set.seed(13)
  res <- suppressWarnings(
    spls_tune(
      spls_input,
      ncomp = 2,
      keepX = c(5, 8, 10),
      keepY = c(5, 8, 10),
      validation = "Mfold",
      folds = 5,
      nrepeat = 1
    )
  )
  set.seed(13)
  check <- suppressWarnings(
    mixOmics::tune.spls(
      spls_input[[1]],
      spls_input[[2]],
      ncomp = 2,
      mode = "regression",
      test.keepX = c(5, 8, 10),
      test.keepY = c(5, 8, 10),
      validation = "Mfold",
      folds = 5,
      nrepeat = 1
    )
  )

  expect_equal(
    attr(res, "datasets_name"),
    c("rnaseq", "metabolome")
  )

  res$call <- "removed from test"
  attr(res, "datasets_name") <- NULL
  check$call <- "removed from test"
  expect_identical(res, check)
})
