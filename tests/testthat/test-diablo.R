test_that("diablo_run works", {
  multiomics_set <- test_get_multidataset()[, 1:3]
  diablo_input <- get_input_mixomics_supervised(multiomics_set, "pheno_group")
  design_matrix <- .diablo_generate_design_matrix("weighted_full", names(multiomics_set))

  ## For a PLS-DA run
  expect_no_error(
    diablo_run(diablo_input, design = "weighted_full", ncomp = 2)
  )

  res <- diablo_run(diablo_input, design = "weighted_full", ncomp = 2)
  res_plsda <- mixOmics::block.plsda(
    diablo_input[names(multiomics_set)],
    diablo_input[["Y"]],
    ncomp = 2,
    design = design_matrix
  )

  res$call <- "removed from test"
  res_plsda$call <- "removed from test"
  expect_identical(res, res_plsda)

  ## For a sPLS run
  keepx_list <- list(
    c(10, 10),
    c(8, 8),
    c(6, 6)
  ) |>
    rlang::set_names(names(multiomics_set))

  expect_no_error(
    diablo_run(diablo_input, design = "weighted_full", ncomp = 2, keepX = keepx_list)
  )

  res <- diablo_run(diablo_input, design = "weighted_full", ncomp = 2, keepX = keepx_list)
  res_splsda <- mixOmics::block.splsda(
    diablo_input[names(multiomics_set)],
    diablo_input[["Y"]],
    ncomp = 2,
    design = design_matrix,
    keepX = keepx_list
  )

  res$call <- "removed from test"
  res_splsda$call <- "removed from test"
  expect_identical(res, res_splsda)
})
