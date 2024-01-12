test_that("so2pls_crossval_o2m_adjR2 works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  expect_error(so2pls_crossval_o2m_adjR2("TEST"), "'omicspls_input' should be a list of length 2 produced by get_input_omicspls().")
  expect_error(so2pls_crossval_o2m_adjR2(list("test")), "'omicspls_input' should be a list of length 2 produced by get_input_omicspls().")

  expect_error(so2pls_crossval_o2m_adjR2(omicspls_input, a = 1:2, ax = 0:1, ay = 0:1, nr_folds = 2), NA)


  set.seed(1)
  cv_adj_res <- suppressMessages(OmicsPLS::crossval_o2m_adjR2(omicspls_input[[1]],
                                                              omicspls_input[[2]],
                                                              a = 1:2,
                                                              ax = 0:1,
                                                              ay = 0:1,
                                                              nr_folds = 2
  ))

  set.seed(1)
  res <- suppressMessages(so2pls_crossval_o2m_adjR2(omicspls_input,
                                                    a = 1:2,
                                                    ax = 0:1,
                                                    ay = 0:1,
                                                    nr_folds = 2
  ))

  expect_equal(res, cv_adj_res, ignore_attr = TRUE)
  expect_equal(attr(res, "datasets_name"), names(omicspls_input))
})


test_that("so2pls_plot_cv_adj works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  cv_adj_res <- suppressMessages(OmicsPLS::crossval_o2m_adjR2(omicspls_input[[1]],
                                                              omicspls_input[[2]],
                                                              a = 1:2,
                                                              ax = 0:1,
                                                              ay = 0:1,
                                                              nr_folds = 2
  ))

  expect_error(so2pls_plot_cv_adj("TEST"),
               "Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.",
               fixed = TRUE
  )

  expect_error(so2pls_plot_cv_adj(data.frame(x = 1:2, y = 1:2)),
               "Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.",
               fixed = TRUE
  )

  expect_s3_class(so2pls_plot_cv_adj(cv_adj_res), "ggplot")
})

test_that("so2pls_print_cv_adj works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  cv_adj_res <- suppressMessages(OmicsPLS::crossval_o2m_adjR2(omicspls_input[[1]],
                                                              omicspls_input[[2]],
                                                              a = 1:2,
                                                              ax = 0:1,
                                                              ay = 0:1,
                                                              nr_folds = 2
  ))

  expect_error(so2pls_print_cv_adj("TEST"),
               "Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.",
               fixed = TRUE
  )

  expect_error(so2pls_print_cv_adj(data.frame(x = 1:2, y = 1:2)),
               "Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.",
               fixed = TRUE
  )

  res <- so2pls_print_cv_adj(cv_adj_res)
  expect_s3_class(res, "tbl_df")
  expect_equal(ncol(res), 4)
})


test_that("so2pls_get_optim_ncomp_adj works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  cv_adj_res <- suppressMessages(OmicsPLS::crossval_o2m_adjR2(omicspls_input[[1]],
                                                              omicspls_input[[2]],
                                                              a = 1:2,
                                                              ax = 0:1,
                                                              ay = 0:1,
                                                              nr_folds = 2
  ))

  expect_error(so2pls_get_optim_ncomp_adj("TEST"),
               "Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.",
               fixed = TRUE
  )

  expect_error(so2pls_get_optim_ncomp_adj(data.frame(x = 1:2, y = 1:2)),
               "Input should be a data-frame with columns 'MSE', 'n', 'nx' and 'ny'. Make sure input is the output from the OmicsPLS::crossval_o2m_adjR2() function.",
               fixed = TRUE
  )

  res <- so2pls_get_optim_ncomp_adj(cv_adj_res)

  expect_equal(length(res), 3)
  expect_equal(names(res), c("n", "nx", "ny"))
})

test_that("so2pls_crossval_o2m works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  set.seed(1)
  cv_adj_res <- suppressMessages(OmicsPLS::crossval_o2m_adjR2(omicspls_input[[1]],
                                                              omicspls_input[[2]],
                                                              a = 1:2,
                                                              ax = 0:1,
                                                              ay = 0:1,
                                                              nr_folds = 2
  ))

  set.seed(1)
  cv_res <- suppressMessages(OmicsPLS::crossval_o2m(omicspls_input[[1]],
                                                    omicspls_input[[2]],
                                                    a = 2:3,
                                                    ax = 2:3,
                                                    ay = 2:3,
                                                    nr_folds = 2
  ))

  expect_error(so2pls_crossval_o2m("TEST"), "'omicspls_input' should be a list of length 2 produced by get_input_omicspls().")
  expect_error(so2pls_crossval_o2m(list("test")), "'omicspls_input' should be a list of length 2 produced by get_input_omicspls().")

  expect_error(so2pls_crossval_o2m(omicspls_input, a = 1:2, ax = 0:1, ay = 0:1, nr_folds = 2), NA)

  set.seed(1)
  res1 <- suppressMessages(so2pls_crossval_o2m(omicspls_input,
                                               a = 2:3,
                                               ax = 2:3,
                                               ay = 2:3,
                                               nr_folds = 2
  ))

  expect_s3_class(res1, "cvo2m")
  expect_equal(
    dimnames(res1$Original),
    list(
      c("ax=2", "ax=3"),
      c("ay=2", "ay=3"),
      c("a=2", "a=3")
    )
  )
  ## need to remove the time element whic
  expect_equal(res1[c("Original", "Sorted", "kcv")], cv_res[c("Original", "Sorted", "kcv")], ignore_attr = TRUE)

  res2 <- suppressMessages(so2pls_crossval_o2m(omicspls_input,
                                               cv_adj_res = cv_adj_res,
                                               nr_folds = 2
  ))

  expect_s3_class(res2, "cvo2m")
  expect_equal(
    dimnames(res2$Original),
    list(
      c("ax=0", "ax=1"),
      c("ay=0", "ay=1"),
      c("a=1", "a=2", "a=3")
    )
  )
})

test_that("so2pls_plot_cv works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  cv_res <- suppressMessages(OmicsPLS::crossval_o2m(omicspls_input[[1]],
                                                    omicspls_input[[2]],
                                                    a = 1:2,
                                                    ax = 0:1,
                                                    ay = 0:1,
                                                    nr_folds = 2
  ))

  expect_error(so2pls_plot_cv("TEST"),
               "Expecting a cvo2m object. Make sure input is the output from the so2pls_crossval_o2m() or OmicsPLS::crossval_o2m() function.",
               fixed = TRUE
  )

  expect_error(so2pls_plot_cv(cv_res), NA)
  expect_s3_class(so2pls_plot_cv(cv_res), "ggplot")

  cv_res2 <- suppressMessages(so2pls_crossval_o2m(omicspls_input,
                                                  a = 1:2,
                                                  ax = 0:1,
                                                  ay = 0:1,
                                                  nr_folds = 2
  ))

  expect_error(so2pls_plot_cv(cv_res2), NA)
  expect_s3_class(so2pls_plot_cv(cv_res2), "ggplot")
})

test_that("so2pls_get_optim_ncomp works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  cv_res <- suppressMessages(OmicsPLS::crossval_o2m(omicspls_input[[1]],
                                                    omicspls_input[[2]],
                                                    a = 1:2,
                                                    ax = 0:1,
                                                    ay = 0:1,
                                                    nr_folds = 2
  ))

  expect_error(so2pls_get_optim_ncomp("TEST"),
               "Expecting a cvo2m object. Make sure input is the output from the so2pls_crossval_o2m() or OmicsPLS::crossval_o2m() function.",
               fixed = TRUE
  )

  res <- so2pls_get_optim_ncomp(cv_res)

  expect_equal(length(res), 3)
  expect_equal(names(res), c("n", "nx", "ny"))
})

test_that("so2pls_crossval_sparsity works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  set.seed(1)
  ref <- suppressMessages(OmicsPLS::crossval_sparsity(omicspls_input[[1]],
                                                      omicspls_input[[2]],
                                                      n = 2,
                                                      nx = 1,
                                                      ny = 1,
                                                      nr_folds = 2,
                                                      keepx_seq = c(5, 10),
                                                      keepy_seq = c(5, 10)
  ))

  set.seed(1)
  res <- suppressMessages(so2pls_crossval_sparsity(omicspls_input,
                                                   n = 2,
                                                   nx = 1,
                                                   ny = 1,
                                                   nr_folds = 2,
                                                   keepx_seq = c(5, 10),
                                                   keepy_seq = c(5, 10)
  ))

  expect_equal(res$Best, ref$Best)
  expect_equal(length(res$Covs), 2)
  expect_equal(length(res$SEcov), 2)
  expect_equal(attr(res, "datasets_name"), names(omicspls_input))
})


test_that("so2pls_get_optim_keep works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  cv_res <- suppressMessages(so2pls_crossval_sparsity(omicspls_input,
                                                      n = 3,
                                                      nx = 1,
                                                      ny = 1,
                                                      nr_folds = 2,
                                                      keepx_seq = c(5, 10, 15),
                                                      keepy_seq = c(1, 2, 3)
  ))

  expect_equal(
    sapply(so2pls_get_optim_keep(cv_res), length),
    c("keepx" = 3, "keepy" = 3)
  )
  expect_equal(
    sapply(so2pls_get_optim_keep(cv_res, use_1sd_rule = FALSE), length),
    c("keepx" = 3, "keepy" = 3)
  )
})


test_that("so2pls_print_cv_sparsity works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("rnaseq", "metabolome")) ## Use the same dataset to get more common samples

  cv_res <- suppressMessages(so2pls_crossval_sparsity(omicspls_input,
                                                      n = 3,
                                                      nx = 1,
                                                      ny = 1,
                                                      nr_folds = 2,
                                                      keepx_seq = c(5, 10, 15),
                                                      keepy_seq = c(1, 2, 3)
  ))

  cv_res_list <- so2pls_get_optim_keep(cv_res)

  expect_error(so2pls_print_cv_sparsity("TEST"), "'cv_res_optim' should be a named list of length 2 as returned by so2pls_get_optim_keep() function.", fixed = TRUE)

  res <- so2pls_print_cv_sparsity(cv_res_list)

  expect_equal(
    res |>
      dplyr::filter(dataset == "rnaseq") |>
      tidyr::pivot_longer(
        cols = starts_with("Joint"),
        names_to = "joint_component",
        values_to = "n_features"
      ) |>
      dplyr::arrange(joint_component) |>
      dplyr::pull(n_features),
    unname(cv_res_list[["keepx"]])
  )

  expect_equal(
    res |>
      dplyr::filter(dataset == "metabolome") |>
      tidyr::pivot_longer(
        cols = starts_with("Joint"),
        names_to = "joint_component",
        values_to = "n_features"
      ) |>
      dplyr::arrange(joint_component) |>
      dplyr::pull(n_features),
    unname(cv_res_list[["keepy"]])
  )
})


test_that("so2pls_plot_cv_sparsity works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  cv_res <- suppressMessages(so2pls_crossval_sparsity(omicspls_input,
                                                      n = 3,
                                                      nx = 1,
                                                      ny = 1,
                                                      nr_folds = 2,
                                                      keepx_seq = c(5, 10, 15),
                                                      keepy_seq = c(1, 2, 3)
  ))

  expect_s3_class(so2pls_plot_cv_sparsity(cv_res), "ggplot")
})


test_that("so2pls_o2m works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("metabolome", "metabolome")) ## Use the same dataset to get more common samples

  set.seed(1)
  ## Crossval adj
  cv_adj_res <- suppressMessages(OmicsPLS::crossval_o2m_adjR2(omicspls_input[[1]],
                                                              omicspls_input[[2]],
                                                              a = 1:2,
                                                              ax = 0:1,
                                                              ay = 0:1,
                                                              nr_folds = 2
  ))
  cv_adj_optim <- so2pls_get_optim_ncomp_adj(cv_adj_res)


  ## Crossval
  cv_res <- suppressMessages(OmicsPLS::crossval_o2m(omicspls_input[[1]],
                                                    omicspls_input[[2]],
                                                    a = 2:3,
                                                    ax = 2:3,
                                                    ay = 2:3,
                                                    nr_folds = 2
  ))
  cv_optim <- so2pls_get_optim_ncomp(cv_res)

  # Crossval sparsity - there are singularity issues so using dummy data
  cv_sparsity_res <- suppressMessages(OmicsPLS::crossval_sparsity(omicspls_input[[1]],
                                                                  omicspls_input[[2]],
                                                                  n = cv_optim[["n"]],
                                                                  nx = cv_optim[["nx"]],
                                                                  ny = cv_optim[["ny"]],
                                                                  nr_folds = 2,
                                                                  keepx_seq = c(8, 10),
                                                                  keepy_seq = c(8, 10)
  ))
  cv_sparsity_optim <- so2pls_get_optim_keep(cv_sparsity_res)

  ## Test on input
  expect_error(so2pls_o2m("test"), "'omicspls_input' should be a list of length 2 produced by get_input_omicspls().")
  expect_error(so2pls_o2m(omicspls_input, cv_res = "TEST"), "'cv_res' argument should be a named vector of length 3 with names: 'n', 'nx', 'ny'.")
  expect_error(so2pls_o2m(omicspls_input, cv_res = cv_optim), NA)
  expect_error(so2pls_o2m(omicspls_input, cv_res = cv_adj_optim), NA)
  expect_error(so2pls_o2m(omicspls_input, cv_res = cv_optim, sparsity_res = "TEST"), "'sparsity_res' argument should be a named list of length 2 with names: 'keepx', 'keepy'.")
  expect_error(so2pls_o2m(omicspls_input), "Need to provide either a cross-validation result through 'cv_res' argument or a value for 'n' argument.")
  expect_error(so2pls_o2m(omicspls_input, n = 2), "Need to provide either a cross-validation result through 'cv_res' argument or a value for 'nx' argument.")
  expect_error(so2pls_o2m(omicspls_input, n = 2, nx = 1), "Need to provide either a cross-validation result through 'cv_res' argument or a value for 'ny' argument.")
  expect_error(so2pls_o2m(omicspls_input, n = 2, nx = 1, ny = 1), NA)
  expect_error(so2pls_o2m(omicspls_input, cv_res = cv_optim, sparse = TRUE), "'sparse' = TRUE: need to provide either a sparsity cross-validation result through 'sparsity_res' argument or a value for 'keepx' argument.")
  expect_error(so2pls_o2m(omicspls_input, cv_res = cv_optim, sparse = TRUE, keepx = 2), "'sparse' = TRUE: need to provide either a sparsity cross-validation result through 'sparsity_res' argument or a value for 'keepy' argument.")
  expect_error(so2pls_o2m(omicspls_input, cv_res = cv_optim, sparse = TRUE, keepx = 2, keepy = 2), NA)
  expect_error(so2pls_o2m(omicspls_input, cv_res = cv_optim, sparsity_res = cv_sparsity_optim), NA)



  o2m_res <- suppressMessages(OmicsPLS::o2m(omicspls_input[[1]],
                                            omicspls_input[[2]],
                                            n = cv_optim[["n"]],
                                            nx = cv_optim[["nx"]],
                                            ny = cv_optim[["ny"]],
                                            sparse = TRUE,
                                            keepx = cv_sparsity_optim[["keepx"]],
                                            keepy = cv_sparsity_optim[["keepy"]]
  ))


  res <- suppressMessages(so2pls_o2m(
    omicspls_input,
    cv_optim,
    cv_sparsity_optim
  ))

  ## for comparison
  o2m_res$flags$call <- o2m_res$flags$time <- NULL
  res$flags$call <- res$flags$time <- NULL

  expect_equal(res, o2m_res, ignore_attr = TRUE)
})

test_that("so2pls_get_components works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("rnaseq", "metabolome"))

  so2pls_res <- so2pls_o2m(
    omicspls_input,
    n = 3,
    nx = 2,
    ny = 1
  )

  expect_equal(
    get_output_so2pls(so2pls_res, use_average_dimensions = TRUE) |>
      so2pls_get_components(),
    list(
      "joint" = paste0("joint component ", 1:3),
      "specific" = list(
        "rnaseq" = paste0("rnaseq specific component ", 1:2),
        "metabolome" = "metabolome specific component 1"
      )
    )
  )

  expect_equal(
    get_output_so2pls(so2pls_res, use_average_dimensions = FALSE) |>
      so2pls_get_components(),
    list(
      "joint" = paste0(c("rnaseq", "metabolome"), " joint component ", rep(1:3, each = 2)),
      "specific" = list(
        "rnaseq" = paste0("rnaseq specific component ", 1:2),
        "metabolome" = "metabolome specific component 1"
      )
    )
  )

  so2pls_res <- so2pls_o2m(
    omicspls_input,
    n = 3,
    nx = 2,
    ny = 0
  )

  expect_equal(
    get_output_so2pls(so2pls_res, use_average_dimensions = TRUE) |>
      so2pls_get_components(),
    list(
      "joint" = paste0("joint component ", 1:3),
      "specific" = list(
        "rnaseq" = paste0("rnaseq specific component ", 1:2)
      )
    )
  )

  so2pls_res <- so2pls_o2m(
    omicspls_input,
    n = 3,
    nx = 0,
    ny = 0
  )
  res <- get_output_so2pls(so2pls_res, use_average_dimensions = TRUE) |>
    so2pls_get_components()
  expect_length(res[["specific"]], 0)
})

test_that("so2pls_plot_summary works", {
  o2m_res <- test_get_so2pls_run()
  multiomics_set <- test_get_multidataset()

  expect_s3_class(so2pls_plot_summary(o2m_res), "patchwork")
  expect_s3_class(so2pls_plot_summary(o2m_res, datasets = "rnaseq"), "patchwork")

  attr(o2m_res, "datasets_name") <- NULL

  expect_s3_class(so2pls_plot_summary(o2m_res), "patchwork")
  expect_s3_class(so2pls_plot_summary(o2m_res, datasets = "X"), "patchwork")
})

test_that("so2pls_plot_joint_components_coefficients works", {
  multiomics_set <- test_get_multidataset()
  omicspls_input <- get_input_omicspls(multiomics_set, c("rnaseq", "metabolome"))

  o2m_res <- suppressMessages(so2pls_o2m(omicspls_input,
                                         n = 2,
                                         nx = 1,
                                         ny = 1,
                                         sparse = TRUE,
                                         keepx = c(5, 5),
                                         keepy = c(5, 5)
  ))

  expect_error(so2pls_plot_joint_components_coefficients("TEST"), "Expecting a o2m object. Make sure input is the output from the OmicsPLS::o2m() function.", fixed = TRUE)
  expect_error(so2pls_plot_joint_components_coefficients(o2m_res), NA)
  expect_error(so2pls_plot_joint_components_coefficients(o2m_res, "rnaseq"), NA)
  expect_error(so2pls_plot_joint_components_coefficients(o2m_res, "X"), "'datasets' argument: 'X' are not existing datasets. Possible values are: 'rnaseq', 'metabolome'.")
})

test_that("so2pls_plot_samples_specific_components works", {
  output_so2pls <- test_get_so2pls_run()

  expect_length(so2pls_plot_samples_specific_components(output_so2pls), 2)
  expect_s3_class(so2pls_plot_samples_specific_components(output_so2pls, dataset = "rnaseq"), "ggplot")
  expect_error(
    so2pls_plot_samples_specific_components(output_so2pls, dataset = "TEST"),
    "'dataset' argument: 'TEST' is not a valid dataset name. Valid names are:.+"
  )
})
