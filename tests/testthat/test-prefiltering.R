test_that("feature_preselection_splsda_factory works", {
  tar_res <- feature_preselection_splsda_factory(
    mo_data,
    group = "pheno_group",
    to_keep_ns = list("snps+A" = 50, "metabolome" = 30)
  )

  ## Testing that we get a list of target stems
  expect_type(tar_res, "list")
  expect_s3_class(tar_res[[1]], "tar_stem")

  ## Testing targets name
  expect_equal(
    sapply(tar_res, function(x) {
      x$settings$name
    }),
    c("splsda_spec", "individual_splsda_input", "individual_splsda_perf", "individual_splsda_run", "filtered_set_slpsda")
  )

  ## Have to make this test separate and using strings because of a weird
  ## whitespace error I cannot fix
  expect_equal(
    tar_res[[1]]$command$expr |> deparse() |> paste0(collapse = ""),
    "expression(tar_group(dplyr::group_by(tibble::tibble(dsn = c(\"snps+A\", \"metabolome\"), tkn = list(50, 30), tkp = NULL), dsn)))"
  )

  # Testing targets command
  expect_equal(
    lapply(tar_res[-1], function(x) {
      x$command$expr
    }),
    list(
      expression(get_input_splsda(mo_data, splsda_spec$dsn, "pheno_group", NULL)),
      expression(perf_splsda(individual_splsda_input)),
      expression(run_splsda(individual_splsda_input, perf_res = individual_splsda_perf,
                            to_keep_n = splsda_spec$tkn, to_keep_prop = splsda_spec$tkp)),
      expression(get_filtered_dataset_splsda(mo_data, individual_splsda_run))
    )
  )

  ## Testing targets pattern
  expect_equal(
    sapply(tar_res, function(x) {
      x$settings$pattern
    }),
    list(
      NULL,
      expression(map(splsda_spec)),
      expression(map(individual_splsda_input)),
      expression(map(individual_splsda_input, individual_splsda_perf, splsda_spec)),
      NULL
    )
  )

  ## Testing iteration method (for dynamic branching it should be a list)
  expect_equal(
    sapply(tar_res, function(x) {
      x$settings$iteration
    }),
    c("group", "list", "list", "list", "vector")
  )

  ## Testing targets format
  expect_equal(
    sapply(tar_res, function(x) {
      x$settings$format
    }),
    c("rds", "rds", "rds", "rds", "rds")
  )

  ## Testing using the multilevel option
  tar_res2 <- feature_preselection_splsda_factory(
    mo_data,
    group = "pheno_group",
    to_keep_ns = list("snps+A" = 50, "metabolome" = 30),
    multilevel = "rep"
  )
  expect_equal(
    sapply(tar_res2[-2], function(x) {
      x$command$expr
    }),
    sapply(tar_res[-2], function(x) {
      x$command$expr
    })
  )
  expect_equal(
    tar_res2[[2]]$command$expr,
    expression(get_input_splsda(mo_data, splsda_spec$dsn, "pheno_group", "rep"))
  )

  ## Testing custom perf arguments
  tar_res3 <- feature_preselection_splsda_factory(
    mo_data,
    group = "pheno_group",
    to_keep_ns = list("snps+A" = 50, "metabolome" = 30),
    ncomp_max = 10,
    folds = 3
  )
  expect_equal(
    sapply(tar_res3[-3], function(x) {
      x$command$expr
    }),
    sapply(tar_res[-3], function(x) {
      x$command$expr
    })
  )
  expect_equal(
    tar_res3[[3]]$command$expr,
    expression(perf_splsda(individual_splsda_input, ncomp_max = 10, folds = 3))
  )

  ## Testing custom complete target name + custom prefix for target names
  tar_res4 <- feature_preselection_splsda_factory(
    mo_data,
    group = "pheno_group",
    to_keep_ns = list("snps+A" = 50, "metabolome" = 30),
    target_name_prefix = "test_",
    filtered_set_target_name = "final_test"
  )
  expect_equal(
    sapply(tar_res4, function(x) {
      x$settings$name
    }),
    c("test_splsda_spec", "test_individual_splsda_input", "test_individual_splsda_perf",
      "test_individual_splsda_run", "final_test")
  )
})
