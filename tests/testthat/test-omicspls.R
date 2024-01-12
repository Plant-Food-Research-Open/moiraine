test_that("get_input_omicspls works", {
  multiomics_set <- test_get_multidataset()

  ## Testing input parameters
  expect_error(get_input_omicspls("TEST"), "Expecting MultiDataSet object")
  expect_error(get_input_omicspls(multiomics_set, "TEST"),
               "'TEST' datasets are not present in mo_data. Possible dataset names are: 'snps+A', 'rnaseq', 'metabolome', 'phenotypes'.",
               fixed = TRUE
  )
  expect_error(get_input_omicspls(multiomics_set), "In 'datasets' argument: expecting a vector of length 2 (as OmicsPLS integrates 2 datasets).", fixed = TRUE)

  ds_list <- get_datasets(multiomics_set)[c("snps+A", "metabolome")]
  res1 <- get_input_omicspls(multiomics_set, c("snps+A", "metabolome"))
  res2 <- get_input_omicspls(multiomics_set, c("snps+A", "metabolome"), scale_data = TRUE)

  ## Testing output
  expect_type(res1, "list")
  expect_equal(names(res1), names(ds_list))
  expect_equal(sapply(res1, ncol), sapply(ds_list, nrow))
  expect_equal(unname(sapply(res1, nrow)), rep(5, 2))

  expect_equal(
    res1[["snps+A"]],
    scale(t(ds_list[["snps+A"]][, MultiDataSet::commonIds(multiomics_set[, c("snps+A", "metabolome")])]), center = TRUE, scale = FALSE)
  )
  expect_equal(
    res2[["snps+A"]],
    scale(t(ds_list[["snps+A"]][, MultiDataSet::commonIds(multiomics_set[, c("snps+A", "metabolome")])]), center = TRUE, scale = TRUE)
  )
})
