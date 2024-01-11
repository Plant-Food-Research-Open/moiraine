test_that("is_equal_or_null works", {

  expect_equal(is_equal_or_null(NULL, 3), TRUE)
  expect_equal(is_equal_or_null(3, 3), TRUE)
  expect_equal(is_equal_or_null(5, 3), FALSE)

  expect_equal(is_equal_or_null("test", "test"), TRUE)
  expect_equal(is_equal_or_null("something", "test"), FALSE)
})

test_that("hclust_matrix_rows works", {
  mat <- matrix(rnorm(12), nrow = 4, ncol = 3)
  res <- hclust_matrix_rows(mat)

  expect_s3_class(res, "dendrogram")
  expect_equal(length(stats::order.dendrogram(res)), 4)
})

test_that(".make_var_list works", {
  ## Expected input
  test <- 3
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC")),
    NA
  )
  test <- letters[1:4]
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC")),
    NA
  )


  test <- list("omicsA" = 1, "omicsB" = 2)
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC")),
    "'test' argument should either be a vector or a named list of length 3"
  )
  test <- list(1, 2, 3)
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC")),
    "'test' argument should be named; names should be: 'omicsA', 'omicsB', 'omicsC'."
  )
  test <- list("omicsA" = 1, "omicsB" = 2, "TEST" = 3)
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC")),
    "'test' argument: 'TEST' are not existing datasets or do not match the 'datasets' argument. Possible values are: 'omicsA', 'omicsB', 'omicsC'."
  )
  test <- list("omicsA" = 1, "omicsB" = 2, "omicsC" = 3)
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC")),
    NA
  )

  test <- 1:2
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC"), fixed_length = 1),
    "'test' argument values should have length 1."
  )
  test <- list("omicsA" = 1:2, "omicsB" = 2, "omicsC" = 3)
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC"), fixed_length = 1),
    "'test' argument values should have length 1."
  )
  test <- 1:3
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC"), fixed_length = 3),
    NA
  )
  test <- NULL
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC"), fixed_length = 3),
    NA
  )
  test <- list("omicsA" = 1, "omicsB" = 2, "omicsC" = 3)
  expect_error(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC"), fixed_length = 1),
    NA
  )

  ## Test output
  test <- NULL
  expect_equal(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC")),
    list("omicsA" = NULL, "omicsB" = NULL, "omicsC" = NULL)
  )
  test <- 3
  expect_equal(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC")),
    list("omicsA" = 3, "omicsB" = 3, "omicsC" = 3)
  )
  test <- letters[1:4]
  expect_equal(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC")),
    list("omicsA" = letters[1:4], "omicsB" = letters[1:4], "omicsC" = letters[1:4])
  )
  test <- list("omicsA" = 1, "omicsB" = 2, "omicsC" = 3)
  expect_equal(
    .make_var_list(test, c("omicsA", "omicsB", "omicsC")),
    list("omicsA" = 1, "omicsB" = 2, "omicsC" = 3)
  )
})

test_that(".check_input_var_fmetadata works", {
  multiomics_set <- test_get_multidataset()

  test <- list("snps+A" = "position", "rnaseq" = "TEST")
  expect_error(
    .check_input_var_fmetadata(test, multiomics_set),
    "'test' argument: 'TEST' is not a column in the features metadata for the rnaseq dataset. Possible values are:.+"
  )
  test <- list("snps+A" = "position", "rnaseq" = "chromosome")
  expect_error(
    .check_input_var_fmetadata(test, multiomics_set),
    NA
  )
  test <- list("snps+A" = NULL, "rnaseq" = "chromosome")
  expect_error(
    .check_input_var_fmetadata(test, multiomics_set),
    NA
  )
})


test_that(".check_input_var_smetadata works", {
  multiomics_set <- test_get_multidataset()

  test <- list("snps+A" = "pheno_group", "rnaseq" = "TEST")
  expect_error(
    .check_input_var_smetadata(test, multiomics_set),
    "'test' argument: 'TEST' is not a column in the samples metadata for the rnaseq dataset. Possible values are:.+"
  )
  test <- list("snps+A" = "pheno_group", "rnaseq" = "time")
  expect_error(
    .check_input_var_smetadata(test, multiomics_set),
    NA
  )
  test <- list("snps+A" = NULL, "rnaseq" = "time")
  expect_error(
    .check_input_var_smetadata(test, multiomics_set),
    NA
  )
})

test_that(".check_input_var_smetadata_common works", {
  multiomics_set <- test_get_multidataset()

  test <- "TEST"
  expect_error(
    .check_input_var_smetadata_common(test, multiomics_set),
    "'test' argument: 'TEST' is not a column in the samples metadata of any dataset. Possible values are:.+"
  )
  test <- "time"
  expect_error(
    .check_input_var_smetadata_common(test, multiomics_set),
    NA
  )
  test <- NULL
  expect_error(
    .check_input_var_smetadata_common(test, multiomics_set),
    NA
  )
})
