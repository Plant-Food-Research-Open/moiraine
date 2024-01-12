test_that("plot_correlation_matrix_full works", {
  mat <- matrix(runif(16), nrow = 4, dimnames = list(paste0("v", 1:4), paste0("v", 1:4)))

  expect_s3_class(
    plot_correlation_matrix_full(mat),
    "ggplot"
  )
})
