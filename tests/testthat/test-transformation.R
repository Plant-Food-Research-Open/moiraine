test_that("transform_vsn works", {
  if (!requireNamespace("vsn", quietly = TRUE)) {
    expect_error(transform_vsn(), "Package \"vsn\" must be installed to use this function.")
    skip("Package \"vsn\" must be installed to use this function.")
  }

  dat <- matrix(rpois(50 * 10, 50), nrow = 50, ncol = 10, dimnames = list(paste0("feature_", 1:50), paste0("sample_", 1:10)))

  ## Testing matrix output
  t1 <- suppressMessages(transform_vsn(dat, return_matrix_only = TRUE))

  expect_true(is.matrix(t1))
  expect_equal(rownames(t1), rownames(dat))
  expect_equal(colnames(t1), colnames(dat))

  ## Testing list output
  t2 <- suppressMessages(transform_vsn(dat, return_matrix_only = FALSE))

  expect_true(is.list(t2))
  expect_equal(names(t2), c("transformed_data", "info_transformation", "transformation"))
  expect_equal(t2$transformed_data, t1)
  expect_true(is.null(t2$info_transformation))
  expect_equal(t2$transformation, "vsn")
})

test_that("transform_vst works", {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    expect_error(transform_vst(), "Package \"DESeq2\" must be installed to use this function.")
    skip("Package \"DESeq2\" must be installed to use this function.")
  }

  dat <- matrix(rpois(50 * 10, 50), nrow = 50, ncol = 10, dimnames = list(paste0("feature_", 1:50), paste0("sample_", 1:10)))

  ## Testing matrix output
  t1 <- suppressMessages(transform_vst(dat, return_matrix_only = TRUE))

  expect_true(is.matrix(t1))
  expect_equal(rownames(t1), rownames(dat))
  expect_equal(colnames(t1), colnames(dat))

  ## Testing list output
  t2 <- suppressMessages(transform_vst(dat, return_matrix_only = FALSE))

  expect_true(is.list(t2))
  expect_equal(names(t2), c("transformed_data", "info_transformation", "transformation"))
  expect_equal(t2$transformed_data, t1)
  expect_s4_class(t2$info_transformation, "DESeqTransform")
  expect_equal(t2$transformation, "vst-deseq2")
})


test_that("transform_bestNormalise_auto works", {
  if (!requireNamespace("bestNormalize", quietly = TRUE)) {
    expect_error(transform_bestNormalise_auto(), "Package \"bestNormalize\" must be installed to use this function.")
    skip("Package \"bestNormalize\" must be installed to use this function.")
  }

  dat <- matrix(rnorm(10 * 50), nrow = 10, ncol = 50, dimnames = list(paste0("feature_", 1:10), paste0("sample_", 1:50)))

  ## Testing matrix output
  set.seed(1)
  t1 <- transform_bestNormalise_auto(dat, return_matrix_only = TRUE)

  expect_true(is.matrix(t1))
  expect_equal(rownames(t1), rownames(dat))
  expect_equal(colnames(t1), colnames(dat))

  ## Testing list output
  set.seed(1)
  t2 <- transform_bestNormalise_auto(dat, return_matrix_only = FALSE)

  expect_true(is.list(t2))
  expect_equal(names(t2), c("transformed_data", "info_transformation", "transformation"))
  expect_equal(t2$transformed_data, t1)
  expect_true(is.list(t2$info_transformation))
  temp <- rep("bestNormalize", 10)
  names(temp) <- paste0("feature_", 1:10)
  expect_equal(sapply(t2$info_transformation, class), temp)
  expect_equal(t2$transformation, "best-normalize-auto")
})

test_that("transform_bestNormalise_manual works", {
  if (!requireNamespace("bestNormalize", quietly = TRUE)) {
    expect_error(transform_bestNormalise_auto(), "Package \"bestNormalize\" must be installed to use this function.")
    skip("Package \"bestNormalize\" must be installed to use this function.")
  }

  dat <- matrix(rnorm(10 * 50), nrow = 10, ncol = 50, dimnames = list(paste0("feature_", 1:10), paste0("sample_", 1:50)))

  expect_error(transform_bestNormalise_manual(dat, "TEST"),
               "'method' argument: 'TEST' not a valid method. Possible methods are: 'arcsinh_x', 'boxcox', 'log_x', 'sqrt_x', 'yeojohnson', 'center_scale', 'exp_x', 'orderNorm'.",
               fixed = TRUE
  )

  ## Testing matrix output
  t1 <- transform_bestNormalise_manual(dat, "center_scale", return_matrix_only = TRUE)

  expect_true(is.matrix(t1))
  expect_equal(rownames(t1), rownames(dat))
  expect_equal(colnames(t1), colnames(dat))

  ## Testing list output
  t2 <- transform_bestNormalise_manual(dat, "center_scale", return_matrix_only = FALSE)

  expect_true(is.list(t2))
  expect_equal(names(t2), c("transformed_data", "info_transformation", "transformation"))
  expect_equal(t2$transformed_data, t1)
  expect_true(is.list(t2$info_transformation))
  temp <- rep("center_scale", 10)
  names(temp) <- paste0("feature_", 1:10)
  expect_equal(sapply(t2$info_transformation, class)[1, ], temp)
  expect_equal(t2$transformation, "best-normalize-manual-center_scale")
})

test_that("transform_logx works", {

  mat1 <- matrix(c(0, 1, 2, 3), nrow = 2)

  expect_equal(
    transform_logx(mat1, return_matrix_only = TRUE, base = 2),
    log2(mat1 + 0.5)
  )
  expect_equal(
    transform_logx(mat1, return_matrix_only = TRUE, base = 10, pre_log_function = \(x){x + 1}),
    log10(mat1 + 1)
  )
  expect_warning(
    transform_logx(mat1, return_matrix_only = TRUE, base = 10, pre_log_function = NULL),
    "The matrix contains zero values; log-transformation will yield `-Inf`."
  )

  expect_equal(
    transform_logx(mat1, return_matrix_only = FALSE, base = 2),
    list(
      transformed_data = log2(mat1 + 0.5),
      info_transformation = list(log_base = 2, pre_log_function = offset_half_min),
      transformation = "log2"
    )
  )
  expect_equal(
    transform_logx(mat1, return_matrix_only = FALSE, base = 10, pre_log_function = \(x){x + 1}),
    list(
      transformed_data = log10(mat1 + 1),
      info_transformation = list(log_base = 10, pre_log_function = \(x){x + 1}),
      transformation = "log10"
    )
  )
})

test_that("offset_half_min works", {
  mat1 <- matrix(c(0, 1, 2, 3), nrow = 2)
  mat2 <- matrix(c(2, 1, 2, 3), nrow = 2)

  expect_equal(
    offset_half_min(mat1),
    matrix(c(0.5, 1.5, 2.5, 3.5), nrow = 2)
  )

  expect_equal(
    offset_half_min(mat2),
    mat2
  )
})

test_that(".get_transformed_matrix works", {
  if (!requireNamespace("bestNormalize", quietly = TRUE)) {
    skip("Package \"bestNormalize\" must be installed.")
  }

  dat <- matrix(rnorm(10 * 50), nrow = 10, ncol = 50, dimnames = list(paste0("feature_", 1:10), paste0("sample_", 1:50)))

  transf1 <- transform_bestNormalise_manual(dat, "center_scale", return_matrix_only = TRUE)
  transf2 <- transform_bestNormalise_manual(dat, "center_scale", return_matrix_only = FALSE)

  expect_equal(.get_transformed_matrix(transf1, TRUE), transf1)
  expect_equal(.get_transformed_matrix(transf2, FALSE), transf1)
})

test_that("transform_dataset works", {
  multiomics_set <- test_get_multidataset()

  ## Testing input
  expect_error(transform_dataset("TEST"), "Expecting MultiDataSet object")
  expect_error(
    transform_dataset(multiomics_set, "TEST"),
    "'TEST' datasets are not present in mo_data. Possible dataset names are:.+"
  )
  expect_error(
    transform_dataset(multiomics_set, "snps+A", "TEST"),
    "'transformation' argument: 'TEST' is not a recognised transformation. Possible values are: 'vsn', 'vst-deseq2', 'best-normalize-auto', 'best-normalize-manual'."
  )
  expect_error(transform_dataset(multiomics_set, "snps+A", "best-normalize-manual"), "'method' argument should be provided for 'best-normalize-manual' transformation.")

  res <- transform_dataset(
    multiomics_set,
    "rnaseq",
    "best-normalize-manual",
    method = "center_scale",
    return_multidataset = TRUE
  )
  expect_s4_class(res, "MultiDataSet")
  expect_equal(names(res), c("snps+A", "rnaseq", "metabolome", "phenotypes"))

  # expect_message(transform_dataset(multiomics_set, dataset = "snps+A", transformation = "best-normalize-auto", method = "center_scale"))
})
