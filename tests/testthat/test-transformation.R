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
    transform_logx(mat1, return_matrix_only = TRUE, log_base = 2),
    log2(zero_to_half_min(mat1))
  )
  expect_equal(
    transform_logx(mat1, return_matrix_only = TRUE, log_base = 10, pre_log_function = \(x){x + 1}),
    log10(mat1 + 1)
  )
  expect_warning(
    transform_logx(mat1, return_matrix_only = TRUE, log_base = 10, pre_log_function = \(x){x}),
    "The matrix contains zero values; log-transformation will yield `-Inf`."
  )

  expect_equal(
    transform_logx(mat1, return_matrix_only = FALSE, log_base = 2),
    list(
      transformed_data = log2(zero_to_half_min(mat1)),
      info_transformation = list(log_base = 2, pre_log_function = zero_to_half_min),
      transformation = "log2"
    )
  )
  expect_equal(
    transform_logx(mat1, return_matrix_only = FALSE, log_base = 10, pre_log_function = \(x){x + 1}),
    list(
      transformed_data = log10(mat1 + 1),
      info_transformation = list(log_base = 10, pre_log_function = \(x){x + 1}),
      transformation = "log10"
    )
  )
})

test_that("zero_to_half_min works", {
  mat1 <- matrix(c(0, 1, 2, 3), nrow = 2)
  mat2 <- matrix(c(2, 1, 2, 3), nrow = 2)

  expect_equal(
    zero_to_half_min(mat1),
    matrix(c(0.5, 1, 2, 3), nrow = 2)
  )

  expect_equal(
    zero_to_half_min(mat2),
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
    "'transformation' argument: 'TEST' is not a recognised transformation. Possible values are: 'vsn', 'vst-deseq2', 'logx', 'best-normalize-auto', 'best-normalize-manual'."
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

test_that("transformation_datasets_factory works - default", {
  tar_res <- transformation_datasets_factory(
    mo_set,
    c("rnaseq" = "vst", "metabolome" = "vsn")
  )

  expect_type(tar_res, "list")
  expect_equal(
    purrr::map_chr(tar_res, \(x) x$settings$name),
    c("transformations_spec", "transformations_runs_list", "transformed_set")
  )
  expect_s3_class(tar_res[[1]], "tar_stem")
  expect_s3_class(tar_res[[2]], "tar_pattern")
  expect_s3_class(tar_res[[3]], "tar_stem")

  expect_equal(
    tar_res[[1]]$command$expr |> test_clean_expression(),
    expression(tibble::tibble(
          dsn = c("rnaseq", "metabolome"),
          transf = c("vst", "vsn"),
          meth = list(NULL, NULL),
          log_b = list(NULL, NULL),
          prelog_f = list(NULL, NULL)
        ) |>
          dplyr::group_by(dsn) |>
          targets::tar_group()) |> test_clean_expression()
  )
  expect_equal(
    tar_res[[2]]$command$expr |> test_clean_expression(),
    expression(transform_dataset(
          mo_set,
          dataset = transformations_spec$dsn,
          transformation = transformations_spec$transf,
          return_matrix_only = FALSE,
          method = transformations_spec$meth[[1]],
          log_base = transformations_spec$log_b[[1]],
          pre_log_function = transformations_spec$prelog_f[[1]]
        )) |> test_clean_expression()
  )
  expect_equal(
    tar_res[[3]]$command$expr,
    str2expression("get_transformed_data(mo_set, transformations_runs_list)")
  )

  ## Adding prefix
  tar_res <- transformation_datasets_factory(
    mo_set,
    c("rnaseq" = "vst", "metabolome" = "vsn"),
    target_name_prefix = "TEST_"
  )
  expect_type(tar_res, "list")
  expect_equal(
    purrr::map_chr(tar_res, \(x) x$settings$name),
    c("TEST_transformations_spec", "TEST_transformations_runs_list", "TEST_transformed_set")
  )
  expect_s3_class(tar_res[[1]], "tar_stem")
  expect_s3_class(tar_res[[2]], "tar_pattern")
  expect_s3_class(tar_res[[3]], "tar_stem")

  ## Changing final set name
  tar_res <- transformation_datasets_factory(
    mo_set,
    c("rnaseq" = "vst", "metabolome" = "vsn"),
    transformed_data_name = "TEST"
  )
  expect_type(tar_res, "list")
  expect_equal(
    purrr::map_chr(tar_res, \(x) x$settings$name),
    c("transformations_spec", "transformations_runs_list", "TEST")
  )
  expect_s3_class(tar_res[[1]], "tar_stem")
  expect_s3_class(tar_res[[2]], "tar_pattern")
  expect_s3_class(tar_res[[3]], "tar_stem")
})

test_that("transformation_datasets_factory works - logx", {
  tar_res <- transformation_datasets_factory(
    mo_set,
    c("rnaseq" = "vst", "metabolome" = "logx")
  )

  expect_type(tar_res, "list")
  expect_equal(tar_res[[1]]$settings$name, "transformations_spec")
  expect_s3_class(tar_res[[1]], "tar_stem")

  expect_equal(
    tar_res[[1]]$command$expr |> test_clean_expression(),
    expression(
      tibble::tibble(
          dsn = c("rnaseq", "metabolome"),
          transf = c("vst", "logx"),
          meth = list(NULL, NULL),
          log_b = list(NULL, 2),
          prelog_f = list(NULL, function(mat) {
            if (!any(mat == 0)) {
              return(mat)
            }

            min_val <- min(mat[mat != 0])
            mat[mat == 0] <- min_val / 2

            return(mat)
          })
        ) |>
          dplyr::group_by(dsn) |>
          targets::tar_group()) |> test_clean_expression()
  )

  tar_res <- transformation_datasets_factory(
    mo_set,
    c("rnaseq" = "logx", "metabolome" = "logx"),
    log_bases = list(rnaseq = 10, metabolome = 2),
    pre_log_functions = list(
      rnaseq = \(x) x + 0.5,
      metabolome = \(x) x + 1
    )
  )

  expect_equal(
    tar_res[[1]]$command$expr |> test_clean_expression(),
    expression(tibble::tibble(
      dsn = c("rnaseq", "metabolome"),
      transf = c("logx", "logx"),
      meth = list(NULL, NULL),
      log_b =  list(10, 2),
      prelog_f = list(
        \(x) x + 0.5,
        \(x) x + 1
      )
    ) |>
      dplyr::group_by(dsn) |>
      targets::tar_group()) |> test_clean_expression()
  )
})
