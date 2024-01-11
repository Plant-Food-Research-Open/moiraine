test_that("create_omics_set works - basics", {
  x <- test_get_data_list()[["A"]]

  ## Error if wrong omics type given
  expect_error(
    create_omics_set(x, omics_type = "TEST"),
    "`omics_type` must be one of \"phenomics\", \"genomics\", \"transcriptomics\", or \"metabolomics\", not \"TEST\"."
  )

  ## Testing output type
  expect_s4_class(create_omics_set(x, omics_type = "genomics"), "SnpSet")
  expect_s4_class(create_omics_set(x, omics_type = "transcriptomics"), "ExpressionSet")
  expect_s4_class(create_omics_set(x, omics_type = "metabolomics"), "MetabolomeSet")
  expect_s4_class(create_omics_set(x, omics_type = "phenomics"), "PhenotypeSet")

  ## Testing empty metadata
  for (i in c("genomics", "transcriptomics", "metabolomics", "phenomics")) {
    t1 <- create_omics_set(x, omics_type = i)

    t1_fdata <- t1@featureData@data
    expect_equal(colnames(t1_fdata), "feature_id")
    expect_equal(rownames(t1_fdata), paste0("featureA_", 1:30))

    expect_equal(dim(t1@phenoData@data), c(15, 0))
  }
})

test_that("create_omics_set works - with fmetadata", {
  data_list <- test_get_data_list()
  fmeta_list <- test_get_fmeta_list()

  ## Testing when features metadata exactly matches features in dataset
  expect_warning(
    create_omics_set(data_list[["A"]], omics_type = "genomics", features_metadata = fmeta_list[["A"]]),
    regexp = NA
  )

  t1 <- create_omics_set(data_list[["A"]], omics_type = "genomics", features_metadata = fmeta_list[["A"]])
  t1_fmeta <- t1@featureData

  expect_s4_class(t1_fmeta, "AnnotatedDataFrame")
  expect_equal(colnames(t1_fmeta), c("feature_id", "chromosome", "position"))
  expect_equal(rownames(t1_fmeta), paste0("featureA_", 1:30))

  ## Testing when features metadata is missing some features
  expect_warning(
    create_omics_set(data_list[["C"]], omics_type = "metabolomics", features_metadata = fmeta_list[["C"]]),
    "10 features are not present in feature metadata."
  )

  t1 <- suppressWarnings(
    create_omics_set(data_list[["C"]], omics_type = "metabolomics", features_metadata = fmeta_list[["C"]])
    )
  t1_fmeta <- t1@featureData

  expect_equal(dim(t1_fmeta@data), c(40, 3))
  expect_equal(rownames(t1_fmeta), paste0("featureC_", 1:40))
  expect_true(all(is.na(t1_fmeta@data[31:40, -1])))
  expect_equal(t1_fmeta@data$feature_id, rownames(t1_fmeta@data))

  ## Testing when features metadata has too many features
  expect_warning(
    create_omics_set(data_list[["D"]], omics_type = "phenomics", features_metadata = fmeta_list[["D"]]),
    "5 features in feature metadata not in dataset, will be removed from metadata."
  )

  t1 <- suppressWarnings(
    create_omics_set(data_list[["D"]], omics_type = "phenomics", features_metadata = fmeta_list[["D"]])
  )
  t1_fmeta <- t1@featureData

  expect_equal(dim(t1_fmeta@data), c(40, 2))
  expect_equal(rownames(t1_fmeta), paste0("featureD_", 1:40))
})

test_that("create_omics_set works - with smetadata", {
  data_list <- test_get_data_list()
  smeta_list <- test_get_smeta_list()

  ## Testing when samples metadata exactly matches samples in dataset
  expect_warning(
    create_omics_set(data_list[["A"]], omics_type = "genomics", samples_metadata = smeta_list[["A"]]),
    regexp = NA
  )

  t1 <- create_omics_set(data_list[["A"]], omics_type = "genomics", samples_metadata = smeta_list[["A"]])
  t1_smeta <- t1@phenoData

  expect_s4_class(t1_smeta, "AnnotatedDataFrame")
  expect_equal(colnames(t1_smeta), c("id", "name", "pheno_group", "time"))
  expect_equal(rownames(t1_smeta), paste0("sample_", 1:15))

  ## Testing when samples metadata is missing some samples
  expect_warning(
    create_omics_set(data_list[["D"]], omics_type = "phenomics", samples_metadata = smeta_list[["D"]]),
    "5 samples are not present in samples metadata."
  )

  t1 <- suppressWarnings(create_omics_set(data_list[["D"]], omics_type = "phenomics", samples_metadata = smeta_list[["D"]]))
  t1_smeta <- t1@phenoData

  expect_equal(dim(t1_smeta@data), c(15, 4))
  expect_equal(rownames(t1_smeta), paste0("sample_", 16:30))
  expect_true(all(is.na(t1_smeta@data[25:30, -1])))
  expect_equal(t1_smeta@data$id, rownames(t1_smeta@data))

  ## Testing when sample metadata has too many samples
  expect_warning(
    create_omics_set(data_list[["C"]], omics_type = "metabolomics", samples_metadata = smeta_list[["C"]]),
    "5 samples in samples metadata not in dataset, will be removed from metadata."
  )

  t1 <- suppressWarnings(create_omics_set(data_list[["C"]], omics_type = "metabolomics", samples_metadata = smeta_list[["C"]]))
  t1_smeta <- t1@phenoData

  expect_equal(dim(t1_smeta@data), c(15, 4))
  expect_equal(rownames(t1_smeta), paste0("sample_", 11:25))
})


test_that("create_omics_set_factory works - no metadata", {
  tar_res <- create_omics_set_factory(
    datasets = c(data_omicsA, data_omicsB, data_omicsC, data_omicsD),
    omics_types = c("genomics", "transcriptomics", "metabolomics", "phenomics")
  )

  ## Testing that we get a list of target stems
  expect_type(tar_res, "list")
  expect_equal(sapply(tar_res, length), c(set = 4))
  expect_s3_class(tar_res$set[[1]], "tar_stem")

  ## Testing targets name
  expect_equal(
    sapply(tar_res$set, function(x) {
      x$settings$name
    }),
    paste0("set_omics", LETTERS[1:4])
  )

  ## Testing targets command
  expect_equal(
    sapply(tar_res$set, function(x) {
      x$command$expr
    }),
    sapply(
      paste0(
        "create_omics_set(dataset = data_omics", LETTERS[1:4],
        ", omics_type = ", c("\"genomics\"", "\"transcriptomics\"", "\"metabolomics\"", "\"phenomics\""),
        ", features_metadata = NULL, samples_metadata = NULL)"
      ),
      str2expression,
      USE.NAMES = FALSE
    )
  )

  ## Testing targets format
  expect_equal(
    sapply(tar_res$set, function(x) {
      x$settings$format
    }),
    rep("rds", 4)
  )

  ## Testing error when datasets and omics_types have different length
  expect_error(
    create_omics_set_factory(datasets = c(data_omicsA, data_omicsB), omics_types = c("genomics")),
    "'datasets' and 'omics_types' vectors must have the same length."
  )

  ## Testing change in target suffixes
  tar_res <- create_omics_set_factory(
    datasets = c(data_omicsA, data_omicsB, data_omicsC),
    omics_types = c("genomics", "transcriptomics", "metabolomics"),
    target_name_suffixes = c("D1", "D2", "D3")
  )

  expect_equal(
    sapply(tar_res$set, function(x) {
      x$settings$name
    }),
    paste0("set_D", 1:3)
  )

  ## Testing error when datasets and target_namex_suffixes have different length
  expect_error(
    create_omics_set_factory(datasets = c(data_omicsA), omics_types = c("genomics"), target_name_suffixes = c("t1", "t2")),
    "'datasets' and 'target_name_suffixes' vectors must have the same length."
  )
})

test_that("create_omics_set_factory works - with metadata", {
  ## Testing feature metadata
  tar_res <- create_omics_set_factory(
    datasets = c(data_omicsA, data_omicsB, data_omicsC),
    omics_types = c("genomics", "transcriptomics", "metabolomics"),
    features_metadatas = c(fmeta_A, NULL, fmetaC)
  )

  expect_equal(
    sapply(tar_res$set, function(x) {
      x$command$expr
    }),
    sapply(
      paste0(
        "create_omics_set(dataset = data_omics", LETTERS[1:3],
        ", omics_type = ", c("\"genomics\"", "\"transcriptomics\"", "\"metabolomics\""),
        ", features_metadata = ", c("fmeta_A", "NULL", "fmetaC"), ", samples_metadata = NULL)"
      ),
      str2expression,
      USE.NAMES = FALSE
    )
  )

  expect_error(
    create_omics_set_factory(datasets = c(data_omicsA), omics_types = c("genomics"), features_metadatas = c(fmetA, fmetB)),
    "'datasets' and 'features_metadatas' vectors must have the same length."
  )

  ## Testing sample metadata
  tar_res <- create_omics_set_factory(
    datasets = c(data_omicsA, data_omicsB, data_omicsC),
    omics_types = c("genomics", "transcriptomics", "metabolomics"),
    samples_metadatas = c(smeta_A, NULL, smetaC)
  )

  expect_equal(
    sapply(tar_res$set, function(x) {
      x$command$expr
    }),
    sapply(
      paste0(
        "create_omics_set(dataset = data_omics", LETTERS[1:3],
        ", omics_type = ", c("\"genomics\"", "\"transcriptomics\"", "\"metabolomics\""),
        ", features_metadata = NULL, samples_metadata = ", c("smeta_A", "NULL", "smetaC"), ")"
      ),
      str2expression,
      USE.NAMES = FALSE
    )
  )

  expect_error(
    create_omics_set_factory(datasets = c(data_omicsA), omics_types = c("genomics"), samples_metadatas = c(smetA, smetB)),
    "'datasets' and 'samples_metadatas' vectors must have the same length."
  )

  ## Testing feature and sample metadata
  tar_res <- create_omics_set_factory(
    datasets = c(data_omicsA, data_omicsB, data_omicsC),
    omics_types = c("genomics", "transcriptomics", "metabolomics"),
    features_metadatas = c(NULL, fmeta_B, fmeta_C),
    samples_metadatas = c(smeta_A, smeta_B, NULL)
  )

  expect_equal(
    sapply(tar_res$set, function(x) {
      x$command$expr
    }),
    sapply(
      paste0(
        "create_omics_set(dataset = data_omics", LETTERS[1:3],
        ", omics_type = ", c("\"genomics\"", "\"transcriptomics\"", "\"metabolomics\""),
        ", features_metadata = ", c("NULL", "fmeta_B", "fmeta_C"),
        ", samples_metadata = ", c("smeta_A", "smeta_B", "NULL"), ")"
      ),
      str2expression,
      USE.NAMES = FALSE
    )
  )
})


test_that("create_multiomics_set works", {
  omics_sets <- test_get_omics_list()

  ## Testing correct input
  expect_error(create_multiomics_set(list()), "sets_list list is empty.")
  expect_error(
    create_multiomics_set(list("TEST")),
    "Elements in sets_list must be SnpSet, ExpressionSet, MetabolomeSet, PhenotypeSet objects."
  )

  expect_warning(create_multiomics_set(omics_sets[1]), NA)
  expect_warning(create_multiomics_set(omics_sets), NA)

  ## Testing one set of each data type - no input names
  t1 <- create_multiomics_set(omics_sets)

  expect_s4_class(t1, "MultiDataSet")
  expect_equal(names(t1), c("snps", "rnaseq", "metabolome", "phenotypes"))
  expect_s4_class(t1[["snps"]], "SnpSet")
  expect_s4_class(t1[["rnaseq"]], "ExpressionSet")
  expect_s4_class(t1[["metabolome"]], "MetabolomeSet")
  expect_s4_class(t1[["phenotypes"]], "PhenotypeSet")
  expect_equal(n_features(t1), c("snps" = 30, "rnaseq" = 35, "metabolome" = 40, "phenotypes" = 40))
  expect_equal(n_samples(t1), c("snps" = 15, "rnaseq" = 15, "metabolome" = 15, "phenotypes" = 15))

  ## Testing two sets of the same data type - no input name
  t2 <- create_multiomics_set(omics_sets[c(1:3, 1)], show_warnings = FALSE)
  expect_equal(names(t2), c("snps+1", "rnaseq", "metabolome", "snps+2"))

  ## Testing one set of each data type - input names
  expect_error(
    create_multiomics_set(omics_sets, LETTERS[1:3], show_warnings = FALSE),
    "dataset_names vector must have same length as sets_list list."
  )
  t3 <- create_multiomics_set(omics_sets, LETTERS[1:4], show_warnings = FALSE)
  expect_equal(names(t3), c("snps+A", "rnaseq+B", "metabolome+C", "phenotypes+D"))

  t4 <- create_multiomics_set(omics_sets, LETTERS[c(1, 1:3)], show_warnings = FALSE)
  expect_equal(names(t4), c("snps+A", "rnaseq+A", "metabolome+B", "phenotypes+C"))

  t5 <- create_multiomics_set(omics_sets, c("", LETTERS[1:3]), show_warnings = FALSE)
  expect_equal(names(t5), c("snps", "rnaseq+A", "metabolome+B", "phenotypes+C"))

  ## Testing two sets of the same data type - input names
  t6 <- create_multiomics_set(omics_sets[c(1:3, 1)], LETTERS[1:4], show_warnings = FALSE)
  expect_equal(names(t6), c("snps+A", "rnaseq+B", "metabolome+C", "snps+D"))

  expect_error(
    create_multiomics_set(omics_sets[c(1:3, 1)], LETTERS[c(1:3, 1)], show_warnings = FALSE),
    "Dataset names for objects of a same type must be unique."
  )

  ## Testing with contradictions in the info
  smeta_list <- test_get_smeta_list()
  smeta_list[["C"]]["sample_11", "pheno_group"] <- "group1"
  omics_sets <- test_get_omics_list(smeta_list = smeta_list)

  expect_error(
    create_multiomics_set(omics_sets),
    "Conflicting information in samples metadata for samples sample_11."
  )
})
