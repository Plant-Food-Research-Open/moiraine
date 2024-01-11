test_that("import_dataset_csv works - features as rows", {
  data_omicsA_file <- test_path("fixtures", "data_omicsA.csv")

  ## Testing correct behaviour
  read_omicsA <- import_dataset_csv(
    data_omicsA_file,
    col_id = "Feature",
    show_col_types = FALSE ## to avoid messages when running tests
  )

  expect_true(is.matrix(read_omicsA))
  expect_equal(dim(read_omicsA), c(30, 15))
  expect_equal(rownames(read_omicsA), paste0("featureA_", 1:30))
  expect_equal(colnames(read_omicsA), paste0("sample_", 1:15))

  ## Testing incorrect behaviours - wrong column name
  expect_error(
    import_dataset_csv(data_omicsA_file, col_id = "Test", show_col_types = FALSE),
    "'Test' is not a column in the csv file. Please specify a valid column name for argument col_id."
  )
})

test_that("import_dataset_csv works - features as columns", {
  data_omicsB_file <- test_path("fixtures", "data_omicsB.csv")
  read_omicsB <- import_dataset_csv(
    data_omicsB_file,
    col_id = "Sample",
    features_as_rows = FALSE,
    show_col_types = FALSE
  )

  expect_equal(dim(read_omicsB), c(35, 15))
  expect_equal(rownames(read_omicsB), paste0("featureB_", 1:35))
  expect_equal(colnames(read_omicsB), paste0("sample_", 6:20))
})

test_that("import_fmetadata_csv works", {
  if (!requireNamespace("stringi", quietly = TRUE)) {
    expect_error(
      suppressMessages(import_fmetadata_csv(fmeta_omicsC_file, "Feature")),
      "Package \"stringi\" must be installed to use this function."
    )
    skip("Package \"stringi\" must be installed to use this function.")
  }

  if (!requireNamespace("textclean", quietly = TRUE)) {
    expect_error(
      suppressMessages(import_fmetadata_csv(fmeta_omicsC_file, "Feature")),
      "Package \"textclean\" must be installed to use this function."
    )
    skip("Package \"textclean\" must be installed to use this function.")
  }

  fmeta_omicsC_file <- test_path("fixtures", "fmeta_omicsC.csv")
  df_fmeta_omicsC <- readr::read_csv(fmeta_omicsC_file, show_col_types = FALSE)

  ## Testing correct behaviour
  read_fmeta_omicsC <- import_fmetadata_csv(
    fmeta_omicsC_file,
    "Feature",
    show_col_types = FALSE
  )

  expect_s3_class(read_fmeta_omicsC, "data.frame")
  expect_equal(dim(read_fmeta_omicsC), c(30, 3))
  expect_equal(rownames(read_fmeta_omicsC), paste0("featureC_", 1:30))
  expect_equal(
    colnames(read_fmeta_omicsC),
    c("feature_id", "name", "retention_time")
  )

  expect_equal(read_fmeta_omicsC$name[1], "one beta-beta")


  ## Testing incorrect behaviours - wrong column name
  expect_error(
    import_fmetadata_csv(fmeta_omicsC_file, "Test", show_col_types = FALSE),
    "'Test' is not a column in the csv file. Please specify a valid column name for argument col_id."
  )
})

test_that("import_fmetadata_gff works - gff as input", {
  gff_file <- test_path("fixtures", "genome_annotation.gff")

  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    expect_error(
      suppressMessages(import_fmetadata_gff(gff_file, "genes")),
      "Package \"GenomicFeatures\" must be installed to use this function."
    )
    skip("Package \"GenomicFeatures\" must be installed to use this function.")
  }

  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    expect_error(
      suppressMessages(import_fmetadata_gff(gff_file, "genes")),
      "Package \"rtracklayer\" must be installed to use this function."
    )
    skip("Package \"rtracklayer\" must be installed to use this function.")
  }

  if (!requireNamespace("stringi", quietly = TRUE)) {
    expect_error(
      suppressMessages(import_fmetadata_gff(gff_file, "genes")),
      "Package \"stringi\" must be installed to use this function."
    )
    skip("Package \"stringi\" must be installed to use this function.")
  }

  if (!requireNamespace("textclean", quietly = TRUE)) {
    expect_error(
      suppressMessages(import_fmetadata_gff(gff_file, "genes")),
      "Package \"textclean\" must be installed to use this function."
    )
    skip("Package \"textclean\" must be installed to use this function.")
  }

  expect_error(
    suppressMessages(import_fmetadata_gff(gff_file, "TEST")),
    "feature_type argument should be 'genes' or 'transcripts'."
  )

  ## Testing the extraction of genes annotation
  annot <- suppressMessages(import_fmetadata_gff(gff_file, "genes"))

  expect_s3_class(annot, "data.frame")
  expect_equal(dim(annot), c(7, 6))
  expect_equal(rownames(annot), annot$feature_id)
  expect_true(all(stringr::str_detect(rownames(annot), "PGSC0003DMG40")))
  expect_equal(names(annot)[1], "feature_id")
  expect_false("gene_id" %in% names(annot))

  ## Testing the extraction of additional fields
  annot_withname <- suppressMessages(
    import_fmetadata_gff(gff_file, "genes", "name")
  )

  expect_equal(dim(annot_withname), c(7, 7))
  expect_equal(names(annot_withname)[7], "name")

  expect_error(
    suppressMessages(import_fmetadata_gff(gff_file, "genes", "TEST"))
  )

  ## Testing the extraction of transcripts annotation
  annot_tr <- suppressMessages(import_fmetadata_gff(gff_file, "transcripts"))

  expect_equal(dim(annot_tr), c(7, 6))
  expect_equal(rownames(annot_tr), annot_tr$feature_id)
  expect_true(all(stringr::str_detect(rownames(annot_tr), "PGSC0003DMT40")))
  expect_equal(names(annot)[1], "feature_id")
  expect_false("tx_name" %in% names(annot))
  expect_false("tx_id" %in% names(annot))

  ## Testing the extraction of additional fields
  annot_tr_withfields <- suppressMessages(
    import_fmetadata_gff(gff_file, "transcripts", c("Parent", "Source_id"))
  )
  expect_equal(dim(annot_tr_withfields), c(7, 8))
  expect_equal(names(annot_tr_withfields)[7:8], c("Parent", "Source_id"))
})

test_that("import_fmetadata_gff works - gtf as input", {
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    skip("Package \"GenomicFeatures\" must be installed to use this function.")
  }

  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    skip("Package \"rtracklayer\" must be installed to use this function.")
  }


  if (!requireNamespace("stringi", quietly = TRUE)) {
    skip("Package \"stringi\" must be installed to use this function.")
  }

  if (!requireNamespace("textclean", quietly = TRUE)) {
    skip("Package \"textclean\" must be installed to use this function.")
  }

  gtf_file <- test_path("fixtures", "genome_annotation.gtf")

  ## Testing the extraction of genes annotation
  annot <- suppressMessages(import_fmetadata_gff(gtf_file, "genes"))

  expect_s3_class(annot, "data.frame")
  expect_equal(dim(annot), c(19, 6))
  expect_equal(rownames(annot), annot$feature_id)
  expect_true(all(stringr::str_detect(rownames(annot), "PGSC0003DMG40")))
  expect_equal(names(annot)[1], "feature_id")
  expect_false("gene_id" %in% names(annot))

  ## Testing the extraction of additional fields
  annot_withname <- suppressMessages(
    import_fmetadata_gff(gtf_file, "genes", "gene_name")
  )

  expect_equal(dim(annot_withname), c(19, 7))
  expect_equal(names(annot_withname)[7], "gene_name")

  expect_error(
    suppressMessages(import_fmetadata_gff(gtf_file, "genes", "TEST"))
  )

  ## Testing the extraction of transcripts annotation
  annot_tr <- suppressMessages(import_fmetadata_gff(gtf_file, "transcripts"))

  expect_equal(dim(annot_tr), c(19, 6))
  expect_equal(rownames(annot_tr), annot_tr$feature_id)
  expect_true(all(stringr::str_detect(rownames(annot_tr), "PGSC0003DMT40")))
  expect_equal(names(annot)[1], "feature_id")
  expect_false("tx_name" %in% names(annot))
  expect_false("tx_id" %in% names(annot))

  ## Testing the extraction of additional fields
  annot_tr_withfields <- suppressMessages(
    import_fmetadata_gff(gtf_file, "transcripts", c("gene_id", "gene_name"))
  )

  expect_equal(dim(annot_tr_withfields), c(19, 8))
  expect_equal(names(annot_tr_withfields)[7:8], c("gene_id", "gene_name"))
})

test_that("import_smetadata_csv works", {
  smeta_omicsA_file <- test_path("fixtures", "smeta_omicsA.csv")
  read_smeta_omicsA <- import_smetadata_csv(
    smeta_omicsA_file,
    "Sample",
    show_col_types = FALSE
  )

  expect_s3_class(read_smeta_omicsA, "data.frame")
  expect_equal(dim(read_smeta_omicsA), c(15, 4))
  expect_equal(rownames(read_smeta_omicsA), paste0("sample_", 1:15))
  expect_equal(
    colnames(read_smeta_omicsA),
    c("id", "name", "pheno_group", "time")
  )

  expect_error(
    import_smetadata_csv(smeta_omicsA_file, "Test", show_col_types = FALSE),
    "'Test' is not a column in the csv file. Please specify a valid column name for argument col_id."
  )
})

test_that("import_dataset_csv_factory works", {
  ## Here we're not looking for an existing file "data/genotype_data.csv"
  ## it's just to test the targets created, not reading data
  tar_res <- import_dataset_csv_factory(
    c(
      "data/genotype_data.csv",
      "data/rnaseq_data.csv"
    ),
    col_ids = c("Marker", "Sample"),
    features_as_rowss = c(TRUE, FALSE),
    target_name_suffixes = c("geno", "transcripto")
  )

  ## Testing that we get a list of target stems
  expect_type(tar_res, "list")
  expect_equal(names(tar_res), c("dataset_file", "data"))
  expect_equal(sapply(tar_res, length), c(dataset_file = 2, data = 2))
  expect_s3_class(tar_res$dataset_file[[1]], "tar_stem")

  ## Testing targets name
  expect_equal(
    sapply(tar_res$dataset_file, function(x) {
      x$settings$name
    }),
    c("dataset_file_geno", "dataset_file_transcripto")
  )
  expect_equal(
    sapply(tar_res$data, function(x) {
      x$settings$name
    }),
    c("data_geno", "data_transcripto")
  )

  ## Testing targets command
  expect_equal(
    sapply(tar_res$dataset_file, function(x) {
      x$command$expr
    }),
    c(expression("data/genotype_data.csv"), expression("data/rnaseq_data.csv"))
  )
  expect_equal(
    sapply(tar_res$data, function(x) {
      x$command$expr
    }),
    c(
      expression(import_dataset_csv(dataset_file_geno, col_id = "Marker", features_as_rows = TRUE)),
      expression(import_dataset_csv(dataset_file_transcripto, col_id = "Sample", features_as_rows = FALSE))
    )
  )

  ## Testing targets format
  expect_equal(
    sapply(tar_res$dataset_file, function(x) {
      x$settings$format
    }),
    c("file", "file")
  )
  expect_equal(
    sapply(tar_res$data, function(x) {
      x$settings$format
    }),
    c("rds", "rds")
  )
})

test_that("import_fmetadata_csv_factory works", {
  tar_res <- import_fmetadata_csv_factory(
    c(
      "data/genotype_fmetadata.csv",
      "data/rnaseq_fmetadata.csv"
    ),
    col_ids = c("Marker", "Info"),
    target_name_suffixes = c("geno", "transcripto")
  )


  ## Testing that we get a list of target stems
  expect_type(tar_res, "list")
  expect_equal(names(tar_res), c("fmetadata_file", "fmetadata"))
  expect_equal(sapply(tar_res, length), c(fmetadata_file = 2, fmetadata = 2))
  expect_s3_class(tar_res$fmetadata[[1]], "tar_stem")

  ## Testing targets name
  expect_equal(
    sapply(tar_res$fmetadata_file, function(x) {
      x$settings$name
    }),
    c("fmetadata_file_geno", "fmetadata_file_transcripto")
  )
  expect_equal(
    sapply(tar_res$fmetadata, function(x) {
      x$settings$name
    }),
    c("fmetadata_geno", "fmetadata_transcripto")
  )

  ## Testing targets command
  expect_equal(
    sapply(tar_res$fmetadata_file, function(x) {
      x$command$expr
    }),
    c(
      expression("data/genotype_fmetadata.csv"),
      expression("data/rnaseq_fmetadata.csv")
    )
  )
  expect_equal(
    sapply(tar_res$fmetadata, function(x) {
      x$command$expr
    }),
    c(
      expression(
        import_fmetadata_csv(fmetadata_file_geno, col_id = "Marker")
      ),
      expression(
        import_fmetadata_csv(fmetadata_file_transcripto, col_id = "Info")
      )
    )
  )

  ## Testing targets format
  expect_equal(
    sapply(tar_res$fmetadata_file, function(x) {
      x$settings$format
    }),
    c("file", "file")
  )
  expect_equal(
    sapply(tar_res$fmetadata, function(x) {
      x$settings$format
    }),
    c("rds", "rds")
  )
})

test_that("import_fmetadata_gff_factory works", {
  tar_res <- import_fmetadata_gff_factory(
    c(
      "data/annotation.gff",
      "data/annotationv2.gtf"
    ),
    feature_types = c("genes", "transcripts"),
    add_fieldss = list(
      c("gene_name", "gene_custom_ID"),
      c("transcript_name")
    ),
    target_name_suffixes = c("geno", "transcripto")
  )

  ## Testing that we get a list of target stems
  expect_type(tar_res, "list")
  expect_equal(names(tar_res), c("fmetadata_file", "fmetadata"))
  expect_equal(sapply(tar_res, length), c(fmetadata_file = 2, fmetadata = 2))
  expect_s3_class(tar_res$fmetadata[[1]], "tar_stem")

  ## Testing targets name
  expect_equal(
    sapply(tar_res$fmetadata_file, function(x) {
      x$settings$name
    }),
    c("fmetadata_file_geno", "fmetadata_file_transcripto")
  )
  expect_equal(
    sapply(tar_res$fmetadata, function(x) {
      x$settings$name
    }),
    c("fmetadata_geno", "fmetadata_transcripto")
  )

  ## Testing targets command
  expect_equal(
    sapply(tar_res$fmetadata_file, function(x) {
      x$command$expr
    }),
    c(expression("data/annotation.gff"), expression("data/annotationv2.gtf"))
  )
  expect_equal(
    sapply(tar_res$fmetadata, function(x) {
      x$command$expr
    }),
    c(
      expression(
        import_fmetadata_gff(
          fmetadata_file_geno,
          feature_type = "genes",
          add_fields = c("gene_name", "gene_custom_ID")
        )
      ),
      expression(
        import_fmetadata_gff(
          fmetadata_file_transcripto,
          feature_type = "transcripts",
          add_fields = "transcript_name")
      )
    )
  )

  ## Testing targets format
  expect_equal(
    sapply(tar_res$fmetadata_file, function(x) {
      x$settings$format
    }),
    c("file", "file")
  )
  expect_equal(
    sapply(tar_res$fmetadata, function(x) {
      x$settings$format
    }),
    c("rds", "rds")
  )
})

test_that("import_smetadata_csv_factory works", {
  tar_res <- import_smetadata_csv_factory(
    c(
      "data/genotype_smetadata.csv",
      "data/rnaseq_smetadata.csv"
    ),
    col_ids = c("Sample", "SampleIDs"),
    target_name_suffixes = c("geno", "transcripto")
  )

  ## Testing that we get a list of target stems
  expect_type(tar_res, "list")
  expect_equal(names(tar_res), c("smetadata_file", "smetadata"))
  expect_equal(sapply(tar_res, length), c(smetadata_file = 2, smetadata = 2))
  expect_s3_class(tar_res$smetadata[[1]], "tar_stem")

  ## Testing targets name
  expect_equal(
    sapply(tar_res$smetadata_file, function(x) {
      x$settings$name
    }),
    c("smetadata_file_geno", "smetadata_file_transcripto")
  )
  expect_equal(
    sapply(tar_res$smetadata, function(x) {
      x$settings$name
    }),
    c("smetadata_geno", "smetadata_transcripto")
  )

  ## Testing targets command
  expect_equal(
    sapply(tar_res$smetadata_file, function(x) {
      x$command$expr
    }),
    c(
      expression("data/genotype_smetadata.csv"),
      expression("data/rnaseq_smetadata.csv")
    )
  )
  expect_equal(
    sapply(tar_res$smetadata, function(x) {
      x$command$expr
    }),
    c(
      expression(
        import_smetadata_csv(smetadata_file_geno, col_id = "Sample")
      ),
      expression(
        import_smetadata_csv(smetadata_file_transcripto, col_id = "SampleIDs")
      )
    )
  )

  ## Testing targets format
  expect_equal(
    sapply(tar_res$smetadata_file, function(x) {
      x$settings$format
    }),
    c("file", "file")
  )
  expect_equal(
    sapply(tar_res$smetadata, function(x) {
      x$settings$format
    }),
    c("rds", "rds")
  )
})
