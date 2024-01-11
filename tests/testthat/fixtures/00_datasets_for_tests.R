here::i_am("tests/testthat/fixtures/00_datasets_for_tests.R")

library(here)
library(tidyverse)

library(mixOmics)
library(MOFA2)
library(OmicsPLS)

#######################################
## Create omics datasets for testing ##
#######################################
set.seed(1)

n_omics <- 4 ## number of omics datasets to generate

omics_names <- LETTERS[1:4] ## omics ID
samples_omics <- list(1:15, 6:20, 11:25, 16:30) ## sample IDs per omics
features_omics <- list(1:30, 1:35, 1:40, 1:40) ## feature IDs per omics

omics_data <- lapply(1:n_omics, function(i) {
  nr <- length(features_omics[[i]]) ## number of features
  nc <- length(samples_omics[[i]]) ## number of samples
  matrix(runif(nr * nc, 0, 100),
    nrow = nr,
    ncol = nc,
    dimnames = list(
      paste0("feature", omics_names[i], "_", features_omics[[i]]),
      paste0("sample_", samples_omics[[i]])
    )
  )
})

## Omics 1 contains 50 missing values
omics_data[[1]][sample(seq_len(prod(dim(omics_data[[1]]))), 50, replace = FALSE)] <- NA
sum(apply(omics_data[[1]], 1, function(x) {
  any(is.na(x))
}))
sum(apply(omics_data[[1]], 2, function(x) {
  any(is.na(x))
}))

## Omics 2 has features as columns rather than rows
omics_data[[2]] <- t(omics_data[[2]])

## Saving datasets
omics_rownames <- c("Feature", "Sample", "Feature", "Feature")
for (i in 1:n_omics) {
  omics_data[[i]] |>
    as_tibble(rownames = omics_rownames[i]) |>
    write_csv(
      here(paste0("tests/testthat/fixtures/data_omics", omics_names[i], ".csv"))
    )
}

#########################################
## Create feature metadata for testing ##
#########################################

## For omics 3, the feature metadata information will have less features than in
## the dataset. For omics 4, the feature metadata information will have more
## features than in the dataset
fmeta_omics <- features_omics
fmeta_omics[[3]] <- fmeta_omics[[3]][1:30]
fmeta_omics[[4]] <- c(fmeta_omics[[4]], 41:45)

fmeta_cols_list <- list(
  list(
    chromosome = function(n) {
      paste0("ch", sample(1:3, n, replace = TRUE))
    },
    position = function(n) {
      sample(1:100, n, replace = TRUE)
    }
  ),
  list(
    chromosome = function(n) {
      paste0("ch", sample(1:3, n, replace = TRUE))
    },
    start = function(n) {
      sample(1:100, n, replace = TRUE)
    },
    end = function(n) {
      sample(1:100, n, replace = TRUE)
    },
    name = function(n) {
      sample(fruit, n, replace = TRUE)
    }
  ),
  list(
    name = function(n) {
      sample(fruit, n, replace = TRUE)
    },
    retention_time = function(n) {
      runif(n, 1, 100)
    }
  ),
  list(unit = function(n) {
    sample(LETTERS, n, replace = TRUE)
  })
)

omics_fmetadata <- lapply(1:n_omics, function(i) {
  res <- purrr::map(
    fmeta_cols_list[[i]],
    ~ .x(length(fmeta_omics[[i]]))
  ) |>
    as.data.frame()
  rownames(res) <- paste0("feature", omics_names[i], "_", fmeta_omics[[i]])

  return(res)
})

omics_fmetadata[[2]][, "end"] <- omics_fmetadata[[2]][, "start"] + 20

## For omics 3, add a non-ASCII character in the table
omics_fmetadata[[3]]$name[1] <- "one \u03b2-\u03b2"

## Saving datasets
fmeta_rownames <- rep("Feature", n_omics)
for (i in 1:n_omics) {
  omics_fmetadata[[i]] |>
    as_tibble(rownames = fmeta_rownames[i]) |>
    write_csv(
      here(
        paste0("tests/testthat/fixtures/fmeta_omics", omics_names[i], ".csv")
      )
    )
}


################################################################
## Create small genome annotation files (GFF/GTF) for testing ##
################################################################
## We'll subset 100 lines from each file
annot_source_file <- paste0(
  "/workspace/hrpoab/Potato_data/Reference/",
  "PGSC_DM_V403_fixed_representative_genes"
)

annot_source_ext <- c(".gff", ".gtf")

annot_source_ext |>
  set_names() |>
  purrr::map(\(x) paste0(annot_source_file, x)) |>
  ## Read template file
  purrr::map(function(x) {
    template_conn <- file(x, open = "r")
    template_lines <- readLines(template_conn)
    close(template_conn)
    return(template_lines)
  }) |>
  ## Write the first 104 lines into test fixture (for GFF files, first 4 lines
  ## are header)
  purrr::imap(function(x, y) {
    fixture_conn <- file(
      here(paste0("tests/testthat/fixtures/genome_annotation", y)),
      open = "w"
    )
    writeLines(x[1:104], fixture_conn)
    close(fixture_conn)
  })


#########################################
## Create sample metadata for testing ##
#########################################

## For omics 3, the sample metadata information will have less samples than in
## the dataset For omics 4, the sample metadata information will have more
## samples than in the dataset
smeta_omics <- samples_omics
smeta_omics[[3]] <- c(smeta_omics[[3]], 31:35)
smeta_omics[[4]] <- smeta_omics[[4]][1:10]

id_samples <- unique(unlist(smeta_omics))

smetadata_df <- data.frame(
  name = paste0("Sample ", id_samples),
  pheno_group = sample(
    paste0("group", 1:2),
    length(id_samples),
    replace = TRUE
  ),
  time = seq_along(id_samples)
)
rownames(smetadata_df) <- paste0("sample_", id_samples)

omics_smetadata <- lapply(1:n_omics, function(i) {
  smetadata_df[paste0("sample_", smeta_omics[[i]]), ]
})

## Saving datasets
smeta_rownames <- rep("Sample", 4)
for (i in 1:n_omics) {
  omics_smetadata[[i]] |>
    as_tibble(rownames = smeta_rownames[i]) |>
    write_csv(
      here(
        paste0(
          "tests/testthat/fixtures/smeta_omics",
          omics_names[i],
          ".csv"
        )
      )
    )
}


############################################
## Create MultiDataSet object for testing ##
############################################

## To use the package functions
library(devtools)
load_all()

data_list <- test_get_data_list()
fmeta_list <- test_get_fmeta_list()
smeta_list <- test_get_smeta_list()

omics_sets <- purrr::map2(
  names(data_list),
  c("genomics", "transcriptomics", "metabolomics", "phenomics"),
  ~ suppressWarnings(
    create_omics_set(
      data_list[[.x]],
      omics_type = .y,
      features_metadata = fmeta_list[[.x]],
      samples_metadata = smeta_list[[.x]]
    )
  )
)

multiomics_set <- create_multiomics_set(
  omics_sets,
  datasets_names = c("A", "", "", "")
)
saveRDS(
  multiomics_set,
  file = here("tests/testthat/fixtures/multiomics_set.rds")
)

# multiomics_set <- readRDS(here("tests/testthat/fixtures/multiomics_set.rds"))

##############################################
## Run PCA on multiomics set for testing ##
##############################################
pca_res <- run_pca(multiomics_set, "rnaseq")
saveRDS(pca_res, file = here("tests/testthat/fixtures/pca_res.rds"))

##############################################
## Run sPLS-DA on multiomics set for testing ##
##############################################
splsda_res <- run_splsda(
  multiomics_set,
  to_keep_n = 10,
  dataset_name = "rnaseq",
  group = "pheno_group",
  ncomp = 2,
  multilevel = NULL
)
saveRDS(splsda_res, file = here("tests/testthat/fixtures/splsda_res.rds"))


##############################################
## Run sPLS on multiomics set for testing ##
##############################################

spls_input <- get_input_spls(
  multiomics_set,
  "regression",
  c("rnaseq", "metabolome")
)

spls_run <- spls_run(
  spls_input,
  ncomp = 2,
  keepX = c(7, 7),
  keepY = c(5, 5)
)
saveRDS(spls_run, file = here("tests/testthat/fixtures/spls_res.rds"))

##############################################
## Run DIABLO on multiomics set for testing ##
##############################################

diablo_input <- get_input_mixomics_supervised(
  multiomics_set[, 1:3],
  "pheno_group"
)

diablo_run <- block.splsda(
  diablo_input[setdiff(names(diablo_input), "Y")],
  diablo_input$Y
)
saveRDS(diablo_run, file = here("tests/testthat/fixtures/diablo_res.rds"))

############################################
## Run MOFA on multiomics set for testing ##
############################################

mofa_input <- get_input_mofa(multiomics_set[, 1:3])
mofa_res <- MOFA2::run_mofa(mofa_input)
saveRDS(mofa_res, file = here("tests/testthat/fixtures/mofa_res.rds"))


##############################################
## Run sO2PLS on multiomics set for testing ##
##############################################

omicspls_input <- get_input_omicspls(multiomics_set, c("rnaseq", "metabolome"))

so2pls_res <- so2pls_o2m(
  omicspls_input,
  n = 2,
  nx = 1,
  ny = 1,
  sparse = TRUE,
  keepx = c(3, 3),
  keepy = c(5, 5)
)

saveRDS(so2pls_res, file = here("tests/testthat/fixtures/so2pls_res.rds"))
