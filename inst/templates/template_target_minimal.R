library(targets)
library(tarchetypes)
library(moiraine)

## add any package you need to use in the pipeline here
tar_option_set(
  packages = c(
    "moiraine",
    "MOFA2",
    "mixOmics",
    "readr",
    "tibble",
    "tidyr",
    "dplyr",
    "ggplot2",
    "patchwork"
  )
)

list(

  # Data import ----------------------------------------------------------------

  ## Data import using a target factory
  import_dataset_csv_factory(
    files = c(
      system.file("extdata/genomics_dataset.csv", package = "moiraine"),
      system.file("extdata/transcriptomics_dataset.csv", package = "moiraine"),
      system.file("extdata/metabolomics_dataset.csv", package = "moiraine")
    ),
    col_ids = c("marker", "gene_id", "sample_id"),
    features_as_rowss = c(TRUE, TRUE, FALSE),
    target_name_suffixes = c("geno", "transcripto", "metabo")
  ),

  ## Genomics features metadata file
  tar_target(
    fmetadata_file_geno,
    system.file("extdata/genomics_features_info.csv", package = "moiraine"),
    format = "file"
  ),

  ## Genomics features metadata import
  tar_target(
    fmetadata_geno,
    import_fmetadata_csv(
      fmetadata_file_geno,
      col_id = "marker",
      col_types = c("chromosome" = "c")
    )
  ),


  ## Metabolomics features metadata import
  import_fmetadata_csv_factory(
    files = c(
      system.file("extdata/metabolomics_features_info.csv", package = "moiraine")
    ),
    col_ids = c("feature_id"),
    target_name_suffixes = c("metabo")
  ),

  ## Transcriptomics features metadata import
  import_fmetadata_gff_factory(
    files = system.file("extdata/bos_taurus_gene_model.gff3", package = "moiraine"),
    feature_types = "genes",
    add_fieldss = c("Name", "description"),
    target_name_suffixes = "transcripto"
  ),

  ## Samples metadata import
  import_smetadata_csv_factory(
    files = system.file("extdata/samples_info.csv", package = "moiraine"),
    col_ids = "animal_id",
    target_name_suffixes = "all"
  ),

  ## Creating omics sets for each dataset
  create_omics_set_factory(
    datasets = c(data_geno, data_transcripto, data_metabo),
    omics_types = c("genomics", "transcriptomics", "metabolomics"),
    features_metadatas = c(fmetadata_geno, fmetadata_transcripto, fmetadata_metabo),
    samples_metadatas = c(smetadata_all, smetadata_all, smetadata_all)
  ),

  ## Creating the MultiDataSet object
  tar_target(
    mo_set,
    create_multiomics_set(
      list(set_geno,
           set_transcripto,
           set_metabo)
    )
  ),

  # Inspecting datasets --------------------------------------------------------

  ## Creating a density plot for each dataset
  tar_target(
    density_plots,
    plot_density_data(
      mo_set,
      combined = FALSE,
      scales = "free"
    )
  ),

  ## Plotting the relationship between features mean and standard deviation
  ## for each dataset
  tar_target(
    mean_sd_plots,
    plot_meansd_data(mo_set)
  ),

  ## Assessing missing values
  tar_target(
    n_missing_values,
    check_missing_values(mo_set)
  )
)
