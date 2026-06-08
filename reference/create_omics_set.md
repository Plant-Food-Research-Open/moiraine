# Create a Biobase set object to store omics data

Creates a Biobase object to store an omics dataset and associated
samples and features metadata.

## Usage

``` r
create_omics_set(
  dataset,
  omics_type = c("phenomics", "genomics", "transcriptomics", "metabolomics"),
  features_metadata = NULL,
  samples_metadata = NULL
)
```

## Arguments

- dataset:

  Matrix, the omics dataset in matrix form with features as rows and
  samples as columns.

- omics_type:

  Character, which type of omics data is being stored? Possible values
  are `'genomics'`, `'transcriptomics'`, `'metabolomics'` and
  `'phenomics'`. Use `'phenomics'` for any other omics.

- features_metadata:

  Data.frame, a feature annotation data-frame with features as rows and
  information about the features as columns. The number of rows and row
  names must match those of the `dataset` matrix.

- samples_metadata:

  Data.frame, a samples information data-frame with samples as rows and
  information about the samples as columns. The number of rows and row
  names must match the number of columns and column names of the
  `dataset` matrix.

## Value

An object derived from Biobase::eSet:

- if `omics_type = 'genomics'`: a Biobase::SnpSet object;

- if `omics_type = 'transcriptomics'`: a Biobase::ExpressionSet object.

- if `omics_type = 'metabolomics'`: a
  [MetabolomeSet](https://plant-food-research-open.github.io/moiraine/reference/MetabolomeSet-class.md)
  object.

- if `omics_type = 'phenomics'` a
  [PhenotypeSet](https://plant-food-research-open.github.io/moiraine/reference/PhenotypeSet-class.md)
  object.

## Examples

``` r
if (FALSE) { # \dontrun{
data_geno <- import_dataset_csv(
"genotype_dataset.csv",
  col_id = "Marker",
  features_as_rows = TRUE
)
geno_info_features <- import_fmetadata_csv(
"genotype_features_info.csv",
  col_id = "Marker"
)
samples_information <- import_smetadata_csv(
"samples_information.csv",
  col_id = "Sample"
)
create_omics_set(
  dataset = data_geno,
  omics_type = "genomics",
  features_metadata = geno_info_features,
  samples_metadata = samples_information
)
} # }
```
