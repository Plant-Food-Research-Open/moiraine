# Target factory for omics sets creation

Creates a list of targets that generate omics sets from targets
containing datasets, features and samples metadata.

## Usage

``` r
create_omics_set_factory(
  datasets,
  omics_types,
  features_metadatas = NULL,
  samples_metadatas = NULL,
  target_name_suffixes = NULL
)
```

## Arguments

- datasets:

  Vector of symbols, the names of the targets containing the omics
  datasets.

- omics_types:

  Character vector, which type of omics data is being stored for each
  dataset? Possible values are `'genomics'`, `'transcriptomics'`,
  `'metabolomics'` and `'phenomics'`. Use `'phenomics'` for any other
  omics. Use `'phenomics'` for any other omics.

- features_metadatas:

  Vector of symbols, the names of the targets containing the features
  metadata data-frame associated with each omics dataset. Use `NULL` if
  no feature metadata exists for a dataset.

- samples_metadatas:

  Vector of symbols, the names of the targets containing the samples
  metadata data-frame associated with each omics dataset. Use `NULL` if
  no samples metadata exists for a dataset.

- target_name_suffixes:

  Character vector, a suffix to add to the name of the targets created
  by this target factory for each dataset. If none provided, the
  suffixes will be extracted from the `datasets` argument. Default value
  is NULL.

## Value

A list of target objects, with three datasets provided, and
`target_name_suffixes = c("geno", "transcripto", "metabo")`, the
following targets will be returned: `set_geno`, `set_transcripto` and
`set_metabo`.

## Examples

``` r
if (FALSE) { # \dontrun{
## in the _targets.R
library(moiraine)
library(targets)

list(
  ## targets to import the different datasets

  ## Example where genomics dataset has no features metadata information
  ## Will generate the following targets: set_geno, set_transcripto
  create_omics_set_factory(
    datasets = c(data_geno, data_transcripto),
    omics_types = c("genomics", "transcriptomics"),
    features_metadata = c(NULL, fmeta_transcripto),
    samples_metadata = c(smeta_geno, smeta_transcripto)
  )
)
} # }
```
