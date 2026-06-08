# Target factory for csv samples metadata import

Creates a list of targets that track each file and import the samples
metadata from each csv file.

## Usage

``` r
import_smetadata_csv_factory(files, col_ids, target_name_suffixes)
```

## Arguments

- files:

  Character vector, vector of paths to the samples metadata csv files.

- col_ids:

  Character vector, the name of the column in each file that contains
  the ID of the rows (i.e. the sample IDs).

- target_name_suffixes:

  Character vector, a suffix to add to the name of the targets created
  by this target factory for each dataset.

## Value

A list of target objects. For example, with two files to import and
`target_name_suffixes = c("geno", "transcripto")`, the factory returns
the following targets:

- `smetadata_file_geno` and `smetadata_file_transcripto`: targets
  tracking the genomics and transcriptomics samples metadata files,
  respectively.

- `smetadata_geno` and `smetadata_transcripto`: targets that import the
  genomics and transcriptomics samples metadata datasets, respectively.

## Examples

``` r
if (FALSE) { # \dontrun{
## in the _targets.R
library(moiraine)

list(
  import_smetadata_csv_factory(
    c(
      "data/genotype_smetadata.csv",
      "data/rnaseq_smetadata.csv"
    ),
    col_ids = c("Sample", "SampleIDs"),
    target_name_suffixes = c("geno", "transcripto")
  )
)
} # }
```
