# Target factory for csv features metadata import

Creates a list of targets that track each file and import the features
metadata from each csv file.

## Usage

``` r
import_fmetadata_csv_factory(files, col_ids, target_name_suffixes)
```

## Arguments

- files:

  Character vector, vector of paths to the features metadata csv files.

- col_ids:

  Character vector, the name of the column in each file that contains
  the features ID.

- target_name_suffixes:

  Character vector, a suffix to add to the name of the targets created
  by this target factory for each dataset.

## Value

A list of target objects. For example, with two files to import and
`target_name_suffixes = c("geno", "transcripto")`, the factory returns
the following targets:

- `fmetadata_file_geno` and `fmetadata_file_transcripto`: targets
  tracking the genomics and transcriptomics features metadata files,
  respectively.

- `fmetadata_geno` and `fmetadata_transcripto`: targets that import the
  genomics and transcriptomics features metadata dataset.

## Examples

``` r
if (FALSE) { # \dontrun{
## in the _targets.R
library(moiraine)

list(
  import_fmetadata_csv_factory(
    c(
      "data/genotype_fmetadata.csv",
      "data/rnaseq_fmetadata.csv"
    ),
    col_ids = c("Marker", "Info"),
    target_name_suffixes = c("geno", "transcripto")
  )
)
} # }
```
