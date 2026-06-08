# Target factory for csv datasets import

Creates a list of targets that track each file and import the dataset
from each csv file.

## Usage

``` r
import_dataset_csv_factory(
  files,
  col_ids,
  features_as_rowss,
  target_name_suffixes
)
```

## Arguments

- files:

  Character vector, vector of paths to the dataset csv files.

- col_ids:

  Character vector, the name of the column in each file that contains
  the ID of the rows (i.e. feature IDs if value in `features_as_rowss`
  is `TRUE` for the corresponding dataset, or sample IDs if value in
  `features_as_rowss` is `FALSE`).

- features_as_rowss:

  Logical vector, do the rows in each file correspond to features?

- target_name_suffixes:

  Character vector, a suffix to add to the name of the targets created
  by this target factory for each dataset.

## Value

A list of target objects. For example, with two files to import and
`target_name_suffixes = c("geno", "transcripto")`, the factory returns
the following targets:

- `dataset_file_geno` and `dataset_file_transcripto`: targets tracking
  the genomics dataset file and the transcriptomics dataset file,
  respectively.

- `data_geno` and `data_transcripto`: targets that import the genomics
  and transcriptomics dataset, respectively.

## Examples

``` r
if (FALSE) { # \dontrun{
## in the _targets.R
library(moiraine)

list(
  import_dataset_csv_factory(
    c(
      "data/genotype_data.csv",
      "data/rnaseq_data.csv"
    ),
    col_ids = c("Marker", "Sample"),
    features_as_rows = c(TRUE, FALSE),
    target_name_suffixes = c("geno", "transcripto")
  )
)
} # }
```
