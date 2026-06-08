# Import a dataset from a csv file

Reads a csv file and returns a matrix in which the rows corresponds to
features (e.g. markers, genes, phenotypes...) and the columns correspond
to samples/observations.

## Usage

``` r
import_dataset_csv(file, col_id, features_as_rows = TRUE, ...)
```

## Arguments

- file:

  Character, path to the dataset csv file.

- col_id:

  Character, the name of the column in the file that contains the ID of
  the rows (i.e. feature IDs if `features_as_rows` is `TRUE`, or sample
  IDs if `features_as_rows` is `FALSE`).

- features_as_rows:

  Logical, do the rows in the file correspond to features? Default value
  is `TRUE`, i.e. the file contains features as rows and samples as
  columns.

- ...:

  Further arguments passed to
  [`readr::read_csv()`](https://readr.tidyverse.org/reference/read_delim.html).

## Value

A matrix with the samples as columns and the features as rows. Feature
IDs are used as row names and sample IDs as column names.

## Examples

``` r
if (FALSE) { # \dontrun{
data_geno <- import_dataset_csv(
  "genotype_dataset.csv",
  col_id = "Marker",
  features_as_rows = TRUE
)
data_pheno <- import_dataset_csv(
  "phenotype_dataset.csv",
  col_id = "Sample",
  features_as_rows = FALSE
)
} # }
```
