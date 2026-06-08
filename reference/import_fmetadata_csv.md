# Import feature metadata from a csv file

Reads a csv file and returns a dataframe in which the rows correspond to
features (e.g. markers, genes, phenotypes...) and columns correspond to
information about the features. Non-ASCII characters are replaced with
ASCII equivalents (using the stringi and textclean packages).

## Usage

``` r
import_fmetadata_csv(file, col_id, ...)
```

## Arguments

- file:

  Character, path to the dataset csv file.

- col_id:

  Character, the name of the column in the file that contains the
  feature IDs.

- ...:

  Further arguments passed to
  [`readr::read_csv()`](https://readr.tidyverse.org/reference/read_delim.html).

## Value

A data-frame with the features as rows and features information as
columns. Feature IDs are used as row names.

## Examples

``` r
if (FALSE) { # \dontrun{
geno_info_features <- import_fmetadata_csv(
  "genotype_features_info.csv",
  col_id = "Marker"
)
} # }
```
