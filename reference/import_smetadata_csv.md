# Import samples metadata from a csv file

Reads a csv file and returns a dataframe in which the rows correspond to
features (e.g. markers, genes, phenotypes...) and columns correspond to
information about the features.

## Usage

``` r
import_smetadata_csv(file, col_id, ...)
```

## Arguments

- file:

  Character, path to the dataset csv file.

- col_id:

  Character, the name of the column in the file that contains the ID of
  the rows (i.e. sample IDs).

- ...:

  Further arguments passed to
  [`readr::read_csv()`](https://readr.tidyverse.org/reference/read_delim.html).

## Value

A data-frame with the samples as rows and the samples properties as
columns. Sample IDs are used as rownames.

## Examples

``` r
if (FALSE) { # \dontrun{
samples_information <- import_smetadata_csv(
  "samples_information.csv",
  col_id = "Sample"
)
} # }
```
