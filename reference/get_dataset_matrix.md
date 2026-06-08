# Get multi-omics dataset as matrix

Extracts an omics dataset as a matrix of measurements from a
MultiDataSet object.

## Usage

``` r
get_dataset_matrix(mo_data, dataset_name, keep_dataset_name = FALSE)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- dataset_name:

  Character, name of the omics dataset to extract.

- keep_dataset_name:

  Logical, should the dataset name be stored in the `'dataset_name'`
  attribute of the resulting matrix? Default value is `FALSE`.

## Value

A matrix of measurements with features as rows and samples as columns.
The name of the dataset is stored in the `'dataset_name'` attribute if
`keep_dataset_name` is `TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
## mo_data is a MultiDataSet object with a dataset called "rnaseq"
mat <- get_dataset_matrix(mo_data, "rnaseq", keep_dataset_name = TRUE)
## with keep_dataset_name = TRUE, can recover dataset name as follows:
attr(mat, "dataset_name")
} # }
```
