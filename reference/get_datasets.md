# Get multi-omics measurement datasets

Returns the multi-omics datasets as a list of matrices from a
MultiDataSet object.

## Usage

``` r
get_datasets(mo_data)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

## Value

A named list of matrices, each with features as rows and samples as
columns.
