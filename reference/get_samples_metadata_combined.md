# Get combined samples metadata data-frame from MultiDataSet

Extracts the samples metadata data-frame (phenoData field) from each
dataset for a MultiDataSet object and combine them into one dataframe.

## Usage

``` r
get_samples_metadata_combined(mo_data, only_common_cols = TRUE)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- only_common_cols:

  Logical, whether to retain only common columns. If `TRUE` (default
  value), only retain the columns that are present in the samples
  metadata of all datasets. If `FALSE`, retain all columns from each
  datasets' sample metadata.

## Value

A data-frame of samples metadata.
