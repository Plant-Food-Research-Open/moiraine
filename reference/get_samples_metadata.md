# Get samples metadata dataframes from MultiDataSet

Extracts the samples metadata data-frame (phenoData field) from each
dataset for a MultiDataSet object.

## Usage

``` r
get_samples_metadata(mo_data)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

## Value

A named list of data-frames, one per dataset in the `mo_data` object.
