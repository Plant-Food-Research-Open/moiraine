# Get features metadata dataframes from MultiDataSet

Extracts the features metadata dataframe (featureData field) from each
dataset for a MultiDataSet object.

## Usage

``` r
get_features_metadata(mo_data)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

## Value

A named list of data-frames, one per dataset in the `mo_data` object.
