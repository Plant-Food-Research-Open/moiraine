# Adding data-frame to features metadata

Adds information from a data-frame to the features metadata of a
MultiDataSet object.

## Usage

``` r
add_features_metadata(mo_data, df)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- df:

  A tibble or data-frame of features information, with at least columns
  `feature_id` (giving the feature IDs) and `dataset` (giving the name
  of the dataset to which the features belong), and one row per feature
  ID. Can contain info about features from different datasets.

## Value

A MultiDataSet object, with info from `df` adding to its corresponding
features metadata.
