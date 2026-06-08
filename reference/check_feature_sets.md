# Checks features assignment to sets

Checks proportion of features from a multi-omics dataset that are
assigned to feature sets.

## Usage

``` r
check_feature_sets(feature_sets, mo_data, datasets = names(mo_data))
```

## Arguments

- feature_sets:

  Named list, where each element corresponds to a feature set, and
  contains a vector of features ID of all features belonging to that
  set.

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- datasets:

  Character vector, the names of the datasets for which features
  assignment should be checked. By default, all datasets in `mo_data`
  are considered.

## Value

A tibble giving for each dataset the number and fraction of features
that are assigned to at least one feature set. The `message` column is
meant to facilitate reporting.
