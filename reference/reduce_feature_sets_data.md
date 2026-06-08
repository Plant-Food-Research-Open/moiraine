# Reduce feature sets to match multi-omics dataset

Removes from a list of feature sets all features that are not present in
the multi-omics dataset.

## Usage

``` r
reduce_feature_sets_data(feature_sets, mo_data, datasets = names(mo_data))
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

A feature sets list, i.e. named list where each element corresponds to a
feature set, containing the ID of all features that belong to this set
and that are present in the multi-omics dataset.
