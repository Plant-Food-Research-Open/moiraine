# Makes list of feature sets from features metadata

Creates a list of feature sets from features metadata.

## Usage

``` r
make_feature_sets_from_fm(mo_data, col_names, combine_omics_sets = FALSE)
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- col_names:

  Named list of character, with one element per dataset for which
  feature sets should be generated. The name of each element should
  correspond to the name of a dataset, and the value should be a column
  name from the feature metadata of the corresponding dataset to use as
  set ID.

- combine_omics_sets:

  Logical, can sets contain features from different omics datasets? If
  `FALSE` (the default), feature sets are created for each omics
  separately. If two sets from different omics have the same ID, they
  will be made unique.

## Value

a named list, where each element corresponds to a set, and contains a
vector of features ID of all features belonging to that set.
