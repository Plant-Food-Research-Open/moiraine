# Get feature IDs from MultiDataSet

Extract the list of feature IDs from each dataset in a MultiDataSet
object.

## Usage

``` r
get_features(mo_data)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

## Value

A named list, with one element per dataset, and each element is a
character vector of feature IDs.
