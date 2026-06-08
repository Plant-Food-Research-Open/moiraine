# Subset a MultiDataSet object by feature

Subsets a MultiDataSet object based on a list of feature IDs provided.

## Usage

``` r
subset_features(mo_data, features_id)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- features_id:

  Character vector, a vector of feature IDs (from across the datasets)
  to select. Also accepts lists (e.g. list with a vector of feature IDs
  per dataset).

## Value

A
[MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
object with only features specified.

## Examples

``` r
if (FALSE) { # \dontrun{
## works with a vector of feature IDs:
subset_features(mo_data, c("featureA", "featureB", "featureC"))

## or with a list of feature IDs (typically one per dataset, but doesn't
## have to be):
subset_features(
  mo_data,
  list(
    c("omics1_featureA", "omics1_featureB", "omics1_featureC"),
    c("omics2_featureA", "omics2_featureB"),
    c("omics3_featureA", "omics3_featureB", "omics3_featureC"),
  )
)
} # }
```
