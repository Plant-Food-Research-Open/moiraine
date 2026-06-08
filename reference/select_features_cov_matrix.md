# Select features based on Coefficient of Variation from matrix

Computes the Coefficient of Variation (COV) for each feature in an omics
dataset from a MultiDataSet object, and select features with the highest
COV values.

## Usage

``` r
select_features_cov_matrix(
  mat,
  to_keep_n = NULL,
  to_keep_prop = NULL,
  with_ties = TRUE
)
```

## Arguments

- mat:

  Matrix of omics measurement, with features as rows and samples as
  columns.

- to_keep_n:

  Integer, the number of features to retain in the dataset. Should be
  less than the number of features in the dataset. If `NULL` or `NA`,
  `to_keep_prop` will be used instead.

- to_keep_prop:

  Numeric, the proportion of features to retain in the dataset. Will be
  ignored if `to_keep_n` is supplied. Value should be \> 0 and \< 1.

- with_ties:

  Should ties be kept together? If `TRUE`, may return more features than
  requested. Default value is `TRUE`.

## Value

A tibble with columns `feature_id`, `cov` and `selected` (logical,
indicates whether the feature is selected based on its COV value). In
addition, the name of the dataset filtered is stored in the return
object attribute `dataset_name` (which can be accessed via
`attr(res, "dataset_name")`).
