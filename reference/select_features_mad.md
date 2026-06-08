# Select features based on Median Absolute Deviation from MultiDataSet

Computes the Median Absolute Deviation (MAD) for each feature in an
omics dataset from a MultiDataSet object, and select features with the
highest MAD values. This is a wrapper function around the
[`get_dataset_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/get_dataset_matrix.md)
and
[`select_features_mad_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/select_features_mad_matrix.md)
functions.

## Usage

``` r
select_features_mad(
  mo_data,
  dataset_name,
  to_keep_n = NULL,
  to_keep_prop = NULL,
  with_ties = TRUE
)
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- dataset_name:

  Character, name of the omics dataset on which to apply feature
  pre-selection.

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

A tibble with columns `feature_id`, `mad` and `selected` (logical,
indicates whether the feature is selected based on its MAD value). In
addition, the name of the dataset filtered is stored in the return
object attribute `dataset_name` (which can be accessed via
`attr(res, "dataset_name")`).
