# Extract top features

Extracts the features with highest contribution to the latent dimensions
constructed by an integration method. Can retain a specific number of
top contributing features for each dataset and latent dimension, or all
features above a minimum importance score.

## Usage

``` r
get_top_features(
  method_output,
  n_features = 10,
  min_importance = NULL,
  latent_dimensions = NULL,
  datasets = NULL,
  mo_data = NULL
)
```

## Arguments

- method_output:

  Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- n_features:

  Integer, the number of features to extract for each latent dimension
  and dataset. Ignored if `min_importance` is set. Default value is
  `10`. Will include all ties.

- min_importance:

  Numeric value between 0 and 1, minimum importance score used to select
  features. Default value is `NULL`, i.e. the top `n_features` features
  are selected instead.

- latent_dimensions:

  Character vector of latent dimensions name. Default value is `NULL`
  (top contributing features will be returned for all latent
  dimensions).

- datasets:

  Character vector of datasets name. Default value is `NULL` (top
  contributing features will be returned for all datasets).

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

## Value

A tibble containing one row per feature and latent dimension, giving the
weight and importance score of the feature for the corresponding latent
dimension. If `mo_data` is supplied, information about the features from
the features metadata will be added to the resulting table.
