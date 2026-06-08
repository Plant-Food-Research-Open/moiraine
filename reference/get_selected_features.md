# Extract selected features

Extracts selected features from the output of an integration method.
Only features with a non-null weight for at least one latent dimension
will be returned. If a `MultiDataSet` object is supplied, information
about the features from the features metadata will be added to the
resulting table.

## Usage

``` r
get_selected_features(
  method_output,
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
