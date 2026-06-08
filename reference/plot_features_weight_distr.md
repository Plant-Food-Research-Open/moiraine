# Plots features weight distribution

Plots the features weight or importance per dataset and latent
dimension.

## Usage

``` r
plot_features_weight_distr(
  method_output,
  latent_dimensions = NULL,
  datasets = NULL,
  features_metric = c("signed_importance", "weight", "importance"),
  top_n = 0,
  mo_data = NULL,
  label_cols = NULL,
  truncate = NULL,
  text_size = 2.5
)
```

## Arguments

- method_output:

  Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- latent_dimensions:

  Character vector giving the latent dimensions to display. Default
  value is `NULL`, i.e. all latent dimensions will be shown.

- datasets:

  Character vector giving the datasets to display. Default value is
  `NULL`, i.e. all datasets will be shown.

- features_metric:

  Character, which attribute should be plotted: can be
  `'signed_importance'` (i.e. importance value but with the weight
  sign), `'importance'` or `'weight'`. Default value is
  `'signed_importance'`.

- top_n:

  Integer, number of top features (in terms of importance) for which the
  label should be shown. Default value is `0`.

- mo_data:

  A `MultiDataSet` object. Only used if `label_cols` is not `NULL`.

- label_cols:

  Character or named list of character, giving for each dataset the name
  of the column in the corresponding features metadata to use as label.
  If one value, will be used for all datasets. If list, the names must
  correspond to the names of the datasets in `mo_data`. If a dataset is
  missing from the list or no value is provided, feature IDs will be
  used as labels. Alternatively, use `feature_id` to get the feature IDs
  as labels.

- truncate:

  Integer, width to which the labels should be truncated (to avoid
  issues with very long labels in plots). If `NULL` (default value), no
  truncation will be performed.

- text_size:

  Numeric, size of the feature labels.

## Value

A
[`patchwork`](https://patchwork.data-imaginist.com/reference/patchwork-package.html)
of plots.
