# Plots top features importance

Plots the top features importance per dataset and latent dimension.

## Usage

``` r
plot_top_features(
  method_output,
  latent_dimensions = NULL,
  group_latent_dims = TRUE,
  datasets = NULL,
  n_features = 20,
  mo_data = NULL,
  label_cols = NULL,
  truncate = NULL,
  nrow = 1
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

- group_latent_dims:

  Logical, for integrations methods that construct datasets- specific
  versions of each latent dimension, should these be grouped? e.g. when
  DIABLO constructs a snps- and rnaseq version of component 1, should
  the two be grouped as "Component 1"? Default value is `TRUE`.

- datasets:

  Character vector giving the datasets to display. Default value is
  `NULL`, i.e. all datasets will be shown.

- n_features:

  Integer, number of top features to display per dataset and latent
  dimension.

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

- nrow:

  Integer, number of rows over which the dataset panels should be
  plotted for each latent dimensions.

## Value

A
[`patchwork`](https://patchwork.data-imaginist.com/reference/patchwork-package.html)
of plots.
