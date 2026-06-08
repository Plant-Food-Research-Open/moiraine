# Plots features weight as a scatterplot

Plots the features weight from a pair of latent dimensions from one or
two dimension reduction analysis as a scatterplot.

## Usage

``` r
plot_features_weight_pair(
  method_output,
  latent_dimensions,
  datasets = NULL,
  features_metric = c("signed_importance", "weight", "importance"),
  include_missing_features = TRUE,
  top_n = 5,
  metric = "geometric",
  label_cols = NULL,
  mo_data = NULL,
  truncate = NULL,
  ncol = NULL,
  label_size = 3
)
```

## Arguments

- method_output:

  A single Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function, or a list of two integration method outputs.

- latent_dimensions:

  Character vector of length 2 (if `method_output` is a single output
  object), or named list of length 2 (if `method_output` is a list of
  two output objects). In the first case, gives the name of the two
  latent dimensions that should be represented. In the second case, the
  names of the list should correspond to the names of the methods, and
  the values should each be a character giving for the corresponding
  method the name of the latent dimension to display.

- datasets:

  Character vector, names of the datasets for which features weight
  should be plotted. Default value is `NULL`, i.e. all relevant datasets
  are shown.

- features_metric:

  Character, the features metric that should be plotted on the y-axis.
  Should be one of `'signed_importance'` (default value), `'weight'` or
  `'importance'`.

- include_missing_features:

  Logical, whether to show the features that were in the input of one
  method but not the other, when comparing the results of two different
  integration methods. Default value is `TRUE`.

- top_n:

  Integer, number of top features (according to the consensus importance
  metric) to highlight in the plot. Default value is `5`.

- metric:

  Character, one of the metrics to use to compute the consensus score.
  Can be one of `'min'`, `'max'`, `'average'`, `'product'`, `'l2'` (for
  L2-norm), `'geometric'` (for geometric mean) or `'harmonic'` (for
  harmonic mean). Default value is `'geometric'`. Names must match those
  of `output_list`.

- label_cols:

  Character or named list of character, giving for each dataset the name
  of the column in the corresponding features metadata to use as label.
  If one value, will be used for all datasets. If list, the names must
  correspond to the names of the datasets in `mo_data`. If a dataset is
  missing from the list or no value is provided, feature IDs will be
  used as labels. Alternatively, use `feature_id` to get the feature IDs
  as labels.

- mo_data:

  A `MultiDataSet` object. Only used if `label_cols` is not `NULL`.

- truncate:

  Integer, width to which the labels should be truncated (to avoid
  issues with very long labels in plots). If `NULL` (default value), no
  truncation will be performed.

- ncol:

  Integer, number of columns of datasets in the combined plot. Default
  value is `NULL`, i.e. will be picked automatically.

- label_size:

  Integer, size of the features label. Default value is `3`.

## Value

A ggplot.
