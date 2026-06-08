# Plots features weight in/not in a set

Plots the distribution of features weight from an integration method,
depending on whether the features belong to a feature set of interest.

## Usage

``` r
plot_features_weight_set(
  method_output,
  feature_set,
  set_name = "set",
  features_metric = c("signed_importance", "weight", "importance"),
  add_missing_features = FALSE,
  mo_data = NULL,
  datasets = NULL,
  latent_dimensions = NULL,
  point_alpha = 0.5,
  add_boxplot = TRUE,
  scales = "free_x"
)
```

## Arguments

- method_output:

  Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- feature_set:

  Character vector, features ID belonging to the features set of
  interest.

- set_name:

  Character, name of the set. Default value is `'set'`.

- features_metric:

  Character, the features metric that should be plotted on the y-axis.
  Should be one of `'signed_importance'` (default value), `'weight'` or
  `'importance'`.

- add_missing_features:

  Logical, whether features that are in a multi-omics dataset (provided
  through the `mo_data` argument) but don't have a weight in the
  integration results (e.g. because they were not selected in the
  pre-processing step) should be added in the results. If `TRUE`
  (default value), they will be added with a weight and importance of 0.

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object. If `add_missing_features` is true, all features in the
  multi-omics dataset with no weight in the integration method result
  will be added with a weight and importance of 0.

- datasets:

  Character vector, name of the datasets for which the features
  importance should be plotted. If `NULL` (default value), all datasets
  will be considered.

- latent_dimensions:

  Character vector, the latent dimensions to represent in the plot. If
  `NULL` (default value), all latent dimensions will be represented.

- point_alpha:

  Numeric between 0 and 1, the opacity of the points in the plot (with 1
  = fully opaque, and 0 = fully transparent). Default value is `0.5`.

- add_boxplot:

  Logical, should a boxplot be drawn on top of the points for
  categorical covariates? Default value is `TRUE`.

- scales:

  Character, value to use for the `scales` argument of
  [`ggplot2::facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html).
  Default value is `'free_x'`.

## Value

a ggplot.
