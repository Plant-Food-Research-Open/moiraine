# Plots features weight against covariate

Plots the features weight or importance from the result of an
integration method against a covariate from the features metadata.

## Usage

``` r
plot_features_weight_covariate(
  method_output,
  mo_data,
  covariate,
  features_metric = c("signed_importance", "weight", "importance"),
  remove_null_weight = FALSE,
  latent_dimensions = NULL,
  colour_by = NULL,
  shape_by = NULL,
  point_alpha = 0.5,
  add_se = TRUE,
  add_boxplot = TRUE,
  scales = "free_x"
)
```

## Arguments

- method_output:

  Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- mo_data:

  A `MultiDataSet` object (will be used to extract samples information).

- covariate:

  Character or named list of character, giving for each dataset the name
  of the column in the corresponding features metadata to use as x-axis
  in the plot. If one value, will be used for all datasets. If list, the
  names must correspond to the names of the datasets in `mo_data`. If a
  dataset is not present in this list, will be excluded from the plot.

- features_metric:

  Character, the features metric that should be plotted on the y-axis.
  Should be one of `'signed_importance'` (default value), `'weight'` or
  `'importance'`.

- remove_null_weight:

  Logical, should features with null weight/importance be removed from
  the plot? Default value is `FALSE`.

- latent_dimensions:

  Character vector giving the latent dimensions to display. Default
  value is `NULL`, i.e. all latent dimensions will be shown.

- colour_by:

  Character or named list of character, giving for each dataset the name
  of column in the corresponding feature metadata to use to colour the
  features in the plot. If one value, will be used for all datasets. If
  list, the names must correspond to the names of the datasets in
  `covariate`. Default value is `NULL`.

- shape_by:

  Character or named list of character, giving for each dataset the name
  of column in the corresponding feature metadata to use as shape for
  the features in the plot. If one value, will be used for all datasets.
  If list, the names must correspond to the names of the datasets in
  `covariate`. Default value is `NULL`.

- point_alpha:

  Numeric between 0 and 1, the opacity of the points in the plot (with 1
  = fully opaque, and 0 = fully transparent). Default value is `0.5`.

- add_se:

  Logical, should a confidence interval be drawn around the smoothing
  curves for numerical covariates? Default value is `TRUE`.

- add_boxplot:

  Logical, should a boxplot be drawn on top of the points for
  categorical covariates? Default value is `TRUE`.

- scales:

  Character, value to use for the `scales` argument of
  [`ggplot2::facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html).
  Default value is `'free_x'`.

## Details

If the covariate is numeric, the function creates a scatter plot, with a
loess curve to summarise the trend between the covariate and the
features weight. If `colour_by` is used, and the corresponding variable
is numeric, the loess curve will not take into account this variable. If
instead the `colour_by` variable is a character or factor, a loess curve
will be fitted separately for each category.

If the covariate is not numeric, the function creates a violin/boxplot.
If `colour_by` is used, and the corresponding variable is numeric, the
violins and boxplots will not take into account this variable. If
instead the `colour_by` variable is a character or factor, a separate
violin and boxplot will be drawn for each category.
