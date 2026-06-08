# Plots sample scores against covariate

Plots the samples score from the result of an integration method against
a covariate from the samples metadata.

## Usage

``` r
plot_samples_score_covariate(
  method_output,
  mo_data,
  covariate,
  latent_dimensions = NULL,
  colour_by = NULL,
  shape_by = NULL,
  point_alpha = 1,
  add_se = TRUE,
  add_boxplot = TRUE,
  ncol = NULL
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

  Character, name of column in one of the samples metadata tables from
  `mo_data` to use as x-axis in the plot.

- latent_dimensions:

  Character vector giving the latent dimensions to display. Default
  value is `NULL`, i.e. all latent dimensions will be shown.

- colour_by:

  Character, name of column in one of the samples metadata tables from
  `mo_data` to use to colour the samples in the plot. Default value is
  `NULL`.

- shape_by:

  Character, name of column in one of the samples metadata tables from
  `mo_data` to use as shape for the samples in the plot. Default value
  is `NULL`.

- point_alpha:

  Numeric between 0 and 1, the opacity of the points in the plot (with 1
  = fully opaque, and 0 = fully transparent). Default value is `1`.

- add_se:

  Logical, should a confidence interval be drawn around the smoothing
  curves for numerical covariates? Default value is `TRUE`.

- add_boxplot:

  Logical, should a boxplot be drawn on top of the points for
  categorical covariates? Default value is `TRUE`.

- ncol:

  Integer, number of columns in the faceted plot. Default value is
  `NULL`.

## Value

a ggplot.

## Details

If the covariate is numeric, the function creates a scatter plot, with a
loess curve to summarise the trend between the covariate and the samples
score. If `colour_by` is used, and the corresponding variable is
numeric, the loess curve will not take into account this variable. If
instead the `colour_by` variable is a character or factor, a loess curve
will be fitted separately for each category.

If the covariate is not numeric, the function creates a violin/boxplot.
If `colour_by` is used, and the corresponding variable is numeric, the
violins and boxplots will not take into account this variable. If
instead the `colour_by` variable is a character or factor, a separate
violin and boxplot will be drawn for each category.
