# Plot of variance explained

Displays the percentage of variance explained by each latent dimension
in the output of a dimension reduction method for each dataset analysed.

## Usage

``` r
plot_variance_explained(
  method_output,
  datasets = NULL,
  latent_dimensions = NULL,
  ncol = NULL,
  free_y_axis = FALSE,
  cumulative = FALSE
)
```

## Arguments

- method_output:

  Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- datasets:

  Character vector giving the datasets to display. Default value is
  `NULL`, i.e. all datasets will be shown.

- latent_dimensions:

  Character vector giving the latent dimensions to display. Default
  value is `NULL`, i.e. all latent dimensions will be shown.

- ncol:

  Integer, number of columns in the faceted plot. Default value is
  `NULL`.

- free_y_axis:

  Logical, whether the y-axis (representing the percentage of variance)
  should have the same range for all datasets. Default value is `FALSE`.

- cumulative:

  Logical, whether the cumulative percentage of variance explained
  should be plotted. Default is `FALSE`.

## Value

A ggplot.
