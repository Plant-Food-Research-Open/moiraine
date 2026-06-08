# Plots sample scores as scatterplot matrix

Plots the samples score from a dimension reduction analysis as a matrix
of scatterplots. If there is only one latent dimension, will be plotted
as boxplot instead.

## Usage

``` r
plot_samples_score(
  method_output,
  latent_dimensions = NULL,
  mo_data = NULL,
  colour_upper = NULL,
  colour_diag = colour_upper,
  colour_lower = colour_upper,
  shape_upper = NULL,
  shape_lower = shape_upper,
  scale_colour_upper = NULL,
  scale_colour_diag = NULL,
  scale_colour_lower = NULL,
  scale_shape_upper = NULL,
  scale_shape_lower = NULL,
  title = NULL,
  point_size = 1.5
)
```

## Arguments

- method_output:

  Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- latent_dimensions:

  Character vector giving the latent dimensions to display. If `NULL`
  (default value), all latent dimensions will be shown.

- mo_data:

  A `MultiDataSet` object (will be used to extract samples information).

- colour_upper:

  Character, name of column in one of the samples metadata tables from
  `mo_data` to use for colouring observations in the upper triangle
  plots. Default value is `NULL`.

- colour_diag:

  Character, name of column in one of the samples metadata tables from
  `mo_data` to use for colouring observations in the diagonal plots. By
  default, will follow `colour_upper`.

- colour_lower:

  Character, name of column in one of the samples metadata tables from
  `mo_data` to use for colouring observations in the lower triangle
  plots. By default, will follow `colour_upper`.

- shape_upper:

  Character, name of column in one of the samples metadata tables from
  `mo_data` to use for shaping observations in the upper triangle plots.
  Default value is `NULL`.

- shape_lower:

  Character, name of column in one of the samples metadata tables from
  `mo_data` to use for shaping observations in the lower triangle plots.
  By default, will follow `shape_upper`.

- scale_colour_upper:

  ggplot2 colour scale to use for the upper triangle plots. Default
  value is `NULL` (if `colour_upper` is not `NULL`, will use ggplot2
  default colour scales).

- scale_colour_diag:

  ggplot2 colour scale to use for the diagonal plots. If `NULL`
  (default), the colour scale used for the upper triangle plots will be
  used if `colour_diag` is equal to `colour_upper`; or the colour scale
  used for the lower triangle plots will be used if `colour_diag` is
  equal to `colour_lower`.

- scale_colour_lower:

  ggplot2 colour scale to use for the lower triangle plots. If `NULL`
  (default), the colour scale used for the upper triangle plots will be
  used.

- scale_shape_upper:

  ggplot2 shape scale to use for the upper triangle plots. Default value
  is `NULL` (if `shape_upper` is not `NULL`, will use ggplot2 default
  shape scale).

- scale_shape_lower:

  ggplot2 shape scale to use for the lower triangle plots. If `NULL`
  (default), the shape scale used for the upper triangle plots will be
  used.

- title:

  Character, title of the plot. If `NULL` (default value), the method
  name from `method_output` will be used to construct the plot title.

- point_size:

  Numeric, the size of points (in pt) in the plot. Default value is 1.5.

## Value

A `ggmatrix` plot.

## Examples

``` r
if (FALSE) { # \dontrun{
## Let's say we've already prepared a MultiDataSet mo_data, in which the
## datasets have samples metadata with columns treatment (discrete),
## weeks (continuous), tissue_type (discrete), disease_score (continuous).

library(ggplot2)

pca_res <- run_pca(mo_data, "metabolome")
output_pca <- get_output_pca(output_pca)
pcs <- paste0("Principal component ", 1:4)

# Simple matrix of scatterplot to visualised PCs two by two
plot_samples_score(
  output_pca,
  pcs
)

# Colouring points according to weeks
plot_samples_score(
  output_pca,
  pcs,
  colour_upper = "weeks"
)
# Adding a custom colour palette
plot_samples_score(
  output_pca,
  pcs,
  colour_upper = "weeks",
  scale_colour_upper = scale_colour_viridis_c()
)

# Adding the treatment as shape
plot_samples_score(
  output_pca,
  pcs,
  colour_upper = "weeks",
  shape_upper = "treatment"
)

# Using the lower triangle of the plots to display disease score
# Again can pass custom colour scale through scale_colour_lower
plot_samples_score(
  output_pca,
  pcs,
  colour_upper = "weeks",
  shape_upper = "treatment",
  colour_lower = "disease_score"
)

# By default the diagonal plots follow the colour of the upper plots,
# but can follow the lower plots instead
plot_samples_score(
  output_pca,
  pcs,
  colour_upper = "weeks",
  shape_upper = "treatment",
  colour_lower = "disease_score",
  colour_diag = "disease_score"
)

# or diagonal can show a different variable
plot_samples_score(
  output_pca,
  pcs,
  colour_upper = "weeks",
  shape_upper = "treatment",
  colour_lower = "tissue_type"
)

# also the lower plots can have a different shape than the upper plots
plot_samples_score(
  output_pca,
  pcs,
  colour_upper = "weeks",
  shape_upper = "treatment",
  colour_lower = "disease_score",
  shape_lower = "tissue_type"
)
} # }
```
