# ggpairs plot with custom colours

Creates a ggpairs plot (see
[`GGally::ggpairs()`](https://ggobi.github.io/ggally/reference/ggpairs.html))
where the colours and shapes can differ between the upper triangle,
lower triangle and diagonal plots.

## Usage

``` r
ggpairs_custom(
  toplot,
  vars,
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

- toplot:

  Tibble in wide format, with observations as rows and variables as
  columns.

- vars:

  Character vector, names of the columns in `toplot` that correspond to
  variables to be used for the plot matrix.

- colour_upper:

  Character, name of column in `toplot` to use for colouring
  observations in the upper triangle plots. Default value is `NULL`.

- colour_diag:

  Character, name of column in `toplot` to use for colouring
  observations in the diagonal plots. By default, will follow
  `colour_upper`.

- colour_lower:

  Character, name of column in `toplot` to use for colouring
  observations in the lower triangle plots. By default, will follow
  `colour_upper`.

- shape_upper:

  Character, name of column in `toplot` to use for shaping observations
  in the upper triangle plots. Default value is `NULL`.

- shape_lower:

  Character, name of column in `toplot` to use for shaping observations
  in the lower triangle plots. By default, will follow `shape_upper`.

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

  Character, title of the plot. Default value is `NULL` (no title added
  to the plot).

- point_size:

  Numeric, the size of points (in pt) in the plot. Default value is 1.5.

## Value

A `ggmatrix` plot.

## Examples

``` r
if (FALSE) { # \dontrun{
library(palmerpenguins)
library(ggplot2)
data("penguins")

vars <- c(
  "bill_length_mm",
  "bill_depth_mm",
  "flipper_length_mm"
)

toplot <- penguins |>
  dplyr::filter(!is.na(bill_length_mm))

# simple scatterplots of the variables
ggpairs_custom(toplot, vars)

# colouring points by species, using custom colour palette
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1")
)

# now adding the sex variable as shape of the points
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1"),
  shape_upper = "sex"
)

# using the lower plots to show the island as colour
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1"),
  shape_upper = "sex",
  colour_lower = "island",
  scale_colour_lower = scale_colour_viridis_d()
)

# showing species and sex in upper plots, body mass and island
# in lower plots
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1"),
  shape_upper = "sex",
  shape_lower = "island",
  colour_lower = "body_mass_g",
  scale_colour_lower = scale_colour_viridis_c(option = "plasma")
)

# same as above, but the diagonal plots show density per year
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1"),
  shape_upper = "sex",
  shape_lower = "island",
  colour_lower = "body_mass_g",
  scale_colour_lower = scale_colour_viridis_c(option = "plasma"),
  colour_diag = "year"
)

# common legend if the diagonal follows the colour of the lower plots
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1"),
  colour_lower = "island",
  scale_colour_lower = scale_colour_brewer(palette = "Accent"),
  colour_diag = "island"
)
} # }
```
