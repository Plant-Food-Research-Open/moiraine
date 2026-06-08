# Wrapper to create plot

Wrapper around
[`plot_x_continuous()`](https://plant-food-research-open.github.io/moiraine/reference/plot_x_continuous.md)
and
[`plot_x_discrete()`](https://plant-food-research-open.github.io/moiraine/reference/plot_x_discrete.md),
will choose which one to use depending on whether the x-axis variable is
continuous or discrete.

## Usage

``` r
plot_x_wrapper(
  toplot,
  x,
  y,
  colour,
  shape,
  point_alpha = 1,
  add_se = TRUE,
  add_boxplot = TRUE,
  facet_wrap = NULL,
  ncol_wrap = NULL,
  facet_grid = NULL,
  scales_facet = "free_y"
)
```

## Arguments

- toplot:

  Tibble with data to plot.

- x:

  Character, name of the column in `toplot` to use as x-axis.

- y:

  Character, name of the column in `toplot` to use as y-axis.

- colour:

  Character, name of the column in `toplot` to use as colour.

- shape:

  Character, name of the column in `toplot` to use as shape.

- point_alpha:

  Numeric between 0 and 1, the opacity of the points in the plot (with 1
  = fully opaque, and 0 = fully transparent). Default value is `1`.

- add_se:

  Logical, should a confidence interval be drawn around the smoothing
  curve for scatterplots? Default value is `TRUE`.

- add_boxplot:

  Logical, should a boxplot be drawn on top of the points for violin
  plots? Default value is `TRUE`.

- facet_wrap:

  Character, name of the column in `toplot` to use for faceting (using
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)).
  Default is `NULL`.

- ncol_wrap:

  Integer, number of columns in the faceted plot if using `facet_wrap`.
  Default value is `NULL`.

- facet_grid:

  Character vector of length 2, name of the columns in `toplot` to use
  for row (first element) and column (second element) faceting (using
  [`ggplot2::facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html)).
  Will be ignored if `facet_wrap` is not `NULL`. Default is `NULL`.

- scales_facet:

  Character, value to use for the `scales` argument of
  [`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  or
  [`ggplot2::facet_grid()`](https://ggplot2.tidyverse.org/reference/facet_grid.html).
  Default value is `'free_y'`.

## Value

A ggplot.
