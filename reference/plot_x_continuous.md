# Scatter plot function

Creates a scatter plot

## Usage

``` r
plot_x_continuous(toplot, x, y, colour, shape, point_alpha = 1, add_se = TRUE)
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

## Value

A ggplot.

## Details

The function adds a loess curve to summarise the trend between the
covariate and the samples score. If `colour_by` is used, and the
corresponding variable is numeric, the loess curve will not take into
account this variable. If instead the `colour_by` variable is a
character or factor, a loess curve will be fitted separately for each
category.
