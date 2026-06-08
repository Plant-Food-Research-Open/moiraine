# Violin plot function function

Creates a violin plot

## Usage

``` r
plot_x_discrete(
  toplot,
  x,
  y,
  colour,
  shape,
  point_alpha = 1,
  add_boxplot = TRUE
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

- add_boxplot:

  Logical, should a boxplot be drawn on top of the points for violin
  plots? Default value is `TRUE`.

## Value

A ggplot.

## Details

If `colour_by` is used, and the corresponding variable is numeric, the
violins and boxplots will not take into account this variable. If
instead the `colour_by` variable is a character or factor, a separate
violin and boxplot will be drawn for each category.
