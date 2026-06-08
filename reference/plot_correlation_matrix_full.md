# Plots a full correlation matrix (corrplot-style)

Generates a plot of a correlation matrix in the style of the corrplot
package.

## Usage

``` r
plot_correlation_matrix_full(
  mat,
  rows_title = NULL,
  cols_title = NULL,
  title = NULL,
  show_cor = TRUE,
  min_show_cor = 0.2,
  round_cor = 2
)
```

## Arguments

- mat:

  Correlation matrix to plot.

- rows_title:

  Character, title for rows. Default value is `NULL`.

- cols_title:

  Character, title for cols. Default value is `NULL`.

- title:

  Character, title of the plot. Default value is `NULL`.

- show_cor:

  Logical, should the correlation values be added to the plot? Default
  value is `TRUE`.

- min_show_cor:

  Numeric, the minimum value below which correlation coefficients values
  are not added to the plot (i.e. only a circle will appear for these
  values but no text). Ignored if `show_cor` is `FALSE`. Default value
  is 0.2.

- round_cor:

  Integer, how many decimal places to show for the correlation
  coefficients. Ignored if `show_cor` is `FALSE`. Default value is 2.

## Value

a ggplot.
