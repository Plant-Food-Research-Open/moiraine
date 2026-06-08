# Correlation plot between latent components

Plots the correlation between either the samples score or the features
weight of the latent components obtained from two different integration
methods.

## Usage

``` r
comparison_plot_correlation(
  output_list,
  by = "both",
  latent_dimensions = NULL,
  include_missing_features = FALSE,
  show_cor = TRUE,
  min_show_cor = 0.2,
  round_cor = 2
)
```

## Arguments

- output_list:

  List of length 2 of integration methods output, each generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function. If named, the names will be used to annotate the plot. See
  details.

- by:

  Character, should the correlation be calculated based on samples score
  ( `by = 'samples'`) or on features weight (`by = 'features'`), or both
  (`by = 'both'`, i.e. two matrices will be plotted). Default value is
  `'both'`.

- latent_dimensions:

  Named list, where each element is a character vector giving the latent
  dimensions to retain in the corresponding element of `output_list`.
  Names must match those of `output_list`. Can be used to filter latent
  dimensions only in certain elements from output_list (see examples).
  If `NULL` (default value), all latent dimensions will be used.

- include_missing_features:

  Logical, see
  [`get_features_weight_correlation`](https://plant-food-research-open.github.io/moiraine/reference/get_features_weight_correlation.md)
  for details. Default value is `FALSE`.

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

A ggplot.

## Details

If `output_list` is unnamed, the different elements in the list will be
differentiated by the name of the method used to produce them (e.g.
DIABLO, sO2PLS, etc). In order to compare different results from a same
integration method (e.g. DIABLO applied to the full vs pre-filtered
data), it is possible to assign names to the elements of `output_list`.
These names will be used in place of the method name in the plot to
identify where the latent dimensions come from.
