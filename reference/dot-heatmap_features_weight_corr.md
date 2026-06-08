# Heatmap of features weight correlation

Constructs a lower triangle heatmap of features weight correlation
between the latent dimensions constructed by several integration
methods.

## Usage

``` r
.heatmap_features_weight_corr(
  output_list,
  include_missing_features = FALSE,
  legend_ncol = 1,
  legend_position = "bottom",
  legend_title = "correlation"
)
```

## Arguments

- output_list:

  List of integration methods output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- include_missing_features:

  Logical, see
  [`get_features_weight_correlation()`](https://plant-food-research-open.github.io/moiraine/reference/get_features_weight_correlation.md)
  for details. Default value is `FALSE`.

- legend_ncol:

  Integer, number of columns in the legend. Default value is `1`.

- legend_position:

  Character, position of the legend. Should be one of `"bottom"`
  (default), `"top"`, `"left"` or `"right`.

- legend_title:

  Character, name to give to the heatmap colour legend. Intended to use
  to shorten the legend title if it goes out of frame.

## Value

a
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
(lower triangle only).
