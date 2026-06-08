# Heatmap of samples score correlation

Constructs a upper triangle heatmap of samples score correlation between
the latent dimensions constructed by several integration methods.

## Usage

``` r
.heatmap_samples_score_corr(output_list, hclust_fw)
```

## Arguments

- output_list:

  List of integration methods output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- hclust_fw:

  Dendrogram of latent dimensions according to the features weight
  correlation (obtained with `.heatmap_features_weight_corr`).

## Value

a
[`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)
(upper triangle only).
