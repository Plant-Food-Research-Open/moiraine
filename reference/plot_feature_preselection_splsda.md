# Diagnostics plots for sPLS-DA-based feature preselection

Displays the PLS-DA classification performance across different number
of latent components for each prefiltered dataset. The classification
error rates are computed with different measures (column facets) and
different distance metrics (colours). A vertical grey bar represents for
each dataset the number of latent components selected for the feature
preselection step. In addition, a circle highlights the measure and
distance metric used to select the number of latent component.

## Usage

``` r
plot_feature_preselection_splsda(
  perf_splsda_res,
  measure = NULL,
  distance = NULL
)
```

## Arguments

- perf_splsda_res:

  A list with the result from the
  [perf_splsda](https://plant-food-research-open.github.io/moiraine/reference/perf_splsda.md)
  for each dataset to be filtered.

- measure:

  Which measure(s) should be displayed? Can be one of `"BER"` or
  `"overall"`. If NULL, all measures will be displayed. Default value is
  `NULL`.

- distance:

  Which measure(s) should be displayed? Can be one of `"max.dist"`,
  `"centroids.dist"` or `"mahalanobis.dist"`. If NULL, all measures will
  be displayed. Default value is `NULL`.

## Value

A ggplot.
