# Diagnostics plots for COV-based feature preselection

Displays the COV distribution across all features in the original (i.e.
non-filtered) datasets, with a vertical red line showing the cut-off
used by the preselection function.

## Usage

``` r
plot_feature_preselection_cov(cov_list)
```

## Arguments

- cov_list:

  A list with the result from the COV calculation for each dataset to be
  filtered, as returned by the
  [`select_features_cov`](https://plant-food-research-open.github.io/moiraine/reference/select_features_cov.md)
  function.

## Value

A ggplot.
