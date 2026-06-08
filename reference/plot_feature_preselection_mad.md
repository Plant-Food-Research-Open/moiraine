# Diagnostics plots for MAD-based feature preselection

Displays the MAD distribution across all features in the original (i.e.
non-filtered) datasets, with a vertical red line showing the cut-off
used by the preselection function.

## Usage

``` r
plot_feature_preselection_mad(mad_list)
```

## Arguments

- mad_list:

  A list with the result from the MAD calculation for each dataset to be
  filtered, as returned by the
  [`select_features_mad`](https://plant-food-research-open.github.io/moiraine/reference/select_features_mad.md)
  function.

## Value

A ggplot.
