# Plot sparsity cross-validation results for sO2PLS

Plots the results of a sparsity cross-validation for an sO2PLS run (from
the `OmicsPLS` package).

## Usage

``` r
so2pls_plot_cv_sparsity(cv_res)
```

## Arguments

- cv_res:

  List, result from a call to the
  [`so2pls_crossval_sparsity`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_crossval_sparsity.md).

## Value

A ggplot.

## Details

The produced plot has one facet for each joint component. The x-axis
corresponds to the number of features retained from the X dataset to
construct the joint component, and the y-axis to the number of features
retained from the Y dataset to construct the joint component. The colour
of each point in the ith facet represents the average covariance
obtained between the joint ith components of the two datasets over the
cross-validation folds. The size of the points' shadow correspond to the
covariance standard error across the cross-validation folds. For each
joint component, the setting yielding the maximum average covariance is
highlighted in orange, and the one yielding the highest average
covariance under the 1-SD rule in red.
