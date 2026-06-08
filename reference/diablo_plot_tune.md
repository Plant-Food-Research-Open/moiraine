# Plots DIABLO tune results

Displays the error rate of a DIABLO run cross-validation to estimate the
optimal number of features to retain from each dataset (`keepX`).

## Usage

``` r
diablo_plot_tune(tune_res)
```

## Arguments

- tune_res:

  The cross-validation results, computed with
  [`diablo_tune()`](https://plant-food-research-open.github.io/moiraine/reference/diablo_tune.md).

## Value

A `ggplot2` object.
