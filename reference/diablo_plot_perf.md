# Plots DIABLO perf results

Displays the error rate of a DIABLO run cross-validation to estimate the
optimal number of components (`ncomp`)

## Usage

``` r
diablo_plot_perf(perf_res)
```

## Arguments

- perf_res:

  The cross-validation results, computed with
  [`perf`](https://rdrr.io/pkg/mixOmics/man/perf.html).

## Value

A `ggplot2` object.
