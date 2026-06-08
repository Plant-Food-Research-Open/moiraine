# Plots sO2PLS joint components samples scores

Plots the samples scores for the average joint components from an sO2PLS
run (from the `OmicsPLS` package).

## Usage

``` r
so2pls_plot_samples_joint_components(so2pls_res, ...)
```

## Arguments

- so2pls_res:

  The output from the [`o2m`](https://rdrr.io/pkg/OmicsPLS/man/o2m.html)
  function.

- ...:

  Further arguments passed to
  [`plot_samples_score()`](https://plant-food-research-open.github.io/moiraine/reference/plot_samples_score.md).

## Value

A ggmatrix plot.
