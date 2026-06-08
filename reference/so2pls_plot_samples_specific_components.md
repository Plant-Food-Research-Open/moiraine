# Plots sO2PLS specific components samples scores

Plots the samples scores for the datasets specific components from an
sO2PLS run (from the `OmicsPLS` package).

## Usage

``` r
so2pls_plot_samples_specific_components(so2pls_res, dataset = NULL, ...)
```

## Arguments

- so2pls_res:

  The output from the [`o2m`](https://rdrr.io/pkg/OmicsPLS/man/o2m.html)
  function.

- dataset:

  Character, the name of the dataset for which the specific components
  should be plotted. Default is `NULL`, i.e. the specific components of
  both datasets are plotted.

- ...:

  Further arguments passed to
  [`plot_samples_score()`](https://plant-food-research-open.github.io/moiraine/reference/plot_samples_score.md).

## Value

A list of ggmatrix plots (one per dataset), or one plot if `dataset` was
used to specify a dataset.
