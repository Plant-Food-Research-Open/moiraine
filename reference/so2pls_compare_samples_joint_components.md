# Compares sO2PLS samples joint component scores between the two datasets

Plots a comparison of the samples joint component scores obtained for
the two datasets in an sO2PLS run (from the `OmicsPLS` package).

## Usage

``` r
so2pls_compare_samples_joint_components(so2pls_res, components = NULL, ...)
```

## Arguments

- so2pls_res:

  The output from the [`o2m`](https://rdrr.io/pkg/OmicsPLS/man/o2m.html)
  function.

- components:

  Optional, an integer vector with the joint components that should be
  plotted. Default is `NULL`, i.e. all joint components are represented.

- ...:

  Further arguments passed to
  [`plot_samples_score_pair`](https://plant-food-research-open.github.io/moiraine/reference/plot_samples_score_pair.md).

## Value

A [patchwork](https://patchwork.data-imaginist.com/index.html) of
ggplots.
