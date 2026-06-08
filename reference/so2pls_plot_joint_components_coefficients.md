# Plots sO2PLS contributions between datasets joint components

Plots the regression coefficients that link the joint components of the
two datasets, from an SO2PLS run (from the `OmicsPLS` package).

## Usage

``` r
so2pls_plot_joint_components_coefficients(so2pls_res, datasets = NULL)
```

## Arguments

- so2pls_res:

  The output from the [`o2m`](https://rdrr.io/pkg/OmicsPLS/man/o2m.html)
  function.

- datasets:

  Optional, a character vector with the names of the datasets that
  should be plotted. Default is `NULL`, i.e. both datasets are
  considered.

## Value

A [patchwork](https://patchwork.data-imaginist.com/index.html) of
ggplots.
