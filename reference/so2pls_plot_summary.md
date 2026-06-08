# Plot summary of sO2PLS run

Plots a summary of variation from an sO2PLS run (from the `OmicsPLS`
package).

## Usage

``` r
so2pls_plot_summary(so2pls_res, datasets = NULL)
```

## Arguments

- so2pls_res:

  The output from the [`o2m`](https://rdrr.io/pkg/OmicsPLS/man/o2m.html)
  function.

- datasets:

  Optional, a character vector with the names of the datasets for which
  selected features should be extracted. Default is `NULL`, i.e. both
  datasets are considered.

## Value

A [patchwork](https://patchwork.data-imaginist.com/index.html) of
ggplots.
