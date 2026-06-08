# Screeplot sO2PLS run

Plots the percentage of variation explained by each latent component in
an sO2PLS (from the `OmicsPLS` package).

## Usage

``` r
so2pls_screeplot(so2pls_res, datasets = NULL)
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

## Details

Note that the plots are set up so that it is possible to add a custom
colour palette to get different colours for each dataset (see example).

## Examples

``` r
if (FALSE) { # \dontrun{
## by default, same colour used for both datasets (cannot find a way to fix that cleanly)
so2pls_screeplot(so2pls_final_res)

## Add a colour palette to get different colour for each dataset
so2pls_screeplot(so2pls_final_res) &
 scale_fill_brewer(palette = "Set1", drop = F)
} # }
```
