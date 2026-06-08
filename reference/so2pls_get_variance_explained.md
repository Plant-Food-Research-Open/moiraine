# Percentage of variance explained for sO2PLS

Generates a table giving the percentage of variance explained by each
component from an sO2PLS in the corresponding dataset.

## Usage

``` r
so2pls_get_variance_explained(so2pls_res)
```

## Arguments

- so2pls_res:

  The output from the [`o2m`](https://rdrr.io/pkg/OmicsPLS/man/o2m.html)
  function.

## Value

A tibble with columns `latent_dimension`, `dataset` and `prop_var_expl`.
