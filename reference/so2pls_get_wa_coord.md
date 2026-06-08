# Computes average sample coordinates for sO2PLS joint components

Computes the average sample coordinates for sO2PLS joint components
across the two datasets.

## Usage

``` r
so2pls_get_wa_coord(so2pls_res)
```

## Arguments

- so2pls_res:

  The output from the [`o2m`](https://rdrr.io/pkg/OmicsPLS/man/o2m.html)
  function.

## Value

A matrix of samples coordinates, with samples in rows and joint
components in columns.
