# Computes average sample coordinates for sPLS components

Computes the average sample coordinates for sPLS components across the
two datasets.

## Usage

``` r
spls_get_wa_coord(spls_res)
```

## Arguments

- spls_res:

  The output from the
  [`spls_run()`](https://plant-food-research-open.github.io/moiraine/reference/spls_run.md)
  or [`mixOmics::spls()`](https://rdrr.io/pkg/mixOmics/man/spls.html)
  function.

## Value

A matrix of samples coordinates, with samples in rows and joint
components in columns.
