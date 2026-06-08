# Get weighted average coordinates

Computes the samples coordinates in the weighted average latent
components space from a DIABLO result object.

## Usage

``` r
diablo_get_wa_coord(diablo_res)
```

## Arguments

- diablo_res:

  The output from
  [`block.splsda`](https://rdrr.io/pkg/mixOmics/man/block.splsda.html)
  or
  [`diablo_run`](https://plant-food-research-open.github.io/moiraine/reference/diablo_run.md).

## Value

A matrix with one row per sample and one column per latent component.
