# Get parameters from DIABLO run

Extracts the `ncomp` and `keepX` parameters from a DIABLO run and format
them into a table.

## Usage

``` r
diablo_get_params(diablo_res)
```

## Arguments

- diablo_res:

  The output from
  [`block.splsda`](https://rdrr.io/pkg/mixOmics/man/block.splsda.html)
  or
  [`diablo_run`](https://plant-food-research-open.github.io/moiraine/reference/diablo_run.md).

## Value

A tibble.
