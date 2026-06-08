# Get parameters from sPLS run

Extracts the `ncomp`, `keepX` and `keepY` parameters from a sPLS run and
format them into a table.

## Usage

``` r
spls_get_params(spls_res)
```

## Arguments

- spls_res:

  The output from [`spls`](https://rdrr.io/pkg/mixOmics/man/spls.html)
  or
  [`spls_run`](https://plant-food-research-open.github.io/moiraine/reference/spls_run.md).

## Value

A tibble.
