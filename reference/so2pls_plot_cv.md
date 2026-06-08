# Plots cross-validation results for sO2PLS

Plots the results of a cross-validation for an sO2PLS run (from the
`OmicsPLS` package)

## Usage

``` r
so2pls_plot_cv(cv_res, nb_col = NULL)
```

## Arguments

- cv_res:

  A `cvo2m` object, output from the
  [`crossval_o2m`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m.html)
  function.

- nb_col:

  Integer, the number of columns to use for the faceted plot. Default
  value is `NULL` (the number of columns will be chosen automatically).

## Value

A ggplot.
