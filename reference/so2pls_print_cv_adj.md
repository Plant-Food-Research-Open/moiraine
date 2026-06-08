# Print adjusted cross-validation results for sO2PLS

Prints the results of an adjusted cross-validation for an sO2PLS run
(from the `OmicsPLS` package)

## Usage

``` r
so2pls_print_cv_adj(cv_res)
```

## Arguments

- cv_res:

  Data-frame, output from the
  [`crossval_o2m_adjR2`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m_adjR2.html)
  function.

## Value

A tibble.
