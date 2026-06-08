# Plot adjusted cross-validation results for sO2PLS

Plots the results of an adjusted cross-validation for an sO2PLS run
(from the `OmicsPLS` package)

## Usage

``` r
so2pls_plot_cv_adj(cv_res, with_labels = TRUE)
```

## Arguments

- cv_res:

  Data-frame, output from the
  [`crossval_o2m_adjR2`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m_adjR2.html)
  function.

- with_labels:

  Boolean, whether the optimal values for `nx` and `ny` for each value
  of `n` should be displayed. Default value is `TRUE`.

## Value

A ggplot
