# Extract optimal number of components from adjusted cross-validation results for sO2PLS

Extracts the optimal number of components (joint and dataset-specific)
estimated via adjusted cross-validation results for sO2PLS.

## Usage

``` r
so2pls_get_optim_ncomp_adj(cv_res)
```

## Arguments

- cv_res:

  Data-frame, output from the
  [`crossval_o2m_adjR2`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m_adjR2.html)
  function.

## Value

A vector with three integer values:

- `n`: optimal number of joint components

- `nx`: optimal number of specific components for dataset X (first
  dataset)

- `ny`: optimal number of specific components for dataset Y (second
  dataset)
