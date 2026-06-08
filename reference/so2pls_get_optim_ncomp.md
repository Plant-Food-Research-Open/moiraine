# Extract optimal number of components from cross-validation results for sO2PLS

Extracts the optimal number of components (joint and dataset-specific)
estimated via cross-validation results for sO2PLS.

## Usage

``` r
so2pls_get_optim_ncomp(cv_res)
```

## Arguments

- cv_res:

  A `cvo2m` object, output from the
  [`crossval_o2m`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m.html)
  function.

## Value

A vector with three integer values:

- `n`: optimal number of joint components

- `nx`: optimal number of specific components for dataset X (first
  dataset)

- `ny`: optimal number of specific components for dataset Y (second
  dataset)
