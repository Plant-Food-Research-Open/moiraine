# Print sparsity cross-validation results for sO2PLS

Prints the results of a sparsity cross-validation for an sO2PLS run
(from the `OmicsPLS` package)

## Usage

``` r
so2pls_print_cv_sparsity(cv_res_optim)
```

## Arguments

- cv_res_optim:

  Named list, output from the
  [`so2pls_get_optim_keep`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_optim_keep.md)
  function.

## Value

A tibble, giving for each dataset (`dataset` column) and joint component
(other columns) the optimal number of features to retain, as well as the
total number of features per dataset to retain.
