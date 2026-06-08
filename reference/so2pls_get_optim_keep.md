# Extract optimal number of features to keep from cross-validation results for sO2PLS

Extracts the optimal number of features to retain from datasets `X` and
`Y` for the joint components.

## Usage

``` r
so2pls_get_optim_keep(cv_res, use_1sd_rule = TRUE)
```

## Arguments

- cv_res:

  List, result from a call to the
  [`so2pls_crossval_sparsity()`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_crossval_sparsity.md)
  or
  [`OmicsPLS::crossval_sparsity()`](https://rdrr.io/pkg/OmicsPLS/man/crossval_sparsity.html).

- use_1sd_rule:

  Boolean, should the 1 standard deviation rule be used when selecting
  the optimal number of features to retain? See Details.

## Value

A list with elements `keepx` and `keepy`, each a vector of length equal
to the number of joint components, where the ith element giving the
number of features to retain from dataset `X` (`keepx`) or `Y` (`keepy`)
for the i-th joint component.

## Details

The 1-SD rule means that we are retaining the smallest number of
features yielding an average covariance that is within 1SD of the
maximum covariance obtained.
