# Wrapper for OmicsPLS::o2m function

Wrapper function around the
[`o2m`](https://rdrr.io/pkg/OmicsPLS/man/o2m.html) function. The main
purpose of this wrapper is to add to the result the names of the
datasets to facilitate plotting.

## Usage

``` r
so2pls_o2m(
  omicspls_input,
  cv_res = NULL,
  sparsity_res = NULL,
  n = NULL,
  nx = NULL,
  ny = NULL,
  sparse = FALSE,
  keepx = NULL,
  keepy = NULL,
  ...
)
```

## Arguments

- omicspls_input:

  A named list of length 2, produced by
  [`get_input_omicspls`](https://plant-food-research-open.github.io/moiraine/reference/get_input_omicspls.md).

- cv_res:

  Named integer vector of length 3, with names `n`, `nx`, `ny`. Should
  be obtained with
  [`so2pls_get_optim_ncomp_adj`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_optim_ncomp_adj.md)
  or
  [`so2pls_get_optim_ncomp`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_optim_ncomp.md).

- sparsity_res:

  Named list of length 2, with names `keepx` and `keepy`. Should be
  obtained with
  [`so2pls_get_optim_keep`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_get_optim_keep.md).

- n:

  Positive integer, number of joint components to compute. Ignored if
  `cv_res` is not `NULL`.

- nx:

  Positive integer, number of specific components to compute for the
  first dataset. Ignored if `cv_res` is not `NULL`.

- ny:

  Positive integer, number of specific components to compute for the
  second dataset. Ignored if `cv_res` is not `NULL`.

- sparse:

  Logical, should feature selection be performed? Default value is
  `FALSE`. If `sparsity_res` is not `NULL`, will be set to `TRUE`.

- keepx:

  Integer or integer vector of length `n`, number of features from the
  first dataset to retain for each joint component. Ignored if
  `sparsity_res` is not `NULL`.

- keepy:

  Integer or integer vector of length `n`, number of features from the
  second dataset to retain for each joint component. Ignored if
  `sparsity_res` is not `NULL`.

- ...:

  Other arguments passed to
  [`o2m`](https://rdrr.io/pkg/OmicsPLS/man/o2m.html).

## Value

A list (see [`o2m`](https://rdrr.io/pkg/OmicsPLS/man/o2m.html)).
