# Wrapper for OmicsPLS::crossval_o2m_adjR2 function

Wrapper function around the
[`crossval_o2m_adjR2`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m_adjR2.html)
function. The main purpose of this wrapper is to add to the result the
names of the datasets to facilitate plotting.

## Usage

``` r
so2pls_crossval_o2m_adjR2(
  omicspls_input,
  a = 1:5,
  ax = seq(0, 10, by = 2),
  ay = seq(0, 10, by = 2),
  nr_folds = 10,
  ...
)
```

## Arguments

- omicspls_input:

  A named list of length 2, produced by
  [`get_input_omicspls`](https://plant-food-research-open.github.io/moiraine/reference/get_input_omicspls.md).

- a:

  Vector of positive integers, number of joint components to test.

- ax:

  Vector of non-negative integers, number of specific components to test
  for the first dataset.

- ay:

  Vector of non-negative integers, number of specific components to test
  for the second dataset.

- nr_folds:

  Positive integer, number of folds to use for the cross-validation.
  Default value is `10`.

- ...:

  Further arguments passed to the
  [`crossval_o2m_adjR2`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m_adjR2.html)
  function.

## Value

A data-frame with four columns: `MSE`, `n`, `nx` and `ny.` Each row
corresponds to an element in `a`.
