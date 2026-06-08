# Wrapper for OmicsPLS::crossval_o2m function

Wrapper function around the
[`crossval_o2m`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m.html)
function. The main purpose of this wrapper is to add to the result the
names of the datasets to facilitate plotting. If the result of a
previous call to
[`crossval_o2m_adjR2`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m_adjR2.html)
or
[`so2pls_crossval_o2m_adjR2`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_crossval_o2m_adjR2.md)
is provided, will be used to set values to test for `a`, `ax` and `ay`.

## Usage

``` r
so2pls_crossval_o2m(
  omicspls_input,
  cv_adj_res = NULL,
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

- cv_adj_res:

  Data-frame returned by
  [`crossval_o2m_adjR2`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m_adjR2.html)
  or
  [`so2pls_crossval_o2m_adjR2`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_crossval_o2m_adjR2.md).
  Default value is `NULL`.

- a:

  Vector of positive integers, number of joint components to test.
  Ignored if `cv_adj_res` is not `NULL`.

- ax:

  Vector of non-negative integers, number of specific components to test
  for the first dataset. Ignored if `cv_adj_res` is not `NULL`.

- ay:

  Vector of non-negative integers, number of specific components to test
  for the second dataset. Ignored if `cv_adj_res` is not `NULL`.

- nr_folds:

  Positive integer, number of folds to use for the cross-validation.
  Default value is `10`.

- ...:

  Further arguments passed to the
  [`crossval_o2m_adjR2`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m_adjR2.html)
  function.

## Value

A list of class "cvo2m" with the original and sorted Prediction errors
and the number of folds used.

## Details

If the result of a previous call to
[`crossval_o2m_adjR2`](https://rdrr.io/pkg/OmicsPLS/man/crossval_o2m_adjR2.html)
or
[`so2pls_crossval_o2m_adjR2`](https://plant-food-research-open.github.io/moiraine/reference/so2pls_crossval_o2m_adjR2.md)
is provided through the `cv_adj_res` parameter, the optimal values for
`n`, `nx` and `ny` are extracted from it, and the values of `a`, `ax`
and `ay` are set as follows:

- `a` = max(`n` - 1, 1):(`n` + 1)

- `ax` = max(`nx` - 1, 0):(`nx` + 1)

- `ay` = max(`ny` - 1, 0):(`ny` + 1)
