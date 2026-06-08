# Perform cross-validation to find the optimal number of features/groups to keep for each joint component for sO2PLS

Computes the optimal number of features/groups to keep for each joint
component for an sO2PLS run. Directly copied from the
[`OmicsPLS::crossval_sparsity()`](https://rdrr.io/pkg/OmicsPLS/man/crossval_sparsity.html)
function, but improved the output for plotting purposes.

## Usage

``` r
so2pls_crossval_sparsity(
  omicspls_input,
  n,
  nx,
  ny,
  nr_folds = 10,
  keepx_seq = NULL,
  keepy_seq = NULL,
  groupx = NULL,
  groupy = NULL,
  tol = 1e-10,
  max_iterations = 100,
  seed = NULL
)
```

## Arguments

- omicspls_input:

  A named list of length 2, produced by
  [`get_input_omicspls()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_omicspls.md).

- n:

  Integer, number of joint PLS components. Must be positive.

- nx:

  Integer, number of orthogonal components in `X`. Negative values are
  interpreted as 0.

- ny:

  Integer, number of orthogonal components in `Y`. Negative values are
  interpreted as 0.

- nr_folds:

  integer, number of folds for the cross-validation. Default value is
  10.

- keepx_seq:

  Numeric vector, how many features/groups to keep for cross-validation
  in each of the joint components of `X`. Sparsity of each joint
  component will be selected sequentially.

- keepy_seq:

  Numeric vector, how many features/groups to keep for cross-validation
  in each of the joint components of `Y`. Sparsity of each joint
  component will be selected sequentially.

- groupx:

  Character vector, group name of each `X`-feature. Its length must be
  equal to the number of features in `X`. The order of the group names
  must corresponds to the order of the features. If `NULL`, no groups
  are considered. Default value is `NULL`.

- groupy:

  Character vector, group name of each `Y`-feature. Its length must be
  equal to the number of features in `Y`. The order of the group names
  must corresponds to the order of the features. If `NULL`, no groups
  are considered. Default value is `NULL`.

- tol:

  Numeric, threshold for which the NIPALS method is deemed converged.
  Must be positive. Default value is `1e-10`.

- max_iterations:

  Integer, maximum number of iterations for the NIPALS method.

- seed:

  Integer, seed to use. Default is `NULL`, i.e. no seed is set inside
  the function.

## Value

A list with the following elements:

- `Best`: a vector giving for each join component the number of features
  to keep from `X` and `Y` that yield the highest covariance between the
  joint components of `X` and `Y` (elements `x1`, `y1`, `x2`, `y2`,
  etc), and the number of features to keep from `X` and `Y` yielding the
  highest covariance under the 1 standard error rule (elements `x_1sd1`,
  `y_1sd1`, `x_1sd2`, `y_1sd2`, etc).

- `Covs`: a list, with as many elements as number of joint components
  (`n`). Each element is a matrix giving the average covariance between
  the joint components of `X` and `Y` obtained across the folds, for
  each tested values of `keepx` (columns) and of `keepy` (rows).

- `SEcov`: a list, with as many elements as number of joint components
  (`n`). Each element is a matrix giving the standard error of the
  covariance between the joint components of `X` and `Y` obtained across
  the folds, for each tested values of `keepx` (columns) and of `keepy`
  (rows).
