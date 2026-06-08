# Applies a log-x transformation to matrix

Applies a log-x transformation (by default log2) through the
[`log()`](https://rdrr.io/r/base/Log.html) function.

## Usage

``` r
transform_logx(
  mat,
  return_matrix_only = FALSE,
  log_base = 2,
  pre_log_function = zero_to_half_min
)
```

## Arguments

- mat:

  Numeric matrix.

- return_matrix_only:

  Logical, should only the transformed matrix be returned? If `TRUE`,
  the function will return a matrix. If `FALSE`, the function instead
  returns a list with the transformed data and potentially other
  information relevant to the transformation. Default value is `FALSE`.

- log_base:

  Numeric, the base with respect to which logarithms are computed.

- pre_log_function:

  Function that will be applied to the matrix before the log
  transformation (e.g. to apply an offset to the values to avoid issues
  with zeros). Default value is the
  [`zero_to_half_min()`](https://plant-food-research-open.github.io/moiraine/reference/zero_to_half_min.md)
  function.

## Value

Depending on the `return_matrix_only`, either a matrix of transformed
data, or a list with the following elements:

- `transformed_data`: matrix of the transformed data;

- `info_transformation`: a list with the log base used and the function
  applied prior to log-transformation.
