# Applies the bestNormalize function to rows of a matrix

Applies an appropriate normalisation method to each feature (row) in a
matrix, via the
[`bestNormalize`](https://petersonr.github.io/bestNormalize/reference/bestNormalize.html)
function from the `bestNormalize` package.

## Usage

``` r
transform_bestNormalise_auto(mat, return_matrix_only = FALSE, ...)
```

## Arguments

- mat:

  Numeric matrix.

- return_matrix_only:

  Logical, should only the transformed matrix be returned? If `TRUE`,
  the function will return a matrix. If `FALSE`, the function instead
  returns a list with the transformed data and potentially other
  information relevant to the transformation. Default value is `FALSE`.

- ...:

  Further arguments passed to the
  [`bestNormalize`](https://petersonr.github.io/bestNormalize/reference/bestNormalize.html)
  function.

## Value

Depending on the `return_matrix_only`, either a matrix of transformed
data, or a list with the following elements:

- `transformed_data`: matrix of the transformed data;

- `info_transformation`: A named list with one element per feature
  (row), giving details of the transformation applied to the feature
  (see output of
  [`bestNormalize`](https://petersonr.github.io/bestNormalize/reference/bestNormalize.html)).
