# Applies Variance Stabilising Normalisation (vsn) to matrix

Applies the Variance Stabilising Normalisation performed by the `vsn`
package via the [`justvsn`](https://rdrr.io/pkg/vsn/man/justvsn.html)
function.

## Usage

``` r
transform_vsn(mat, return_matrix_only = FALSE, ...)
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

  Further arguments passed to
  [`vsn2`](https://rdrr.io/pkg/vsn/man/vsn2.html).

## Value

Depending on the `return_matrix_only`, either a matrix of transformed
data, or a list with the following elements:

- `transformed_data`: matrix of the transformed data;

- `info_transformation`: NULL.
