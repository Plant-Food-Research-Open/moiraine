# Applies a normalisation method from bestNormalize to rows of a matrix

Applies an chosen normalisation method to each feature (row) in a
matrix, via the `bestNormalize` package.

## Usage

``` r
transform_bestNormalise_manual(mat, method, return_matrix_only = FALSE, ...)
```

## Arguments

- mat:

  Numeric matrix.

- method:

  Character, name of the normalisation method to apply. Possible values
  are `"arcsinh_x"`, `"boxcox"`, `"center_scale"`, `"exp_x"`, `"log_x"`,
  `"orderNorm"`, `"sqrt_x"`, `"yeojohnson"`. See Details.

- return_matrix_only:

  Logical, should only the transformed matrix be returned? If `TRUE`,
  the function will return a matrix. If `FALSE`, the function instead
  returns a list with the transformed data and potentially other
  information relevant to the transformation. Default value is `FALSE`.

- ...:

  Further arguments passed to the `method` function from the
  `bestNormalize` package.

## Value

Depending on the `return_matrix_only`, either a matrix of transformed
data, or a list with the following elements:

- `transformed_data`: matrix of the transformed data;

- `info_transformation`: A named list with one element per feature
  (row), giving details of the transformation applied to the feature
  (see output for the `bestNormalize` function corresponding to
  `method`).

## Details

Applies a normalisation method implemented in the `bestNormalize`
package. The `method` argument corresponds to the function from the
`bestNormalize` package that will be applied to the rows of the matrix.
See the [vignette from the `bestNormalize`
package](https://cran.r-project.org/web/packages/bestNormalize/vignettes/bestNormalize.html)
for more information about the transformations.
