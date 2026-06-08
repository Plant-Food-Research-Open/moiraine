# Applies Variance Stabilising Transformation (DESeq2) to matrix

Applies the Variance Stabilising Transformation (VST) performed by the
`DESeq2` package via the
[`varianceStabilizingTransformation`](https://rdrr.io/pkg/DESeq2/man/varianceStabilizingTransformation.html)
function. Includes a size factor normalisation prior to the VST. Only
applies to a matrix of count.

## Usage

``` r
transform_vst(mat, return_matrix_only = FALSE)
```

## Arguments

- mat:

  Numeric matrix, must contain integers only.

- return_matrix_only:

  Logical, should only the transformed matrix be returned? If `TRUE`,
  the function will return a matrix. If `FALSE`, the function instead
  returns a list with the transformed data and potentially other
  information relevant to the transformation. Default value is `FALSE`.

## Value

Depending on the `return_matrix_only`, either a matrix of transformed
data, or a list with the following elements:

- `transformed_data`: matrix of the transformed data;

- `info_transformation`: A
  [`DESeqTransform`](https://rdrr.io/pkg/DESeq2/man/DESeqTransform.html)
  object, with details about the transformation.
