# Applies a transformation to a dataset from a MultiDataSet object

Applies a transformation to a dataset from a `MultiDataSet` object.
Implemented transformations are: Variance Stabilising Normalisation
(from the `vsn` package), Variance Stabilising Transformation (from the
`DESeq2` package - only for count data), and appropriate feature-wise
normalisation through the `BestNormalise` package.

## Usage

``` r
transform_dataset(
  mo_data,
  dataset,
  transformation,
  return_multidataset = FALSE,
  return_matrix_only = FALSE,
  verbose = TRUE,
  log_base = 2,
  pre_log_function = zero_to_half_min,
  method,
  ...
)
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- dataset:

  Character, name of the dataset to transform.

- transformation:

  Character, transformation to be applied. Possible values are: `vsn`,
  `vst-deseq2`, `logx` `best-normalize-auto` or `best-normalize-manual`.
  See `Details`.

- return_multidataset:

  Logical, should a `MultiDataSet` object with the original data
  replaced by the transformed data returned? If `FALSE`, the output of
  the function depends on `return_matrix_only`. Default value is
  `FALSE`.

- return_matrix_only:

  Logical, should only the transformed matrix be returned? If `TRUE`,
  the function will return a matrix. If `FALSE`, the function instead
  returns a list with the transformed data as well as other information
  relevant to the transformation. Ignored if `return_multidataset` is
  `TRUE`. Default value is `FALSE`.

- verbose:

  Logical, should information about the transformation be printed?
  Default value is `TRUE`.

- log_base:

  Numeric, the base with respect to which logarithms are computed.
  Default value is `2`. Only used if `transformation = 'logx'`.

- pre_log_function:

  Function that will be applied to the matrix before the log
  transformation (e.g. to apply an offset to the values to avoid issues
  with zeros). Default value is the
  [`zero_to_half_min()`](https://plant-food-research-open.github.io/moiraine/reference/zero_to_half_min.md)
  function. Only used if `transformation = 'logx'`.

- method:

  Character, if `transformation = 'best-normalize-manual'`, which
  normalisation method should be applied. See possible values in
  [`transform_bestNormalise_manual()`](https://plant-food-research-open.github.io/moiraine/reference/transform_bestNormalise_manual.md).
  Ignored for other transformations.

- ...:

  Further arguments passed to the
  [`bestNormalize::bestNormalize()`](https://petersonr.github.io/bestNormalize/reference/bestNormalize.html)
  function or the `method` function from the `bestNormalize` package.

## Value

- if `return_multidataset = TRUE`: a
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object, in which the original data for the transformed dataset has
  been replaced.

- if `return_multidataset = FALSE` and `return_matrix_only = TRUE`: a
  matrix with the transformed data.

- if `return_multidataset = FALSE` and `return_matrix_only = FALSE`: a
  list with two elements, `transformed_data` containing a matrix of
  transformed data, and `info_transformation` containing information
  about the transformation (depends on the transformation applied).

## Details

Currently implemented transformations and recommendations based on
dataset type:

- `vsn`: Variance Stabilising normalisation, implemented in the
  [`vsn::justvsn()`](https://rdrr.io/pkg/vsn/man/justvsn.html) function
  from the `vsn` package. This method was originally developed for
  microarray intensities. This transformation is recommended for
  microarray, metabolome, chemical or other intensity-based datasets. In
  practice, applies the
  [`transform_vsn()`](https://plant-food-research-open.github.io/moiraine/reference/transform_vsn.md)
  function.

- `vst-deseq2`: Variance Stabilising Transformation, implemented in the
  [`DESeq2::varianceStabilizingTransformation()`](https://rdrr.io/pkg/DESeq2/man/varianceStabilizingTransformation.html)
  function from the `DESeq2` package. This method is applicable to count
  data only. This transformation is recommended for RNAseq or similar
  count-based datasets. In practice, applies the
  [`transform_vst()`](https://plant-food-research-open.github.io/moiraine/reference/transform_vst.md)
  function.

- `logx`: log-transformation (default to log2, but base can be
  specified). In practice, applies the
  [`transform_logx()`](https://plant-food-research-open.github.io/moiraine/reference/transform_logx.md)
  function.

- `best-normalize-auto`: most appropriate normalisation method
  automatically selected from a number of options, implemented in the
  [`bestNormalize::bestNormalize()`](https://petersonr.github.io/bestNormalize/reference/bestNormalize.html)
  function from the `bestNormalize` package. This transformation is
  recommended for phenotypes that are each measured on different scales
  (since the transformation method selected will potentially be
  different across the features), preferably with a reasonable number of
  features (less than 100) to avoid large computation times. In
  practice, applies the
  [`transform_bestNormalise_auto()`](https://plant-food-research-open.github.io/moiraine/reference/transform_bestNormalise_auto.md)
  function.

- `best-normalize-manual`: performs the same transformation (specified
  through the `method` argument) to each feature of a dataset. This
  transformation is recommended for phenotypes data in which the
  different phenotypes are measured on the same scale. The different
  normalisation methods are:

  - `"arcsinh_x"`: data is transformed as `log(x + sqrt(x^2 + 1))`;

  - `"boxcox"`: Box Cox transformation;

  - `"center_scale"`: data is centered and scaled;

  - `"exp_x"`: data is transformed as `exp(x)`;

  - `"log_x"`: data is transformed as `log_b(x+a)` (`a` and `b` either
    selected automatically per variable or passed as arguments);

  - `"orderNorm"`: Ordered Quantile technique;

  - `"sqrt_x"`: data transformed as `sqrt(x + a)` (`a` selected
    automatically per variable or passed as argument),

  - `"yeojohnson"`: Yeo-Johnson transformation.
