# Target factory for datasets transformation

Create a list of targets to apply some transformation methods to one or
more datasets in a `MultiDataSet` object.

## Usage

``` r
transformation_datasets_factory(
  mo_data_target,
  transformations,
  return_matrix_only = FALSE,
  target_name_prefix = "",
  transformed_data_name = NULL,
  log_bases = 2,
  pre_log_functions = zero_to_half_min,
  methods,
  ...
)
```

## Arguments

- mo_data_target:

  Symbol, the name of the target containing the `MultiDataSet` object.

- transformations:

  Named character vector, name of each element is the name of a dataset
  to transform, corresponding element gives the type of transformation
  to apply to the dataset (e.g.
  `c(rnaseq = 'vst-deseq2', phenotypes = 'best-normalize-auto')`). See
  Details for a list of available transformations. If
  `'best-normalize-auto'` is selected, need to provide the `methods`
  argument as well.

- return_matrix_only:

  Logical, should only the transformed matrix be returned for each
  transformation? If `TRUE`, only transformed matrices will be stored.
  If `FALSE`, instead for each transformation, a list with the
  transformed data and potentially other information relevant to the
  transformation will be saved. Default value is `FALSE`.

- target_name_prefix:

  Character, a prefix to add to the name of the targets created by this
  target factory. Default value is `""`.

- transformed_data_name:

  Character, the name of the target containing the `MultiDataSet` with
  transformed data to be created. If `NULL`, will be selected
  automatically. Default value is `NULL`.

- log_bases:

  Numeric or named numeric list, gives for each dataset for which the
  `'logx'` transformation is selected the log base to use. If one value,
  will be used for all concerned datasets. Otherwise, can specify a
  different log-base for each concerned dataset by passing a named list.

- pre_log_functions:

  Function or named list of functions, gives for each dataset for which
  the \`'logx“ transformation is selected the function that will be
  applied to the matrix before the log transformation (e.g. to apply an
  offset to the values to avoid issues with zeros). Default value is the
  [`zero_to_half_min()`](https://plant-food-research-open.github.io/moiraine/reference/zero_to_half_min.md)
  function. If one value, will be used for all concerned datasets.
  Otherwise, can specify a different log-base for each concerned dataset
  by passing a named list.

- methods:

  Character or named character list, gives for each dataset for which
  the `'best-normalize-manual'` transformation is selected the
  normalisation method that should be applied. See possible values in
  Details. If one value, will be used for all concerned datasets.
  Otherwise, can specify a different method for each concerned dataset
  by passing a named list.

- ...:

  Further arguments passed to the
  [`transform_dataset`](https://plant-food-research-open.github.io/moiraine/reference/transform_dataset.md)
  function or the `method` function from the `bestNormalize` package.
  Only relevant for `'best-normalize-XX'` transformations.

## Value

A list of target objects. With `target_name_prefix = ""` and
`transformed_data_name = NULL`, the following targets are created:

- `transformations_spec`: generates a grouped tibble where each row
  corresponds to one dataset to be tranformed, with the columns
  specifying each dataset name and the transformation to apply.

- `transformations_runs_list`: a dynamic branching target that runs the
  [`transform_dataset()`](https://plant-food-research-open.github.io/moiraine/reference/transform_dataset.md)
  function on each dataset. Returns a list.

- `transformed_set`: a target that returns the `MultiDataSet` object
  with the original data replaced by the transformed data.

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

## Examples

``` r
if (FALSE) { # \dontrun{
## in the _targets.R
library(moiraine)

list(
  ## add code here to load the different datasets

  ## the following target creates a MultiDataSet object from previously
  ## created omics sets (geno_set, trans_set, etc)
  tar_target(
    mo_set,
    create_multiomics_set(geno_set, trans_set, metabo_set, pheno_set)
  ),

  ## Example 1
  transformation_datasets_factory(mo_set,
    c(
      rnaseq = "vst-deseq2",
      metabolome = "vsn",
      phenotypes = "best-normalize-auto"
    ),
    return_matrix_only = FALSE,
    transformed_data_name = "mo_set_transformed"
  ),

  ## Example 2 - with a log2 transformation for both datasets
  transformation_datasets_factory(
    mo_set_complete,
    c(
      "rnaseq" = "logx",
      "metabolome" = "logx"
    ),
    log_bases = 2,
    pre_log_functions = zero_to_half_min
  ),

  ## Example 3 - with different log bases for each dataset and a different
  ## preprocessing function to be run before applying the log
  transformation_datasets_factory(
    mo_set_complete,
    c(
      "rnaseq" = "logx",
      "metabolome" = "logx"
    ),
    log_bases = list(rnaseq = 10, metabolome = 2),
    pre_log_functions = list(
      rnaseq = \(x) x + 0.5,
      metabolome = zero_to_half_min
     )
  )
)
} # }
```
