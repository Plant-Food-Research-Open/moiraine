# Run PCA on matrix

Runs a Principal Component Analysis on an omics matrix, using the
[`pcaMethods::pca()`](https://rdrr.io/pkg/pcaMethods/man/pca.html)
function.

## Usage

``` r
run_pca_matrix(mat, n_pcs = 10, scale = "none", center = TRUE, method = NULL)
```

## Arguments

- mat:

  Matrix of omics measurement, with features as rows and samples as
  columns.

- n_pcs:

  numeric, number of Principal Components to compute. Default value is
  10.

- scale:

  character, type of scaling that should be applied to the dataset
  before running the PCA. Should be one of `'none'`, `'pareto'`,
  `'vector'`, `'uv'` (see
  [`pcaMethods::pca()`](https://rdrr.io/pkg/pcaMethods/man/pca.html)).
  Default value is `'none'`.

- center:

  boolean, should the dataset be centred prior to running the PCA?
  Default value is `TRUE`.

- method:

  character, type of PCA that should be applied to the dataset. See
  [`pcaMethods::listPcaMethods()`](https://rdrr.io/pkg/pcaMethods/man/listPcaMethods.html).
  for a list of available methods. Default value is `'svd'` for datasets
  with no missing value, and `'nipals'` for datasets with missing
  values.

## Value

A [pcaMethods::pcaRes](https://rdrr.io/pkg/pcaMethods/man/pcaRes.html)
object containing the result from the PCA analysis. The attribute
`dataset_name` specifies the name of the dataset analysed.
