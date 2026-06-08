# Run PCA on MultiDataSet

Runs a Principal Component Analysis on an omics dataset from a
MultiDataSet object. This is a wrapper function around the
[`get_dataset_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/get_dataset_matrix.md)
and
[`run_pca_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/run_pca_matrix.md)
functions.

## Usage

``` r
run_pca(
  mo_data,
  dataset_name,
  n_pcs = 10,
  scale = "none",
  center = TRUE,
  method = NULL
)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- dataset_name:

  Character, name of the omics dataset on which a PCA should be run.

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

## Details

To facilitate the use of dynamic branching with the `targets` package,
the `dataset_name` attribute of the resulting object is set as the value
of the `dataset_name` parameter, and can be accessed via
`attr(res_pca, "dataset_name")`.
