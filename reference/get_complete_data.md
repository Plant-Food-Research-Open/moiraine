# Get MultiDataSet object with imputed values

Replace missing values with imputed values for each dataset of a
MultiDataSet object, based on the results of a Principal Component
Analysis applied to the corresponding dataset.

## Usage

``` r
get_complete_data(mo_data, pca_result)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- pca_result:

  A list in which each element is the result of a PCA run on a different
  dataset, computed with the
  [`run_pca()`](https://plant-food-research-open.github.io/moiraine/reference/run_pca.md)
  function.

## Value

A
[MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
object, for which the assay of each dataset is the imputed dataset.

## Details

Uses the
[`pcaMethods::completeObs()`](https://rdrr.io/pkg/pcaMethods/man/completeObs-nniRes-method.html)
function to impute missing values.
