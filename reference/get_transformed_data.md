# Get MultiDataSet with transformed data

Replace the original datasets with transformed datasets in a
MultiDataSet object from the results of transformations applied to the
datasets.

## Usage

``` r
get_transformed_data(mo_data, transformation_result)
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- transformation_result:

  A list in which each element is the result of a transformation applied
  to a different dataset, computed with the
  [`transform_dataset`](https://plant-food-research-open.github.io/moiraine/reference/transform_dataset.md)
  function.

## Value

A
[`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
object, for which the assay of each dataset is the imputed dataset.
