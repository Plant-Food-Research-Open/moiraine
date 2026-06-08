# Get table with transformation applied to each dataset

From the results of transformations on datasets, generates a table
giving for each dataset the transformation that was applied to it.

## Usage

``` r
get_table_transformations(
  transformation_result,
  best_normalize_details = FALSE
)
```

## Arguments

- transformation_result:

  A list in which each element is the result of a transformation applied
  to a different dataset, computed with the
  [transform_dataset](https://plant-food-research-open.github.io/moiraine/reference/transform_dataset.md)
  function.

- best_normalize_details:

  Logical, should information about the transformations selected by
  bestNormalize for each feature be displayed? Default value is `FALSE`.

## Value

A tibble with columns `'Dataset'` and `'Transformation'`. If
`best_normalize_details = TRUE`, an additional column `'Details'` lists
the chsoen transformation applied to each feature of the corresponding
dataset for a bestNormalize transformation.
