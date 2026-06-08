# Evaluate feature selection against features label

Compares the selection of features with different feature labels (e.g.
result of a DE analysis) for each latent dimension.

## Usage

``` r
evaluate_feature_selection_table(
  method_output,
  mo_data,
  col_names,
  latent_dimensions = NULL
)
```

## Arguments

- method_output:

  Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- col_names:

  Named character vector, giving for each dataset the name of the column
  in the features metadata table that contains the features label. If a
  dataset is not present in this vector, will be excluded from the
  resulting table.

- latent_dimensions:

  Character vector, the latent dimensions to include in the resulting
  table. If `NULL` (default value), all latent dimensions will be
  represented.

## Value

A tibble, with for each dataset and latent dimension the number of
selected and non-selected features per feature label.
