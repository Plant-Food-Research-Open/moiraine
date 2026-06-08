# Per-dataset density plot for MultiDataSet object

Displays a density plot of the values in each dataset of a MultiDataSet
object.

## Usage

``` r
plot_density_data(
  mo_data,
  datasets = names(mo_data),
  combined = TRUE,
  scales = "fixed"
)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- datasets:

  Character vector, names of the datasets to include in the plot. By
  default, all datasets are included.

- combined:

  Logical, should the different datasets be represented in the same
  plot? If `FALSE` (default value), each dataset will be represented in
  its own subplot. Default value is `TRUE`.

- scales:

  Character, how should the axes be plotted if `combined = FALSE`. Can
  be either `'fixed'`, i.e. the same limits will be applied to the axes
  of each subplot; or `'free'`, i.e. the axis limits will be adapted to
  each subplot. Ignored if `combined = TRUE`. Default value is `'fixed`.

## Value

A ggplot.
