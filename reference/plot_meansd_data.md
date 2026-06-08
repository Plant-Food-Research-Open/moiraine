# Per-dataset mean-sd trend plots for MultiDataSet object

Displays for each dataset in a MultiDataSet object the trend between
features mean and standard deviation across all samples.

## Usage

``` r
plot_meansd_data(
  mo_data,
  datasets = names(mo_data),
  by_rank = FALSE,
  colour_log10 = TRUE
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

- by_rank:

  Logical, should the x-axis display the rank of the features (ordered
  by mean) rather than the features mean? Default value is `FALSE`, i.e.
  the x axis represents the mean of the features.

- colour_log10:

  Should the colour legend be on the log10 scale? Default value is
  `TRUE`.

## Value

A ggplot.
