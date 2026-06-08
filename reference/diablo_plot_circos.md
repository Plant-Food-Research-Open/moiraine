# Plots DIABLO circos plot

Displays the DIABLO circos plot, but uses available feature metadata to
display feature names.

## Usage

``` r
diablo_plot_circos(diablo_res, mo_data, label_cols, truncate = NULL, ...)
```

## Arguments

- diablo_res:

  The output from
  [`block.splsda`](https://rdrr.io/pkg/mixOmics/man/block.splsda.html)
  or
  [`diablo_run`](https://plant-food-research-open.github.io/moiraine/reference/diablo_run.md).

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- label_cols:

  Character or named list of character, giving for each dataset the name
  of the column in the corresponding features metadata to use as label.
  If one value, will be used for all datasets. If list, the names must
  correspond to the names of the datasets in `mo_data`. Default value is
  `feature_id`, i.e. by default the ID of the features will be used as
  label.

- truncate:

  Integer, width to which the labels should be truncated (to avoid
  issues with very long labels in plots). If `NULL` (default value), no
  truncation will be performed.

- ...:

  Additional arguments passed to
  [`circosPlot`](https://rdrr.io/pkg/mixOmics/man/circosPlot.html).
