# Plots DIABLO features correlation circle

Displays the DIABLO correlation circle plot, but uses available feature
metadata to display feature names.

## Usage

``` r
diablo_plot_var(
  diablo_res,
  mo_data,
  label_cols = "feature_id",
  truncate = NULL,
  ...
)
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
  [`plotVar`](https://rdrr.io/pkg/mixOmics/man/plotVar.html).

## Value

A plot (see [`plotVar`](https://rdrr.io/pkg/mixOmics/man/plotVar.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
# Use the default features ID for the plot
diablo_plot_var(
  diablo_final_run,
  mo_data,
  "feature_id",
  overlap = FALSE,
  cex = rep(2, 3),
  comp = 1:2
)

# Using a different column from the feature metadata of each omics dataset
diablo_plot_var(
  diablo_final_run,
  mo_presel_supervised,
  c(
    "snps" = "feature_id",
    "rnaseq" = "gene_name",
    "metabolome" = "comp_name"
  ),
  overlap = FALSE,
  cex = rep(2, 3),
  comp = 1:2
)
} # }
```
