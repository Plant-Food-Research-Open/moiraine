# Plots sPLS features correlation circle

Displays the sPLS correlation circle plot, but uses available feature
metadata to display feature names.

## Usage

``` r
spls_plot_var(
  spls_res,
  mo_data,
  label_cols = "feature_id",
  truncate = NULL,
  ...
)
```

## Arguments

- spls_res:

  The output from
  [`mixOmics::spls()`](https://rdrr.io/pkg/mixOmics/man/spls.html) or
  [`spls_run()`](https://plant-food-research-open.github.io/moiraine/reference/spls_run.md).

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
  [`mixOmics::plotVar()`](https://rdrr.io/pkg/mixOmics/man/plotVar.html).

## Value

A plot (see
[`mixOmics::plotVar()`](https://rdrr.io/pkg/mixOmics/man/plotVar.html)).

## Examples

``` r
if (FALSE) { # \dontrun{
# Use the default features ID for the plot
spls_plot_var(
  spls_final_run,
  mo_data,
  "feature_id",
  overlap = FALSE,
  cex = c(3, 3),
  comp = 1:2
)

# Using a different column from the feature metadata of each omics dataset
spls_plot_var(
  spls_final_run,
  mo_presel_supervised,
  c(
    "rnaseq" = "Name",
    "metabolome" = "name"
  ),
  overlap = FALSE,
  cex = c(3, 3),
  comp = 1:2
)
} # }
```
