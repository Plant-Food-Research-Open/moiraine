# Get filtered MultiDataSet object based on variability measure

Selects most highly variable features from omics datasets based on
features' variability (e.g. MAD or COV).

## Usage

``` r
get_filtered_dataset_variability(mo_data, var_list)
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- var_list:

  A list with the result from the MAD or COV calculation for each
  dataset to be filtered, as returned by the
  [`select_features_mad`](https://plant-food-research-open.github.io/moiraine/reference/select_features_mad.md)
  or
  [`select_features_cov`](https://plant-food-research-open.github.io/moiraine/reference/select_features_cov.md)
  function.

## Value

A
[`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Goal: keep 20% of features in dataset1, and 50% of features in dataset2
to_keep_prop <- c("dataset1" = 0.2, "dataset_2" = 0.5)

# 1) compute MAD values and select features for dataset1 and dataset2
mad_list <- lapply(names(to_keep_prop), function(i) {
  select_features_mad(mo_data, i, to_keep_prop[i])
})

# 2) Get the filtered dataset
mo_data_filtered <- get_filtered_dataset_variability(mo_data, mad_list)
} # }
```
