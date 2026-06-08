# Get filtered MultiDataSet object based on sPLS-DA runs

Selects features most associated with the phenotype of interest from
omics datasets based on results from sPLS-DA applied to the
corresponding omics datasets.

## Usage

``` r
get_filtered_dataset_splsda(mo_data, splsda_res_list)
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- splsda_res_list:

  A list with the result from a sPLS-DA run for each dataset to be
  filtered, as returned by the
  [`run_splsda`](https://plant-food-research-open.github.io/moiraine/reference/run_splsda.md)
  function.

## Value

A
[`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
object.

## Details

Note that the sPLS-DA method can select the same feature for several
latent components, so the number of features retained for a dataset
might be less than the number specified in the `to_keep` argument.

## Examples

``` r
if (FALSE) { # \dontrun{
# Goal: keep 20% of features in dataset1, and 50% of features in dataset2
# outcome_group is the outcome of interest in the samples metadata
to_keep_prop <- c("dataset1" = 0.2, "dataset_2" = 0.5)

# 1) assess optimal number of latent components for dataset1 and dataset2
splsda_perf_runs <- lapply(names(to_keep_prop), function(i) {
  perf_splsda(mo_data, i, "outcome_group")
})

# 2) run sPLS-DA with optimal number of latent components for dataset1 and dataset2
splsda_runs <- lapply(splsda_perf_runs, function(x) {
  run_splsda(mo_data, x, to_keep_prop = to_keep_prop[attr(x, "dataset_name")])
})

# 3) Get the filtered dataset
mo_data_filtered <- get_filtered_dataset_splsda(mo_data, splsda_runs)
} # }
```
