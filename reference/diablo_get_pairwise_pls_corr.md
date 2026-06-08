# Get pairwise correlations from PLS run

Computes the correlation matrix between datasets based on the first
component of a PLS run.

## Usage

``` r
diablo_get_pairwise_pls_corr(pairwise_pls_result)
```

## Arguments

- pairwise_pls_result:

  List containing the results from the pairwise PLS runs, computed by
  the
  [`run_pairwise_pls`](https://plant-food-research-open.github.io/moiraine/reference/run_pairwise_pls.md)
  function.

## Value

A matrix of correlation coefficients between datasets.

## Details

The correlation coefficient between two datasets is computed as the
correlation coefficient between the first component of each dataset
obtained from a Projection to Latent Structure (PLS) run, from the
`mixOmics` package.

## Examples

``` r
if (FALSE) { # \dontrun{
## get mixomics input from MultiDataSet object
mixomics_data <- get_input_mixomics_supervised(mo_set, "outcome_group")

## Get the list of dataset names
datasets <- setdiff(names(mixomics_data), "Y")

## Get all possible pairwise combinations of dataset names
ds_pairs <- utils::combn(datasets, 2)

## run PLS for each pair of datasets
pls_res_list <- lapply(1:ncol(ds_pairs), function(i) {
  run_pairwise_pls(mo_set, ds_pairs[, i])
})

## extract the pairwise correlation matrix
diablo_get_pairwise_pls_corr(pls_res_list)
} # }
```
