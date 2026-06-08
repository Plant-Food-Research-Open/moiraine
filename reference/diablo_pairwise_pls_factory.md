# Target factory for pairwise PLS and design matrix estimation for DIABLO run

Creates a list of targets to perform a PLS run on each pair of datasets,
and uses the results to assess the correlation between datasets and
create a design matrix for the DIABLO algorithm.

## Usage

``` r
diablo_pairwise_pls_factory(
  mixomics_data,
  ...,
  threshold = 0.8,
  low_val = 0.1,
  high_val = 1,
  y_val = 1,
  target_name_prefix = ""
)
```

## Arguments

- mixomics_data:

  A `mixOmics` input object created with
  [`get_input_mixomics_supervised`](https://plant-food-research-open.github.io/moiraine/reference/get_input_mixomics_supervised.md).

- ...:

  Additional parameters to be passed to the
  [`run_pairwise_pls`](https://plant-food-research-open.github.io/moiraine/reference/run_pairwise_pls.md)
  function.

- threshold:

  Numeric, correlation value above which datasets are considered as
  highly correlated (see Details). Default value is 0.8.

- low_val:

  Numeric, value in the design matrix for datasets that are not highly
  correlated. Default value is 0.1.

- high_val:

  Numeric, value in the design matrix for datasets that are highly
  correlated. Default value is 1.

- y_val:

  Numeric, value in the design matrix between datasets and the outcome
  (Y). Default value is 1.

- target_name_prefix:

  Character, prefix to add to the name of the targets created by the
  factory. Default value is `""`.

## Value

A list of targets. For example, with `target_name_prefix = ""`, the
following targets are created:

- `diablo_pairs_datasets`: a target that generates a list of all
  possible pairs of dataset names.

- `diablo_pls_runs_list`: a dynamic branching target that runs the PLS
  algorithm on each possible pair of datasets. The target returns a list
  with the PLS results for each pair of datasets.

- `diablo_pls_correlation_matrix`: a target that computes from the PLS
  results list a correlation matrix between the datasets.

- `diablo_design_matrix`: a target that constructs from the datasets
  correlation matrix a design matrix to use for the DIABLO algorithm.

## Examples

``` r
if (FALSE) { # \dontrun{
## in the _targets.R file
library(moiraine)

list(
  ## code to import the datasets, etc
  ## mo_set is the target containing the MultiDataSet object
  tar_target(
    mixomics_input,
    get_input_mixomics_supervised(mo_set, "outcome_group")
  ),
  diablo_pairwise_pls_factory(mixomics_input)
)
} # }
```
