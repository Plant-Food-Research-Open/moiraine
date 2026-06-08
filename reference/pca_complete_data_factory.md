# Target factory for PCA run and missing values imputation on each omics dataset

Creates a list of targets that perform a PCA run for each omics dataset
from a `MultiDataSet` object using dynamic branching, and imputes the
missing values in those datasets using the results of the PCA runs.

## Usage

``` r
pca_complete_data_factory(
  mo_data_target,
  dataset_names = NULL,
  target_name_prefix = "",
  complete_data_name = NULL,
  ...
)
```

## Arguments

- mo_data_target:

  Symbol, the name of the target containing the `MultiDataSet` object.

- dataset_names:

  Character vector, the names of the datasets on which a PCA should be
  run. If `NULL`, a PCA will be run on all datasets. Default value is
  `NULL`.

- target_name_prefix:

  Character, prefix to add to the name of the targets created by the
  factory. Default value is `""`.

- complete_data_name:

  Character, the name of the target containing the `MultiDataSet` with
  missing data imputed to be created. If `NULL`, will be selected
  automatically. Default value is `NULL`.

- ...:

  Further arguments passed to the
  [`run_pca_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/run_pca_matrix.md)
  function.

## Value

A List of targets. If `target_name_prefix = ""` and
`complete_data_name = NULL`, the following targets are created:

- `dataset_names_pca`: target containing a character vector that gives
  the names of the datasets on which a PCA should be run.

- `dataset_mats_pca`: a dynamic branching target that applies the
  [`get_dataset_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/get_dataset_matrix.md)
  function to each dataset specified in `dataset_names`. The results are
  saved in a list. Note that because it is using dynamic branching, the
  names of the list are not meaningful. Rather, use
  `sapply(pca_pca_runs_listruns_list, attr, "dataset_name")` to assess
  which element of the list corresponds to which omics dataset.

- `pca_runs_list`: a dynamic branching target that applies the
  [`run_pca_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/run_pca_matrix.md)
  function to each matrix in `dataset_mats_pca`. The results are saved
  in a list. Note that because it is using dynamic branching, the names
  of the list are not meaningful. Rather, use
  `sapply(pca_runs_list, attr, "dataset_name")` to assess which element
  of the list corresponds to which omics dataset.

- `complete_set`: a target that returns a `MultiDataSet` in which
  missing values have been imputed.

## Examples

``` r
if (FALSE) { # \dontrun{
## in the _targets.R
library(moiraine)
list(
  # ... code for importing datasets etc

  ## mo_set is the target containing the MultiDataSet object
  ## Example 1: running a PCA on all datasets
  run_pca_factory(mo_set),

  ## Example 2: running a PCA on 'rnaseq' and 'metabolome' datasets
  run_pca_factory(
    mo_set,
    c("rnaseq", "metabolome"),
    complete_data_name = "mo_data_complete"
  )
)
} # }
```
