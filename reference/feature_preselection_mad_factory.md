# Target factory for feature preselection based on Median Absolute Deviation

Creates a list of targets to perform feature preselection on datasets
from a `MultiDataSet` object by retaining features with the highest
Median Absolute Deviation (MAD).

## Usage

``` r
feature_preselection_mad_factory(
  mo_data_target,
  to_keep_ns,
  to_keep_props = NULL,
  with_ties = TRUE,
  target_name_prefix = "",
  filtered_set_target_name = NULL
)
```

## Arguments

- mo_data_target:

  Symbol, the name of the target containing the `MultiDataSet` object.

- to_keep_ns:

  Named integer vector, the number of feature to retain in each dataset
  to be prefiltered (names should correspond to a dataset name). Value
  should be less than the number of features in the corresponding
  dataset. Set to `NULL` in order to use `to_keep_props` instead.

- to_keep_props:

  Named numeric vector, the proportion of features to retain in each
  dataset to be prefiltered (names should correspond to a dataset name).
  Value should be \> 0 and \< 1. Will be ignored if `to_keep_ns` is not
  `NULL`.

- with_ties:

  Should ties be kept together? If `TRUE`, may return more features than
  requested. Default value is `TRUE`.

- target_name_prefix:

  Character, a prefix to add to the name of the targets created by this
  target factory. Default value is `""`.

- filtered_set_target_name:

  Character, the name of the final target containing the filtered
  `MultiDataSet` object. If NULL, a name will automatically be supplied.
  Default value is `NULL`.

## Value

A list of target objects. With `target_name_prefix = ""` and
`filtered_set_target_name = NULL`, the following targets are created:

- `mad_spec`: a target that generates a grouped tibble where each row
  corresponds to one dataset to be filtered, with the columns specifying
  each dataset name, and associated values from `to_keep_ns`,
  `to_keep_props` and `with_ties`.

- `mad_mat`: a dynamic branching target that run the
  [`get_dataset_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/get_dataset_matrix.md)
  function for each dataset.

- `individual_mad_values`: a dynamic branching target that runs the
  [`select_features_mad_matrix()`](https://plant-food-research-open.github.io/moiraine/reference/select_features_mad_matrix.md)
  function for each dataset.

- `filtered_set_mad`: a target to retain from the original
  `MultiDataSet` object only features selected based on their MAD
  values.

## Examples

``` r
if (FALSE) { # \dontrun{
## in the _targets.R
library(moiraine)

list(
  ## add code here to load the different datasets

  ## the following target creates a MultiDataSet object from previously
  ## created omics sets (geno_set, trans_set, etc)
  tar_target(
    mo_set,
    create_multiomics_set(geno_set, trans_set, metabo_set, pheno_set)
  ),
  feature_preselection_mad_factory(
    mo_set,
    to_keep_ns = c("rnaseq" = 1000, "metabolome" = 500),
    filtered_set_target_name = "mo_set_filtered"
  ),

  ## Another example using to_keep_props
  feature_preselection_mad_factory(
    mo_set,
    to_keep_ns = NULL,
    to_keep_props = c("rnaseq" = 0.3, "metabolome" = 0.5),
    filtered_set_target_name = "mo_set_filtered"
  )
)
} # }
```
