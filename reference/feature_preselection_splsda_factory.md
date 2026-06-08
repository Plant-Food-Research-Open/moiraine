# Target factory for feature preselection based on sPLS-DA

Creates a list of targets to perform feature preselection on datasets
from a `MultiDataSet` object with sPLS-DA (from the `mixOmics` package).

## Usage

``` r
feature_preselection_splsda_factory(
  mo_data_target,
  group,
  to_keep_ns,
  to_keep_props = NULL,
  target_name_prefix = "",
  filtered_set_target_name = NULL,
  multilevel = NULL,
  seed_perf = NULL,
  seed_run = NULL,
  ...
)
```

## Arguments

- mo_data_target:

  Symbol, the name of the target containing the `MultiDataSet` object.

- group:

  Character, the column name in the samples information data-frame to
  use as samples group.

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

- target_name_prefix:

  Character, a prefix to add to the name of the targets created by this
  target factory. Default value is `""`.

- filtered_set_target_name:

  Character, the name of the final target containing the filtered
  `MultiDataSet` object. If NULL, a name will automatically be supplied.
  Default value is `NULL`.

- multilevel:

  Character vector of length 1 or 3 to be used as information about
  repeated measurements. See
  [`get_input_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_splsda.md)
  for details. Default value is `NULL` (no repeated measurements).

- seed_perf:

  Named integer vector, the seed to use for the
  [`perf_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/perf_splsda.md)
  function for each dataset. The length and names should match those of
  `to_keep_ns` or `to_keep_props`. If not named, the values will be used
  in order of the datasets in `to_keep_ns` or `to_keep_props`. Default
  value is `NULL`, i.e. no seed is set.

- seed_run:

  Named integer vector, the seed to use for the
  [`run_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/run_splsda.md)
  function for each dataset. The length and names should match those of
  `to_keep_ns` or `to_keep_props`. If not named, the values will be used
  in order of the datasets in `to_keep_ns` or `to_keep_props`. Default
  value is `NULL`, i.e. no seed is set.

- ...:

  Further arguments passed to the
  [perf_splsda](https://plant-food-research-open.github.io/moiraine/reference/perf_splsda.md)
  function.

## Value

A list of target objects. With `target_name_prefix = ""` and
`filtered_set_target_name = NULL`, the following targets are created:

- `splsda_spec`: generates a grouped tibble where each row corresponds
  to one dataset to be filtered, with the columns specifying each
  dataset name, and associated values from `to_keep_ns` and
  `to_keep_props`.

  - `individual_splsda_input`: a dynamic branching target that runs the
    [`get_input_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_splsda.md)
    function for each dataset.

- `individual_splsda_perf`: a dynamic branching target that runs the
  [`perf_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/perf_splsda.md)
  function for each dataset.

- `individual_splsda_run`: a dynamic branching target that runs the
  [`run_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/run_splsda.md)
  function for each dataset, using the results from
  `individual_splsda_perf` to guide the number of latent components to
  construct.

- `filtered_set_slpsda`: a target to retain from the original
  `MultiDataSet` object only features selected in each sPLS-DA run.

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
  feature_preselection_splsda_factory(
    mo_set,
    group = "outcome_group",
    to_keep_ns = c("rnaseq" = 1000, "metabolome" = 500),
    filtered_set_target_name = "mo_set_filtered",
    folds = 10 ## example of an argument passed to perf_splsda
  ),

  ## Another example using to_keep_props
  feature_preselection_splsda_factory(
    mo_set,
    group = "outcome_group",
    to_keep_ns = NULL,
    to_keep_props = c("rnaseq" = 0.3, "metabolome" = 0.5),
    filtered_set_target_name = "mo_set_filtered",
    folds = 10 ## example of an argument passed to perf_splsda
  )
)
} # }
```
