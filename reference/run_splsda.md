# Performs sPLS-DA on omics dataset from MultiDataSet object

Performs a sPLS-DA (implemented in the `mixOmics`) package on a omics
dataset from a MultiDataSet object. This is intended for feature
preselection in the omics dataset (see
[`get_filtered_dataset_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/get_filtered_dataset_splsda.md)).

## Usage

``` r
run_splsda(
  splsda_input,
  perf_res,
  to_keep_n = NULL,
  to_keep_prop = NULL,
  ncomp = NULL,
  seed = NULL
)
```

## Arguments

- splsda_input:

  Input for the sPLS-DA functions from mixOmics, created with
  [`get_input_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/get_input_splsda.md).

- perf_res:

  Result of the
  [`perf_splsda()`](https://plant-food-research-open.github.io/moiraine/reference/perf_splsda.md)
  function. If not supplied, sPLS-DA will be run on dataset specified by
  argument `dataset_name` with number of latent components specified by
  argument `comp`.

- to_keep_n:

  Integer, the number of features to retain in the dataset. Should be
  less than the number of features in the dataset. If `NULL` or `NA`,
  `to_keep_prop` will be used instead.

- to_keep_prop:

  Numeric, the proportion of features to retain in the dataset. Will be
  ignored if `to_keep_n` is supplied. Value should be \> 0 and \< 1.

- ncomp:

  Integer, number of latent components to construct. Ignored if
  `perf_res` is supplied. Default value is `NULL`.

- seed:

  Integer, seed to use. Default is `NULL`, i.e. no seed is set inside
  the function.

## Value

A list as per the output of the
[`mixOmics::splsda()`](https://rdrr.io/pkg/mixOmics/man/splsda.html)
function.

## Details

This function uses the
[`mixOmics::plsda()`](https://rdrr.io/pkg/mixOmics/man/plsda.html)
function from the `mixOmics` package. Note that the sPLS-DA method can
select the same feature for several latent components, so the number of
features retained for a dataset might be less than the number specified
in the `to_keep_n` argument.
