# Generate input data for MOFA2 package

Creates an object that can be used as input for the MOFA or MEFISTO
analysis implemented in the MOFA2 package. It contains the omics
datasets as well as features and samples metadata. Should not be called
directly; instead use
[get_input_mofa](https://plant-food-research-open.github.io/moiraine/reference/get_input_mofa.md)
or
[get_input_mefisto](https://plant-food-research-open.github.io/moiraine/reference/get_input_mefisto.md).

## Usage

``` r
get_input_mofa2(
  mo_data,
  datasets,
  covariates,
  groups,
  options_list,
  only_common_samples
)
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- datasets:

  Character vector, the names of the datasets from `mo_data` to include
  in the analysis.

- covariates:

  Character or character vector of length 2, the column name(s) in the
  samples metadata data-frames to use as continuous covariates. if NULL,
  creates an input object for MOFA. If not null, creates an input object
  for MEFISTO.

- groups:

  Character, the column name in the samples metadata data-frames to use
  as group.

- options_list:

  A named list. Should contain at most 3 or 4 elements (depending on
  whether the input is for MOFA or MEFISTO), named 'data_options',
  'model_options', 'training_options' and 'mefisto_options' (latter only
  for MEFISTO input).

- only_common_samples:

  Logical, whether only the samples present in all datasets should be
  returned. Default value is `FALSE`.

## Value

A [`MOFA`](https://rdrr.io/pkg/MOFA2/man/MOFA.html) object.
