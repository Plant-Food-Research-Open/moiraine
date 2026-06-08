# Generate MEFISTO input data

Creates an object that can be used as input for the MEFISTO analysis
implemented in the MOFA2 package. It contains the omics datasets as well
as features and samples metadata.

## Usage

``` r
get_input_mefisto(
  mo_data,
  covariates,
  datasets = names(mo_data),
  groups = NULL,
  options_list = NULL,
  only_common_samples = FALSE
)
```

## Arguments

- mo_data:

  A
  [`MultiDataSet-class`](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- covariates:

  Character or character vector of length 2, the column name(s) in the
  samples metadata data-frames to use as continuous covariates.

- datasets:

  Character vector, the names of the datasets from `mo_data` to include
  in the analysis.

- groups:

  Character, the column name in the samples metadata data-frames to use
  as groups (use
  [`get_samples_metadata`](https://plant-food-research-open.github.io/moiraine/reference/get_samples_metadata.md)
  to view the samples metadata data-frame for each dataset).

- options_list:

  A named list. Should contain at most 4 elements, named 'data_options',
  'model_options', 'training_options' and 'mefisto_options'. Provide
  respectively the data, model, training and mefisto options to apply
  for the MEFISTO run. See
  [`get_default_data_options`](https://rdrr.io/pkg/MOFA2/man/get_default_data_options.html),
  [`get_default_model_options`](https://rdrr.io/pkg/MOFA2/man/get_default_model_options.html),
  [`get_default_training_options`](https://rdrr.io/pkg/MOFA2/man/get_default_training_options.html)
  and
  [`get_default_mefisto_options`](https://rdrr.io/pkg/MOFA2/man/get_default_mefisto_options.html).

- only_common_samples:

  Logical, whether only the samples present in all datasets should be
  returned. Default value is `FALSE`.

## Value

A [`MOFA`](https://rdrr.io/pkg/MOFA2/man/MOFA.html) object.
