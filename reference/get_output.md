# Extract output of integration method in standard format

Extract samples score and features weight from the result of an
integration method. The `get_output()` function provides a wrapper
around the methods' specific `get_output_*()` functions.

## Usage

``` r
get_output(method_output, use_average_dimensions = TRUE)

get_output_pca(method_output)

get_output_splsda(method_output)

get_output_spls(method_output, use_average_dimensions = TRUE)

get_output_diablo(method_output, use_average_dimensions = TRUE)

get_output_mofa2(method_output)

get_output_so2pls(method_output, use_average_dimensions = TRUE)
```

## Arguments

- method_output:

  The output of an integration method.

- use_average_dimensions:

  Logical, should the (weighted) average of the samples scores for each
  latent dimension across the datasets be used? If `FALSE`, a separate
  set of sample scores will be returned for each dataset for each of the
  latent dimensions. Only applies to sPLS, DIABLO and sO2PLS results.
  Default value is `TRUE`.

## Value

An S3 object of class `output_dimension_reduction`, i.e. a named list,
with the following elements:

- `features_weight`: tibble of features weight (loadings) for each
  latent dimension, with columns `feature_id`, `dataset`,
  `latent_dimension`, `weight` (unscaled feature weight for the
  corresponding latent dimension), `importance` (which corresponds to
  the scaled absolute weight, i.e. 1 = feature has the maximum absolute
  weight for the corresponding latent dimension and dataset, 0 = the
  feature was not selected for the corresponding latent dimension)

- `samples_score`: tibble of samples score for each latent component,
  with columns `sample_id`, `latent_dimension`, `score` (unscaled
  samples score for the corresponding latent dimension)

- `variance_explained`: tibble of the fraction of variance explained by
  each latent component for the relevant datasets.
