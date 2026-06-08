# Filter latent dimensions in list

Filters latent dimensions by name in a list of outputs from integration
methods.

## Usage

``` r
.filter_output_dimensions_list(
  output_list,
  latent_dimensions,
  all_present = FALSE,
  fixed_length = NULL
)
```

## Arguments

- output_list:

  List of integration method outputs each generated via one of the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- latent_dimensions:

  Named list, where each element is a character vector giving the latent
  dimensions to retain in the corresponding element of `output_list`.
  Names must match those of `output_list`.

- all_present:

  Logical, whether there should be one element in `latent_dimensions`
  for each element of `output_list`. If `TRUE`, an error will be
  returned if the length and names of `output_list` and
  `latent_dimensions` do not match. Default value is `FALSE`.

- fixed_length:

  Integer, expected length of each element of `latent_dimensions`. If
  `NULL` (default value), the length of elements in `latent_dimensions`
  can vary.

## Value

A list of output similar to `output_list`, but the samples score table
or features weight table have been filtered.
