# Filter latent dimensions

Filters latent dimensions by name in the output of an integration
method.

## Usage

``` r
.filter_output_dimensions(
  method_output,
  latent_dimensions,
  fixed_length = NULL,
  method_name = attr(method_output, "method")
)
```

## Arguments

- method_output:

  Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- latent_dimensions:

  Character vector giving the latent dimensions to retain in the
  method's output.

- fixed_length:

  Integer, expected length of `latent_dimensions`. If `NULL` (default
  value), the length of `latent_dimensions` will not be checked.

- method_name:

  Character, name of the method to use in the error message.

## Value

Similar to `method_output`, but the samples score table or features
weight table have been filtered.
