# Filter datasets

Filters datasets by name in the output of an integration method.

## Usage

``` r
.filter_output_datasets(
  method_output,
  datasets,
  fixed_length = NULL,
  method_name = attr(method_output, "method")
)
```

## Arguments

- method_output:

  Integration method output generated via the
  [`get_output()`](https://plant-food-research-open.github.io/moiraine/reference/get_output.md)
  function.

- datasets:

  Character vector giving the datasets to retain in the features weight
  table of the method's output.

- fixed_length:

  Integer, expected length of `datasets`. If `NULL` (default value), the
  length of `datasets` will not be checked.

- method_name:

  Character, name of the method to use in the error message.

## Value

Similar to `method_output`, but the features weight table has been
filtered.
