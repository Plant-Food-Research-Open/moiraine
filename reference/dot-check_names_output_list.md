# Check names of output list

Checks that the names of a list of outputs from several integration
methods are unique. If not named, the name of the method will be used as
name.

## Usage

``` r
.check_names_output_list(output_list)
```

## Arguments

- output_list:

  List of integration methods output, each generated via one of the
  `get_output_*()` function.

## Value

output_list (but named if it wasn't).
