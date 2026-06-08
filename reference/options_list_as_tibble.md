# Returns options list as a tibble

Transforms a list of options (or parameters) into a tibble with the name
of the options (parameters) in one column, and their value in a second
column. Vector values are collapsed to span only one column.

## Usage

``` r
options_list_as_tibble(options_list)
```

## Arguments

- options_list:

  A named list, where each element corresponds to one option or
  parameter and the name of the element corresponds to the name of the
  option/parameter.

## Value

a tibble, with the `Parameter` column giving the list of the options or
parameters, and the `Value` column giving the values of the
corresponding option or parameter.
