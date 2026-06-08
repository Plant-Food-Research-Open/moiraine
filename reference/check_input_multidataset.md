# Check a MultiDataSet input

Checks a MultiDataSet object provided as an input. In particular, checks
that

1.  the input object is a `MultiDataSet` object, 2) the datasets stored
    match the datasets named provided (if any). Will restrict the
    MultiDataSet object to only necessary datasets.

## Usage

``` r
check_input_multidataset(x, datasets = NULL)
```

## Arguments

- x:

  A input object that hopefully is a `MultiDataSet` object.

- datasets:

  Character vector of dataset names that should be in `x`. If `NULL`
  (default value), dataset names will not be checked

## Value

The MultiDataSet object restricted to the datasets required (if
`datasets` is not `NULL`.
