# Check for missing values in MultiDataSet

Checks if there are missing values in each omics dataset of a
MultiDataSet object.

## Usage

``` r
check_missing_values(mo_data)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

## Value

Invisible logical vector indicating whether missing values are present
in each dataset.
