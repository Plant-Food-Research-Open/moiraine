# Round values in omics dataset from MultiDataSet object

Rounds the values of a given omics dataset within a MultiDataSet object.
Can also limit the range of possible values.

## Usage

``` r
round_dataset(
  mo_data,
  dataset_name,
  ndecimals = 0,
  min_val = -Inf,
  max_val = Inf
)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- dataset_name:

  Character, the name of the dataset for which the matrix data should be
  changed.

- ndecimals:

  Integer, the number of decimals to keep in the dataset. Default value
  is 0.

- min_val:

  Numeric, the minimum value allowed in the dataset. Values below
  `min_val` will be set to `min_val`.

- max_val:

  Numeric, the maximum value allowed in the dataset. Values above
  `max_val` will be set to `max_val`.

## Value

A MultiDataSet object.

## Examples

``` r
if (FALSE) { # \dontrun{
## Let's imagine that we imputed missing values in the genomics dataset from
## mo_data using NIPALS-PCA. The imputed values are continuous, but the
## dataset contains dosage values for a diploid organism (i.e. values can
## be 0, 1, 2). We'll round the imputed values and make sure they can't be
## negative or higher than 2.
round_dataset(mo_data, "snps", min_val = 0, max_val = 2)
} # }
```
