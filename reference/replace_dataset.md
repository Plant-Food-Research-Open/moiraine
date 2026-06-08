# Replace matrix dataset within a MultiDataSet object

Replaces the matrix for an omics dataset by a new matrix in a
MultiDataSet object.

## Usage

``` r
replace_dataset(mo_data, dataset_name, new_data)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- dataset_name:

  Character, the name of the dataset for which the matrix data should be
  changed.

- new_data:

  Matrix, the new data. Should have features as rows and samples as
  columns. Rownames should match the corresponding feature IDs, colnames
  should match the corresponding sample IDs.

## Value

A MultiDataSet object.
