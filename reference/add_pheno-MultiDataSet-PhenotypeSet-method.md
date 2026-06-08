# Adds a PhenotypeSet to a MultiDataSet object.

Adds a PhenotypeSet to a MultiDataSet object.

## Usage

``` r
# S4 method for class 'MultiDataSet,PhenotypeSet'
add_pheno(object, pheno_set, warnings = TRUE, ...)
```

## Arguments

- object:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- pheno_set:

  A
  [PhenotypeSet](https://plant-food-research-open.github.io/moiraine/reference/PhenotypeSet-class.md)
  object.

- warnings:

  Logical, should warnings be displayed? Default is `TRUE`.

- ...:

  Further arguments passed to the
  [`MultiDataSet::add_eset()`](https://rdrr.io/pkg/MultiDataSet/man/add_eset.html)
  function.

## Value

A new
[MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
object with a slot filled.
