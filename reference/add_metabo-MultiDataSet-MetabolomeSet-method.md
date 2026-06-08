# Adds a MetabolomeSet to a MultiDataSet object.

Adds a MetabolomeSet to a MultiDataSet object.

## Usage

``` r
# S4 method for class 'MultiDataSet,MetabolomeSet'
add_metabo(object, met_set, warnings = TRUE, ...)
```

## Arguments

- object:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- met_set:

  A
  [MetabolomeSet](https://plant-food-research-open.github.io/moiraine/reference/MetabolomeSet-class.md)
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
