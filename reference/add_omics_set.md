# Adds an omics set to a MultiDataSet object

Adds a omics set to an existing MultiDataSet object.

## Usage

``` r
add_omics_set(mo_data, omics_set, ds_name, ...)
```

## Arguments

- mo_data:

  A
  [MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
  object.

- omics_set:

  A Biobase::eSet object, created via
  [`create_omics_set()`](https://plant-food-research-open.github.io/moiraine/reference/create_omics_set.md).
  Currently accepted objects: Biobase::SnpSet, Biobase::ExpressionSet,
  [MetabolomeSet](https://plant-food-research-open.github.io/moiraine/reference/MetabolomeSet-class.md),
  [PhenotypeSet](https://plant-food-research-open.github.io/moiraine/reference/PhenotypeSet-class.md).

- ds_name:

  Character, name of the dataset (will be used as suffix for the name of
  the dataset in the resulting MultiDataSet object).

- ...:

  Further arguments passed to
  `[MultiDataSet::add_snps()], [MultiDataSet::add_rnaseq()], [add_metabo()] or[add_pheno()] (depending on `omics_set\`
  class).

## Value

A
[MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
object, the `mo_data` with `omics_set` as an additional dataset.

## Examples

``` r
if (FALSE) { # \dontrun{
add_omics_set(mo_data, omics_set, "exp1")
} # }
```
