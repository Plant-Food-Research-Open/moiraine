# Create a MultiDataSet object to store multi-omics data

Creates a MultiDataSet object from a list of Biobase Set objects to
store the different omics sets.

## Usage

``` r
create_multiomics_set(sets_list, datasets_names = NULL, show_warnings = TRUE)
```

## Arguments

- sets_list:

  List of Biobase::eSet objects, created via
  [`create_omics_set()`](https://plant-food-research-open.github.io/moiraine/reference/create_omics_set.md).
  Currently accepted objects: Biobase::SnpSet, Biobase::ExpressionSet,
  [MetabolomeSet](https://plant-food-research-open.github.io/moiraine/reference/MetabolomeSet-class.md),
  [PhenotypeSet](https://plant-food-research-open.github.io/moiraine/reference/PhenotypeSet-class.md).
  Note that feature IDs must be unique between the sets, i.e. a same ID
  cannot be used in several omics sets.

- datasets_names:

  Optional, vector of character, name for each Set object. Will be
  appended to the data type in the resulting object. If the `sets_list`
  list contains several objects of the same data type (e.g. several
  SnpSets), their names must be unique. If "" is provided, no name will
  be appended to the data type for the corresponding dataset.

- show_warnings:

  Logical, should warnings be displayed when adding a set to the
  MultiDataSet object? Default value is `TRUE`.

## Value

[MultiDataSet::MultiDataSet](https://rdrr.io/pkg/MultiDataSet/man/MultiDataSet-class.html)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
## set_geno, set_transcripto and set_metabo are all Set objects
## Generating a MultiDataSet object with standard name
create_multiomics_set(
  list(set_geno, set_transcripto, set_metabo)
)

## Adding custom names for genomics and metabolomics datasets
## but not for the transcriptomics dataset
create_multiomics_set(
  list(set_geno, set_transcripto, set_metabo),
  datasets_names = c("genome1", "", "lcms")
)
} # }
```
