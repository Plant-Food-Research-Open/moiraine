# Pairwise PLS datasets comparison

Runs a Projection to Latent Structure (PLS) analysis on a pair of omics
datasets, as per the mixOmics package.

## Usage

``` r
run_pairwise_pls(mixomics_data, datasets_name, ..., verbose = TRUE)
```

## Arguments

- mixomics_data:

  A `mixOmics` input object created with
  [`get_input_mixomics_supervised`](https://plant-food-research-open.github.io/moiraine/reference/get_input_mixomics_supervised.md).

- datasets_name:

  Character vector of length 2, the names of the two omics datasets to
  analyse.

- ...:

  Additional parameters to be passed to the
  [`pls`](https://rdrr.io/pkg/mixOmics/man/pls.html) function.

- verbose:

  Logical, should details be printed during execution? Default value is
  `TRUE`.

## Value

A named list; each element is an object of class
[`pls`](https://rdrr.io/pkg/mixOmics/man/pls.html), which provides the
result of the PLS run. The name of the datasets analysed is stored as a
character vector in the `datasets_name` attribute.

## Details

Note only one latent component is computed as only the first latent
component will be used to assess the correlation between datasets.

## Examples

``` r
if (FALSE) { # \dontrun{
run_pairwise_pls(mo_set, c("rnaseq", "metabolome"))
} # }
```
