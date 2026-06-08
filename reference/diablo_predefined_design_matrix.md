# Generate a design matrix for DIABLO

Generates a design matrix for DIABLO, following a predesigned pattern as
recommended by the mixOmics authors.

## Usage

``` r
diablo_predefined_design_matrix(
  datasets_name,
  design_matrix = c("null", "weighted_full", "full")
)
```

## Arguments

- datasets_name:

  Character vector, the names of the datasets to integrate. Should
  include the value `"Y"` to represent the samples outcome groups.

- design_matrix:

  Character, the type of design matrix to generate. Should be one of
  `"null"`, `"weighted_full"` or `"full"`.

## Value

A matrix with as many rows and columns as the length of `datasets_name`,
and filled with either 0 (`design_matrix = "null"`), 0.1
(`design_matrix = "weighted_full"`) or 1 (`design_matrix = "full"`).
Values in the diagonal are set to 0, and values in the `"Y"` row and
columns are set to 1.
