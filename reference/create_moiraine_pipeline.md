# Creates a target script file from template

Creates a target script file form a template multi-omics integration
pipeline. This function only creates the script file but does not
execute it.

## Usage

``` r
create_moiraine_pipeline(file = "_targets.R", overwrite = FALSE)
```

## Arguments

- file:

  Name (and path) of the file to be created. Should end with `.R`.
  Default value (recommended) is `"_targets.R"` (in the current
  directory).

- overwrite:

  Logical, should existing file be overwritten?

## Value

The file name (invisibly).
