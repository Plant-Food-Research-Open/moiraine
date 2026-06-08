# Target factory for GFF/GTF features metadata import

Creates a list of targets that track each file and import the features
metadata from each GFF/GTF file.

## Usage

``` r
import_fmetadata_gff_factory(
  files,
  feature_types,
  add_fieldss,
  target_name_suffixes
)
```

## Arguments

- files:

  Character vector, vector of paths to the samples metadata GFF or GTF
  files.

- feature_types:

  Character vector, the type of features to extract from each annotation
  file. Currently supports `'genes'` and `'transcripts'`.

- add_fieldss:

  List, where each element is a character vector of field names in each
  GFF/GTF file to extract that are not imported by default. If a
  character vector is provided, will be used for all files to read in.

- target_name_suffixes:

  Character vector, a suffix to add to the name of the targets created
  by this target factory for each dataset.

## Value

A list of target objects. For example, with two files to import and
`target_name_suffixes = c("geno", "transcripto")`, the factory returns
the following targets:

- `fmetadata_file_geno` and `fmetadata_file_transcripto`: targets
  tracking the genomics and transcriptomics annotation files,
  respectively.

- `fmetadata_geno` and `fmetadata_transcripto`: targets that import the
  genomics and transcriptomics features metadata datasets, respectively.

## Examples

``` r
if (FALSE) { # \dontrun{
## in the _targets.R
library(moiraine)

list(
  import_fmetadata_gff_factory(
    c(
      "data/annotation.gff",
      "data/annotationv2.gtf"
    ),
    feature_types = c("genes", "transcripts"),
    add_fieldss = list(
      c("gene_name", "gene_custom_ID"),
      c("transcript_name")
    ),
    target_name_suffixes = c("geno", "transcripto")
  )
)
} # }
```
