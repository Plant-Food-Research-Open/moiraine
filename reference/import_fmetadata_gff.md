# Import features metadata from a GFF/GTF file

Reads a GFF or GTF annotation file and returns a dataframe in which the
rows correspond to features (e.g. genes or transcripts) and columns
correspond to information about the features. Non-ASCII characters are
replaced with ASCII equivalents (using the stringi and textclean
packages).

## Usage

``` r
import_fmetadata_gff(file, feature_type, add_fields = NULL)
```

## Arguments

- file:

  Character, path to the dataset GFF or GTF file.

- feature_type:

  Character, the type of feature to extract from the annotation file.
  Currently supports `'genes'` and `'transcripts'`.

- add_fields:

  Character vector, fields in the GFF/GTF file to extract that are not
  imported by default (only use if you've run the function once and
  realised that some fields are not extracted by the function).

## Value

A data-frame with the features as rows and features information as
columns. Feature IDs are used as row names.

## Examples

``` r
if (FALSE) { # \dontrun{
import_fmetadata_gff(
  "bos_taurus_gene_model.gff3",
  "genes",
  add_fields = c("name", "description")
)
} # }
```
