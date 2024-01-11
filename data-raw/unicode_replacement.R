## ASCII symbols to replace non-ASCII symbols in molecules' name
library(readr)

unicode_replacement <- read_tsv(
  "data-raw/unicode_symbols_replacement.csv",
  quote = ""
)

usethis::use_data(unicode_replacement, overwrite = TRUE, internal = TRUE)
