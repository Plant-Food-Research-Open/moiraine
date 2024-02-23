#' Creates a target script file from template
#'
#' Creates a target script file form a template multi-omics integration
#' pipeline. This function only creates the script file but does not execute it.
#'
#' @param file Name (and path) of the file to be created. Should end with `.R`.
#'   Default value (recommended) is `"_targets.R"` (in the current directory).
#' @param overwrite Logical, should existing file be overwritten?
#' @returns The file name (invisibly).
#' @export
create_moiraine_pipeline <- function(file = "_targets.R", overwrite = FALSE) {
  if (!stringr::str_detect(file, "\\.R$")) {
    stop("'file' argument should end with '.R'.")
  }

  ## Check that the path to file is valid
  dir_name <- dirname(file)
  if (!file.exists(dir_name)) {
    stop("Directory '", dir_name, "' does not exist.")
  }

  ## Check whether the file already exists
  if (file.exists(file) & !overwrite) {
    stop("File '",
         file,
         "' already exists. Use 'overwrite = TRUE' to overwrite it.")
  }

  ## Get the template targets file
  targets_template_file <- system.file(
    "templates",
    "template_target_minimal.R",
    package = "moiraine"
  )

  file.copy(targets_template_file, file, overwrite = overwrite)

  message("File '", file, "' created.\n")

  return(invisible(file))
}
