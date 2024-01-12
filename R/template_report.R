#' Creates Rmd report from template
#'
#' Creates a Rmarkdown report to present the results from the integration analysis.
#' This function only creates a `.Rmd` file and does not knit the document.
#'
#' @param file Name (and path) of the file to be created. Should end with `.Rmd`
#' if `use_quarto` is `FALSE`, or `.qmd` if `use_quarto` is `TRUE`.
#' @param add_sections Character vector, names of the sections to include in the
#' report. Possible values are 'spls', 'so2pls', 'mofa', 'diablo' and 'comparison'.
#' By default, all sections are included.
#' @param overwrite Logical, should existing file be overwritten?
#' @param target_project Character, name of the current targets project (i.e. value to be used for
#' the `TAR_PROJECT` environment variable). If none provided, will be read from the `TAR_PROJECT`
#' environment variable or set to `"main"` if the former is not set.
#' @param use_quarto Boolean, whether to use a Quarto report. Default value is `FALSE`, i.e.
#' uses a Rmarkdown report.
#' @return Invisible character, the path and name of the generated `.Rmd` (or `.qmd`) file.
#' @export
create_report <- function(file,
                          add_sections = c("spls", "so2pls", "mofa", "diablo", "comparison"),
                          overwrite = FALSE,
                          target_project = Sys.getenv("TAR_PROJECT", "main"),
                          use_quarto = FALSE) {
  if (use_quarto) {
    if (!stringr::str_detect(file, "\\.qmd$")) stop("'file' argument should end with '.qmd'.")
    prefix <- "quarto_"
    ext <- ".qmd"
  } else {
    if (!stringr::str_detect(file, "\\.Rmd$")) stop("'file' argument should end with '.Rmd'.")
    prefix <- ""
    ext <- ".Rmd"
  }


  ## Check that the path to file is valid
  dir_name <- dirname(file)
  if (!file.exists(dir_name)) stop("Directory '", dir_name, "' does not exist.")

  ## Check whether the file already exists
  if (file.exists(file) & !overwrite) stop("File '", file, "' already exists. Use 'overwrite = TRUE' to overwrite it.")

  ## Copy the content of the template Rmd report
  report_template_file <- system.file("templates", paste0(prefix, "report_template", ext), package = "moiraine")
  lines_text <- readLines(report_template_file)

  to_replace <- c(
    "_REPORT_FILE_PATH_" = file,
    "_TARGETS_PROJECT_" = target_project
  )

  out_lines <- lines_text

  for (i in names(to_replace)) {
    out_lines <- gsub(
      i,
      to_replace[i],
      out_lines
    )
  }

  .check_names(
    add_sections,
    c("spls", "so2pls", "mofa", "diablo", "comparison"),
    "In 'add_sections' argument, the following section names are not recognised: '_W_'. Possible section names are: '_C_'."
  )

  for (i in add_sections) {
    template_file <- system.file("templates", paste0(prefix, "fragment_report_", i, "_template", ext), package = "moiraine")
    lines_template <- readLines(template_file)
    out_lines <- c(out_lines, lines_template)
  }

  writeLines(out_lines, file)
  message("File '", file, "' created.\n")

  return(invisible(file))
}

#' Creates a target script file from template
#'
#' Creates a target script file form a template multi-omics integration pipeline.
#' This function only creates the script file but does not execute it.
#'
#' @param file Name (and path) of the file to be created. Should end with `.R`.
#' Default value (recommended) is `"_targets.R"` (in the current directory).
#' @param overwrite Logical, should existing file be overwritten?
#' @return The file name (invisibly).
#' @export
create_targets_pipeline <- function(file = "_targets.R", overwrite = FALSE) {
  if (!stringr::str_detect(file, "\\.R$")) stop("'file' argument should end with '.R'.")

  ## Check that the path to file is valid
  dir_name <- dirname(file)
  if (!file.exists(dir_name)) stop("Directory '", dir_name, "' does not exist.")

  ## Check whether the file already exists
  if (file.exists(file) & !overwrite) stop("File '", file, "' already exists. Use 'overwrite = TRUE' to overwrite it.")

  ## Get the template targets file
  targets_template_file <- system.file("templates", "template_target.R", package = "moiraine")

  file.copy(targets_template_file, file, overwrite = overwrite)

  message("File '", file, "' created.\n")

  return(invisible(file))
}

#' Make Quarto report template from Rmd template
#'
#' Generates the Quarto report templates from corresponding Rmd report
#' templates.
#'
#' @noRd
.make_quarto_template <- function() {
  ## Constructing template for report header and first part
  main_rmd <- here::here("inst/templates/report_template.Rmd")

  main_rmd_con <- file(main_rmd, "r")
  main_rmd_lines <- readLines(main_rmd_con)
  close(main_rmd_con)

  yaml_header_end <- which(main_rmd_lines == "---")[2]
  yaml_header <- main_rmd_lines[seq_len(yaml_header_end)]
  content <- main_rmd_lines[-seq_len(yaml_header_end)]

  yaml_header[yaml_header == "date: '`r format(Sys.Date(), \"%B %d, %Y\")`'"] <- "date: today"
  yaml_header[yaml_header == "output:"] <- "format:"
  yaml_header[yaml_header == "  html_document:"] <- "  html:"
  yaml_header <- yaml_header[!(yaml_header %in% c("     toc_float: true", "     theme: flatly"))]

  tmp <- which(yaml_header == "  html:")
  yaml_header <- c(
    yaml_header[seq_len(tmp)],
    "     toc-location: right",
    "     embed-resources: true",
    yaml_header[-seq_len(tmp)]
  )

  content <- content |>
    stringr::str_remove_all(" \\{\\.tabset \\.tabset-pills\\}") |>
    stringr::str_remove_all(" \\{\\.active\\}")

  main_qmd <- here::here("inst/templates/quarto_report_template.qmd")
  main_qmd_con <- file(main_qmd, "w")
  writeLines(c(yaml_header, content), con = main_qmd_con)
  close(main_qmd_con)

  dir("inst/templates/", pattern = "fragment.+Rmd") |>
    purrr::walk(
      function(.x) {
        file_con <- file(here::here("inst/templates", .x), "r")
        content <- readLines(file_con)
        close(file_con)

        content <- content |>
          stringr::str_remove_all(" \\{\\.tabset \\.tabset-pills\\}") |>
          stringr::str_remove_all(" \\{\\.active\\}")

        new_file <- paste0(
          "quarto_",
          stringr::str_replace(.x, "\\.Rmd$", ".qmd")
        )
        new_file_con <- file(here::here("inst/templates", new_file), "w")
        writeLines(content, con = new_file_con)
        close(new_file_con)
      }
    )
}
