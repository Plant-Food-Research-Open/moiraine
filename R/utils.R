#' Check a character vector against reference vector
#'
#' Check values in a character vector against values in a reference vector.
#' Generates a useful error message listing the discrepancies if needed.
#'
#' @param names_to_check Character vector to be checked.
#' @param correct_names Character vector used as reference.
#' @param message Character, the error message to display in case of
#'   discrepancy.
#' @param wrong_names_code Character, the code used in `message` to be replaced
#'   with the values from `names_to_check` that are not in `correct_names`.
#' @param correct_names_code Character, the code used in `message` to be
#'   replaced with the values from `correct_names`.
#'
#' @noRd
.check_names <- function(names_to_check,
                         correct_names,
                         message,
                         wrong_names_code = "_W_",
                         correct_names_code = "_C_") {
  wrong_names <- setdiff(names_to_check, correct_names)

  if (length(wrong_names)) {
    msg <- stringr::str_replace(
      message,
      wrong_names_code,
      paste0(wrong_names, collapse = "', '")
    )
    msg <- stringr::str_replace(
      msg,
      correct_names_code,
      paste0(correct_names, collapse = "', '")
    )
    stop(msg, call. = FALSE)
  }

  invisible(NULL)
}
