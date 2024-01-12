#' Generate omicsPLS input data
#'
#' Creates an object that can be used as input for the omicsPLS package.
#' It contains the omics datasets restricted to common samples. Each dataset
#' is feature-centred.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @param datasets Character vector of length 2, the names of the datasets from \code{mo_data} to include in the analysis.
#' @param scale_data Boolean, should the datasets be scaled? Default value is `FALSE`.
#' @return A list, in which each element corresponds to one omics dataset, with samples as rows and features as columns.
#' @export
get_input_omicspls <- function(mo_data, datasets = names(mo_data), scale_data = FALSE) {
  mo_data <- check_input_multidataset(mo_data, datasets)
  if (length(datasets) != 2) stop("In 'datasets' argument: expecting a vector of length 2 (as OmicsPLS integrates 2 datasets).")

  common_samples <- MultiDataSet::commonIds(mo_data)

  ds_list <- lapply(get_datasets(mo_data), function(x) {
    scale(t(x[, common_samples, drop = FALSE]), center = TRUE, scale = scale_data)
  })

  return(ds_list)
}
