#' Number of features in each dataset of MultiDataSet object
#'
#' Gives the number of features in each dataset of a MultiDataSet object.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @return A named vector, where each element is the number of features in the
#' corresponding dataset.
#' @export
n_features <- function(mo_data) {
  return(vapply(mo_data@featureData, function(y) nrow(y[[1]]), numeric(1)))
}

#' Number of samples in each dataset of MultiDataSet object
#'
#' Gives the number of samples in each dataset of a MultiDataSet object.
#'
#' @param mo_data A \code{\link[MultiDataSet]{MultiDataSet-class}} object.
#' @return A named vector, where each element is the number of sample in the
#' corresponding dataset.
#' @export
n_samples <- function(mo_data) {
  return(vapply(mo_data@phenoData, function(y) nrow(y[[1]]), numeric(1)))
}
