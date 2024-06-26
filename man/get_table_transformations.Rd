% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transformation.R
\name{get_table_transformations}
\alias{get_table_transformations}
\title{Get table with transformation applied to each dataset}
\usage{
get_table_transformations(
  transformation_result,
  best_normalize_details = FALSE
)
}
\arguments{
\item{transformation_result}{A list in which each element is the result of a
transformation applied to a different dataset, computed with the
\link{transform_dataset} function.}

\item{best_normalize_details}{Logical, should information about the
transformations selected by bestNormalize for each feature be displayed?
Default value is \code{FALSE}.}
}
\value{
A tibble with columns \code{'Dataset'} and \code{'Transformation'}. If
\code{best_normalize_details = TRUE}, an additional column \code{'Details'} lists the
chsoen transformation applied to each feature of the corresponding dataset
for a bestNormalize transformation.
}
\description{
From the results of transformations on datasets, generates a table giving for
each dataset the transformation that was applied to it.
}
