% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output_utils.R
\name{.filter_output_dimensions_list}
\alias{.filter_output_dimensions_list}
\title{Filter latent dimensions in list}
\usage{
.filter_output_dimensions_list(
  output_list,
  latent_dimensions,
  all_present = FALSE,
  fixed_length = NULL
)
}
\arguments{
\item{output_list}{List of integration method outputs each generated via one of the
\code{get_output()} function.}

\item{latent_dimensions}{Named list, where each element is a character vector
giving the latent dimensions to retain in the corresponding element of \code{output_list}.
Names must match those of \code{output_list}.}

\item{all_present}{Logical, whether there should be one element in \code{latent_dimensions}
for each element of \code{output_list}. If \code{TRUE}, an error will be
returned if the length and names of \code{output_list} and \code{latent_dimensions} do
not match. Default value is \code{FALSE}.}

\item{fixed_length}{Integer, expected length of each element of \code{latent_dimensions}.
If \code{NULL} (default value), the length of elements in \code{latent_dimensions} can vary.}
}
\value{
A list of output similar to \code{output_list}, but the samples score table
or features weight table have been filtered.
}
\description{
Filters latent dimensions by name in a list of outputs from integration methods.
}
