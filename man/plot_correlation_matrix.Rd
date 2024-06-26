% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_plot.R
\name{plot_correlation_matrix}
\alias{plot_correlation_matrix}
\title{Plot correlation matrix}
\usage{
plot_correlation_matrix(cormat, ...)
}
\arguments{
\item{cormat}{A correlation matrix.}

\item{...}{Other arguments to be passed to the \code{\link[corrplot:corrplot]{corrplot::corrplot()}}
function.}
}
\value{
A correlation plot
}
\description{
Plots a correlation matrix using the corrplot package.
}
