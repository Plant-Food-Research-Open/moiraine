% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods_comparison.R
\name{get_features_weight_correlation}
\alias{get_features_weight_correlation}
\title{Get features weight correlation}
\usage{
get_features_weight_correlation(output_list, include_missing_features = FALSE)
}
\arguments{
\item{output_list}{List of integration methods output, each generated via the
\code{get_output()} function. If named, the names will be added at the beginning of
each latent dimension' label. If unnamed, the name of the integration method
will be used instead.}

\item{include_missing_features}{Logical, whether features missing in some of the
output should be included in the calculation (see Details). Default value is \code{FALSE}.}
}
\value{
A correlation matrix.
}
\description{
Constructs the correlation matrix between the features weight of the latent
dimensions obtained with different integration methods.
}
\details{
If \code{include_missing_features} is \code{FALSE} (default behaviour), and some
features are present in the output of one integration method but not
the other (e.g. because a different pre-filtering was applied to the input
data of the two methods), these features will be ignored. This does not mean
that features that were selected by one method but not the other are discarded;
in that case the feature will be assigned a weight of 0 for the method that did
not select it. This is the recommended behaviour, should only be changed in
specific scenarios (e.g. to check whether using all features in a dataset vs
doing a variance-based preselection affect which features are deemed most important).
If \code{include_missing_features} is \code{TRUE}, missing features will
be assigned a weight of 0.
}
