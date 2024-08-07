% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_plot.R
\name{ggpairs_custom}
\alias{ggpairs_custom}
\title{ggpairs plot with custom colours}
\usage{
ggpairs_custom(
  toplot,
  vars,
  colour_upper = NULL,
  colour_diag = colour_upper,
  colour_lower = colour_upper,
  shape_upper = NULL,
  shape_lower = shape_upper,
  scale_colour_upper = NULL,
  scale_colour_diag = NULL,
  scale_colour_lower = NULL,
  scale_shape_upper = NULL,
  scale_shape_lower = NULL,
  title = NULL,
  point_size = 1.5
)
}
\arguments{
\item{toplot}{Tibble in wide format, with observations as rows and variables
as columns.}

\item{vars}{Character vector, names of the columns in \code{toplot} that
correspond to variables to be used for the plot matrix.}

\item{colour_upper}{Character, name of column in \code{toplot} to use for
colouring observations in the upper triangle plots. Default value is
\code{NULL}.}

\item{colour_diag}{Character, name of column in \code{toplot} to use for colouring
observations in the diagonal plots. By default, will follow \code{colour_upper}.}

\item{colour_lower}{Character, name of column in \code{toplot} to use for
colouring observations in the lower triangle plots. By default, will follow
\code{colour_upper}.}

\item{shape_upper}{Character, name of column in \code{toplot} to use for shaping
observations in the upper triangle plots. Default value is \code{NULL}.}

\item{shape_lower}{Character, name of column in \code{toplot} to use for shaping
observations in the lower triangle plots. By default, will follow
\code{shape_upper}.}

\item{scale_colour_upper}{ggplot2 colour scale to use for the upper triangle
plots. Default value is \code{NULL} (if \code{colour_upper} is not \code{NULL}, will use
ggplot2 default colour scales).}

\item{scale_colour_diag}{ggplot2 colour scale to use for the diagonal plots.
If \code{NULL} (default), the colour scale used for the upper triangle plots
will be used if \code{colour_diag} is equal to \code{colour_upper}; or the colour
scale used for the lower triangle plots will be used if \code{colour_diag} is
equal to \code{colour_lower}.}

\item{scale_colour_lower}{ggplot2 colour scale to use for the lower triangle
plots. If \code{NULL} (default), the colour scale used for the upper triangle
plots will be used.}

\item{scale_shape_upper}{ggplot2 shape scale to use for the upper triangle
plots. Default value is \code{NULL} (if \code{shape_upper} is not \code{NULL}, will use
ggplot2 default shape scale).}

\item{scale_shape_lower}{ggplot2 shape scale to use for the lower triangle
plots. If \code{NULL} (default), the shape scale used for the upper triangle
plots will be used.}

\item{title}{Character, title of the plot. Default value is \code{NULL} (no title
added to the plot).}

\item{point_size}{Numeric, the size of points (in pt) in the plot. Default
value is 1.5.}
}
\value{
A \code{ggmatrix} plot.
}
\description{
Creates a ggpairs plot (see \code{\link[GGally:ggpairs]{GGally::ggpairs()}}) where the colours and shapes
can differ between the upper triangle, lower triangle and diagonal plots.
}
\examples{
\dontrun{
library(palmerpenguins)
library(ggplot2)
data("penguins")

vars <- c(
  "bill_length_mm",
  "bill_depth_mm",
  "flipper_length_mm"
)

toplot <- penguins |>
  dplyr::filter(!is.na(bill_length_mm))

# simple scatterplots of the variables
ggpairs_custom(toplot, vars)

# colouring points by species, using custom colour palette
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1")
)

# now adding the sex variable as shape of the points
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1"),
  shape_upper = "sex"
)

# using the lower plots to show the island as colour
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1"),
  shape_upper = "sex",
  colour_lower = "island",
  scale_colour_lower = scale_colour_viridis_d()
)

# showing species and sex in upper plots, body mass and island
# in lower plots
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1"),
  shape_upper = "sex",
  shape_lower = "island",
  colour_lower = "body_mass_g",
  scale_colour_lower = scale_colour_viridis_c(option = "plasma")
)

# same as above, but the diagonal plots show density per year
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1"),
  shape_upper = "sex",
  shape_lower = "island",
  colour_lower = "body_mass_g",
  scale_colour_lower = scale_colour_viridis_c(option = "plasma"),
  colour_diag = "year"
)

# common legend if the diagonal follows the colour of the lower plots
ggpairs_custom(
  toplot,
  vars,
  colour_upper = "species",
  scale_colour_upper = scale_colour_brewer(palette = "Set1"),
  colour_lower = "island",
  scale_colour_lower = scale_colour_brewer(palette = "Accent"),
  colour_diag = "island"
)
}
}
