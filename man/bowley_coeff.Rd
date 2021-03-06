% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic.R
\name{bowley_coeff}
\alias{bowley_coeff}
\title{Bowley's coefficient}
\usage{
bowley_coeff(x)
}
\arguments{
\item{x}{a numerical vector}

\item{na.rm}{a logical value indicating whether NA values should be stripped
before the computation proceeds.}
}
\description{
\code{bowley_coeff} computes the Bowley's Coefficient.
}
\details{
Bowley's coefficient assess the symmetry of the data.
}
\examples{
bowley_coeff(runif(1000))
}
