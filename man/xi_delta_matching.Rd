% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{xi_delta_matching}
\alias{xi_delta_matching}
\title{Matching the spike and slab intersection.}
\usage{
xi_delta_matching(xi, delta, n)
}
\arguments{
\item{xi}{A positive penalty parameter.}

\item{delta}{A positive number. Taken to be the inverse of the intended
precision, e.g., \code{delta = qnorm(.975, lower.tail = TRUE)}.}

\item{n}{A positive integer. Number of observations.}
}
\value{
A numerical value.
}
\description{
The function \code{xi_delta_matching} is used to match the intersection
of the spike and slab prior distributions to a predefined number
\code{delta}, as a function of a penalty parameter \code{xi}.
}
\details{
The root of \code{xi_delta_matching} is searched by
\code{\link{set_spike_and_slab}}. The root value \code{xi} for a given
\code{delta} and \code{n} is used to set the spike variance.
}
