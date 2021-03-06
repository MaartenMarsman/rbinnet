% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maximum_pseudolikelihood_functions.R
\name{multidimensional_update}
\alias{multidimensional_update}
\title{Multidimensional Newton Raphson pseudoposterior updates.}
\usage{
multidimensional_update(x, sigma, index, suff_stat, prior_var = 1)
}
\arguments{
\item{x}{An \code{n} by \code{p} matrix containing binary coded
responses (i.e., coded \code{0,1}) for \code{n} independent observations
on \code{p} variables in the network or graph.}

\item{sigma}{A \code{p} by \code{p} numeric matrix with pairwise association
estimates on the off-diagonal elements and threshold estimates on the
diagonal elements.}

\item{index}{An indexing matrix used to convert the matrix \code{sigma} to a
vector \code{eta} during optimization. See \code{\link{indexing}}.}

\item{suff_stat}{A \code{p} by \code{p} matrix of sufficient statistics.}

\item{prior_var}{The variance of the prior distribution stipulated on the
Ising model parameters. Currently a normal distribution is used for all
model parameters, with a mean equal to zero and a variance equal to
\code{prior_var}. Defaults to \code{1}.}
}
\value{
A \code{p} by \code{p} numeric matrix with updates pairwise
  association estimates on the off-diagonal elements and updated threshold
  estimates on the diagonal elements.
}
\description{
The function \code{multidimensional_update} updates all pseudoposterior
values all at once.
}
\details{
The function \code{multidimensional_update} is used by
\code{optimize_pseudoposterior} to optimize the pseudoposterior parameters.
}
