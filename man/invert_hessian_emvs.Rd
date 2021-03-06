% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emvs_functions.R
\name{invert_hessian_emvs}
\alias{invert_hessian_emvs}
\title{Invert the log joint pseudoposterior Hessian matrix.}
\usage{
invert_hessian_emvs(
  x,
  sigma,
  gamma,
  index,
  prior_var_intercepts = 1,
  spike_var,
  slab_var
)
}
\arguments{
\item{x}{An \code{n} by \code{p} matrix containing binary coded
responses (i.e., coded \code{0,1}) for \code{n} independent observations
on \code{p} variables in the network or graph.}

\item{sigma}{A \code{p} by \code{p} numeric matrix with pairwise association
estimates on the off-diagonal elements and threshold estimates on the
diagonal elements.}

\item{gamma}{A \code{p} by \code{p} numeric matrix with expected inclusion
probabilities for the pairwise association parameters (i.e., edges).}

\item{index}{An indexing matrix used to convert the matrix \code{sigma} to a
vector \code{eta} during optimization. See \code{\link{indexing}}.}

\item{prior_var_intercepts}{The variance of the prior distribution
stipulated on the Ising model's threshold parameters. Currently a normal
distribution is used for all threshold parameters, with a mean equal to
zero and a variance equal to \code{prior_var_intercepts}. Defaults to
\code{1}.}

\item{spike_var, slab_var}{The \code{p} by \code{p} matrices of variances
that are used in the specification of the spike and slab prior
distributions that are stipulated on the association parameters.}
}
\value{
A \code{p(p-1)/2} by \code{p(p-1)/2} numeric matrix.
}
\description{
The function \code{invert_hessian_emvs} computes the inverse Hessian of the
log of the joint pseudoposterior distribution.
}
\details{
The function \code{invert_hessian_emvs} is used by
\code{multidimensional_update_emvs} for the multidimensional Newton Raphson
step.
}
