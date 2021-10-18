#' One dimensional Newton Raphson pseudoposterior updates.
#'
#' The function \code{onedimensional_update} updates all pseudoposterior values
#' one at a time.
#'
#' The function \code{onedimensional_update} is used by
#' \code{\link{optimize_pseudoposterior}} to provide reasonable starting values. After
#' which it runs \code{\link{multidimensional_update}}.
#'
#' @inheritParams compute_log_pseudoposterior
#'
#' @param suff_stat A \code{p} by \code{p} matrix of sufficient statistics.
#'
#' @return A \code{p} by \code{p} numeric matrix with updates pairwise
#'   association estimates on the off-diagonal elements and updated threshold
#'   estimates on the diagonal elements.
onedimensional_update <- function(x, sigma, suff_stat, prior_var = 1) {
  p <- ncol(x)

  # update main effects -------------------------------------------------------
  for(s in 1:p) {
    tmp <- sigma[s, s] + x[, -s] %*% sigma[s, -s]
    prob <- expit(tmp)

    first_derivative <- suff_stat[s, s] - sum(prob)
    first_derivative <- first_derivative - sigma[s, s] / prior_var

    second_derivative <- -sum(prob * (1 - prob)) - 1 / prior_var

    sigma[s, s] <- sigma[s, s] - first_derivative / second_derivative
  }

  # update interaction effects ------------------------------------------------
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      tmp <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
      prob_s <- expit(tmp)

      tmp <- sigma[t, t] + x[, -t] %*% sigma[-t, t]
      prob_t <- expit(tmp)

      first_derivative <- 2 * suff_stat[s, t] -  x[, t] %*% prob_s
      first_derivative <- first_derivative - x[, s] %*% prob_t
      first_derivative <- first_derivative - sigma[s, t] / prior_var

      tmp <- prob_s * (1 - prob_s)
      tmp2 <- prob_t * (1 - prob_t)
      second_dervative <- -x[, t] %*% tmp - x[, s] %*% tmp2
      second_dervative <- second_derivative - 1 / prior_var

      sigma[s, t] <- sigma[s, t] - first_derivative / second_dervative
      sigma[t, s] <- sigma[s, t]
    }
  }

  return(sigma)
}

#' Multidimensional Newton Raphson pseudoposterior updates.
#'
#' The function \code{multidimensional_update} updates all pseudoposterior
#' values all at once.
#'
#' The function \code{multidimensional_update} is used by
#' \code{optimize_pseudoposterior} to optimize the pseudoposterior parameters.
#'
#' @inheritParams onedimensional_update
#'
#' @param index An indexing matrix used to convert the matrix \code{sigma} to a
#' vector \code{eta} during optimization. See \code{\link{indexing}}.
#'
#' @return A \code{p} by \code{p} numeric matrix with updates pairwise
#'   association estimates on the off-diagonal elements and updated threshold
#'   estimates on the diagonal elements.
multidimensional_update <- function(x, sigma, index, suff_stat, prior_var = 1) {
  n <- nrow(x)
  p <- ncol(x)

  # compute vector of first-order derivatives (main effects) ------------------
  prob <- matrix(0, nrow = n, ncol = p)
  derivatives <- vector(length = p * (p + 1) / 2)
  for(s in 1:p) {
    phi <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
    prob[, s] <- expit(phi)
    derivatives[s] <- suff_stat[s, s] - sum(prob[, s])
    derivatives[s] <- derivatives[s] - sigma[s, s] / prior_var
  }

  # compute vector of first-order derivatives (interaction effects) -----------
  for(s in 1:(p - 1)) {
    for(t in (s + 1): p) {
      row <- p + index[index[, 2] == s & index[, 3] == t, 1]
      derivatives[row] <- 2 * suff_stat[s, t] - x[, t] %*% prob[, s]
      derivatives[row] <- derivatives[row] - x[, s] %*% prob[, t]
      derivatives[row] <- derivatives[row] - sigma[s, t] / prior_var
    }
  }

  # compute the inverse Hessian -----------------------------------------------
  inv_hessian <- invert_hessian(sigma = sigma,
                                index = index,
                                x = x,
                                prior_var = prior_var)

  # convert to eta values -----------------------------------------------------
  eta <- vector(length = p * (p + 1) / 2)
  eta[1:p] <- diag(sigma)
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row <- index[index[,2] == s & index[, 3] == t, 1] + p
      eta[row] <- sigma[s, t]
    }
  }

  # Newton - Raphson update ---------------------------------------------------
  eta <- eta - inv_hessian %*% derivatives

  # revert to sigma values ----------------------------------------------------
  diag(sigma) <- eta[1:p]
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      row <- p + index[index[,2] == s & index[, 3] == t, 1]
      sigma[s, t] <- eta[row]
      sigma[t, s] <- eta[row]
    }
  }
  return(sigma)
}

#' Invert the log pseudoposterior Hessian matrix.
#'
#' The function \code{invert_hessian} computes the inverse Hessian of the log
#' pseudoposterior distribution.
#'
#' The function \code{invert_hessian} is used by
#' \code{multidimensional_update} for the multidimensional Newton Raphson step,
#' and by \code{\link{compute_standard_deviation}} to compute the
#' pseudoposterior standard deviations.
#'
#' @inheritParams multidimensional_update
#'
#' @return A \code{p(p-1)/2} by \code{p(p-1)/2} numeric matrix.
invert_hessian <- function(x, sigma, index, prior_var = 1) {
  # the Hessian matrix is build up as
  #           (main_hessian    , cross_hessian)
  #           (cross_hessian^T , int_hessian)
  # main_hessian contains partial derivatives w.r.t the main effects.
  # int_hessian contains partial derivatives w.r.t the interaction effects.

  p <- ncol(x)
  n <- nrow(x)

  # precomputation ------------------------------------------------------------
  pq <- matrix(0, nrow = n, ncol = p)
  for(s in 1:p) {
    phi <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
    pq[, s] <- nexpit(phi)
  }

  # compute int_hessian -------------------------------------------------------
  int_hessian <- matrix(0, nrow = p * (p - 1) / 2, ncol = p * (p-1) / 2)
  for(row in 1:(p * (p - 1) / 2)) {
    s <- index[row, 2]
    t <- index[row, 3]

    # partial derivative w.r.t. sigma[s, t]; s < t ----------------------------
    int_hessian[row, row] <- -x[, t] %*% pq[, s] - x[, s] %*% pq[, t]
    int_hessian[row, row] <- int_hessian[row, row] - 1 / prior_var

    # partial derivative w.r.t. sigma[s, t] & sigma [s, r]; s < t -------------
    I <- which(index[, 2] == s & index[, 3] != t)
    for(i in I) {
      r <- index[i, 3]
      column <- index[i, 1]
      int_hessian[row, column] <- -(x[, t] * x[, r]) %*% pq[, s]
    }

    # partial derivative w.r.t. sigma[s, t] & sigma [r, s]; s < t -------------
    I <- which(index[, 3] == s)
    for(i in I) {
      r <- index[i, 2]
      column <- index[i, 1]
      int_hessian[row, column] <- -(x[, t] * x[, r]) %*% pq[, s]
    }

    # partial derivative w.r.t. sigma[s, t] & sigma [k, t]; s < t -------------
    I <- which(index[, 3] == t & index[, 2] != s)
    for(i in I) {
      k <- index[i, 2]
      column <- index[i, 1]
      int_hessian[row, column] <- -(x[, s] * x[, k]) %*% pq[, t]
    }

    # partial derivative w.r.t. sigma[s, t] & sigma [t, k]; s < t -------------
    I <- which(index[, 2] == t)
    for(i in I) {
      k <- index[i, 3]
      column <- index[i, 1]
      int_hessian[row, column] <- -(x[, s] * x[, k]) %*% pq[, t]
    }
  }

  # compute cross_hessian -----------------------------------------------------
  cross_hessian <- matrix(0, nrow = p, ncol = p * (p - 1) / 2)
  for(s in 1:p) {
    I <- which(index[, 2] == s)
    for(i in I) {
      t <- index[i, 3]
      column <- index[i, 1]
      cross_hessian[s, column] <- - x[, t] %*% pq[, s]
    }

    I <- which(index[, 3] == s)
    for(i in I) {
      t <- index[i, 2]
      column <- index[i, 1]
      cross_hessian[s, column] <- - x[, t] %*% pq[, s]
    }
  }

  # compute inv_hessian -------------------------------------------------------
  inv_hessian <- matrix(data = 0,
                        nrow = p * (p + 1) / 2,
                        ncol = p * (p + 1) / 2)

  # indices for main effects and interactions
  index_main <- 1:p
  index_int <- (1:(p * (p + 1) / 2))[-index_main]

  # inverse of main for Schur complement --------------------------------------
  inv_main <- matrix(0, nrow = p, ncol = p)
  for(s in 1:p)
    inv_main[s, s] <- 1 / (-sum(pq[, s]) - 1 / prior_var)

  # use Schur complement for inverting ----------------------------------------
  inv_hessian[index_int, index_int] <- int_hessian -
    t(cross_hessian) %*% inv_main %*% cross_hessian
  inv_hessian[index_int, index_int] <- solve(inv_hessian[index_int, index_int])
  inv_hessian[index_main, index_int] <-
    -inv_main %*% cross_hessian %*% inv_hessian[index_int, index_int]
  inv_hessian[index_int, index_main] <- t(inv_hessian[index_main, index_int])
  inv_hessian[index_main, index_main] <- inv_main -
    inv_hessian[index_main, index_int] %*% t(cross_hessian) %*% inv_main

  return(inv_hessian)
}

#' Compute the proportional log pseudoposterior density value.
#'
#' The function \code{\link{compute_log_pseudoposterior}} computes the proportional
#' log pseudoposterior density if \code{prior_var} is finite, and computes the
#' log pseudolikelihood if \code{prior_var} is set to \code{Inf}.
#'
#' The pseudolikelihood is the product of all full-conditional distributions of
#' the full Ising model. Used by \code{\link{optimize_pseudoposterior}} to check for
#' convergence of the Newton Raphson procedure.
#'
#' @param x An \code{n} by \code{p} matrix containing binary coded
#'   responses (i.e., coded \code{0,1}) for \code{n} independent observations
#'   on \code{p} variables in the network or graph.
#'
#' @param prior_var The variance of the prior distribution stipulated on the
#'   Ising model parameters. Currently a normal distribution is used for all
#'   model parameters, with a mean equal to zero and a variance equal to
#'   \code{prior_var}. Defaults to \code{1}.
#'
#' @param sigma A \code{p} by \code{p} numeric matrix with pairwise association
#'   estimates on the off-diagonal elements and threshold estimates on the
#'   diagonal elements.
#'
#' @return A single numeric value.
compute_log_pseudoposterior <- function(x, sigma, prior_var = 1) {
  # returns log of pseudolikelihood if prior_var = Inf
  # returns log of (proportional) pseudoposterior otherwise

  p <- ncol(x)

  log_pseudoposterior <- 0
  for(s in 1:p) {
    phi <- sigma[s, s] + x[, -s] %*% sigma[-s, s]
    log_pseudoposterior <- log_pseudoposterior + x[, s] %*% phi
    log_pseudoposterior <- log_pseudoposterior - sum(log(1 + exp(phi)))
    if(is.finite(prior_var)) {
      # zero-mean normal prior distribution for the main effects --------------
      log_pseudoposterior <- log_pseudoposterior +
        dnorm(x = sigma[s, s], sd = sqrt(prior_var), log = TRUE)
    }
  }
  if(is.finite(prior_var)) {
    # zero-mean normal prior distribution for the interaction effects ---------
    for(s in 1:(p - 1)) {
      for(t in (s + 1):p) {
        log_pseudoposterior <- log_pseudoposterior +
          dnorm(x = sigma[s, t], sd = sqrt(prior_var), log = TRUE)
      }
    }
  }

  return(log_pseudoposterior)
}

#' Compute the asymptotic pseudoposterior standard deviations.
#'
#' The function \code{compute_standard_deviation} computes the pseudoposterior
#' standard deviations based on an asymptotic normal approximation. If
#' \code{prior_var} is set to \code{Inf}, it returns the standard errors of the
#' pseudolikelihood parameters instead.
#'
#' The function \code{compute_standard_deviation} calls
#' \code{\link{invert_hessian}} to produce the standard deviations
#' (standard errors).
#'
#' @inheritParams compute_log_pseudoposterior
#'
#' @return A \code{p} by \code{p} numeric matrix with standard deviations for
#'   the pairwise association parameters on the off-diagonal elements and for
#'   the threshold paramaters on the diagonal elements.
compute_standard_deviation <- function(x, sigma, prior_var = 1) {
  p <- ncol(x)

  index <- indexing(p)
  inv_hessian <- invert_hessian(sigma = sigma,
                                index = index,
                                x = x,
                                prior_var = prior_var)

  standard_deviation <- matrix(0, nrow = p, ncol = p)
  for(s in 1:(p - 1)) {
    for(t in (s + 1):p) {
      tmp <- diag(inv_hessian)[index[index[, 2] == s & index[, 3] == t, 1] + p]
      standard_deviation[s, t] <- sqrt(-tmp)
      standard_deviation[t, s] <- sqrt(-tmp)
    }
  }
  diag(standard_deviation) <- sqrt(-diag(inv_hessian)[1:p])

  return(standard_deviation)
}

#' Optimizing the Ising model's pseudoposterior distribution
#'
#' The function \code{fit_pseudoposterior} can be used to optimize the Ising
#' model's pseudoposterior distribution. It provides the posterior modes and
#' standard deviations of both the thresholds and pairwise associations. The
#' function produces maximum pseudolikelihood estimates and standard errors
#' whenever \code{prior_var} is set to \code{Inf}.
#'
#' Optimization is performed with respect to all \code{p} full-conditional
#' distributions at once. Standard deviations (standard errors) are based on
#' assumption of asymptotic normality. The function \code{fit_pseudoposterior}
#' calls \code{\link{optimize_pseudoposterior}} and
#' \code{\link{compute_standard_deviation}}.
#'
#' @param x An \code{n} by \code{p} matrix containing binary coded
#'   responses (i.e., coded \code{0,1}) for \code{n} independent observations
#'   on \code{p} variables in the network or graph.
#'
#' @param prior_var The variance of the prior distribution stipulated on the
#'   Ising model parameters. Currently a normal distribution is used for all
#'   model parameters, with a mean equal to zero and a variance equal to
#'   \code{prior_var}. Defaults to \code{1}.
#'
#' @param iteration_max The maximum number of Newton-Raphson iterations used in
#'   optimization. Defaults to \code{1e2}. A warning is issued if procedure has
#'   not converged in \code{iteration_max} iterations.
#'
#' @return A list containing \code{mu} and \code{sigma}, a numeric vector of
#'   length \code{p} containing the threshold estimates and a \code{p} by
#'   \code{p} numeric matrix of pairwise association estimates. Also included
#'   are \code{s.d.mu} and \code{s.d.sigma}, which contain the corresponding
#'   posterior standard deviations.
#'
#' @examples
#' library("IsingSampler")
#' ### Simulate dataset ###
#' # Input:
#' p <- 6 # Number of nodes
#' n <- 1000 # Number of samples
#' # Ising parameters:
#' Graph <- matrix(data = 0, nrow = p, ncol = p)
#' Graph[lower.tri(Graph)] <- runif(n = p * (p - 1) / 2, min = -.5, max = 1)
#' Graph <- Graph + t(Graph)
#' Thresholds <- -rowSums(Graph) / 2
#' # Simulate:
#' Data <- IsingSampler(n = n, graph = Graph, thresholds = Thresholds)
#' ### Fit using fit_pseudoposterior ###
#' Fit_ml <- fit_pseudoposterior(x = Data, prior_var = Inf)
#' Fit_bayes <- fit_pseudoposterior(x = Data, prior_var = 1)
#' # Plot results:
#' library("qgraph")
#' layout(t(1:3))
#' qgraph(Fit_ml$sigma,fade = FALSE)
#' title("ML estimated network")
#' qgraph(Fit_bayes$sigma,fade = FALSE)
#' title("Bayes estimated network")
#' qgraph(Graph,fade = FALSE)
#' title("Original network")
fit_pseudoposterior <- function(x, prior_var = 1, iteration_max = 1e2) {

  if(!exists("x"))
    stop("Data matrix should be provided in x.", call. = FALSE)

  if(length(prior_var) != 1)
    stop("A single prior variance should be specified", call. = FALSE)

  if(prior_var <= 0)
    stop("Prior variance should be positive number.", call. = FALSE)

  estimates <- try(optimize_pseudoposterior(x = x,
                                            prior_var = prior_var,
                                            iteration_max = iteration_max)$sigma, silent = TRUE)
  if(class(estimates) != "try-error") {
    standard_deviation <- try(compute_standard_deviation(x = x,
                                                         sigma = estimates,
                                                         prior_var = prior_var), silent = TRUE)
    if(exists("standard_deviation")) {
      sd.mu <- diag(standard_deviation)
      sd.sigma <- standard_deviation
      diag(sd.sigma) <- 0
    } else {
      warning("Asymptotic standard deviations could not be computed.",
              call. = FALSE)
    }
  } else {
    stop("Estimation failed for unknown reasons. Try ``optimize_pseudoposterior'' directly or contact the author.",
         call. = FALSE)
  }
  mu <- diag(estimates)
  diag(estimates) <- 0

  fit <- list(mu = mu, sigma = estimates, sd.mu = sd.mu, sd.sigma = sd.sigma)
  return(fit)
}

#' Optimization of the Ising model's pseudoposterior distribution
#'
#' The function \code{optimize_pseudoposterior} can be used to optimize the
#' Ising model's pseudoposterior distribution. It provides the posterior modes
#' and of both the thresholds and pairwise associations. The function produces
#' maximum pseudolikelihood estimates whenever \code{prior_var} is set to
#' \code{Inf}.
#'
#' Optimization is performed with respect to all \code{p} full-conditional
#' distributions at once.
#'
#' @inheritParams fit_pseudoposterior
#'
#' @return A \code{p} by \code{p} numeric matrix with pairwise association
#'   estimates on the off-diagonal elements and threshold estimates on the
#'   diagonal elements.
#'
#' @seealso \code{\link{fit_pseudoposterior}}
optimize_pseudoposterior <- function(x, prior_var = 1, iteration_max = 1e2) {
  p <- ncol(x)
  n <- nrow(x)

  # parameter matrix and sufficient statistics --------------------------------
  sigma <- matrix(0, nrow = p, ncol = p)
  suff_stat <- t(x) %*% x

  # compute starting values ---------------------------------------------------
  for(iteration in 1:5) {
    sigma <- onedimensional_update(sigma = sigma,
                                   x = x,
                                   suff_stat = suff_stat,
                                   prior_var = prior_var)
  }

  log_pseudoposterior <- compute_log_pseudoposterior(sigma = sigma,
                                                     x = x,
                                                     prior_var = prior_var)

  # index allows us to convert sigma (matrix) to eta (vector) and revert back -
  index <- indexing (p)

  for(iteration in 1:iteration_max) {
    # update pseudolikelihood parameters --------------------------------------
    sigma <- multidimensional_update(sigma = sigma,
                                     index = index,
                                     x = x,
                                     suff_stat = suff_stat,
                                     prior_var = prior_var)

    # update log of pseudoposterior -------------------------------------------
    log_pseudoposterior_old <- log_pseudoposterior
    log_pseudoposterior <- compute_log_pseudoposterior(sigma = sigma,
                                                       x = x,
                                                       prior_var = prior_var)
    difference <- abs(log_pseudoposterior - log_pseudoposterior_old)
    if(difference < sqrt(.Machine$double.eps)) break

    if(iteration == iteration_max)
      warning(paste("The optimization procedure did not convergence in", iteration_max, "iterations.",
                 sep = " "), call. = FALSE)
  }


  return(list(sigma = sigma))
}
