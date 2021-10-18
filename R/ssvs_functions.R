#' Gibbs structure selection procedure.
#'
#' The function \code{structure_selection_ssvs} explores the space of possible
#' Ising network structures using the Gibbs sampler.
#'
#' The function \code{structure_selection_ssvs} is used by
#' \code{select_structure} for exploring the space of network structures. It
#' makes use of a Bayesian spike and slab approach to edge inclusion, and
#' explores the posterior space of network structures using the Gibbs sampler.
#'
#' @param x An \code{n} by \code{p} matrix containing binary coded
#'   responses (i.e., coded \code{0,1}) for \code{n} independent observations
#'   on \code{p} variables in the network or graph.
#'
#' @param number_iterations Integer. The number of iterations the Gibbs sampler
#'   uses to explore the structure space. Defaults to \code{1e5}.
#'
#' @param number_burnin_iterations Integer. The number of burnin iterations for
#'   the Gibbs sampler. The output for the \code{number_burnin_iterations} is
#'   not recorded. After \code{number_burnin_iterations}, the Gibbs sampler
#'   runs for \code{number_iterations} iterations. Defaults to \code{0}.
#'
#' @param spike_var,slab_var The \code{p} by \code{p} matrices of variances
#'   that are used in the specification of the spike and slab prior
#'   distributions that are stipulated on the association parameters.
#'
#' @param prior_var_intercepts The variance of the prior distribution
#'   stipulated on the Ising model's threshold parameters. Currently a normal
#'   distribution is used for all threshold parameters, with a mean equal to
#'   zero and a variance equal to \code{prior_var_intercepts}. Defaults to
#'   \code{1}.
#'
#' @param theta The prior inclusion probability. The value \code{theta = 0.5},
#'   in combination with \code{hierarchical = FALSE} stipulates a uniform prior
#'   on the space of network structures. Defaults to \code{0.5}.
#'
#' @param hierarchical Logical. If TRUE, a beta prior distribution is
#'  imposed on the prior inclusion probability \code{theta} with
#'  hyperparameters \code{alpha} and \code{beta}. A uniform prior on the
#'  inclusion probability, a beta with \code{alpha = beta = 1}, stipulates a
#'  uniform prior on network structure complexity.
#'
#' @param alpha,beta The hyperparameters of the beta prior distribution
#'   stipulated on the prior inclusion probability \code{theta} if
#'   \code{hierarchical = TRUE}. Default to \code{1}.
#'
#' @param output_samples Logical. If \code{output_samples = TRUE},
#'   \code{structure_selection_ssvs} returns all posterior samples of the Ising
#'   model parameters. If \code{output_samples = FALSE},
#'   \code{structure_selection_ssvs} returns the posterior means instead.
#'
#' @param sigma A \code{p} by \code{p} numeric matrix with pairwise association
#'   estimates on the off-diagonal elements and threshold estimates on the
#'   diagonal elements. Optional. Can be used to specify starting values.
#'
#' @param include A \code{p} by \code{p} binary matrix. If
#'   \code{include[i,j] = 1}, the corresponding inclusion variable
#'   \code{gamma[i, j]} is sampled. Otherwise, if \code{include[i,j] = 0}, the
#'   inclusion variable \code{gamma[i,j]} and the association \code{sigma[i,j]}
#'   are fixed to zero, which means that the edge between nodes \code{i}
#'   and \code{j} is excluded from the structures that are considered.
#'   Used to screen promising edges with \code{\link{screen_edges}}.
#'
#' @param components One of "normal", "t", or "laplace". The type of
#'   distribution that is used in the spike and slab set-up. The variance
#'   of the slab component is matched to \code{slab_var}. Defaults to
#'   \code{components = "normal"}.
#'
#' @param df Used with \code{components = "t"}. The prior degrees of
#'   freedom used in the T distribution. Defaults to \code{df = 5}.
#'
#' @return A list that contains \code{structures}, an \code{S} by
#'   \code{p(p-1)/2} matrix of binary variables that encodes the \code{S}
#'   structures that were visited by the Gibbs sampler. Also included is
#'   \code{posterior_probability},
#'   the \code{S} by \code{1} matrix of posterior probabilities that matches
#'   the rows of \code{structures}. If \code{output_samples = TRUE} the list
#'   contains an \code{number_iterations} by \code{p(p + 1) / 2} matrix
#'   \code{sigma_samples} that contains samples from the posterior distribution
#'   of the Ising model parameters, and if \code{hierarchical = TRUE}, it also
#'   includes an \code{number_iterations} by \code{1} matrix of
#'   \code{theta_samples} which contains samples from the posterior
#'   distribution of the prior inclusion probability \code{theta}. Otherwise,
#'   if \code{output_samples = FALSE}, the list contains a \code{p(p + 1) / 2}
#'   by \code{1} matrix \code{sigma_eap} of posterior means of the Ising model
#'   parameters, and if \code{hierarchical = TRUE}, it also includes
#'   \code{theta_eap}, the posterior mean of the prior inclusion probability.
structure_selection_ssvs <- function(x, number_iterations = 1e5,
                 number_burnin_iterations = 0,
                 spike_var, slab_var, prior_var_intercepts,
                 hierarchical = FALSE, alpha = 1, beta = 1, theta = 0.5,
                 output_samples = FALSE, sigma, include,
                 components, df) {

  if (!requireNamespace("pgdraw", quietly = TRUE)) {
    stop("Package `pgdraw' is needed for simulating Polya-Gamma variates. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop("Package `statmod' is needed for simulating inverse Gaussian variates. Please install it.", call. = FALSE)
  }
  p <- ncol(x)
  n <- nrow(x)
  sufC <- colSums(x)

  if(hasArg("include")) {
    number_edges <- sum(include[lower.tri(include)])
  } else {
    number_edges <- choose(p, 2)
  }

  # prior specification -------------------------------------------------------
  if(components != "normal") {
    r <- spike_var[1, 2] / slab_var[1, 2]
  }

  # starting values -----------------------------------------------------------
  if(!hasArg("sigma")) {
    sigma <- matrix(0, nrow = p, ncol = p)
  }
  if(hasArg("include")) {
    sigma[include == 0] <- 0
    gamma <- include
  } else {
    include <- gamma <- matrix(1, nrow = p, ncol = p)
    diag(gamma) <- 0
  }
  omega <- matrix(0, nrow = n, ncol = p)

  # parameter output (optional) -----------------------------------------------
  edge_names <- par_names <- matrix(0, nrow = p, ncol = p)
  for(s in 1:(p-1)) {
    for(t in (s+1):p) {
      edge_names [t, s] <- paste("(", s, ",", t, ")", sep ="")
      par_names [t, s] <- paste("sigma(", s, ",", t, ")", sep ="")
    }
  }
  diag(par_names) <-  paste("mu(", 1:p, ")", sep ="")
  if(output_samples == TRUE) {
    sigma_samples <- matrix(0, nrow = number_iterations, ncol = p + choose(p, 2))
   colnames(sigma_samples) <- par_names[lower.tri(par_names, diag = TRUE)]
    if(hierarchical == TRUE) {
      theta_samples <- matrix(0, nrow = number_iterations, ncol = 1)
    }
  } else {
    sigma_mean <- sigma[lower.tri(sigma, diag = TRUE)]
    theta_mean <- theta
  }

  for(iteration in (-number_burnin_iterations):number_iterations) {
    for(s in 1:p) {
      # sample polya-gamma variates -------------------------------------------
      chi <- x[, -s] %*% sigma[-s, s]
      omega[, s] <- pgdraw::pgdraw(b = 1, c = sigma[s, s] + chi)

      # sample main effects ---------------------------------------------------
      post_var <- 1 + sum(omega[, s]) * prior_var_intercepts
      post_var <- prior_var_intercepts / post_var

      post_mean <- sufC[s] - n / 2  - omega[, s] %*% chi
      post_mean <- post_var * post_mean

      sigma[s, s] <- rnorm(n = 1, mean = post_mean, sd = sqrt(post_var))
    }

    for(s in 1:(p - 1)) {
      for(t in (s + 1):p) {
        # if edge is not screened ---------------------------------------------
        if(include[s, t] == 1) {
          # sample inclusion variable -----------------------------------------
          tmp1 <- theta / (1 - theta)
          v1 <- slab_var[s, t]
          v0 <- spike_var[s, t]
          tmp2 <- (1/v1 - 1/v0) * sigma[s, t] ^ 2 / 2
          probability <- tmp1 / (tmp1 + sqrt(v1 / v0) * exp(tmp2))

          gamma[s, t] <- rbinom(n = 1, size = 1, prob = probability)
          gamma[t, s] <- gamma[s, t]

          # sample interaction parameter --------------------------------------
          tmp1 <- slab_var[s, t] * gamma[s, t]
          tmp1 <- tmp1 + spike_var[s, t] * (1 - gamma[s, t])
          tmp2 <- omega[, s] %*% x[, t] + omega[, t] %*% x[, s]
          post_var <- tmp1 / (1 + tmp1 * tmp2)

          post_mean <- 2 * x[, s] %*% x[, t] - sufC[s] / 2 - sufC[t] / 2
          tmp1 <- sigma[s, s] + x[, -c(s, t)] %*% sigma[-c(s, t), s]
          tmp2 <- sigma[t, t] + x[, -c(t, s)] %*% sigma[-c(t, s), t]
          post_mean <- post_mean - sum(omega[, s] * x[, t] * tmp1)
          post_mean <- post_mean - sum(omega[, t] * x[, s] * tmp2)
          post_mean <- post_var * post_mean

          sigma[s, t] <- rnorm(n = 1,
                               mean = post_mean,
                               sd = sqrt(post_var))
          sigma[t, s] <- sigma[s, t]

          #sample hierarchical variance ---------------------------------------
          if(components != "normal") {
            if(components == "t") {
              #slab_var[t, s] is fixed at n * Var(sigma[s, t]) ----------------
              a <- df
              b <- slab_var[t, s] * (a - 1)
              rr <- gamma[s, t] + (1 - gamma[s, t]) * r
              phi <- 1 / rgamma(n = 1,
                                shape = a + 0.5,
                                rate = b + sigma[s, t] ^ 2 / (2 * rr))
              #update subdiagonal elements of spike / slab variance matrices --
              spike_var[s, t] <- r * phi
              slab_var[s, t] <- phi
            }
            if(components == "laplace") {
              #slab_var[t, s] is fixed at n * Var(sigma[s, t]) ----------------
              lambda <- sqrt(slab_var[t, s] / 2)
              rr <- gamma[s, t] + (1 - gamma[s, t]) * r

              phi <- 1 / statmod::rinvgauss(n = 1,
                                            mean = sqrt(rr) /
                                              (lambda * abs(sigma[s, t])),
                                            shape = 1 / lambda ^ 2)
              #update subdiagonal elements of spike and slab variance matrices --
              spike_var[s, t] <- r * phi
              slab_var[s, t] <- phi
            }
          }
        }
      }
    }

    # sample prior inclusion probability --------------------------------------
    if(hierarchical == TRUE) {
      theta <- rbeta(n = 1,
                     shape1 = alpha + sum(gamma) / 2,
                     shape2 = beta + number_edges - sum(gamma) / 2)
    }

    # collect output after burnin ---------------------------------------------
    if(iteration > 0) {
      if(output_samples == TRUE) {
        sigma_samples[iteration, ] <- sigma[lower.tri(sigma, diag = TRUE)]
        if(hierarchical == TRUE) {
          theta_samples[iteration] <- theta
        }
      } else  {
        sigma_mean <- (iteration - 1) * sigma_mean / iteration +
          sigma[lower.tri(sigma, diag = TRUE)] / iteration

        theta_mean <- (iteration - 1) * theta_mean / iteration +
          theta / iteration
      }

      if(iteration == 1) {
        #if inclusion[s, t] == 0 then gamma[s, t] == 0 ------------------------
        structures <- matrix(data = 2 * gamma[lower.tri(gamma)] - 1,
                             nrow = 1, ncol = choose(p, 2))

        posterior_probability <- 1
      } else {
        tmp_structure <- matrix(data = 2 * gamma[lower.tri(gamma)] - 1,
                                nrow = 1, ncol = choose(p, 2))
        matching <- structures %*% t(tmp_structure)
        if(sum(matching == choose(p, 2)) == 1) {
          match <- which.max(matching)[1]
          posterior_probability[match] <- posterior_probability[match] + 1
        } else {
          structures <- rbind(structures, tmp_structure)
          posterior_probability <- c(posterior_probability, 1)
        }
      }
    }
  }

  # output --------------------------------------------------------------------
  posterior_probability <- posterior_probability / sum(posterior_probability)
  structures <- structures / 2 + .5
  colnames(structures) <- edge_names[lower.tri(edge_names)]
  if(output_samples == TRUE) {
    if(hierarchical == TRUE) {
      output <- list(structures = structures,
                     posterior_probability = posterior_probability,
                     sigma_samples = sigma_samples, theta = theta_samples)
    } else {
      output <- list(structures = structures,
                     posterior_probability = posterior_probability,
                     sigma_samples = sigma_samples)
    }
  } else {
    sigma_mean <- matrix(sigma_mean, ncol = 1)
    row.names(sigma_mean) <- par_names[lower.tri(par_names, diag = TRUE)]
    if(hierarchical == TRUE) {
      output <- list(structures = structures,
                     posterior_probability = posterior_probability,
                     sigma_eap = sigma_mean, theta_eap = theta_mean)
    } else {
      output <- list(structures = structures,
                     posterior_probability = posterior_probability,
                     sigma_eap = sigma_mean)
    }
  }
  return(output)
}


#' Bayesian structure selection for the Ising model.
#'
#' The function \code{select_structure} explores the space of possible
#' Ising network structures using the Gibbs sampler.
#'
#' The function \code{structure_selection_ssvs} is used by
#' \code{select_structure} for exploring the space of network structures. It
#' makes use of a Bayesian spike and slab approach to edge inclusion, and
#' explores the posterior space of network structures using the Gibbs sampler.
#'
#' @inheritParams structure_selection_ssvs
#'
#' @param spike_var,slab_var The \code{p} by \code{p} matrices of variances
#'   that are used in the specification of the spike and slab prior
#'   distributions that are stipulated on the association parameters. Optional,
#'   if \code{spike_var} and \code{slab_var} are not specified, it makes use of
#'   \code{precision} to specify \code{spike_var} and \code{slab_var} using
#'   \code{\link{set_spike_and_slab}}.
#'
#' @param precision A number between zero and one. The prior precision that is
#'   desired for edge selection. Equal to one minus the desired type-1 error.
#'   Needs to be specified if \code{spike_var} and \code{slab_var} are
#'   unspecified. Defaults to \code{.975}.
#'
#' @return A list that contains \code{structures}, an \code{S} by
#'   \code{p(p-1)/2} matrix of binary variables that encodes the \code{S}
#'   structures that were visited by the Gibbs sampler. Includes edge names
#'   for easy reference. Also included is \code{posterior_probability},
#'   the \code{S} by \code{1} matrix of posterior probabilities that matches
#'   the rows of \code{structures}. If \code{output_samples = TRUE} the list
#'   contains an \code{number_iterations} by \code{p(p + 1) / 2} matrix
#'   \code{sigma_samples} that contains samples from the posterior distribution
#'   of the Ising model parameters, and if \code{hierarchical = TRUE}, it also
#'   includes an \code{number_iterations} by \code{1} matrix of
#'   \code{theta_samples} which contains samples from the posterior
#'   distribution of the prior inclusion probability \code{theta}. Otherwise,
#'   if \code{output_samples = FALSE}, the list contains a \code{p(p + 1) / 2}
#'   by \code{1} matrix \code{sigma_eap} of posterior means of the Ising model
#'   parameters, and if \code{hierarchical = TRUE}, it also includes
#'   \code{theta_eap}, the posterior mean of the prior inclusion probability.
#'
#' @examples
#' \dontrun{
#'   library("IsingSampler")
#'   ### Simulate dataset ###
#'   # Input:
#'   p <- 6 # Number of nodes
#'   n <- 1000 # Number of samples
#'   # Ising parameters:
#'   Graph <- matrix(data = 0, nrow = p, ncol = p)
#'   Graph[lower.tri(Graph)] <- rbinom(n = p * (p - 1) / 2,
#'                                     size = 1,
#'                                     prob = 0.2)
#'   Graph <- Graph * runif(n = p ^ 2, min = 0.5, max = 2)
#'   Graph <- Graph + t(Graph)
#'   Thresholds <- -rowSums(Graph) / 2
#'   # Simulate:
#'   Data <- IsingSampler(n = n, graph = Graph, thresholds = Thresholds)
#'   ### Fit using fit_pseudoposterior ###
#'   selection <- select_structure(x = Data,
#'                                 number_iterations = 1e3,
#'                                 number_burnin_iterations = 0,
#'                                 hierarchical = FALSE)
#'   incl_prob <- matrix(0, nrow = p, ncol = p)
#'   incl_prob[lower.tri(incl_prob, diag = FALSE)] <-
#'     t(selection$structures) %*% selection$posterior_probability
#'   incl_prob <- incl_prob + t(incl_prob)
#'   mps <- matrix(0, nrow = p, ncol = p)
#'   mps[lower.tri(mps, diag = TRUE)] <- selection$sigma_eap
#'   mps <- mps + t(mps)
#'   mps[incl_prob < 0.5] <- 0
#'   # Plot results:
#'   library("qgraph")
#'   layout(t(1:2))
#'   qgraph(mps,fade = FALSE)
#'   title("Median probability structure")
#'   qgraph(Graph,fade = FALSE)
#'   title("Original network")
#' }
select_structure <- function(x, spike_var, slab_var, theta = 0.5, alpha = 1,
                             beta = 1, hierarchical = FALSE,
                             number_iterations = 1e5,
                             number_burnin_iterations = 1e2,
                             prior_var_intercepts = 1, output_samples = FALSE,
                             sigma, include, precision = .975,
                             components = "normal", df = 5) {
  if(!hasArg("x"))
    stop("No data.", call. = FALSE)

  if(!hasArg("spike_var") | !hasArg("slab_var")) {
    spike_and_slab <- try(set_spike_and_slab (x = x,
                                              precision = precision),
                          silent = TRUE)
    if(class(spike_and_slab) == "try-error") {
      stop("Could not set spike and slab parameters.", .call. = FALSE)
    }
    spike_var <- spike_and_slab$spike_var
    slab_var <- spike_and_slab$slab_var
    sigma <- spike_and_slab$sigma_ml
  }

  if(hasArg("sigma") & hasArg("include")) {
    gibbs_sample <- structure_selection_ssvs(number_iterations = number_iterations,
                         number_burnin_iterations =
                           number_burnin_iterations, x = x,
                         spike_var = spike_var, slab_var = slab_var,
                         prior_var_intercepts = prior_var_intercepts,
                         hierarchical = hierarchical, alpha = alpha,
                         beta = beta, theta = theta,
                         output_samples = output_samples, sigma = sigma,
                         include = include, components = components, df = df)
  }

  if(hasArg("sigma") & !hasArg("include")) {
    gibbs_sample <- structure_selection_ssvs(number_iterations = number_iterations,
                         number_burnin_iterations =
                           number_burnin_iterations,
                         x = x, spike_var = spike_var, slab_var = slab_var,
                         prior_var_intercepts = prior_var_intercepts,
                         hierarchical = hierarchical, alpha = alpha,
                         beta = beta, theta = theta,
                         output_samples = output_samples,
                         sigma = sigma, components = components, df = df)
  }

  if(!hasArg("sigma")) {
    gibbs_sample <- structure_selection_ssvs(number_iterations = number_iterations,
                         number_burnin_iterations = number_burnin_iterations,
                         x = x, spike_var = spike_var, slab_var = slab_var,
                         prior_var_intercepts = prior_var_intercepts,
                         hierarchical = hierarchical, alpha = alpha,
                         beta = beta, theta = theta,
                         output_samples = output_samples,
                         components = components, df = df)
  }

  return(gibbs_sample)
}
