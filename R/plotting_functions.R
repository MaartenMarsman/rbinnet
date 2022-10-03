#' Plot posterior likelihood distribution of all structures
#'
#' The function is based in the ggplot2 environment and plots the posterior likelihood
#' of the structures sorted from the most to the least probable. The function can either
#' plot the raw posterior probabilities or transformed as Bayes Factors depicting the 
#' Bayes Factor of the respective structure against the most probable structure. 
#' 
#' Additional arguments can be passed to the function using by adding + and the common 
#' ggplot2 commands.
#'
#' @param output Output of the rbinnet \code{select_structure} function
#' @param as.BF Determines if the raw posterior probabilities are plotted or the transformed Bayes Factor
#'
#' @return Plot of the posterior likelihood of the structures. 


plot_structure_probability <- function(output, as.BF = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  require("ggplot2")
  
  tmp <- sort(output$structure$posterior_probability, decreasing = TRUE) #tmp[1] is most plausible
  if(as.BF){
    BF1s <- tmp / tmp[1] # BF best structure vs. others
    data <- data.frame(structures = 1:length(BF1s), BayesFactor = BF1s)
    ggplot2::ggplot(data, aes(x = structures, y = BayesFactor)) +
      geom_point(size = 4, shape = 1) + 
      scale_y_continuous(trans = "log10") +
      theme_classic()+
      labs(x = "Structures", 
           y = expression(log(BF[1][s])))+
      geom_hline(yintercept = 1/10, linetype = "dashed", size = 1.5)  +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 1.1),
            axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            axis.ticks.length = unit(.25, "cm"))
  } else {
    data <- data.frame(structures = 1:length(tmp), Probs = tmp)
    ggplot2::ggplot(data, aes(x = structures, y = Probs)) +
      geom_point(size = 4, shape = 1) + 
      theme_classic()+
      labs(x = "Structures", 
           y = "Posterior Structure Probability")+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 1.1),
            axis.text = element_text(size = 14), axis.title = element_text(size = 16),
            axis.ticks.length = unit(.25, "cm"))
  }
}

# ---------------------------------------------------------------------------------------------------------------
# Plot most probable structure

#' Most probable structures plot.
#' 
#' The function plots the network of the most probable structures.
#' 
#' This function uses \code{qgraph} for the visualization of the network. 
#' Additional arguments from the \code{qgraph} package can be passed on within this function.  
#' 
#' @param output Output of the rbinnet \code{select_structure} function
#' @param top_n Number of most probable structures plotted.
#' @param ... Additional arguments to be passed on to qgraph
#'
#' @return Network plot of most probable structures.
 
plot_structure <- function(output, top_n = 3, ...){
  
  if (!requireNamespace("qgraph", quietly = TRUE)) {
    stop("Package \"qgraph\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  highest_probs <- tail(sort(output$structure$posterior_probability), top_n)
  most_probable <- match(highest_probs, output$structure$posterior_probability)
  
  for(index in most_probable){
    graph_structure <- vector_to_matrix(output$structure$structures[index, ], output$nodes)
    diag(graph_structure) <- 1
    
    qgraph::qgraph(graph_structure, layout = "circle", theme = "TeamFortress", 
                   title = paste("Posterior Probability = ", output$structure$posterior_probability[index]), 
                   ...) 
  }
}

# ---------------------------------------------------------------------------------------------------------------
# plot for edge inclusion BF

#' Inclusion evidence of edges plot.
#' 
#' The function creates a network showing the evidential strength for inclusion /exclusion. 
#' Red edges indicate evidence for exclusion, blue edges evidence for inclusion, and grey edges 
#' absence of evidence with a Bayes Factor below the evidence threshold.
#' 
#' This function uses \code{qgraph} for the visualization of the network. 
#' Additional arguments from the \code{qgraph} package can be passed on within this function.    
#'
#' @param output Output of the rbinnet \code{select_structure} function
#' @param evidence_thresh numeric larger 1. Specifies the inclusion Bayes Factor above which
#' edges are considered included. Default set to 10. The higher, the more conservative the evidence plot.  
#' @param ... Additional arguments to be passed on to qgraph
#'
#' @return Function returns a network with edges indicating evidential strength for in-/exclusion.

plot_edge_BF <- function(output, evidence_thresh = 10, ...) {
  
  if (!requireNamespace("qgraph", quietly = TRUE)) {
    stop("Package \"qgraph\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  graph <- vector_to_matrix(output$parameters$inc_BF, output$nodes)
  diag(graph) <- 1
  
  # assign a color to each edge (inclusion - blue, exclusion - red, no conclusion - grey)
  graph_color <- graph
  graph_color <-  ifelse(graph < evidence_thresh & graph > 1/evidence_thresh, graph_color <- "#bfbfbf", graph_color <- "#36648b")
  graph_color[graph < (1/evidence_thresh)] <- "#990000"
  
  qgraph::qgraph(matrix(1, ncol = output$nodes, nrow = output$nodes),
                 theme = "TeamFortress", edge.width = 6, 
                 edge.color = graph_color, # specifies the color of the edges, depending on inclusion probability
                 ...
  )
}

# ---------------------------------------------------------------------------------------------------------------
# Plot median probability model

#' Median probability model plot.
#' 
#' The function depicts the median probability structure of the Ising model, eliminating edges 
#' with an inclusion probability below the specified \code{exc_prob}-parameter 
#' (i.e., default set to exc_prob 0.5). Within the network, nodes represent the 
#' variables of the dataset and lines the respective node interaction. The thicker 
#' the line the stronger the interaction between two nodes.
#' 
#' This function uses \code{qgraph} for the visualization of the network. 
#' Additional arguments from the \code{qgraph} package can be passed on within this function. 
#'
#' @param output Output of the rbinnet \code{select_structure} function
#' @param exc_prob Number between 0 and 1. Specifies minimum inclusion probability for edges
#' to be depicted in the network.
#' @param ... Additional arguments to be passed on to qgraph
#'



plot_mpm <- function(output, exc_prob = .5, ...) {
  
  if (!requireNamespace("qgraph", quietly = TRUE)) {
    stop("Package \"qgraph\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  graph <- vector_to_matrix(output$parameters$sigma_eap, output$nodes, diag = F)
  
  # Exclude edges with a inclusion probability lower .5
  inc_probs_m <- vector_to_matrix(output$parameters$inclusion_probabilities, output$nodes)
  graph[inc_probs_m < exc_prob] <- 0
  diag(graph) <- 1
  
  # Plot
  qgraph::qgraph(graph, maximum=1, ...) 
  
}

# ---------------------------------------------------------------------------------------------------------------
# HDI plot 

#' Highest density interval of parameters plot.
#' 
#' The function plots the 95% highest density interval (HDI) of the posterior sample of the Ising
#' parameters. Dots represent the median of the posterior distribution and lines the respective 
#' 95% HDI. 
#' 
#' The function is based in the ggplot2 environment. Additional arguments can be passed on 
#' using the common "+" notation and a respective ggplot command. 
#'
#' @param output Output of the rbinnet \code{select_structure} function
#' @param thresholds binary Determines whether the threshold parameters should be plotted
#'
#' @return Forest plot depicting the median of parameters and their respective 95% HDI

plot_parameter_HDI <- function(output, thresholds = F) {
  
  if(is.null(output$parameters$sigma_samples)){
    stop("Samples of the posterior distribution required. Set \"output_samples = TRUE\".")
  }
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  require("ggplot2")
  
  hdi_intervals <- as.data.frame(apply(output$parameters$sigma_samples, MARGIN = 2, FUN = hdi_interval))
  posterior_medians <- apply(output$parameters$sigma_samples, MARGIN = 2, FUN = median)
  
  posterior <- cbind(colnames(hdi_intervals), data.frame(posterior_medians, row.names = NULL), data.frame(t(hdi_intervals), row.names = NULL))
  colnames(posterior) <- c("parameter", "posterior_medians", "lower", "upper")
  
  if(!thresholds) {
    # Filter interaction parameters (sigma), exclude thresholds
    index <- grep('^s', posterior$parameter)
    posterior <- posterior[index, ]
  }
  
  
  ggplot2::ggplot(data = posterior, aes(x = parameter, y = posterior_medians, ymin = lower, 
                                        ymax = upper)) + 
    geom_pointrange(position=position_dodge(width=c(0.5)), size = .9) + 
    theme_bw() + 
    coord_flip() + 
    ylab("Highest Density Interval of Parameter")+
    xlab("") +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1.3) + 
    theme(axis.text=element_text(size=14), panel.border = element_blank(), 
          axis.line = element_line(colour = "black", size = 1.1), axis.ticks.length=unit(.2, "cm"),
          axis.ticks = element_line(size= .8),
          axis.title.x = element_text(size=16,face="bold"), plot.title = element_text(size = 18, face = "bold"))
}

# ---------------------------------------------------------------------------------------------------------------
# Sigma samples

#' Plot posterior distribution of the parameter samples.
#' 
#' The function plots the density of the posterior distribution of the parameter samples. 
#'
#' @param output Output of the rbinnet \code{select_structure} function
#' @param parameter A vector stating which parameters of the model should be plotted. 
#' Options "all", "sigma", and "mu". Default is set to "all"
#'
#' @return Depicts density plots of parameters' posterior distribution.

plot_parameter_distribution <- function(output, parameter = c("all")){
  
  if(is.null(output$parameters$sigma_samples)){
    stop("Samples of the posterior distribution required. Set \"output_samples = TRUE\".")
  }
  
  sigma_samples <- output$parameters$sigma_samples
  par(mfrow = c(2,2))
  
  if(parameter == "sigma") {
    index_sigma <- grep('^s', colnames(sigma_samples))
    for(index in index_sigma){
      plot(density(sigma_samples[, index]), main = colnames(sigma_samples)[index], 
           xlab = "Parameter Estimate", cex = 4, lwd = 1.2, bty = "n")
    }
    
  } else if (parameter == "mu"){
    index_mu <- grep('^m', colnames(sigma_samples))
    for(index in index_mu){
      plot(density(sigma_samples[, index]), main = colnames(sigma_samples)[index], 
           xlab = "Parameter Estimate", cex = 4, lwd = 1.2, bty = "n")
    }
    
  } else {
    # Plot mu 
    index_mu <- grep('^m', colnames(sigma_samples))
    for(index in index_mu){
      plot(density(sigma_samples[, index]), main = colnames(sigma_samples)[index], 
           xlab = "Parameter Estimate", cex = 4, lwd = 1.2, bty = "n")
    }
    
    # Plot sigma
    index_sigma <- grep('^s', colnames(sigma_samples))
    for(index in index_sigma){
      plot(density(sigma_samples[, index]), main = colnames(sigma_samples)[index], 
           xlab = "Parameter Estimate", cex = 4, lwd = 1.2, bty = "n")
    }
  }
  par(mfrow = c(1,1))
}
