# Plot Structure probabilities
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

plot_structure <- function(output, top_n = 3){
  
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
                   color= c("#f0ae0e"), vsize = 8, repulsion = .9,
                   legend = F, legend.cex = 0.55, 
                   title = paste("Posterior Probability = ", output$structure$posterior_probability[index])) 
  }
}

# plot for edge inclusion BF
plot_edge_BF <- function(output, evidence_thresh = 10) {
  
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
                 theme = "TeamFortress", 
                 color= c("#f0ae0e"),vsize = 14, repulsion = .9, maximum = 1,
                 legend = F,label.cex = 1.2, edge.width = 6, 
                 edge.color = graph_color # specifies the color of the edges
  )
}

# ---------------------------------------------------------------------------------------------------------------
# Plot median probability model

plot_mpm <- function(output, exc_prob = .5) {
  
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
  qgraph::qgraph(graph, theme = "TeamFortress", 
                 color= c("#f0ae0e"), vsize = 14, repulsion = .9,
                 maximum=1, legend = F,  label.cex = 1.2) 
  
}

# ---------------------------------------------------------------------------------------------------------------
# HDI plot 

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
