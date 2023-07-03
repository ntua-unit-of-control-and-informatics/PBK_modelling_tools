library(deSolve)
library(ggplot2)
library(pheatmap)

#---------------------
# Custom AUC function
#---------------------
AUC <- function(x, y){
  individual_auc <- c()
  for (i in 1:(length(x)-1)){
    individual_auc[i] <- (y[i]+y[i+1])*(x[i+1]-x[i])/2
  }
  return(sum(individual_auc))
}

#---------------------------------------------------
# SENSITIVITY ANALYSIS FUNCTION for PBPK Models
#---------------------------------------------------
# The function has the following input variables:
# - model: a function of ODEs ready to be used with the deSolve library
# - thetas: a vector with the values of the parameters, which will be analysed
# - thetas_names: a vector with the names (as characters) of the thetas
# - constant_params: a list with the constant parameters used in the model
# - ranges: the perturbation applied to the parameters
# - targets: a vector with the names of the outputs of interest
# - ode_settings: a list with the settings to be used in the deSolve solver
# - heatmap: binary variable that defines if the heatmap plots will be returned

PBPK_sensitivity <- function(model, thetas, thetas_names, constant_params, ranges, targets,
                             ode_settings, heatmap = FALSE){
  
  if(is.null(model)){
    stop("The ODEs of the PBPK model must be provided as a function compatible
         with the deSolve library")
  }
  if(is.null(thetas)){
    stop("A vector with the values of the parameters for the 
         analysis must be provided")
  }
  
  # The variation dp of each parameter will be equal to "ranges" value 
  if(!is.numeric(ranges) | length(ranges) != 1 | ranges<=0 | ranges>1){
    # "ranges" must be a single numeric value in (0,1]
    stop("For local sensitivity analysis \"ranges\"  should be 
         a single numeric value in (0,1]")
  }else{dp <- ranges}
  
  # The total number of parameters to be analysed. 
  N_thetas <- length(thetas) 
  # The number of target outputs of interest
  N_targets <- length(targets)
  
  # Take all objects from the ode_settings list 
  inits <- ode_settings$inits # Initial conditions of the ODEs
  sample_time <- ode_settings$sample_time # Time points of solution
  # Define the solver to use later
  solver <- ifelse(ode_settings$solver == "default", "bdf", ode_settings$solver)
  
  # Assign the corresponding names to the parameters
  names(thetas) <- thetas_names
  # Merge into a vector the constant parameters and the parameters of the sensitivity test
  params <- c(params_list,thetas)
  
  # Solve the ODEs for the initial values of parameters
  solution_0 <- deSolve::ode(times = sample_time,  func = model, y = inits, parms = params, 
                             method=solver, rtol = 1e-5, atol = 1e-5)
  
  if(sum(!(targets %in% colnames(solution_0))) != 0){
    stop("As \"targets\" should be provided the name(s) of one or more
    of the outputs (compertments) of the PBPK model")
  }
  
  # Calculate AUC of each target for the initial parameters
  AUC_0 <- c()
  # Keep the maximum concentration of each output variable
  Cmax_0 <- c()
  for (i in 1:N_targets) {
    AUC_0[i] <- AUC(solution_0[,"time"],solution_0[,targets[i]])
    Cmax_0[i] <- max(solution_0[,targets[i]])
  }
  
  # Assign the optimal values of thetas to thetas_0 vector
  thetas_0 <- thetas
  # Initialize a vector to store the sensitivity indexes
  SI_auc <- matrix(NA, nrow = N_thetas, ncol = N_targets)
  SI_cmax <- matrix(NA, nrow = N_thetas, ncol = N_targets)
  for (i in 1:N_thetas) {
    # Update the thetas vector with the initial values of parameters
    thetas <- thetas_0
    # Apply the desired perturbation to the i-th parameter
    thetas[i] <- thetas[i]*(1 + dp)
    # Assign names to the parameters
    names(thetas) <- thetas_names
    # Merge into a vector the constant parameters and the parameters of the sensitivity test
    params <- c(params_list,thetas)
    
    # Solve the ode system for given initial values and time
    solution <- deSolve::ode(times = sample_time,  func = model, y = inits, parms = params, 
                             method=solver, rtol = 1e-5, atol = 1e-5)
    
    for (j in 1:N_targets) {
      # Calculate AUC for the target compartment j
      AUC_j <- AUC(solution[,"time"],solution[,targets[j]])
      # Keep the max value for the target compartment j
      Cmax_j <- max(solution[,targets[j]])
      # Calculate sensitivity index of parameter i 
      # Relative Sensitivity of AUC = (dAUC/AUC)/(dp/p)
      SI_auc[i,j] <- ((AUC_j-AUC_0[j])/AUC_0[j])/((thetas[i]- thetas_0[i])/thetas_0[i])
      # Relative Sensitivity of Cmax = (dCmax/Cmax)/(dp/p)
      SI_cmax[i,j] <- ((Cmax_j-Cmax_0[j])/Cmax_0[j])/((thetas[i]- thetas_0[i])/thetas_0[i])
    }
  }
  
  rownames(SI_auc) <- thetas_names
  colnames(SI_auc) <-  paste0("AUC_", targets)
  
  rownames(SI_cmax) <- thetas_names
  colnames(SI_cmax) <- paste0("Cmax_", targets)
  
  if(heatmap){
    
    if(length(targets)==1){
      stop("Provide more than 1 targets in order to create a heatmap
           or turn \"heatmap\" to \"FALSE\".")
    }
    
    
    #########################################################################
    
    heatmap1 <- pheatmap::pheatmap(as.matrix(abs(SI_auc)),
                                   cellwidth = 60, cellheight = 60,
                                   border_color = NA,
                                   fontsize = 20,
                                   fontsize_row = 15, 
                                   fontsize_col = 15,
                                   col = colorRampPalette(RColorBrewer::brewer.pal(8, "YlOrRd"))(25), 
                                   main="AUC Senesitivity Indexes")
    
    AUC_heatmap <- recordPlot(heatmap1)
    
    heatmap2 <- pheatmap::pheatmap(as.matrix(abs(SI_cmax)),
                                   cellwidth = 60, cellheight = 60,
                                   border_color = NA,
                                   fontsize = 20,
                                   fontsize_row = 15, 
                                   fontsize_col = 15,
                                   col = colorRampPalette(RColorBrewer::brewer.pal(8, "YlOrRd"))(25),
                                   main="Cmax Senesitivity Indexes")
    Cmax_heatmap <- recordPlot(heatmap2)
    
    plot_list <- list("AUC_heatmap"=AUC_heatmap, "Cmax_heatmap"=Cmax_heatmap)
  }
  
  # Create a list to return the results 
  
  data_list <- list("Method"="Local Sensitivity", "Targets"=targets,
                    "Change of parameters" = dp,
                    "Parameters" = thetas_names, 
                    "Normalized AUC Sensitivity Coefficients"=data.frame(SI_auc), 
                    "Normalized C_max Sensitivity Coefficients"=data.frame(SI_cmax),
                    "Plot_list"=plot_list)
  
  return(data_list)
}
