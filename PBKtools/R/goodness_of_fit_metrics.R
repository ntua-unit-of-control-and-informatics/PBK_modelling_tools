#' Mean Squared Error
#' 
#' Estimates the Mean Squared Error 
#'
#' @param  observed vector with the observed values
#' @param predicted vector with the predicted values
#' @export
mse_custom <- function(observed, predicted){
  mean((observed - predicted)^2)
}

#' Mean Absolute Percentage Error
#' 
#' Estimates the mean absolute percentage error.
#'
#' @param  observed vector with the observed values
#' @param predicted vector with the predicted values
#' @export
mape <- function(observed, predicted){
  mean(abs(observed-predicted)*100/observed)
}

#' Root Mean Squared Error
#' 
#' Estimates the root mean squared error.
#'
#' @param  observed vector with the observed values
#' @param predicted vector with the predicted values
#' @export
rmse <- function(observed, predicted){
  sqrt(mean((observed-predicted)^2)) 
}

#' Absolute Average-Fold Error  
#' 
#' Estimates the absolute average-fold error.
#'
#' @param  predictions list with the predicted values
#' @param observations list with the observed values
#' @export
AAFE <- function(predictions, observations, times=NULL){
  y_obs <- unlist(observations)
  y_pred <- unlist(predictions)
  # Total number of observations
  N<- length(y_obs)
  log_ratio <- rep(NA, N) 
  for ( i in 1:N){
    log_ratio[i] <- abs(log((y_pred[i]/y_obs[i]), base = 10))
  }
  aafe <- 10^(sum(log_ratio)/N) 
  return(aafe)
}

#' PBK Objective Function
#' 
#' A custom metric for PBK models.
#'
#' @param  observed list with the observed values
#' @param predicted list with the predicted values
#' @export
PBKOF <- function(observed, predicted, comp.names =NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted)){
    stop(" The observations and predictions must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted)){
    stop(" The observations and predictions must have the same compartments")
  }
  Ncomp <- length(observed) # Number of compartments
  I <- rep(NA, Ncomp) # Compartment discrepancy index
  N_obs <- rep(NA, Ncomp) #Number of observations per compartment
  #loop over the compartments
  for (i in 1:Ncomp){
    Et <- 0 #relative error with observations
    St <- 0  #relative error with simulations
    N <- length(observed[[i]]) # number of observations for compartment i
    # Check if observations and predictions have equal length
    if(N != length(predicted[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    N_obs[i] <- N # populate the N_obs vector
    for (j in 1:N){
      # sum of relative squared errors (error = observed - predicted)
      Et <- Et + ( abs(observed[[i]][j] - predicted[[i]][j])  / observed[[i]][j] )  ^2
      St <- St + ( abs(observed[[i]][j] - predicted[[i]][j])  / predicted[[i]][j] )  ^2
    }
    
    # root mean of the square of observed values
    RMEt <- sqrt(Et/N)
    # root mean of the square of simulated values
    RMSt <- sqrt( St/N)
    
    I[i] <- (RMEt + RMSt)/2   
  }
  # Total number of observations
  Ntot <- sum(N_obs)
  # Initialise the consolidated discrepancy index
  Ic <-0
  for (i in 1:Ncomp){
    # Give weight to compartments with more observations (more information)
    Ic <- Ic +  I[i]* N_obs[i]/Ntot
  }
  # Name the list of compartment discrepancy indices
  if ( !is.null(comp.names)){
    names(I) <- comp.names
  }else if (!is.null(names(observed))){
    names(I) <- names(observed)
  } else if (!is.null(names(predicted)) && is.null(comp.names) ){
    names(I) <- names(predicted)
  }
  return(Ic)
  #return(list(Total_index = Ic, Compartment_index= I))
}

#' Weighted Sum of Squared Residuals
#' 
#' Estimates the weighted sum of squared residuals.
#'
#' @param  observed list with the observed values
#' @param predicted list with the predicted values
#' @param weights_values list with the weights
#' @export 
WSSR <- function(observed, predicted, weights_values, comp.names =NULL){
  # Check if the user provided the correct input format
  if (!is.list(observed) || !is.list(predicted) || !is.list(weights_values)){
    stop(" The observations, predictions and weights_values must be lists")
  }
  # Check if the user provided equal length lists
  if (length(observed) != length(predicted) || length(observed) != length(weights_values)){
    stop(" The observations, predictions and weights_values must have the same compartments")
  }
  
  # Define the number of observed outputs
  N_outputs <- length(predicted)
  # Define the number of observations per output
  N_obs <- rep(NA, N_outputs)
  
  # A vector to store the values of the weighted squared sum residuals of each compartment
  outputs_res <- c()
  for (i in 1:N_outputs) { # loop over the observed outputs
    N_obs[i] <- length(observed[[i]])
    
    # Check that all observed, predicted and weights_values vectors have the same length
    if(N_obs[i] != length(predicted[[i]]) || N_obs[i] != length(weights_values[[i]])){
      stop(paste0("Compartment ",i," had different length in the observations and predictions"))
    }
    # The number of observations for output i
    N <- N_obs[i]
    
    # Initiate a variable to estimate the sum of squared residuals for output j
    sq_weighted_res_j <- 0
    for (j in 1:N) { #loop over the experimental points i compartment i
      sq_weighted_res_j <- sq_weighted_res_j + ((observed[[i]][j] - predicted[[i]][j]) / weights_values[[i]][j])^2   
    }
    outputs_res[i] <- sq_weighted_res_j
  }
  
  WSSR_results <- sum(outputs_res)
  
  return(WSSR_results)
}
