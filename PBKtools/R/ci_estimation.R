#' Confidence Intervals
#' This code has been developed to be used with the profile_likelihood and identifiability_analysis
#' functions, in order to estimate the confidence intervals of the parameters, based on their
#' likelihood profiles. The method is an implementation of the methodology described in the
#' publication of Raue et al., 2009 entitled as "Structural and practical
#' identifiability analysis of partially observed dynamical models by exploiting
#' the profile likelihood" (doi: https://doi.org/10.1093/bioinformatics/btp358).
#'
#' @param  to be completed
#' @return to be completed
#' @export
ci_estimation <- function(pl_results, alpha=0.95, df=1, theta_hat, global_optimum,
                          lb, ub){
  # INPUT VARIABLES:
  # pl_results        dataframe with 2 columns containing the theta values 
  #                   and the corresponding chi^2 values   
  # alpha             probability of chi-squared to estimate the quantile (it is
  #                   the "p" variable of qchisq() function)
  # df                degrees of freedom of qchisq()
  # theta_hat         the optimal value of the parameter theta
  # global_optimum    scalar defining the global minimum of the objective function
  # lb                the lower bound of the theta_hat parameters
  # ub                the upper bound of the theta_hat parameters
  
  # As confidence interevals of each theta are considered the borders of the
  # following region: 
  # {theta | chi^2(theta) - chi^2(theta_hat) < Delta_alpha} where
  Delta_alpha = qchisq(alpha, df)
  threshold = Delta_alpha + global_optimum
  
  # Estimate the lower bound of the interval.
  # Check if the threshold was exceeded at the backwards search. 
  if(pl_results[1,2] > threshold){
    # Use those 2 theta_i values that yield minimum distances from the threshold
    # and use interpolation to estimate the lower bound of the parameter. 
    
    # This is the last value of the parameter that yields chi^2 lower than the 
    # threshold
    theta_under_threshold <- pl_results[2,1]
    # This is the chi^2 value that this parameters has
    chi2_under_threshold <- pl_results[2,2]
    # This is the last value of the parameter that yields chi^2 grater than the 
    # threshold
    theta_over_threshold <- pl_results[1,1]
    # This is the chi^2 value that this parameters has
    chi2_over_threshold <- pl_results[1,2]
    
    slope <- (chi2_over_threshold - chi2_under_threshold)/(theta_over_threshold - theta_under_threshold)
    dy <- threshold - chi2_under_threshold
    dx <- dy/slope
    lower_bound <- theta_under_threshold + dx
    
  }else{
    # If threshold not exceeded,
    # take the last 20% of points (at least 3) and perform linear fit
    x <- pl_results[pl_results[,1] <= theta_hat[1],1]
    y <- pl_results[pl_results[,1] <= theta_hat[1],2]
    n_last20 <- max(3, length(which(x - theta_hat[i] > 0.8 * (min(x) - theta_hat[i]))))
    x <- head(x, n_last20)
    y <- head(y, n_last20)
    slope <- sum((x - mean(x))*(y - mean(y)))/sum((x - mean(x))^2)
    
    # If slope < 0, return Inf
    if (slope > 0) lower_bound = -Inf
    
    # Extrapolate until threshold is passed
    dy <- threshold - head(y, 1)
    dx <- dy/slope
    lower_bound <- head(x, 1) + dx
    
    # Test if extrapolation takes the point of passage very far
    # Set to Inf in that case
    if (lower_bound > 10*(max(x) - min(x)) | lower_bound < lb){
      lower_bound <- -Inf
    }else{
      added_row <- as.data.frame(cbind(lower_bound[[1]], threshold[[1]]))
      colnames(added_row) <- colnames(pl_results)
      pl_results <- rbind(added_row, pl_results)
    }
  }
  
  # Check if the threshold was exceeded at the backwards search. 
  if(pl_results[dim(pl_results)[1],2] > threshold){
    # Use those 2 theta_i values that yield minimum distances from the threshold
    # and use interpolation to estimate the lower bound of the parameter. 
    
    # This is the last value of the parameter that yields chi^2 lower than the 
    # threshold
    theta_under_threshold <- pl_results[(dim(pl_results)[1]-1),1]
    # This is the chi^2 value that this parameters has
    chi2_under_threshold <- pl_results[(dim(pl_results)[1]-1),2]
    # This is the last value of the parameter that yields chi^2 grater than the 
    # threshold
    theta_over_threshold <- pl_results[dim(pl_results)[1],1]
    # This is the chi^2 value that this parameters has
    chi2_over_threshold <- pl_results[dim(pl_results)[1],2]
    
    slope <- (chi2_over_threshold - chi2_under_threshold)/(theta_over_threshold - theta_under_threshold)
    dy <- threshold - chi2_under_threshold
    dx <- dy/slope
    upper_bound <- theta_under_threshold + dx
  }else{
    # If threshold not exceeded,
    # take the last 20% of points (at least 3) an perform linear fit
    x <- pl_results[pl_results[,1] >= theta_hat[1],1]
    y <- pl_results[pl_results[,1] >= theta_hat[1],2]
    n_last20 <- max(3, length(which(x - theta_hat[1] > 0.8 * (max(x) - theta_hat[1]))))
    x <- tail(x, n_last20)
    y <- tail(y, n_last20)
    slope <- sum((x - mean(x))*(y - mean(y)))/sum((x - mean(x))^2)
    
    # If slope < 0, return Inf
    if (slope < 0) upper_bound <- Inf
    # Extrapolate until threshold is passed
    dy <- threshold - tail(y, 1)
    dx <- dy/slope
    upper_bound <- tail(x, 1) + dx
    
    # Test if extrapolation takes the point of passage very far
    # Set to Inf in that case
    if (upper_bound > 10*(max(x) - min(x)) | upper_bound > ub){
      upper_bound <- Inf
    }else{
      added_row <- as.data.frame(cbind(upper_bound[[1]], threshold[[1]]))
      colnames(added_row) <- colnames(pl_results)
      pl_results <- rbind(pl_results,added_row)
    }
  }
  
  return(list("Lower_bound" = lower_bound,
              "Upper_bound" = upper_bound)) 
}