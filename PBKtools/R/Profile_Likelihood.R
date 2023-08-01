#' Profile Likelihood
#'
#' This code is an implementation of the Profile Likelihood-based identifiability
#' analysis, as described
#' in the publication of Raue et al., 2009 entitled as "Structural and practical
#' identifiability analysis of partially observed dynamical models by exploiting
#' the profile likelihood" (doi: https://doi.org/10.1093/bioinformatics/btp358).
#' 
#' @param to be completed
#' @return to be completed
#' @export
profile_likelihood <- function(obj_f, 
                               i,
                               thetas,
                               thetas_names, 
                               constant_params = NULL,
                               data_df,
                               errors_df,
                               lb, ub, N_samples,
                               alpha, df, q, global_optimum, 
                               min_step_coef, max_step_coef,
                               break_at_bounds = FALSE,
                               # nlopt settings for the main optimization problem
                               opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
                                           "xtol_rel" = 1e-06, 
                                           "ftol_rel" = 1e-06,
                                           "ftol_abs" = 0.0,
                                           "xtol_abs" = 0.0 ,
                                           "maxeval" = 300,
                                           "print_level" = 1),
                               # nlopt settings for the estimation of theta_step
                               opts_theta_step = list("algorithm" = 'NLOPT_LN_SBPLX',
                                                      "xtol_rel" = 1e-05, 
                                                      "ftol_rel" = 1e-05,
                                                      "ftol_abs" = 0.0,
                                                      "xtol_abs" = 0.0 ,
                                                      "maxeval" = 50,
                                                      "print_level" = 0),
                               create_txt = FALSE){
  
  # INPUT VARIABLES:
  # obj_f             function that returns the value of the objective function
  # i                 integer defining the position of the parameter to estimate the profile-likelihood
  # thetas            vector containing the optimal values of parameters 
  # thetas_names      vector of characters with the names of the parameters
  # constant_params   vector with any extra constant parameters of the model and their names
  # data_df           dataframe containing the data used in the obj_f 
  # errors_df         dataframe with the meassured experimental errors (or SD values) 
  # lb                vector with the lower bounds of the "thetas" parameters
  # ub                vector with the upper bounds of the "thetas" parameters
  # N_samples         integer defining the number of samples taken around the theta optimal
  #                   value (N samples will be taken from each side of the theta)
  # alpha             probability of chi-squared to estimate the quantile (it is
  #                   the "p" variable of qchisq() function)
  # df                degrees of freedom of qchisq()
  # global_optimum    scalar defining the global minimum of the objective function
  # q                 a variable used in the estimation of the adaptive step (see Raue et al., 2009)
  # max_step_coef     coefficient defining the maximum permitted step  
  # min_step_coef     coefficient defining the minimum permitted step
  # break_at_bounds   logical; if TRUE the the sampling of the parameter stops because the bounds were exceeded
  # opts              list with the options selected for the minimization of the objective function
  #                   (check the nloptr package for more details)
  # opts_theta_step   list with the options selected for the estimation of the adaptive step
  #                   (check the nloptr package for more details)
  # create_txt        logical; if TRUE a txt file will be created at the current working directory 
  #                   to save the samples of the parameters and the corresponding values of the objective function
  
  # Estimate the Delta_alpha parameter
  Delta_alpha <- qchisq(alpha, df)
  
  # Take the name of the i-th parameter
  theta_i_name <- names(thetas)[i] # the name of the theta_i parameter
  
  # Function to estimate the theta_step by solving the equation
  # chi^2(theta) - chi^2(theta_hat) - q*Delta_alpha = 0 
  theta_step_estimation <- function(theta_step, theta_last, obj_f, constant_params, index, current_score, q, Delta_alpha, direction){
    i <- index
    x <- theta_last
    if(direction=='forward'){
      x[i] <- x[i] + theta_step
    }else if(direction=='backward'){
      x[i] <- x[i] - theta_step}
    chi2_last <- current_score
    
    chi2 <- obj_f(x = x, constant_theta = NULL, constant_theta_name = NULL, params_names = names(x),
                  constant_params = constant_params,
                  data_df = data_df, errors_df = errors_df)
    
    return(abs(chi2 - chi2_last - q*Delta_alpha))
  }
  
  # Set the threshold. The threshold is estimated as global_optimum + Delta_alpha 
  threshold =  global_optimum + Delta_alpha
  
  if(create_txt) sink(paste0("check_progress_",thetas_names[i], ".txt")) #keep a txt to check the progress while running
  
  # Forward search
  cat("Forward search begins...", "\n")
  # Initiate a counter for the number of iterations
  iter_counter <- 0
  
  # Inititate a parameter to check if the threshold is exceeded in each iteration
  current_score <- global_optimum
  
  # Create a 2-columns dataframe to store the tested theta_i values
  # and the corresponding values of the objective function
  forward_results_df <- data.frame(matrix(NA, nrow=N_samples, ncol = 2))
  colnames(forward_results_df) <- c(theta_i_name, "Likelihood")
  
  # take as constant_theta the parameter that will be fixed to a constant value
  constant_theta <- thetas[i]
  
  # Initiate theta_last and set it equal with the theta_hat
  theta_last <- thetas
  
  # set the current score equal global optimum to estimate the first step
  current_score <- global_optimum 
  while (iter_counter < N_samples & current_score < threshold) {
    iter_counter <- iter_counter + 1
    # Estimation of theta_step 
    theta_step <-  nloptr(x0 = 0.25*abs(thetas[[i]]),
                                  eval_f = theta_step_estimation,
                                  lb = min_step_coef*abs(thetas[[i]]),
                                  ub = max_step_coef*abs(thetas[[i]]),
                                  opts = opts_theta_step,
                                  theta_last = theta_last,
                                  index = i,
                                  current_score = current_score, 
                                  q = q, 
                                  Delta_alpha=Delta_alpha,
                                  direction='forward',
                                  obj_f=obj_f, 
                                  constant_params=constant_params)$solution
    
    # If theta_step exceeds the limits, take the value of the limit
    if(theta_step > max_step_coef*abs(thetas[[i]])){
      theta_step <- max_step_coef*abs(thetas[[i]])
    }else if(theta_step < min_step_coef*abs(thetas[[i]])){
      theta_step <- min_step_coef*abs(thetas[[i]])
    }
    
    # The parameter whose profile likelihood is examined
    constant_theta <- constant_theta + theta_step
    
    #Check if the constant_theta exceeded the corresponding upper boundary
    if(constant_theta > ub[i]){
      f_exit = 3
      break
    } 
    # The name of the constant theta parameter
    constant_theta_name <- theta_i_name
    # Take as x0 all the parameters except the one which will be fixed
    x0 <- theta_last[-i]
    # The names of the rest parameters
    params_names <- names(x0)
    
    # define the lower and upper bounds of the parameters
    set.seed(12312)
    optimization<- nloptr(x0 = x0,
                                  eval_f = obj_f,
                                  lb = lb[-i],
                                  ub = ub[-i],
                                  opts = opts,
                                  constant_theta = constant_theta,
                                  constant_theta_name = constant_theta_name,
                                  params_names = params_names,
                                  constant_params=constant_params,
                                  data_df = data_df,
                                  errors_df = errors_df)
    
    cat("Calculating PL for parameter ", theta_i_name, " and iter = ", iter_counter,". ", "step =", theta_step , constant_theta_name ," = ", constant_theta , "=> Likelihood = ", optimization$objective, "\n")  
    
    forward_results_df[iter_counter,] <- c(constant_theta, optimization$objective)
    
    # update current score
    current_score <- optimization$objective
    
    #update theta_last vector with the new optimized values
    for (k in 1:length(theta_last)) {
      if(k==i){
        theta_last[k] <- constant_theta
      }else if(k>i){
        theta_last[k] <- optimization$solution[k-1]   
      }else{
        theta_last[k] <- optimization$solution[k]} 
    }
    
  }
  if(iter_counter >= N_samples){
    f_exit <- 1
  }else if(current_score > threshold){
    f_exit <- 2
  }
  cat("Forward search ended after ", iter_counter, "iterations.", "\n")
  
  
  
  # Backward search
  cat("Backward search begins...", "\n")
  # Initiate a counter for the number of iterations
  iter_counter <- 0
  # Inititate a parameter to check if the threshold is exceeded in each iteration
  current_score <- global_optimum
  # Create a 2-columns dataframe to store the tested theta_i values
  # and the corresponding values of the objective function
  backward_results_df <- data.frame(matrix(NA, nrow=N_samples, ncol = 2))
  colnames(backward_results_df) <- c(theta_i_name, "Likelihood")
  # take as constant_theta the parameter that will be fixed to a constant value
  constant_theta <- thetas[i]
  # Initiate theta_last and set it equal with the theta_hat
  theta_last <- thetas
  # set the current score equal global optimum to estimate the first step
  current_score <- global_optimum 
  while (iter_counter < N_samples & current_score < threshold) {
    iter_counter <- iter_counter + 1
    # Estimation of theta_step 
    theta_step <-  nloptr(x0 = 0.25*abs(thetas[[i]]),
                                  eval_f = theta_step_estimation,
                                  lb  = min_step_coef*abs(thetas[[i]]),
                                  ub = max_step_coef*abs(thetas[[i]]),
                                  opts = opts_theta_step,
                                  theta_last = theta_last, index = i,
                                  current_score = current_score, q = q, Delta_alpha=Delta_alpha,
                                  direction='backward',
                                  obj_f=obj_f, 
                                  constant_params=constant_params)$solution
    
    # If theta_step exceeds the limits, take the value of the limit
    if(theta_step > max_step_coef*abs(thetas[[i]])){
      theta_step <- max_step_coef*abs(thetas[[i]])
    }else if(theta_step < min_step_coef*abs(thetas[[i]])){
      theta_step <- min_step_coef*abs(thetas[[i]])
    }
    
    # The parameter whose profile likelihood is examined
    constant_theta <- constant_theta - theta_step
    #Check if the constant_theta exceeded the corresponding lower boundary
    if(constant_theta < lb[i]){
      b_exit <- 3
      break
    } 
    # The name of the constant theta parameter
    constant_theta_name <- theta_i_name
    # Take as x0 all the parameters except the one which will be fixed
    x0 <- theta_last[-i]
    # The names of the rest parameters
    params_names <- names(x0)
    
    set.seed(12312)
    optimization<- nloptr(x0 = x0,
                                  eval_f = obj_f,
                                  lb  = lb[-i],
                                  ub = ub[-i],
                                  opts = opts,
                                  constant_theta = constant_theta,
                                  constant_theta_name = constant_theta_name,
                                  params_names = params_names,
                                  constant_params=constant_params,
                                  data_df = data_df,
                                  errors_df = errors_df)
    
    cat("Calculating PL for parameter ", theta_i_name, " and iter = ", iter_counter,". ", "step =", theta_step , constant_theta_name ," = ", constant_theta , "=> Likelihood = ", optimization$objective, "\n")  
    
    backward_results_df[iter_counter,] <- c(constant_theta, optimization$objective)
    
    # update current score
    current_score <- optimization$objective
    
    #update theta_last vector with the new optimized values
    for (k in 1:length(theta_last)) {
      if(k==i){
        theta_last[k] <- constant_theta
      }else if(k>i){
        theta_last[k] <- optimization$solution[k-1]   
      }else{
        theta_last[k] <- optimization$solution[k]} 
    }
  }
  if(iter_counter >= N_samples){
    b_exit <- 1
  }else if(current_score > threshold){
    b_exit <- 2
  }
  cat("Backward search ended after ", iter_counter, "iterations.", "\n")
  
  results_df <- rbind(backward_results_df, forward_results_df, c(thetas[i], global_optimum))
  results_df <- results_df[order(results_df[,1]),]
  results_df <- results_df[complete.cases(results_df),]
  if(create_txt) sink()
  return(list("plik"=results_df, "b_exit"=b_exit, "f_exit"=f_exit))
}