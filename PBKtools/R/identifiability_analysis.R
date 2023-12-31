#' Identifiability Analysis
#'
#'This is a function to implement the framework for identifiability analysis
#'as described in the publication of Raue et al., 2009 entitled as "Structural and practical
#'identifiability analysis of partially observed dynamical models by exploiting
#'the profile likelihood" (doi: https://doi.org/10.1093/bioinformatics/btp358).
#' 
#' @param obj_f             function that returns the value of the objective function
#' @param thetas            vector containing the optimal values of parameters 
#' @param thetas_names      vector of characters with the names of the parameters
#' @param data_df           dataframe containing the data used in the obj_f 
#' @param errors_df         dataframe with the meassured experimental errors (or SD values) 
#' @param lb                vector with the lower bounds of the "thetas" parameters
#' @param ub                vector with the upper bounds of the "thetas" parameters
#' @param N_samples         integer defining the number of samples taken around 
#' the theta optimal value (N samples will be taken from each side of the theta)
#' @param alpha             probability of chi-squared to estimate the quantile (it is the "p" variable of qchisq() function)
#' @param df                degrees of freedom of qchisq()
#' @param q                 a variable used in the estimation of the adaptive step (see Raue et al., 2009)
#' @param global_optimum    scalar defining the global minimum of the objective function
#' @param min_step_coef     coefficient defining the minimum permitted step
#' @param max_step_coef     coefficient defining the maximum permitted step  
#' @param N_cores           integer defining the number of cores for the parallel processing
#' @param constant_params   vector with any extra constant parameters of the model and their names
#' @param break_at_bounds   logical; if TRUE the the sampling of the parameter stops because the bounds were exceeded
#' @param opts              list with the options selected for the minimization 
#' of the objective function (check the nloptr package for more details)
#' @param opts_theta_step   list with the options selected for the estimation 
#' of the adaptive step (check the nloptr package for more details)
#' @param create_txt        logical; if TRUE a txt file will be created at the current working directory 
#' to save the samples of the parameters and the corresponding values of the objective function
#' @export
Identifiability_analysis <- function(obj_f, thetas, thetas_names, data_df, errors_df,
                                     lb, ub,
                                     N_samples = 50,
                                     alpha = 0.95, df = 1,
                                     q = 0.5,
                                     global_optimum, 
                                     min_step_coef = 1e-03, max_step_coef = 0.2,
                                     N_cores,
                                     constant_params = NULL,
                                     exported_to_cluster = NULL,
                                     break_at_bounds = FALSE,
                                     opts = list("algorithm" = "NLOPT_LN_NELDERMEAD",
                                                 "xtol_rel" = 1e-06, 
                                                 "ftol_rel" = 1e-06,
                                                 "ftol_abs" = 0.0,
                                                 "xtol_abs" = 0.0 ,
                                                 "maxeval" = 300,
                                                 "print_level" = 0),
                                     opts_theta_step = list("algorithm" = 'NLOPT_LN_SBPLX',
                                                            "xtol_rel_step" = 1e-05, 
                                                            "ftol_rel_step" = 1e-05,
                                                            "ftol_abs_step" = 0.0,
                                                            "xtol_abs_step" = 0.0 ,
                                                            "maxeval_step" = 50,
                                                            "print_level_step" = 0),
                                     create_txt = TRUE){
  
  # Number of parameters tested in the identifiability analysis
  N_parameters <- length(thetas)
  # prepare the input for parallel processing
  X <- vector("list", N_parameters)
  for(i in 1:N_parameters){
    X[[i]] <- list(index=i, thetas=thetas, thetas_names=thetas_names,
                   constant_params = constant_params,
                   data_df = data_df,
                   errors_df = errors_df,
                   lb=lb, ub=ub, N_samples=N_samples, 
                   alpha=alpha, df=df,
                   global_optimum = global_optimum,
                   q = q,
                   max_step_coef = max_step_coef,
                   min_step_coef = min_step_coef,
                   break_at_bounds = break_at_bounds,
                   opts = opts,
                   opts_theta_step=opts_theta_step,
                   create_txt = create_txt)
  }
  
  # parallel_func is a wrapper of the profile_likelihood function in order to 
  # take all the input variable of the profile_likelihood as a list and un-list
  # the to provide them to profile_likelihood. It just returns the output
  # of the profile_likelihood function.
  parallel_func <- function(X){
    with(as.list(X),{
      profile_likelihood(obj_f, 
                         i=index,
                         thetas,
                         thetas_names, 
                         constant_params,
                         data_df,
                         errors_df,
                         lb, ub, 
                         N_samples,
                         alpha,
                         df,
                         q, 
                         global_optimum, 
                         min_step_coef,
                         max_step_coef,
                         break_at_bounds,
                         opts,
                         opts_theta_step,
                         create_txt)
    })
  }
  start.time <- Sys.time()
  # Set up the cluster.
  cluster <- makeCluster(N_cores)
  # Export to the cluster any function or parameter that the obj_f needs to work.
  clusterExport(cl=cluster, c(names(exported_to_cluster),"obj_f", "profile_likelihood"))
  output <- parLapply(cluster, X, parallel_func)
  # Terminate the cluster.
  stopCluster(cluster)
  total.duration <- Sys.time() - start.time
  
  # Collect all the results of interest and present them in a dataframe.
  results_df <- data.frame(matrix(NA, ncol = 6, nrow = length(thetas)))
  rownames(results_df) <- thetas_names
  colnames(results_df) <- c("Non-Identifiability" , "Optimal", "Lower_Bound", "Upper_Bound", "Exit_code_B", "Exit_code_F")
  for (i in 1:length(thetas)) {
    pl_results <- output[[i]]$'plik'
    confidence_intervals <- ci_estimation(pl_results,  alpha = alpha, df=df,
                                          theta_hat=thetas[[i]], global_optimum=global_optimum,
                                          lb=lb[i], ub=ub[i])
    results_df$Lower_Bound[i] <- confidence_intervals$Lower_bound
    results_df$Upper_Bound[i] <- confidence_intervals$Upper_bound
    results_df$Optimal[i] <- thetas[[i]]
    results_df$Exit_code_B[i] <- output[[i]]$b_exit
    results_df$Exit_code_F[i] <- output[[i]]$f_exit
    if(confidence_intervals$Lower_bound == -Inf & confidence_intervals$Upper_bound == Inf){
      results_df$'Non-Identifiability'[i] <- "Structural"
    }else if(confidence_intervals$Lower_bound == -Inf & confidence_intervals$Upper_bound != Inf |
             confidence_intervals$Lower_bound != -Inf & confidence_intervals$Upper_bound == Inf){
      results_df$'Non-Identifiability'[i] <- "Practical"
    }else{
      results_df$'Non-Identifiability'[i] <- "Identifiable"
    }
  }
  
  return(list("Likelihood_profiles" = output, "Identiafiability_Analysis" = results_df, "Total_duration" = total.duration))
}