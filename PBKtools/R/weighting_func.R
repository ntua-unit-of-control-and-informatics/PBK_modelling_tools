#' Estimating Weights for data points. 
#' 
#' This custom function is designed to assign weight values to data points within 
#' a dataset. Its primary purpose is to identify turning points or points adjacent 
#' to turning points, as well as the first and last points in the data. These 
#' particular points are considered more significant within the dataset, and as a result, 
#' the function assigns higher weights to these points.
#' Additionally, the function applies higher weight values to points that are 
#' farther apart from each other. This weighting approach aims to enhance 
#' the model's performance by prioritizing crucial data points and appropriately 
#' accounting for their relative distances.
#' 
#' @param  df a dataframe whose the 1st column contains the values for the 
#' independent variable (x) and the rest columns contain the data for the dependent variables (y).
#' @param a_first the value of the factor to estimate the weight for the first point
#' of the data
#' @param a_last the value of the factor to estimate the weight for the last point
#' of the data
#' @param a_tp value of the factor to estimate the weight for turning points (or the
#' points that are next to the turning points)
#' 
#' @return A list with the weight values for each dependent variable of the dataset.
#' @export
Weighting.func <- function(df, a_first=2, a_last=2, a_tp=2){
  x_values <- df[,1] # Take the x values 
  y_values <- df[,-1] # Take th y values
  
  # Number of x values
  Nx <- dim(df)[1] 
  # Number of y columns (subtract one column because it's the x values)
  Ny <- dim(y_values)[2]
  
  # Transform the data df to a list of dataframes. Each dataframe will have 
  # the x values and y values of a tissue.
  sub_list <- list()
  
  # Create a list to save the weights
  weights_list <- list()
  for (i in 1:Ny) { # Loop over the y outputs
    
    # Create a 2-column dataframe with x and y values
    sub_list[[i]] <- cbind(x_values, y_values[,i])
    names(sub_list)[i] <- colnames(y_values)[i]
    colnames(sub_list[[i]]) <- colnames(df)[c(1,i+1)]
    
    weights_list <- append(weights_list, NA)
    # Select a df from the sub_list to calulcate its weights
    sub_df <- sub_list[[i]]
    N <- dim(sub_df)[1] # Number of x values of sub_df
    
    # Estimate the weight of the 1st point
    weights_list[[i]][1] <- a_first*(abs(sub_df[1,1] - sub_df[2,1]) + abs(sub_df[1,2] - sub_df[2,2]))/(abs(sub_df[N,1] - sub_df[1,1]) + abs(sub_df[N,2] - sub_df[1,2]))
    # Estimate the weight of the Last point
    weights_list[[i]][N] <- a_last*(abs(sub_df[N,1] - sub_df[N-1,1]) + abs(sub_df[N,2] - sub_df[N-1,2]))/(abs(sub_df[N,1] - sub_df[1,1]) + abs(sub_df[N,2] - sub_df[1,2]))
    
    # next_point (binary variable) indicates if j+1 point is a turning point
    next_point <- 0
    prev_point<- 0
    for (j in 2:(N-2)) {
      if(next_point){ # j is defined as turning point in previous loop
        a <- a_tp
        weights_list[[i]][j] <- (abs(sub_df[j,1] - sub_df[j-1,1]) + abs(sub_df[j,1] - sub_df[j+1,1]) +
                                   a*(abs(sub_df[j,2] - sub_df[j-1,2]) + abs(sub_df[j,2] - sub_df[j+1,2])))/(abs(sub_df[N,1] - sub_df[1,1]) + abs(sub_df[N,2] - sub_df[1,2]))
        
        # Check if the j+1 point is turning point to define next_point value for the next loop
        if((sub_df[j+1,2]<sub_df[j,2] & sub_df[j+1,2]<sub_df[j+2,2]) | (sub_df[j+1,2]>sub_df[j,2] & sub_df[j+1,2]>sub_df[j+2,2])){
          prev_point <- T
          next_point <- T
        }else{
          # prev_point is used to identify if the point of the next loop (which is not a turning point) is after
          # a turning point
          prev_point <- T
          next_point <- F
        }
      }else{ # j is not defined as turning point in previous loop
        # Check if the j+1 point is turning point to define next_point value for the next loop
        # or check if the previous point was a turning point
        if((sub_df[j+1,2]<sub_df[j,2] & sub_df[j+1,2]<sub_df[j+2,2]) | (sub_df[j+1,2]>sub_df[j,2] & sub_df[j+1,2]>sub_df[j+2,2])){
          a <- a_tp
          prev_point <- F
          next_point <- T
        }else{
          a <- ifelse(prev_point, a_tp, 1)
          prev_point <- F
          next_point <- F
        }
        weights_list[[i]][j] <- (abs(sub_df[j,1] - sub_df[j-1,1]) + abs(sub_df[j,1] - sub_df[j+1,1]) +
                                   a*(abs(sub_df[j,2] - sub_df[j-1,2]) + abs(sub_df[j,2] - sub_df[j+1,2])))/(abs(sub_df[N,1] - sub_df[1,1]) + abs(sub_df[N,2] - sub_df[1,2]))
      }
    }
    
    # Define a value and the weight value for the second from the end data point (j=N-1).
    a <- ifelse(next_point | prev_point, a_tp, 1)
    j <- N-1
    weights_list[[i]][j] <- (abs(sub_df[j,1] - sub_df[j-1,1]) + abs(sub_df[j,1] - sub_df[j+1,1]) +
                               a*(abs(sub_df[j,2] - sub_df[j-1,2]) + abs(sub_df[j,2] - sub_df[j+1,2])))/(abs(sub_df[N,1] - sub_df[1,1]) + abs(sub_df[N,2] - sub_df[1,2]))
  }
  names(weights_list) <- colnames(y_values)
  return(weights_list)
}
