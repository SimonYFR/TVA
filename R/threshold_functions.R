#################################################################
#
#------ This file contains the functions related to thresholding 
#------ (pval plots, grids or suggesting cutoffs)
#
#################################################################

#' Plot p-values in the pval one-step elimination
#'
#' Plot the ordered p-values in the one-step elimination. \cr
#' This allows to choose a thresholds according to a targeted support size. \cr
#' This is done by simulating a one-step p-value elimination.
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w is the column name of the weights
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns the plot of the ordered p-values, allowing to see the corresponding size of support for each p-values cutoff in the pval one-step elimination
#' @export
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' fes = c('fes_1')
#' y = 'outcome'
#' w = 'weights'
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' F1 = c(0,1,0,0,0,1,0,1,0,0)
#' Y  = c(5,4,3,5,4,5,4,2,3,2)
#' W  = c(1,1,1,2,1,2,2,1,1,2)
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3, fes_1 = F1, outcome = Y, weights=W)
#' plot_pval_OSE(data,arms,y, fes,w,FALSE)


plot_pval_OSE <- function(data, arms, y, fes=c(), w=NULL, compare_to_zero=FALSE){
  check = check_inputs_integrity(data, arms, y, fes, 1, w, 'pval_OSE', compare_to_zero)
  if (!check$integrity){
    stop(check$message)
  }
  
  prepared_data = prepare_data(data,arms, y, fes,w,compare_to_zero)
  X = prepared_data$X
  data = prepared_data$data
  arms = prepared_data$arms
  variables = prepared_data$variables
  marginals_colnames = prepared_data$marginals_colnames
  thresholds = pval_OSE(X,y,variables,1)$thresholds
  
  thresholds_OSE = thresholds[which(names(thresholds) %in% marginals_colnames)] %>% sort() %>% data.frame()
  thresholds_OSE = thresholds_OSE %>% setNames(.,c('threshold'))
  thresholds_OSE$size_of_support = c(1:nrow(thresholds_OSE))
  
  plot = ggplot2::ggplot(data=thresholds_OSE, ggplot2::aes(x=size_of_support, y=threshold)) +
    ggplot2::geom_line(linetype = "dashed") +
    ggplot2::geom_point() +
    ggplot2::scale_y_continuous(trans='log10') +
    ggplot2::xlab("Marginal support size") +
    ggplot2::ylab("Threshold") +
    ggplot2::theme_bw()
  
  return(plot)
}

#' Plot p-values in the pval multi-step elimination
#'
#' Plot the ordered p-values in the multi-step elimination. \cr
#' This allows to choose a thresholds_OSE according to a targeted support size. \cr
#' This is done by simulating a multi-step p-value elimination.
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w is the column name of the weights
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns the plot of the ordered p-values, allowing to see the corresponding size of support for each p-values cutoff in the pval multi-step elimination
#' @export
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' fes = c('fes_1')
#' y = 'outcome'
#' w = 'weights'
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' F1 = c(0,1,0,0,0,1,0,1,0,0)
#' Y  = c(5,4,3,5,4,5,4,2,3,2)
#' W  = c(1,1,1,2,1,2,2,1,1,2)
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3, fes_1 = F1, outcome = Y, weights=W)
#' plot_pval_MSE(data,arms,y, fes,w,FALSE)

plot_pval_MSE <- function(data,arms,y, fes=c(),w=NULL,compare_to_zero=FALSE){
  check = check_inputs_integrity(data, arms, y, fes, 1, w, 'pval_OSE', compare_to_zero)
  
  if (!check$integrity){
    stop(check$message)
  }
  
  prepared_data = prepare_data(data,arms,y, fes,w,compare_to_zero)
  X = prepared_data$X
  data = prepared_data$data
  arms = prepared_data$arms
  variables = prepared_data$variables
  marginals_colnames = prepared_data$marginals_colnames
  
  pval_MSE = pval_MSE(X,y,variables,0)
  max_pvals = (pval_MSE$pvals)
  max_pvals = max_pvals[names(max_pvals) %in% marginals_colnames]
  thresholds_MSE = data.frame( threshold = cummin(max_pvals), max_pval = max_pvals, size_of_support = rev(c(1:length(max_pvals))))
  
  plot = ggplot2::ggplot(data=thresholds_MSE, ggplot2::aes(x=size_of_support)) +
    ggplot2::geom_line(ggplot2::aes(y = max_pval), color="black", linetype="dashed") +
    ggplot2::geom_line(ggplot2::aes(y = threshold), color="steelblue") +
    ggplot2::geom_point(ggplot2::aes(y = max_pval), color="black")+
    ggplot2::scale_y_continuous(trans='log10') + 
    ggplot2::xlab("Marginal support size") +
    ggplot2::ylab("Maximum p-value") +
    ggplot2::theme_bw()
  
  return(plot)
}

#' Plot betas in the beta one-step elimination
#'
#' Plot the ordered absolute values of betas in the beta one-step elimination. \cr
#' This allows to choose a beta cutoff according to a targeted support size. \cr
#' This is done by simulating a one-step beta elimination.
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w is the column name of the weights
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns the plot of the ordered betas, allowing to see the corresponding size of support for each beta cutoff in the beta one-step elimination
#' @export
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' fes = c('fes_1')
#' y = 'outcome'
#' w = 'weights'
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' F1 = c(0,1,0,0,0,1,0,1,0,0)
#' Y  = c(5,4,3,5,4,5,4,2,3,2)
#' W  = c(1,1,1,2,1,2,2,1,1,2)
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3, fes_1 = F1, outcome = Y, weights=W)
#' plot_beta_OSE(data,arms,y,fes,w,FALSE)


plot_beta_OSE <- function(data,arms,y, fes=c(),w=NULL,compare_to_zero=FALSE){
  check = check_inputs_integrity(data, arms, y, fes, 1, w, 'pval_OSE', compare_to_zero)
  if (!check$integrity){
    stop(check$message)
  }
  
  prepared_data = prepare_data(data,arms,y,fes,w,compare_to_zero)
  X = prepared_data$X
  data = prepared_data$data
  arms = prepared_data$arms
  variables = prepared_data$variables
  marginals_colnames = prepared_data$marginals_colnames
  
  beta = beta_OSE(X,y,variables,0)$beta %>% abs()
  
  beta_OSE = beta[which(names(beta) %in% marginals_colnames)] %>% sort(.,decreasing=TRUE) %>% data.frame()
  beta_OSE = beta_OSE %>% setNames(.,c('beta'))
  beta_OSE$size_of_support = c(1:nrow(beta_OSE))
  
  plot = ggplot2::ggplot(data=beta_OSE, ggplot2::aes(x=size_of_support, y=beta)) +
    ggplot2::geom_line(linetype = "dashed")+
    ggplot2::geom_point() + 
    ggplot2::xlab("Marginal support size") +
    ggplot2::ylab("Beta") +
    ggplot2::theme_bw()
  
  return(plot)
}

#' Plot cutoff vs marginal support size
#'
#' Plot cutoff vs marginal support size for each method. \cr
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes (optional) is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w (optional) is the column name of the weights
#' @param estim_func (optional) is the estimation method (pval_OSE, pval_MSE or beta_OSE), pval_OSE by default
#' @param compare_to_zero (optional) is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns the plot of the ordered betas, allowing to see the corresponding size of support for each beta cutoff in the beta one-step elimination
#' @export
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' fes = c('fes_1')
#' y = 'outcome'
#' w = 'weights'
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' F1 = c(0,1,0,0,0,1,0,1,0,0)
#' Y  = c(5,4,3,5,4,5,4,2,3,2)
#' W  = c(1,1,1,2,1,2,2,1,1,2)
#' estim_func = "pval_MSE"
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3, fes_1 = F1, outcome = Y, weights=W)
#' plot_cutoff_vs_support_size(data=data, arms=arms, y=y, fes=fes, w=w, estim_func=estim_func, compare_to_zero=FALSE)
#' 

plot_cutoff_vs_support_size <- function(data,arms,y, fes=c(),w=NULL, estim_func="pval_OSE", compare_to_zero=FALSE){
  if (estim_func=='beta_OSE'){
    plot = plot_beta_OSE(data,arms,y, fes=fes,w=w,compare_to_zero=compare_to_zero)
  }else if (estim_func=='pval_OSE'){
    plot = plot_pval_OSE(data,arms,y, fes=fes,w=w,compare_to_zero=compare_to_zero)
  }else if (estim_func=='pval_MSE'){
    plot = plot_pval_MSE(data,arms,y, fes=fes,w=w,compare_to_zero=compare_to_zero)
  }else{
    return("estim_func must be either beta_OSE, pval_OSE or pval_MSE")
  }
  return(plot)
}

#' Pval cutoff grid
#'
#' Simulate many pval one-step or multiple-step elimination for different pval cutoffs values.\cr
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w is the column name of the weights
#' @param estim_func is the support estimation function we want to use, either pval_OSE or pval_MSE
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns a dataframe with columns:
#' * 'pval_cutoff': all the p-values
#' * 'marginal_support_size' : the corresponding support size one would get for each of these p-values
#' * 'number_of_pools': the corresponding number of pools one would get for such a support
#' @export
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' fes = c('fes_1')
#' y = 'outcome'
#' w = 'weights'
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' F1 = c(0,1,0,0,0,1,0,1,0,0)
#' Y  = c(5,4,3,5,4,5,4,2,3,2)
#' W  = c(1,1,1,2,1,2,2,1,1,2)
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3, fes_1 = F1, outcome = Y, weights=W)
#' grid_pval(data=data,arms=arms,y=y,fes=fes,w=w,compare_to_zero=FALSE)


grid_pval <- function(data,arms,y,fes=c(),w=NULL,estim_func='pval_OSE', compare_to_zero=FALSE, clusters=NULL){
  
  check = check_inputs_integrity(data, arms, y, fes, 1, w, estim_func, compare_to_zero, clusters)
  
  if (!check$integrity){
    stop(check$message)
  }
  
  prepared_data = prepare_data(data,arms,y, fes,w,compare_to_zero)
  X = prepared_data$X
  data = prepared_data$data
  arms = prepared_data$arms
  variables = prepared_data$variables
  marginals_colnames = prepared_data$marginals_colnames
  
  f = get(estim_func)
  thresholds = f(X,y,variables,0)$thresholds 
  
  m_thresholds = thresholds[marginals_colnames] %>% sort(decreasing = FALSE)
  
  gridval = floor(m_thresholds / 10^(floor(log(m_thresholds, base = 10))-1))/10 * 10^(floor(log(m_thresholds, base = 10))) #round the pvals
  gridval = gridval[1: (length(gridval)/3) %>% ceiling()] %>% unname() %>% unique()
  
  marginal_support_sizes = c()
  number_of_pools = c()
  differ_from_zero = c()
  rsqr = c()
  adj_rsqr = c()
  
  for (threshold in gridval){
    total_support = names(thresholds[which(thresholds<=threshold)]) #compute the total_support
    marginal_support = intersect(total_support,marginals_colnames) %>% sort() #take the marginal support
    
    pooled_data = pool_data(data,arms,marginal_support,compare_to_zero)
    
    pool_ids = paste("pool_id",c(1:max(pooled_data$pool_id)),sep="_")
    fes_support = sort(intersect(total_support,fes))
    
    if (length(marginal_support)>0){
      pooled_ols = get_pooled_ols(pooled_data,y,fes_support,w,pool_ids,clusters) #do pooled OLS
      
      ols_coefs = pooled_ols$coefficients 
      pools_coefs = ols_coefs[grep("pool_id_", ols_coefs %>% names, value = TRUE)] #take pools_coefficients
      two_bests = (pools_coefs %>% sort(.,decreasing=TRUE))[1:2] %>% names() #take the two bests
      two_bests_differ_from_zero = all(pooled_ols$p.value[two_bests] < 0.05)
      ols_rsqr = pooled_ols$r.squared
      ols_adj_rsqr = pooled_ols$adj.r.squared
    }else{
      two_bests_differ_from_zero = FALSE
      ols_rsqr = 0
      ols_adj_rsqr = 0
    }
    
    differ_from_zero = c(differ_from_zero, two_bests_differ_from_zero)
    number_of_pools = c(number_of_pools, pooled_data$pool_id %>% max() +1)
    marginal_support_sizes = c(marginal_support_sizes, length(marginal_support))
    rsqr = c(rsqr, ols_rsqr )
    adj_rsqr = c(adj_rsqr, ols_adj_rsqr)
  }
  
  grid = data.frame( threshold = gridval
                     ,marginal_support_size= marginal_support_sizes
                     ,number_of_pools=number_of_pools 
                     ,differ_from_zero=differ_from_zero
                     ,rsqr = rsqr
                     ,adj_rsqr = adj_rsqr)
  return(grid)
}


#' Lambdas cutoff grid
#'
#' Simulate many beta one-step elimination for different lambdas cutoffs values.\cr
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w is the column name of the weights
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns a dataframe with columns:
#' * 'lambda_cutoff': all the lambdas
#' * 'marginal_support_size' : the corresponding support size one would get for each of these p-values
#' * 'number_of_pools': the corresponding number of pools one would get for such a support
#' @export
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' fes = c('fes_1')
#' y = 'outcome'
#' w = 'weights'
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' F1 = c(0,1,0,0,0,1,0,1,0,0)
#' Y  = c(5,4,3,5,4,5,4,2,3,2)
#' W  = c(1,1,1,2,1,2,2,1,1,2)
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3, fes_1 = F1, outcome = Y, weights=W)
#' grid_lambda(data=data,arms=arms,y=y,fes=fes,w=w,compare_to_zero=FALSE)


grid_lambda <- function(data,arms,y,fes=c(),w=NULL,compare_to_zero=FALSE, clusters=NULL){
  estim_func = "beta_OSE"
  check = check_inputs_integrity(data, arms, y, fes, 1, w, estim_func, compare_to_zero, clusters)
  
  if (!check$integrity){
    stop(check$message)
  }
  
  prepared_data = prepare_data(data,arms,y, fes,w,compare_to_zero)
  X = prepared_data$X
  data = prepared_data$data
  arms = prepared_data$arms
  variables = prepared_data$variables
  marginals_colnames = prepared_data$marginals_colnames
  
  f = get(estim_func)
  
  thresholds = f(X,y,variables,0)$beta %>% abs()
  m_thresholds = thresholds[marginals_colnames] %>% sort(decreasing = TRUE) 
  
  gridval = floor(m_thresholds / 10^(floor(log(m_thresholds, base = 10))-1))/10 * 10^(floor(log(m_thresholds, base = 10))) #round the pvals
  gridval = gridval[1: (length(gridval)/3) %>% ceiling()] %>% unname() %>% unique()
  
  marginal_support_sizes = c()
  number_of_pools = c()
  differ_from_zero = c()
  rsqr = c()
  adj_rsqr = c()
  
  for (threshold in gridval){
    total_support = names(thresholds[which(thresholds>=threshold)]) #compute the total_support
    marginal_support = intersect(total_support,marginals_colnames) %>% sort() #take the marginal support
    
    pooled_data = pool_data(data,arms,marginal_support,compare_to_zero)
    
    pool_ids = paste("pool_id",c(1:max(pooled_data$pool_id)),sep="_")
    fes_support = sort(intersect(total_support,fes))
    
    if (length(marginal_support)>0){
      pooled_ols = get_pooled_ols(pooled_data,y,fes_support,w,pool_ids,clusters) #do pooled OLS
      
      ols_coefs = pooled_ols$coefficients 
      pools_coefs = ols_coefs[grep("pool_id_", ols_coefs %>% names, value = TRUE)] #take pools_coefficients
      two_bests = (pools_coefs %>% sort(.,decreasing=TRUE))[1:2] %>% names() #take the two bests
      two_bests_differ_from_zero = all(pooled_ols$p.value[two_bests] < 0.05)
      ols_rsqr = pooled_ols$r.squared
      ols_adj_rsqr = pooled_ols$adj.r.squared
    }else{
      two_bests_differ_from_zero = FALSE
      ols_rsqr = 0
      ols_adj_rsqr = 0
    }
    
    differ_from_zero = c(differ_from_zero, two_bests_differ_from_zero)
    number_of_pools = c(number_of_pools, pooled_data$pool_id %>% max() +1)
    marginal_support_sizes = c(marginal_support_sizes, length(marginal_support))
    rsqr = c(rsqr, ols_rsqr )
    adj_rsqr = c(adj_rsqr, ols_adj_rsqr)
  }
  
  grid = data.frame( threshold = gridval
                     ,marginal_support_size= marginal_support_sizes
                     ,number_of_pools=number_of_pools 
                     ,differ_from_zero=differ_from_zero
                     ,rsqr = rsqr
                     ,adj_rsqr = adj_rsqr)
  return(grid)
}



#' Cutoff grid
#'
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w is the column name of the weights
#' @param estim_func is the support estimation function we want to use, either pval_OSE or pval_MSE
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns a dataframe with columns:
#' * 'pval_cutoff': all the p-values
#' * 'marginal_support_size' : the corresponding support size one would get for each of these p-values
#' * 'number_of_pools': the corresponding number of pools one would get for such a support
#' @export
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' fes = c('fes_1')
#' y = 'outcome'
#' w = 'weights'
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' F1 = c(0,1,0,0,0,1,0,1,0,0)
#' Y  = c(5,4,3,5,4,5,4,2,3,2)
#' W  = c(1,1,1,2,1,2,2,1,1,2)
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3, fes_1 = F1, outcome = Y, weights=W)
#' grid_pval(data=data,arms=arms,y=y,fes=fes,w=w,compare_to_zero=FALSE)

grid_cutoff <- function(data,arms,y,fes=c(),w=NULL, estim_func='pval_OSE', compare_to_zero=FALSE, clusters=NULL){
  if (estim_func=='beta_OSE'){
    grid = grid_lambda(data=data,arms=arms,y=y,fes=fes,w=w, compare_to_zero=compare_to_zero, clusters=clusters)
  }else if ((estim_func=='pval_OSE')|(estim_func=='pval_MSE')){
    grid = grid_pval(data=data,arms=arms,y=y,fes=fes,w=w, estim_func=estim_func,compare_to_zero=compare_to_zero, clusters=clusters)
  }else{
    return("estim_func must be either beta_OSE, pval_OSE or pval_MSE")
  }
  return(grid)
}


#' Find elbow
#'
#' Find the elbow of a set of points.
#' @param X is a numeric vector
#' @param Y is a numeric vector, same size as X
#' @return returns the value of X corresponding to the elbow of  the curve
#' @export

elbow <- function(X,Y){
  smoothed_Y = smooth.spline(X, Y, spar=0.65)$y
  d1 <- diff(smoothed_Y) / diff(X) # first derivative
  d2 <- diff(d1) / diff(X[-1]) # second derivative
  
  return(X[which(d2==min(d2))+1])
}


#' Suggest one cutoff that could be used the in pval OSE, pval MSE or beta OSE
#'
#' Suggest one cutoff that could be used the in pval OSE, pval MSE or beta OSE.\cr
#' This is done in two steps. First, we take the p-value corresponding to a target, which can be provided by the user, or which corresponds to the elbow of the R-squared vs support size curve. \cr
#' Then, we look for p-values arround this target and check if some of them produce a final pooled OLS with two best effects that are statistically different from 0.
#' @param support_size_target (optional) is the prior we might have on the number of marginals in the support. It is NULL by default. \cr
#' If NULL, the code will take the support size that corresponds to the elbow of the R-squared vs support size curve.
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w is the column name of the weights
#' @param estim_func can be equal to "pval_OSE", "pval_MSE" or "beta_OSE"
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns a pval cutoff suggestion with the according number of pools it will produce
#' @export
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' fes = c('fes_1')
#' y = 'outcome'
#' w = 'weights'
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' F1 = c(0,1,0,0,0,1,0,1,0,0)
#' Y  = c(5,4,3,5,4,5,4,2,3,2)
#' W  = c(1,1,1,2,1,2,2,1,1,2)
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3, fes_1 = F1, outcome = Y, weights=W)
#' suggest_cutoff(data=data,arms=arms,y=y,fes=fes,w=w,compare_to_zero=FALSE)

suggest_cutoff <- function(data,arms,y,support_size_target=NULL, fes=c(),w=NULL,estim_func='pval_OSE',compare_to_zero=FALSE, clusters=NULL){
  check = check_inputs_integrity(data, arms, y, fes, cutoff = NULL, w, estim_func, compare_to_zero, clusters)
  
  if (!check$integrity){
    stop(check$message)
  }

  grid = grid_cutoff(data=data,arms=arms,y=y,fes=fes,w=w,estim_func=estim_func,compare_to_zero=compare_to_zero, clusters=clusters)

  if (is.null(support_size_target)){
    elbow = elbow(grid$marginal_support_size, grid$rsqr)
    cat('Elbow is ',elbow,'\n')
    support_size_target = grid$marginal_support_size[which.min(abs(grid$marginal_support_size - elbow))]
  }else{
    support_size_target = grid$marginal_support_size[which.min(abs(grid$marginal_support_size - support_size_target))]
  }
  
  cat('Target is ',support_size_target,'\n')
  
  optimums = grid[(grid$marginal_support_size %>% dplyr::between(.,round(support_size_target/2),support_size_target*2)) & (grid$differ_from_zero),]
  
  if (nrow((optimums))==0){
    best = grid[grid$marginal_support_size == support_size_target,]
  }else{
    optimums$distance_from_target = abs(optimums$marginal_support_size - support_size_target)
    best = optimums[order(optimums$marginal_support_size <= support_size_target, optimums$distance_from_target, decreasing=FALSE),][1,]
  }
  
  return(best)
  
  
}


#################################################################
# END
#################################################################

