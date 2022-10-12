
#' Winners curse
#'
#' Apply the winners curse algorithm to the final pooled OLS
#' @param pooled_ols the result of the final pooled OLS
#' @param pool_ids the names of the pool id dummy columns which coefficients are in pooled_ols
#' @param alpha is the alpha used in the winners curse algorithm, default to 0.05
#' @param beta is the beta used in the winners curse algorithm, default to NULL\cr
#' if NULL, beta will be computed with beta = alpha / 10
#' @return returns a list with the median unbiased estimate of the best effect and the alpha confidence interval lower and upper bounds
#' @export


winners_curse <- function(pooled_ols,pool_ids,alpha=0.05,beta=NULL){
  
  if (is.null(beta) | beta>alpha){
    beta=alpha/10
  }
  
  cat("Computing the winners curse effect with alpha =",alpha," and beta =",beta,"\n")
  
  pooled_effects <- pooled_ols$coefficients[pool_ids]
  pooled_pval <-pooled_ols$p.value[pool_ids]
  ntreat <- length(pool_ids) + 1
  
  first_pool <- names(which.max(pooled_effects))
  second_pool <- names(which.max(pooled_effects[pooled_effects!=pooled_effects[first_pool]]))
  
  first_scaled_effect <- sqrt(nobs(pooled_ols)) * pooled_effects[first_pool]
  second_scaled_effect <- max(0, sqrt(nobs(pooled_ols)) * pooled_effects[second_pool]) #in case
  
  first_pool_var <- nobs(pooled_ols) * (pooled_ols$std.error[first_pool])^2
  
  hybrid_results_scaled = get_hybrid_estimate(first_scaled_effect, second_scaled_effect, first_pool_var, ntreat, alpha, beta)
  hybrid_results <- (1/sqrt(nobs(pooled_ols))) * hybrid_results_scaled
  return(hybrid_results)
}


#' Winners curse
#'
#' Apply the winners curse algorithm to the final pooled OLS
#' @param first_Y is the scaled best effect
#' @param second Y is the sclaed second best effect
#' @param first_Y_var is the sclaed best effect variation
#' @param ntreat is the number of estimated effects
#' @param alpha is the alpha used in the winners curse algorithm, default to 0.05
#' @param beta is the beta used in the winners curse algorithm
#' @return returns a list with the median unbiased estimate of the best effect and the alpha confidence interval lower and upper bounds
#' @export


get_hybrid_estimate <- function(first_Y, second_Y, first_Y_var, ntreat, alpha, beta) {
  #Projected joint confidence rectangle
  c_beta = qnorm((1-beta/2)^(1/ntreat))
  
  #Hybrid interval
  gamma <- (alpha - beta)/(1- beta)
  
  median = get_hybrid_mu_alpha(first_Y = first_Y, second_Y = second_Y, c_beta = c_beta, first_Y_var = first_Y_var , alpha = 0.50)
  CI_lb = get_hybrid_mu_alpha(first_Y = first_Y, second_Y = second_Y, c_beta = c_beta, first_Y_var = first_Y_var , alpha = gamma/2)
  CI_ub = get_hybrid_mu_alpha(first_Y = first_Y, second_Y = second_Y, c_beta = c_beta, first_Y_var = first_Y_var , alpha = 1 - (gamma/2))
  
  return(c(median = median ,CI_lb=CI_lb, CI_ub=CI_ub))
}

#' Compute the hybrid mu_alpha
#'
#' Solves the equation F_TN^H (Y(theta^hat),mu_alpha,theta^hat,Z_theta^hat) = 1-alpha
#' @param first_Y is the scaled best effect
#' @param second Y is the sclaed second best effect
#' @param c_beta is the coefficient computed for the joint confidence interval of level beta
#' @param first_Y_var is the sclaed best effect variation
#' @param alpha is the alpha used in the winners curse algorithm, default to 0.05
#' @return returns the solution mu_alpha
#' @export

get_hybrid_mu_alpha <- function(first_Y, second_Y, c_beta, first_Y_var, alpha) {
  function_to_solve = function (x) ptmvnorm_extended(q = first_Y, mu = x, sigma = first_Y_var, lb=max(second_Y,x - c_beta*sqrt(first_Y_var)), ub=x + c_beta*sqrt(first_Y_var))- (1-alpha)
  mu_alpha = stats::uniroot(function_to_solve, lower = floor(first_Y - c_beta * sqrt(first_Y_var)), upper = ceiling(first_Y + c_beta * sqrt(first_Y_var)), extendInt = "yes")[1] %>% unname()
  return(mu_alpha[[1]])
}


#' Extended truncated normal
#'
#' Extends the truncated normal functions in case the point is outside bounds
#' @param q the point where we evaluate the function
#' @param mu the mean of the normal
#' @param sigma the standard deviation of the normal
#' @param lb the lower truncation point
#' @param ub the upper truncation point
#' @return returns the extended truncated normal function
#' @export


ptmvnorm_extended <- function(q, mu, sigma, lb, ub) {
  if (q >= ub)  { return(1)  }
  if (q <= lb)  { return(0)  }
  return(TruncatedNormal::ptmvnorm(q, mu, sigma, lb, ub))
  
} 