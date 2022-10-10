
#' pval multi-step elimination
#'
#' Estimate support with p-value multi-step elimination
#' @param X is the regression matrix with the marginals, the fixed effects, the outcome of interest
#' @param y is the outcome of interest
#' @param variables is a vector containing all the variables we want to regress on
#' @param pval_cutoff is the cutoff on p-values
#' @return returns the estimated support
#' @export


pval_MSE <- function(X,y,variables,pval_cutoff){
  #need to check if some columns are always 0, especially the FE
  
  current_variables = variables
  deselect_list <- c()
  deselect_pval <- c()
  current_sp_formula <- as.formula(paste0(y,"~",paste0(c(current_variables,0),collapse = "+")))
  current_model_ols <- estimatr::lm_robust(formula = current_sp_formula, data = X, se_type = "classical")
  current_pvals = current_model_ols$p.value
  current_pvals[is.na(current_pvals)]=1
  current_max_pval <- max(current_pvals)
  n=length(variables)
  i=0
  print("Starting the multiple step elimination procedure")
  while (current_max_pval > pval_cutoff) {
    i=i+1
    cat("\rProgress: ",i," variables eliminated on ", n)
    deselect_name <- names(current_pvals[which.max(current_pvals)])
    deselect_list <- c(deselect_list, deselect_name)
    deselect_pval <- c(deselect_pval, current_pvals[which.max(current_pvals)])
    current_variables <- current_variables[current_variables != deselect_name]
    if (length(current_variables)==0){
      print("pval_cutoff is too strict, no variable survived")
      return(c())
    }
    current_sp_formula <- as.formula(paste0(y,"~",paste0(c(current_variables,0),collapse = "+")))
    current_model_ols <- estimatr::lm_robust(formula = current_sp_formula, data = X,se_type = "classical")
    current_pvals = current_model_ols$p.value
    current_pvals[is.na(current_pvals)]=1
    current_max_pval <- max(current_pvals)
  }
  cat("\n",current_max_pval,"\n")
  support = variables[!(variables %in% deselect_list)]
  
  result = list(support=support, pvals = deselect_pval, eliminated_variables = deselect_list)
  return(result)
  
}

#' From lambda to pval cutoff
#'
#' Compute the equivalent p-value between Puffer_N LASSO lambda and pval one-step elimination
#' @param X is the regression matrix with the marginals, the fixed effects, the outcome of interest
#' @param y is the outcome of interest
#' @return returns the pval equivalent to lambda 
#' @export



lambda_to_pval <- function(X,y,variables,lambda){
  formula <- as.formula(paste0(y,"~",paste0(c(variables,"0"),collapse = "+")))
  model_ols <- estimatr::lm_robust(formula = formula, data = X,  se_type = "classical")
  pval = 2*(1-pnorm(lambda/sqrt(model_ols$res_var)))
  return(pval)
}

#' From pval to lambda cutoff
#'
#' Compute the equivalent lambda between Puffer_N LASSO lambda and pval one-step elimination
#' @param X is the regression matrix with the marginals, the fixed effects, the outcome of interest
#' @param y is the outcome of interest
#' @return returns the lambda equivalent to pval 
#' @export

pval_to_lambda <- function(X,y,variables,pval){
  formula <- as.formula(paste0(y,"~",paste0(c(variables,"0"),collapse = "+")))
  model_ols <- estimatr::lm_robust(formula = formula, data = X,  se_type = "classical")
  lambda = model_ols$res_var %>% sqrt() * qnorm(1-pval/2)
  return(lambda)
}

#' pval one-step elimination
#'
#' Estimate support with p-value one-step elimination
#' @param X is the regression matrix with the marginals, the fixed effects, the outcome of interest
#' @param y is the outcome of interest
#' @param variables is a vector containing all the variables we want to regress on
#' @param pval_cutoff is the cutoff on p-values
#' @return returns the estimated support
#' @export

pval_OSE<- function(X,y,variables,pval_cutoff){
  #need to check if some columns are always 0, especially the FE
  
  lambda = pval_to_lambda(X,y,variables,pval_cutoff)
  
  formula <- as.formula(paste0(y,"~",paste0(c(variables,"0"),collapse = "+")))
  model_ols <- estimatr::lm_robust(formula = formula, data = X,  se_type = "classical")
  support = names(model_ols$p.value[which(model_ols$p.value<=pval_cutoff)])
  
  result = list(support=support, pvals = model_ols$p.value, equiv_lambda = lambda)
  
  return(result)
}

#' Puffer transform
#'
#' Perform the puffer transformation on X and Y
#' @param X a matrix with the observations
#' @param Y a vector with the outcome of interest
#' @return returns FX and FY, where puffer transformation has been applied to X and Y
#' @export


puffer_transform <- function(X,Y) {
  X_svd = svd(X)
  puffer_F = X_svd$u %*% solve(diag(X_svd$d)) %*% t(X_svd$u)
  FX = puffer_F %*% X
  FY = puffer_F %*% Y
  return(list(FX,FY))
}

#' N transform
#'
#' Perform the N transformation on X and Y
#' @param X a matrix with the observations
#' @return returns XN, where N transformation has been applied to X
#' @export

N_transform <- function(X){
  full_N_2 = (solve(t(X) %*% X))
  N = sqrt(full_N_2[row(full_N_2)==col(full_N_2)])
  X_N = t(t(X) * N) 
  return(X_N)
}

#' Puffer_N LASSO
#'
#' Estimate support with Puffer_N LASSO
#' @param X is the regression matrix with the marginals, the fixed effects, the outcome of interest
#' @param y is the outcome of interest
#' @param variables is a vector containing all the variables we want to regress on
#' @param lambda_cutoff is the cutoff on beta values
#' @return returns the estimated support
#' @export


puffer_N_LASSO <- function(X,y,variables,lambda_cutoff){
  #Estimate residual deviation to compute equivalent pval_cutoff
  pval_cutoff = lambda_to_pval(X,y,variables,lambda_cutoff)
  print("Equivalent p-value cutoff to current lambda cutoff is ")
  print(pval_cutoff)
  
  X_matrix = X[,c(variables)] %>% as.matrix()
  Y = X[,y] %>% as.matrix()
  
  X_N = N_transform(X_matrix)
  
  puffer_result = puffer_transform(X_N,Y)
  X_PT = puffer_result[[1]]
  Y_PT = puffer_result[[2]]
  
  lambda_cutoff_glm = lambda_cutoff/nrow(X_PT) #this normalization is needed as glmnet considers a different optimization function
  lasso_model=glmnet::glmnet(X_PT,Y_PT,alpha=1,lambda = lambda_cutoff_glm, intercept=FALSE,standardize=FALSE) #intercept is already in X dataframe and scaling is treated in the main function
  support=names(coef(lasso_model)[,1][which(abs(coef(lasso_model)[,1]) > 0)])
  
  result = list(support=support)
  return (result)
}

#' beta one-step elimination
#'
#' Estimate support with beta one-step elimination
#' @param X is the regression matrix with the marginals, the fixed effects, the outcome of interest
#' @param y is the outcome of interest
#' @param variables is a vector containing all the variables we want to regress on
#' @param lambda_cutoff is the cutoff on beta values
#' @return returns the estimated support
#' @export

beta_OSE <- function(X,y,variables,lambda_cutoff){
  formula <- as.formula(paste0(y,"~",paste0(c(variables,"0"),collapse = "+")))
  model_ols <- estimatr::lm_robust(formula = formula, data = X,  se_type = "classical")
  support = names(model_ols$p.value[which(abs(model_ols$coefficients)>=lambda_cutoff)])
  
  result = list(support=support, beta = model_ols$coefficients)
  return(result)
}

#' Puffer LASSO
#'
#' Estimate support with Puffer LASSO
#' @param X is the regression matrix with the marginals, the fixed effects, the outcome of interest
#' @param y is the outcome of interest
#' @param variables is a vector containing all the variables we want to regress on
#' @param lambda_cutoff is the cutoff on beta values
#' @return returns the estimated support
#' @export



puffer_LASSO <- function(X,y,variables,lambda_cutoff){
  #do we add the intercept column ?
  X_matrix = X[,variables] %>% as.matrix()
  Y = X[,y] %>% as.matrix()
  
  puffer_result = puffer_transform(X_matrix,Y)
  X_PT = puffer_result[[1]]
  Y_PT = puffer_result[[2]]
  
  lambda_cutoff_glm = lambda_cutoff/nrow(X_PT) #this normalization is needed as glmnet considers a different optimization function
  
  lasso_model=glmnet::glmnet(X_PT,Y_PT,alpha=1,lambda=lambda_cutoff_glm,intercept=FALSE,standardize=FALSE)
  support=names(coef(lasso_model)[,1][which(abs(coef(lasso_model)[,1]) > 0)])
  
  result = list(support=support)
  
  return (result)
}








