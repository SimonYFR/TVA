#################################################################
#
#------ This file contains all functions related to pooling and
#------ the main do_TVA wrapper function
#
#################################################################


# Import the pipe function
#' @importFrom magrittr %>%
#' @export

`%>%` = magrittr::`%>%`


#' Check user's inputs integrity
#'
#' Do a bunch of checks on the inputs
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w is the column name of the weights
#' @param cutoff is the cutoff used in the support estimation. This can be empty.
#' @param estim_func is the estimation function we should used. Possible arguments are :\cr
#' 1. 'pval_MSE': a multiple step elimination on p-values\cr
#' 2. 'pval_OSE': a one step elimination on p-values\cr
#' 3. 'puffer_N_LASSO': a LASSO OLS with a Puffer_N transformation\cr
#' 4. 'beta_OSE': a one step elimination on beta values\cr
#' 5. 'puffer_LASSO': a LASSO OLS with a Puffer transformation\cr
#' (2) and (3) should be equivalent, as well as (4) and (5)
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns a message telling if the data's integrity is validated or not, and if not, for what reason
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
#' check_inputs_integrity(data=data, arms=arms, y=y, fes=fes, w=w, cutoff = 0.3, estim_func='pval_OSE', FALSE)

check_inputs_integrity <- function(data, arms, y, fes=c(), cutoff=NULL, w=NULL, estim_func='pval_OSE', compare_to_zero=FALSE, clusters=NULL){
  
  if (!(class(data) == "data.frame")){
    return(list(integrity=FALSE, message="data should be a dataframe"))
  }
  
  if (!is.null(w)){
    if (!((class(w) == "character") & (w %in% names(data)))){
      return(list(integrity=FALSE, message="w should either be NULL or a column name of data"))
    }
  }
  
  if (!is.null(clusters)){
    if (!((class(clusters) == "character") & (clusters %in% names(data)))){
      return(list(integrity=FALSE, message="clusters should either be NULL or a column name of data"))
    }
  }
  
  if (!((class(y) == "character") & (y %in% names(data)))){
    return(list(integrity=FALSE, message="y should be a column name of data"))
  }
  
  if (!(estim_func %in% c('pval_OSE','pval_MSE','beta_OSE','puffer_N_LASSO','puffer_LASSO'))){
    return(list(integrity=FALSE, message="estim_func can only equal 'pval_OSE','pval_MSE','beta_OSE','puffer_N_LASSO', or 'puffer_LASSO' "))
  }
  
  if (!is.null(cutoff)) {
    if (!((class(cutoff) == "numeric") & (cutoff>=0))){
      return(list(integrity=FALSE,message="cutoff should be a positive real number"))
    }
  }
  
  if (!(class(compare_to_zero) == "logical")){
    return(list(integrity=FALSE,message="compare_to_zero should be TRUE or FALSE"))
  }
  
  if (!length(arms)>0){
    return(list(integrity=FALSE,message="arms can't be empty"))
  }
  
  test = TRUE
  for (e in unlist(arms)){
    test = test & (class(e) == "character")
  }
  if (!test){
    return(list(integrity=FALSE,message="arms should contain only strings or list of strings"))
  }
  
  if (!all(intersect(names(data),arms %>% unlist()) == arms %>% unlist() )){
    return(list(integrity=FALSE,message="strings inside arms should be column names of data"))
  }
  
  if (any(duplicated(unlist(arms)))){
    return(list(integrity=FALSE,message="arms can not contained columns names duplicates"))
  }
  
  
  if (unlist(arms) %>% length() == length(arms) ){
    test = TRUE
    for (arm in arms){
      test = test & (all(data[,arm]>=0))
    }
    if (!test){
      return(list(integrity=FALSE,message="arms columns should contain positive (>=0) values"))
    }
  }else{
    if (!all( (data[,unlist(arms)]==0) | (data[,unlist(arms)]==1) )){
      return(list(integrity=FALSE,message="if arms is a list of lists, all the column names should refer to dummy columns with values that are 0 or 1"))
    }
    for (e in arms){
      if (length(e)==0){
        return(list(integrity=FALSE,message="arms can't contain any empty list"))
      }
      else if (length(e)>1){
        if (!all( data[,e] %>% rowSums() < 2 )){
          return(list(integrity=FALSE,message="if arms is a list of lists, all the column names inside a common list must be mutually exclusive indicators. One observation can't have two different dosages on the same arm."))
        }
      }
    }
    if (!names(arms) %>% length() == arms %>% length()){
      return(list(integrity=FALSE,message="if arms is a list of lists, all lists inside arms should be named"))
    }
    if (any(duplicated(names(arms)))){
      return(list(integrity=FALSE,message="if arms is a list of lists, two lists can't have the same name"))
    }
  }
  
  
  test = TRUE
  for (fe in fes){
    test = test & (class(fe) == "character")
  }
  if (!test){
    return(list(integrity=FALSE,message="fes should be a vector of strings"))
  }
  
  if (!all(intersect(names(data),fes) == fes )){
    return(list(integrity=FALSE,message="strings inside fes should be column names of data"))
  }
  
  if (!all(data[,w]>0)){
    return(list(integrity=FALSE,message="column w should have strictly positive values"))
  }
  

  test = TRUE
  test = test & (length(intersect(arms %>% unlist(),fes))==0)
  test = test & (length(intersect(arms %>% unlist(),y))==0)
  test = test & (length(intersect(arms %>% unlist(),w))==0)
  test = test & (length(intersect(y,fes))==0)
  test = test & (length(intersect(w,fes))==0)
  test = test & (length(intersect(y,w))==0)
  
  if (!test){
    return(list(integrity=FALSE,message="arms, y, fes and w should be different column names"))
  }
  
  return(list(integrity=TRUE,message=""))
}


#' From dummy columns to dosages columns
#'
#' Going from the dummy columns to a unique column per arm with the dosages
#' @param data is a dataframe with the data 
#' @param arms is a named list, and each element of the list is a vector containing the dummy columns of one arm dosages, in the order of increasing dosages and without the dosage 0 dummy. The name of this element must be the arm name. 
#' @return returns the dataframe "data" with one column per arm showing the dosage.
#' @export
#' @examples
#' sms_lvl_1 = c(0,0,0,1,0,1)
#' sms_lvl_2 = c(0,1,1,0,0,0)
#' incentive_lvl_1 = c(1,0,1,0,1,0)
#' incentive_lvl_2 = c(0,0,0,1,0,1)
#' information_lvl_1 = c(1,0,0,1,0,1)
#' information_lvl_2 = c(0,0,0,0,1,0)
#' information_lvl_3 = c(0,1,1,0,0,0)
#' arms = list(sms = c('sms_lvl_1','sms_lvl_2'),
#'             incentive = c('incentive_lvl_1','incentive_lvl_2'),
#'             info = c('information_lvl_1','information_lvl_2','information_lvl_3'))
#' df = data.frame(sms_lvl_1, sms_lvl_2, incentive_lvl_1, incentive_lvl_2, information_lvl_1, information_lvl_2, information_lvl_3)
#' create_dosages_from_dummies(data,arms)


create_dosages_from_dummies <- function(data,arms){
  for (k in seq_along(arms)){
    data[,names(arms)[k]] = rowSums(t(t(data[,arms[[k]] ])*seq_along(arms[[k]])))
  }
  return(list(data=data,arms=names(arms)))
}


#' Going from vector to string
#'
#' Create a string representing the vector
#' @param vector is a vector of any size and containing integers
#' @param first_char is a string that appears at the beginning of the vector, equals "c" by default
#' @return returns a string representing the vector in format 'c_A_B_C_D_...' where A,B,C,D,.. are the integers in the vector
#' @export
#' @examples
#' vector_to_string(c(1,2,3,4,5))
#' vector_to_string(c(10))

vector_to_string <- function(vector,first_char="c"){
  return(paste0(c(paste0(c(paste0(first_char,"_"),paste0(vector,collapse='_')),collapse=''),''),collapse=''))
}

#' Going from string to vector
#'
#' Create a vector from its string representation 
#' @param string is a string in the format 'c_A_B_C_D_..' where A,B,C,D,.. are integers
#' @return returns a vector containing the values represented by the string, c(A,B,C,D,...)
#' @export
#' @examples
#' string_to_vector('c_1_2_3_4_5')

string_to_vector <- function(string){
  return(stringr::str_split(substring(string,3,nchar(string)),'_')[[1]] %>% as.integer())
}

#' Going from policy concise string to a full name
#'
#' Create the policy full name from its string representation
#' @param policy is a policy in the format 'p_A_B_C_D_..' where A,B,C,D,... are integers representing the dosage for arm 1,2,3,4,...\cr
#' For example, 'p_1_0_2' means arm 1 has dosage 1, arm 2 is not activated and arm 3 has dosage 2
#' @param arms is a vector containing the arms' names and of the same length as the vector representation of the policy
#' @return returns the full name of the policy with the dosage of each arm
#' @export
#' @examples
#' get_policy_fullname('p_1_0_2', c('financial_incentive','reminder','information'))

get_policy_fullname <- function(policy,arms){
  policy_vector = string_to_vector(policy)
  full_name = paste(paste(arms,policy_vector,sep="="),collapse=" and ")
  return(full_name)
}



#' Create marginals matrix
#'
#' Create the marginals matrix with zeros everywhere
#' @param max_dosage_per_arm is a vector containing the maximum dosage possible for each arm
#' @param n_obs is the number of observations
#' @return returns the empty marginals matrix, in the dataframe format, with all indicators equal to zero
#' @export
#' @examples
#' create_empty_marginals_matrix(c(1,1,3),15)

create_empty_marginals_matrix <- function(max_dosage_per_arm,n_obs){
  cat("Creating an empty marginals matrix","\n")
  grid_list=list()
  for (max_dosage in max_dosage_per_arm){
    grid_list=append(grid_list,list(0:max_dosage))
  }
  grid=expand.grid(grid_list)
  
  marginals_colnames = apply(grid, 1, function(x) vector_to_string(x,first_char = "m"))
  marginals_matrix = setNames(data.frame(matrix(0,ncol = (nrow(grid)) , nrow = n_obs) ),  marginals_colnames)
  
  return(marginals_matrix)
}

#' Fill marginals matrix
#'
#' Fill the empty marginals matrix already created
#' @param marginals_matrix is the marginals matrix in the dataframe format created previously
#' @param data is the data with the dosage on each arm
#' @param arms is a vector containing the column names that represent the dosages of each arm in data
#' @param n_obs is the number of observations, equal to the number of rows in "data"
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns the marginals matrix with the right indicators inside
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3)
#' n_obs = nrow(data)
#' empty_marginals_matrix = create_empty_marginals_matrix(c(1,1,3),n_obs)
#' fill_marginals_matrix(empty_marginals_matrix, data, arms, n_obs,FALSE)

fill_marginals_matrix <- function(marginals_matrix,data,arms,n_obs,compare_to_zero){
  cat("Filling the marginals matrix","\n")
  M=length(arms)
  marginals_colnames = colnames(marginals_matrix)
  
  for (marginal in marginals_colnames){
    marginal_vector = string_to_vector(marginal)
    policy_dominates_marginal = (rowSums((t(t(data[,arms]) - marginal_vector)< 0))==0)
    policy_resembles_marginal = compare_to_zero | (rowSums((t(t(data[,arms]==0) - marginal_vector==0)!=0))==0)
    marginals_matrix[,marginal] = 1*policy_dominates_marginal*policy_resembles_marginal
  }
  
  marginals_matrix = marginals_matrix[, colSums(marginals_matrix != 0) > 0]
  control_marginal_index = which(names(marginals_matrix) %in% c(vector_to_string(rep(0,M), first_char = "m") ))
  if (length(control_marginal_index) > 0){
    marginals_matrix = marginals_matrix[ , -control_marginal_index]
  }
  return(marginals_matrix)
}

#' Weight observations
#'
#' Weight the observations directly in the matrix. \cr
#' Performing OLS on the output matrix is equivalent to performing an OLS on the original matrix and specifying a weight column.
#' @param X is contains all the observations (without the weights)
#' @param W is a vector containing the weights
#' @return returns the weighted observations, where observation i is weighted by sqrt(w_i)
#' @export
#' @examples
#' X = data.frame(outcome = c(2,2,3,2), financial_incentive = c(0,1,0,1), reminder = c(1,1,0,0), information = c(0,1,2,3))
#' W = c(3,4,5,2)
#' weight_observations(X,W)

weight_observations <- function(X,W){
  cat("Weighting the observations","\n")
  rW = sqrt(W)
  X = X * rW
  return (X)
}

#' Prepare the data 
#'
#' Prepare the data for the support estimation
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w is the column name of the weights
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns a list containing :\cr
#' * X: a dataframe containing the marginals matrix, the fixed effects, the outcome and the intercept. X is ready for support estimation\cr
#' * data: the original dataframe with new columns for arm dosages (if dummy arm columns were provided by the user)
#' * arms: a vector containing the columns with the dosages on each arms (no more dummy columns if it was provided by the user). \cr
#' * marginals_colnames: a vector containing all the marginals names (also called the alphas)\cr
#' * variables: the list of variables that should be used in the regression (fixed effects, marginals and the intercept)
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
#' prepare_data(data,arms,y, fes,w,FALSE,FALSE)

prepare_data <- function(data, arms, y, fes, w, compare_to_zero){
  cat("Preparing the data","\n")
  n_obs = nrow(data)
  
  if (is.null(w))  {   W=rep(1,n_obs)   }  else   {   W=data[,w]   }
  
  #check if arms is a list of string or a list of lists of strings
  if (unlist(arms) %>% length() != length(arms) ){ 
    preprocess = create_dosages_from_dummies(data,arms)
    data = preprocess$data
    arms = preprocess$arms
  }
  
  #creating the vector of maximum dosage per arm
  max_dosage_per_arm = sapply(data[,arms], max, na.rm = TRUE)
  
  #Initializing the marginal space matrix
  marginals_matrix = create_empty_marginals_matrix(max_dosage_per_arm,n_obs) %>% fill_marginals_matrix(.,data,arms,n_obs,compare_to_zero)
  marginals_colnames = names(marginals_matrix)
  
  #Creating the X matrix on which we will estimate the alphas
  X = cbind(marginals_matrix, data[,c(fes,y)])
  
  #que faire si intercept existe déjà dans X ?
  X['intercept'] = 1
  
  X = weight_observations(X,W)
  variables = c(fes,marginals_colnames,'intercept')
  
  return(list(X=X,data=data,arms=arms,variables=variables,marginals_colnames=marginals_colnames))
}



#' Pool the data
#'
#' Take the estimated support and pool each observation in its right pool
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param marginal_support_strings is a vector containing strings that represent all the marginals in the support
#' @param compare_to_zero is a boolean.\cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @return returns the dataframe "data" with new columns: \cr
#' * a pool_id column that gives the pool id of the observation's pool 
#' * one dummy column per pool_id, equal to 1 if the observation belongs to this pool id, 0 otherwise 
#' * one column per marginal_j, equal to 1 if the observation i is influenced by the marginal n°j, 0 otherwise. There are as many columns as marginals in the support. 
#' * a "pool_influences" column, with a string of the format "i_x1_x2_x3_.." where xj is equal to 1 if the marginals n°j influences the observation, 0 otherwise. This is basically the definition of the pool the observation belongs to.
#' * a "pool_influences_list" column, with a string that gives all the marginals that influence the observation 
#' @export
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3)
#' marginal_support_strings = c('m_0_1_1', 'm_1_1_2', 'm_1_0_1', 'm_1_0_2')
#' pool_data(data,arms,marginal_support_strings,FALSE)

pool_data <- function(data,arms,marginal_support_strings,compare_to_zero){
  n_obs = nrow(data)
  S = length(marginal_support_strings)
  
  cat("Pooling the raw data with estimated marginal support of size ",S,"\n")

  data$pool_influences = "i"
  data$pool_influences_list = ""
  data$pool_id = 0
  if (S>=1){
    for (i in 1:S){
      a_i_string = marginal_support_strings[i] #get the current marginal in string format
      a_i = string_to_vector(a_i_string) #get the current marginal in vector format
      policy_dominates_a_i = (rowSums((t(t(data[,arms]) - a_i)< 0))==0) 
      policy_resembles_a_i = compare_to_zero | (rowSums((t(t(data[,arms]==0) - a_i==0)!=0))==0) 
      indicators = 1*policy_dominates_a_i*policy_resembles_a_i
      data[,paste("marginal",as.character(i),sep="_")] = indicators
      data$pool_influences = paste(data$pool_influences,indicators,sep='_') 
      
      string_replace = c('0'='','1'=paste(', ',a_i_string))
      data$pool_influences_list = paste(data$pool_influences_list,string_replace[as.character(indicators)] %>% unname(),sep='')
    }
    data$pool_influences_list = gsub("^.{0,3}", "", data$pool_influences_list)
    data$pool_id = as.numeric(as.factor(data$pool_influences))-1 #this gives an id to each pool_influences, -1 ensures that i_0_0_.._0 has id = 0 
    #create dummy columns
    pool_dummy = data.frame(lme4::dummy(data$pool_id))
    pool_ids = paste("pool_id",stringr::str_sub(names(pool_dummy),2),sep="_")
    names(pool_dummy) = pool_ids
    data = cbind(data,pool_dummy)
  }
  return(data)
}

#' Give information for each pool
#'
#' Compute some useful information to understand what's inside each pool
#' @param data is the dataframe containing all our data and a "pool_id" column giving the pool id of each observation
#' @param arms is a vector containing the column names of all the arms
#' @return returns a list containing :\cr
#' * pools_summary: a dataframe containing information on each pool, with columns:
#'     * "pool_id": the pool id
#'     * "n_unique_policies": the number of unique policies inside this pool
#'     * "n_obs": the number of observation in this pool
#'     * "pool_influences": a string of format "i_x1_x2_x3..", where x_i equals 1 if the marginal number i influences this pool, 0 otherwise
#'     * "pool_influences_list": a string that gives all the marginals that influence this pool
#'     * "pool_minimum": the smallest unique policy inside this pool
#'     * "min_arm1", "min_arm2" etc... : one column per arm, giving the minimum dosage for each arm inside this pool
#'     * "pool_minimum_fullname": the full name of the smallest unique policy
#'     * "policy_examples": 5 or less examples of unique policies that are inside this pool
#' 
#' * unique_policy: a dataframe containing information on each unique policy, with columns:
#'     * "policy": a string in format "p_A_B_C_..." where A is the dosage on arm 1, B on arm 2 etc... This is the definition of the unique policy
#'     * "policy_fullname": the full name of the unique policy
#'     * "arm1", "arm2" etc...: the value on each arm of the unique policy
#'     * "pool_id": the id of the pool the unique policy belongs to
#'     * "pool_influences": a string of format "i_x1_x2_x3..", where x_i equals 1 if the marginal number i influences this pool, 0 otherwise
#'     * "pool_influences_list": a string that gives all the marginals that influence this pool
#'     
#' 
#' @export
#' @examples
#' arms = c('financial_incentive','reminder','information')
#' A1 = c(0,0,0,0,0,1,1,1,1,1)
#' A2 = c(1,1,0,0,1,1,0,0,1,1)
#' A3 = c(0,1,2,3,0,3,2,1,0,1)
#' pool_id = c(0,1,0,0,0,2,0,0,0,2)
#' data = data.frame(financial_incentive = A1, reminder = A2, information = A3, pool_id = pool_id)
#' pools_info(data,arms)

pools_info <- function(data,arms){
  cat("Gathering informations on pools","\n")
  unique_policy = data[!duplicated(data[,c(arms,'pool_influences','pool_influences_list')]), ][,c(arms,'pool_influences','pool_id','pool_influences_list')] %>% as.data.frame(row.names = 1:nrow(.)) #taking all the unique policies by pool_id
  unique_policy = unique_policy %>% dplyr::mutate(., policy = apply(unique_policy[,arms], 1, vector_to_string(., first_char = "p")))  %>% dplyr::arrange(., pool_id,policy) #creating policy string column
  unique_policy$policy_fullname = sapply(unique_policy[,'policy'], get_policy_fullname, arms=arms) #create policy full name column
  
  a0 = stats::aggregate(unique_policy$pool_influences, by=list(pool_id=unique_policy$pool_id), FUN=length) %>% setNames(.,c('pool_id','n_unique_policies')) #counting number of unique policies by pool_id
  
  a1 = stats::aggregate(data$pool_influences, by=list(pool_id=data$pool_id), FUN=length) %>% setNames(.,c('pool_id','n_obs')) #counting the number of observations by pool_id
  
  a2 = stats::aggregate(data$pool_influences, by=list(pool_id=data$pool_id), FUN=dplyr::first) %>% setNames(.,c('pool_id','pool_influences')) #taking the pool definition (in terms of influence) by pool_id
  
  a3 = stats::aggregate(data$pool_influences_list, by=list(pool_id=data$pool_id), FUN=dplyr::first) %>% setNames(.,c('pool_id','pool_influences_list')) #taking the pool definition (in terms of influence) by pool_id
  
  a4 = stats::aggregate(data[,arms], by=list(pool_id=data$pool_id), FUN=min) #taking the minimum value by arm by pool_id
  colnames(a4)[-1] <- paste("min", colnames(a4)[-1], sep = "_")
  
  first_5_examples = unique_policy[,c('pool_id','policy')] %>% dplyr::group_by(pool_id) %>% dplyr::slice(1:5) #taking 5 policies examples by pool_id
  last_example = unique_policy[,c('pool_id','policy')] %>% dplyr::group_by(pool_id) %>% dplyr::slice(6:6) %>% dplyr::mutate(.,policy='...') #taking the potential 6th one and overwrite it as "..." to show there are more than 5
  a5 = rbind(first_5_examples, last_example) %>% dplyr::group_by(pool_id) %>%  dplyr::summarize(policy_examples = paste((policy),collapse=", ")) %>% as.data.frame()
  
  pools_summary = merge(a0,a1,by='pool_id') %>% merge(.,a2,by='pool_id') %>% merge(.,a3,by='pool_id') %>% merge(.,a4,by='pool_id') %>% merge(.,a5,by='pool_id')
  pools_summary$pool_minimum = apply(pools_summary[,colnames(a4)[-1]], 1, function(x) vector_to_string(x, first_char = "p"))
  pools_summary$pool_minimum_fullname = sapply(pools_summary[,'pool_minimum'], get_policy_fullname, arms=arms)
  
  pools_summary = pools_summary[c(setdiff(names(pools_summary), 'policy_examples'), 'policy_examples')] #put policy example column at the end
  return(list(pools_summary = pools_summary, unique_policy = unique_policy))
}

#' Perform the final OLS
#'
#' Perform the final OLS on pooled data, regressing on fes and pool_id dummy columns
#' @param data is the dataframe containing all our data and pool_id dummy columns
#' @param fes is a vector containing the column names of all the fixed effects
#' @param y is the column name of the outcome of interest
#' @param w is the column name of the weights
#' @param pool_ids is a vector containing the column names of each pool_id dummy column
#' @param clusters (optional) is the column name that corresponds to the clusters in the data, that should be used in the final pooled OLS. Please refer to estimatr::lm_robust documentation for more information on this parameter.
#' @return returns the ols result
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
#' p1 = c(0,1,0,0,0,0,0,0,0,0)
#' p2 = c(0,0,0,0,0,1,0,0,0,1)
#' pool_ids = c('pool_id_1','pool_id_2')
#' data = data.frame(pool_id_1 = p1, pool_id_2 = p2, financial_incentive = A1, reminder = A2, information = A3, fes_1 = F1, outcome = Y, weights=W)
#' get_pooled_ols(data,y, fes, w, pool_ids)

get_pooled_ols <- function(data,y, fes=c(),w=NULL,pool_ids, clusters){
  cat("Performing the final OLS on pooled data","\n")
  
  pooled_ols_variables = c(pool_ids, fes)
  formula = as.formula(paste0(y,"~",paste0(pooled_ols_variables ,collapse = "+")))
  
  if (is.null(w))         {   W=NULL         }  else   {   W=data[,w]      }
  if (is.null(clusters))  {   Clusters=NULL  }  else   {   Clusters=data[,clusters]    }
  if (is.null(clusters))  {   se_type='HC2'  }  else   {   se_type='CR0'   }
  pooled_ols = estimatr::lm_robust(formula = formula, data = data, weights = W, clusters=Clusters, se_type=se_type)

  return(pooled_ols)
}




############# ------- TVA Main Function ------- #############


#' TVA Main function
#'
#' Perform the whole TVA algorithm
#' @param data is the dataframe containing all our data
#' @param arms is a vector containing the column names of all the arms
#' @param fes (optional) is a vector containing the column names of all the fixed effects, empty by default
#' @param y is the column name of the outcome of interest
#' @param w (optional) is the column name of the weights
#' @param cutoff (optional) is the cutoff used in the support estimation. If empty, TVA will suggest its own cutoff based on the penalty / support size tradeoff.
#' @param estim_func (optional) is the estimation function we should used. Possible arguments are :\cr
#' 1. 'pval_MSE': a multiple step elimination on p-values\cr
#' 2. 'pval_OSE': a one step elimination on p-values\cr
#' 3. 'puffer_N_LASSO': a LASSO OLS with a Puffer_N transformation\cr
#' 4. 'beta_OSE': a one step elimination on beta values\cr
#' 5. 'puffer_LASSO': a LASSO OLS with a Puffer transformation\cr
#' (2) and (3) should be equivalent, as well as (4) and (5)
#' @param compare_to_zero (optional) is a boolean, FALSE by default. \cr
#' If TRUE, the code considers that a policy dominates a marginal if all dosages are greater\cr
#' If FALSE, then they must also have the exact same activated arms (the zeros of the policy vectors are at identical indexes)
#' @param clusters (optional) is the column name that corresponds to the clusters in the data, that should be used in the final pooled OLS. Please refer to estimatr::lm_robust documentation for more information on this parameter.
#' @return returns a list containing:\cr
#' * data: the data with new columns giving pooling information\cr
#' * marginal_support: a dataframe with all the marginals in the support and their according id\cr
#' * pools_summary: a dataframe with information on each pool\cr
#' * unique_policy: a dataframe with all the possible unique policies and their according pool id\cr
#' * fes_support: the intersection between the estimated support and the fixed effects\cr
#' * pooled_ols: the result of the final OLS on the pooled data\cr
#' * winners_effect: the result of the best pooled policy effect, downsized by the winners curse algorithm
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
#' TVA(data,arms,y, fes,w,0.3,'pval_OSE',FALSE,FALSE)

do_TVA <- function(data, arms, y, fes=c(), w = NULL, cutoff = NULL, estim_func = 'pval_OSE', compare_to_zero = FALSE, clusters = NULL){
  
  #check if fake_weights column already exists
  check = check_inputs_integrity(data, arms, y, fes, cutoff, w, estim_func, compare_to_zero, clusters)
  
  if (!check$integrity){
    stop(check$message)
  }
  
  if (is.null(cutoff)) {      #If no cutoff was provided, suggest cutoff
    cat("No cutoff provided. Suggesting cutoff... \n")
    quiet_suggest_pval = purrr::quietly(suggest_pval_cutoff)  #This is to avoid all the additional printing
    cutoff = quiet_suggest_pval(data=data, arms=arms, y=y, fes=fes, w=w, estim_func = estim_func, compare_to_zero=compare_to_zero, clusters=clusters)$result$threshold
    cat("Suggested cutoff is ", cutoff, " (run \"suggest_pval_cutoff\" for details) \n")
  }
  
  
  #Data processing
  prepared_data = prepare_data(data,arms,y, fes,w,compare_to_zero)
  X = prepared_data$X
  data = prepared_data$data
  arms = prepared_data$arms
  variables = prepared_data$variables
  marginals_colnames = prepared_data$marginals_colnames

  
  #Get support
  f = get(estim_func)
  total_support = f(X,y,variables,cutoff)$support
  marginal_support_strings = sort(intersect(total_support,marginals_colnames))
  
  if (marginal_support_strings %>% length() == 0){
    stop("Estimated support is empty, current cutoff does not differentiate any policy with the control")
  }
  cat("Estimated support is:","\n")
  cat(marginal_support_strings,"\n")
  fes_support = sort(intersect(total_support,fes))

  #Pool policies
  data = pool_data(data,arms,marginal_support_strings,compare_to_zero)
  
  #Create marginals ids
  marginal_support = data.frame(marginal = marginal_support_strings)
  marginal_support$marginal_id = as.numeric(as.factor(marginal_support$marginal))
  
  #Give info about pools
  pools_info = pools_info(data,arms)
  pools_summary = pools_info$pools_summary
  unique_policy = pools_info$unique_policy
  
  #Final pooled ols
  pool_ids = paste("pool_id",c(1:max(data$pool_id)),sep="_")
  pooled_ols = get_pooled_ols(data,y,fes_support,w,pool_ids,clusters)

  #Apply winners curse
  winners_effect = winners_curse(pooled_ols,pool_ids,alpha=0.05,beta=0.005)
  
  #Return result
  result =  list(data=data
                 ,marginal_support=marginal_support
                 ,pools_summary=pools_summary
                 ,unique_policy = unique_policy 
                 ,fes_support=fes_support 
                 ,pooled_ols=pooled_ols 
                 ,winners_effect=winners_effect
  )
  
  cat("Returning result","\n")
  return(result)
}



#################################################################
# END
#################################################################





