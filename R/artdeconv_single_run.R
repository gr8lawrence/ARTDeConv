#' ARTdeConv function on a single run
#' 
#' This is an artdeconv function on a single run. 
#' 
#' @param Y the bulk matrix
#' @param Theta_0 \eqn{m_0 \times K_0} reference matrix for the fixed part of the signature matrix
#' @param Theta_it the initial value of the signature matrix
#' @param s_it the initial value of the "cell size" vector
#' @param P_it the initial value of the proportions matrix
#' @param m0 the number of signature genes for well characterized cell types in a tissue
#' @param k0 the number of well characterized cell types in a tissue
#' @param meds a vector of length \eqn{K} of pre-specified medians of cell proportions 
#' @param ranges a vector of length \eqn{K} of pre-specified ranges of cell-type proportions (1/ranges as weights for regularization parameters)
#' @param alpha1 the tuning parameter for regularizing the known part of Theta
#' @param alpha2 the tuning parameter for regularizing the unknown part of Theta
#' @param beta the tuning parameter for regularizing P
#' @param max_iter the maximal number of iterations the algorithm shall perform; default value is \eqn{10^5}
#' @param tol tolerance for convergence indication; default value is \eqn{10^{-6}}
#' 
#' @return A list with following items
#' \itemize{
#' \item{Y_hat the estimated Y from the product of estimated factors}
#' \item{Theta_hat the estimated signature matrix from deconvolution}
#' \item{s_hat the estimated "cell size" vector from deconvolution}
#' \item{P_hat the estimated Theta from deconvolution}
#' }
#' 
#' @export

artdeconv_single_solve <- function(Y, Theta_0, Theta_it, s_it, P_it, m, n, k, m0, k0, meds, ranges, alpha1, alpha2, beta, max_iter = 1e5, tol = 1e-5) {
  
  ## pass the initiate values to "old" values
  Theta_old <- Theta_it
  s_old <- s_it
  P_old <- P_it
  Y_old <- Theta_old %*% diag(s_old) %*% P_old
  obj_old <- artdeconv_obj_fun(Y, Y_old, Theta_old, P_old, m, n, k, m0, k0, Theta_0, meds, ranges, alpha1, alpha2, beta) # calculate the initial objective function value
  
  ## initiate a residual and objective function vector
  res_v <- vector()
  obj_v <- vector()
  Delta <- get_Delta(k, k0, m, m0) # get the Delta matrix
  Delta_c <- 1 - Delta
  n_iter <- 0 # initiate the counter
  delt_obj <- Inf # set the initial delt_obj to inf 
  
  ## start the loop
  for (t in seq_len(max_iter)) {
    if (delt_obj <= tol) break # if delt_Y gets under the tolerance, break out of the loop
    Theta_new <- mu_Theta(Y, s_old, P_old, Theta_0, Delta, Delta_c, Theta_old, m, n, alpha1, alpha2)
    P_new <- mu_P(Y, Theta_new, s_old, meds, P_old, m, n, beta, ranges)
    s_new <- mu_s(Y, Theta_new, s_old, P_new)
    # if (!fixed_s) {
    #   s_new = mu_s(Y, Theta_new, s_old, P_new)
    # } else {
    #   s_new = s_old
    # }
    Y_new <- Theta_new %*% diag(s_new) %*% P_new # calculating the updated Y
    
    ## calculate the updated objective funtion and its change
    obj_v[t] <- obj_new <- artdeconv_obj_fun(Y, Y_new, Theta_new, P_new, m, n, k, m0, k0, Theta_0, meds, ranges, alpha1, alpha2, beta)
    res_v[t] <- norm(Y_new - Y, type = "F")/sqrt(m * n) ## storing all residuals of Y
    delt_obj <- abs(obj_new - obj_old)/obj_old # calculate the change in objective function values
    
    ## passing the updates to the next round
    Y_old <- Y_new
    Theta_old <- Theta_new
    s_old <- s_new
    P_old <- P_new
    obj_old <- obj_new
  }
  
  ## normalize the proportion matrix P so its row sums equal 1
  P_new <- apply(P_new, 2, function(x) x/sum(x))
  return(list(Y_hat = Y_new,
              Theta_hat = Theta_new,
              s_hat = s_new,
              P_hat = P_new,
              res_v = res_v,
              obj_v = obj_v,
              fixed_s = 'size-free',
              weights = ifelse(length(ranges) == 1, 'uniform', 'range-adaptive'),
              n_iter = t - 1))
}

#' ARTdeConv function for a bi-factor model with fixed S 

artdeconv_single_solve_s_fixed <- function(Y, Theta_0, Theta_it, P_it, m, n, k, m0, k0, meds, ranges, alpha1, alpha2, beta, max_iter = 1e5, tol = 1e-5) {
  
  ## pass the initiate values to "old" values
  Theta_old <- Theta_it
  s_old <- rep(1, k)
  P_old <- P_it
  Y_old <- Theta_old %*% P_old
  obj_old <- artdeconv_obj_fun(Y, Y_old, Theta_old, P_old, m, n, k, m0, k0, Theta_0, meds, ranges, alpha1, alpha2, beta) # calculate the initial objective function value
  
  ## initiate a residual and objective function vector
  res_v <- vector()
  obj_v <- vector()
  Delta <- get_Delta(k, k0, m, m0) # get the Delta matrix
  Delta_c <- 1 - Delta
  n_iter <- 0 # initiate the counter
  delt_obj <- Inf # set the initial delt_obj to inf 
  
  ## start the loop
  for (t in seq_len(max_iter)) {
    if (delt_obj <= tol) break # if delt_Y gets under the tolerance, break out of the loop
    Theta_new <- mu_Theta(Y, s_old, P_old, Theta_0, Delta, Delta_c, Theta_old, m, n, alpha1, alpha2)
    P_new <- mu_P(Y, Theta_new, s_old, meds, P_old, m, n, beta, ranges)
    Y_new <- Theta_new %*% P_new # calculating the updated Y
    
    ## calculate the updated objective funtion and its change
    obj_v[t] <- obj_new <- artdeconv_obj_fun(Y, Y_new, Theta_new, P_new, m, n, k, m0, k0, Theta_0, meds, ranges, alpha1, alpha2, beta)
    res_v[t] <- norm(Y_new - Y, type = "F")/sqrt(m * n) ## storing all residuals of Y
    delt_obj <- abs(obj_new - obj_old)/obj_old # calculate the change in objective function values
    
    ## passing the updates to the next round
    Y_old <- Y_new
    Theta_old <- Theta_new
    P_old <- P_new
    obj_old <- obj_new
  }
  
  ## normalize the proportion matrix P so its row sums equal 1
  P_new <- apply(P_new, 2, function(x) x/sum(x))
  return(list(Y_hat = Y_new,
              Theta_hat = Theta_new,
              s_hat = s_old,
              P_hat = P_new,
              res_v = res_v,
              obj_v = obj_v,
              fixed_s = 'size-fixed',
              weights = ifelse(length(ranges) == 1, 'uniform', 'range-adaptive'),
              n_iter = t - 1))
}


#' The Bare Bone Version of The Single Run For Cross Validation

artdeconv_single_solve_for_cv <- function(Y, Theta_0, Theta_it, s_it, P_it, m, n, k, m0, k0, meds, ranges, alpha1, alpha2, beta, max_iter = 1e5, tol = 1e-5) {
  
  ## pass the initial values to "old" variables
  Theta_old <- Theta_it
  s_old <- s_it
  P_old <- P_it
  Y_old <- Theta_old %*% diag(s_old) %*% P_old
  obj_old <- artdeconv_obj_fun(Y, Y_old, Theta_old, P_old, m, n, k, m0, k0, Theta_0, meds, ranges, alpha1, alpha2, beta)
  Delta <- get_Delta(k, k0, m, m0)
  Delta_c <- 1 - Delta
  n_iter <- 0
  delt_obj <- Inf
  
  ## start the loop
  for (t in seq_len(max_iter)) {
    if (delt_obj <= tol) break
    Theta_new <- mu_Theta(Y, s_old, P_old, Theta_0, Delta, Delta_c, Theta_old, m, n, alpha1, alpha2)
    P_new <- mu_P(Y, Theta_new, s_old, meds, P_old, m, n, beta, ranges)
    s_new <- mu_s(Y, Theta_new, s_old, P_new)
    Y_new <- Theta_new %*% diag(s_new) %*% P_new
    obj_new <- artdeconv_obj_fun(Y, Y_new, Theta_new, P_new, m, n, k, m0, k0, Theta_0, meds, ranges, alpha1, alpha2, beta)
    delt_obj <- abs(obj_new - obj_old)/obj_old
    
    ## passing the updates to the next round
    Theta_old <- Theta_new
    s_old <- s_new
    P_old <- P_new
    obj_old <- obj_new
  }
  
  ## normalize the proportion matrix P so its row sums equal 1
  P_new <- apply(P_new, 2, function(x) x/sum(x))
  return(list(Y_hat = Y_new,
              Theta_hat = Theta_new,
              s_hat = s_new,
              P_hat = P_new))
}

#' The Bare Bone Version of The Bi-factor Model (with Fixed S) Single Run For Cross Validation

artdeconv_single_solve_s_fixed_for_cv <- function(Y, Theta_0, Theta_it, P_it, m, n, k, m0, k0, meds, ranges, alpha1, alpha2, beta, max_iter = 1e5, tol = 1e-5) {
  
  ## pass the initial values to "old" variables
  Theta_old <- Theta_it
  s_old <- rep(1, k)
  P_old <- P_it
  Y_old <- Theta_old %*% P_old
  obj_old <- artdeconv_obj_fun(Y, Y_old, Theta_old, P_old, m, n, k, m0, k0, Theta_0, meds, ranges, alpha1, alpha2, beta)
  Delta <- get_Delta(k, k0, m, m0)
  Delta_c <- 1 - Delta
  n_iter <- 0
  delt_obj <- Inf
  
  ## start the loop
  for (t in seq_len(max_iter)) {
    if (delt_obj <= tol) break
    Theta_new <- mu_Theta(Y, s_old, P_old, Theta_0, Delta, Delta_c, Theta_old, m, n, alpha1, alpha2)
    P_new <- mu_P(Y, Theta_new, s_old, meds, P_old, m, n, beta, ranges)
    Y_new <- Theta_new %*% P_new
    obj_new <- artdeconv_obj_fun(Y, Y_new, Theta_new, P_new, m, n, k, m0, k0, Theta_0, meds, ranges, alpha1, alpha2, beta)
    delt_obj <- abs(obj_new - obj_old)/obj_old
    
    ## passing the updates to the next round
    Theta_old <- Theta_new
    P_old <- P_new
    obj_old <- obj_new
  }
  
  ## normalize the proportion matrix P so its row sums equal 1
  P_new <- apply(P_new, 2, function(x) x/sum(x))
  return(list(Y_hat = Y_new,
              Theta_hat = Theta_new,
              s_hat = s_old,
              P_hat = P_new))
}