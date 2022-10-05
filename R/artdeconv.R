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

artdeconv_single_solve <- function(Y, Theta_0, Theta_it, s_it, P_it, m, n, k, m0, k0, meds, ranges, alpha1, alpha2, beta, max_iter = 1e5, tol = 1e-6, fixed_s = FALSE) {
  ## get the initial Y norm
  Y_norm = norm(Y, type = "F")
  ## get the Delta matrix
  Delta = get_Delta(k, k0, m, m0)
  Delta_c = 1 - Delta
  ## initiate the counter
  n_iter = 0
  ## initiate a residual vector
  res_v = vector()
  ## initiate matrices
  Theta_old = Theta_it
  if (fixed_s) {
    message("s is fixed to be 1 for all its values")
    s_old = rep(1, k)
  } else {
    s_old = s_it 
  } 
  P_old = P_it
  Y_old = Theta_old %*% diag(s_old) %*% P_old 
  ## set the initial delt_Y to inf and a vector for all delt_Y values
  delt_Y = Inf
  delt_v = vector()
  ## vectors for saving other norms to see convergence
  P_norms = vector()
  Y_norms = vector()
  ## start the loop
  for (t in seq_len(max_iter)) {
    ## if delt_Y gets under the tolerance, break out of the loop
    if (delt_Y <= tol) break
    Theta_new = mu_Theta(Y, s_old, P_old, Theta_0, Delta, Delta_c, Theta_old, m, n, alpha1, alpha2)
    ## update P 
    P_new = mu_P(Y, Theta_new, s_old, meds, P_old, m, n, beta, ranges)
    ## updating s depending on whether we want s to be fixed
    if (!fixed_s) {
      s_new = mu_s(Y, Theta_new, s_old, P_new)
    } else {
      s_new = s_old
    }
    ## calculating the updated Y
    Y_new = Theta_new %*% diag(s_new) %*% P_new
    ## storing all needed results
    P_norms[t] = norm(P_new, type = "F")
    Y_norms[t] = Y_norm = norm(Y_new, type = "F")
    res_v[t] = res_new = norm(Y_new - Y_old, type = "F")
    ## calculating the convergence criterion
    delt_v[t] = delt_Y = res_new/Y_norm
    ## passing the updates to the next round
    Y_old = Y_new
    Theta_old = Theta_new
    s_old = s_new
    P_old = P_new
  }
  ## normalize the proportion matrix P so its row sums equal 1 
  P_new = apply(P_new, 2, function(x) x/sum(x))
  ## return the results
  return(list(Y_hat = Y_new,
              Theta_hat = Theta_new,
              s_hat = s_new,
              P_hat = P_new,
              res_v = res_v,
              delt_v = delt_v,
              P_norms = P_norms,
              Y_norms = Y_norms,
              fixed_s = ifelse(fixed_s, 'size-fixed', 'size-free'),
              weights = ifelse(length(ranges) == 1, 'uniform', 'range-adaptive'),
              n_iter = t - 1))
}

#' Adaptive regularized tri-factor NMF deconvolution
#' 
#' TODO: write a complete description of the package
#' 
#' @param Y the bulk matrix
#' @param Theta_0 \eqn{m_0 \times K_0} reference matrix for the fixed part of the signature matrix
#' @param m0 the number of signature genes for well characterized cell types in a tissue
#' @param k0 the number of well characterized cell types in a tissue
#' @param meds a vector of length \eqn{K} of pre-specified medians of cell proportions 
#' @param ranges a vector of length \eqn{K} of pre-specified ranges of cell-type proportions (1/`ranges` is used as weights for regularization parameters)
#' @param alpha1 the tuning parameter for regularizing the known part of Theta
#' @param alpha2 the tuning parameter for regularizing the unknown part of Theta
#' @param beta the tuning parameter for regularizing P
#' @param n_start an integer indicating the number of restarts of the algorithm; default value is 10
#' @param parallel a logical value indicating whether or not to run [artdeconv()] in parallel; default is `TRUE` 
#' @param ... other parameters that can be passed to the function (see [artdeconv_single_solve()])
#' 
#' @seealso [artdeconv_single_solve()]
#' 
#' @return A list with following items
#' * Y_hat: the estimated Y from the product of estimated factors
#' * Theta_hat: the estimated signature matrix from deconvolution
#' * s_hat: the estimated "cell size" vector from deconvolution
#' * P_hat the estimated Theta from deconvolution
#' 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar% 
#' @export
artdeconv <- function(Y, Theta_0, m0, k0, meds, ranges, alpha1, alpha2, beta, n_start = 10, parallel = TRUE, ...) {
  ## parameters
  m = nrow(Y)
  n = ncol(Y)
  k = ncol(Theta_0)
  
  ## the parallel version 
  if (parallel) {
    ## run all restarts in parallel
    rl = foreach::foreach(i = seq_len(n_start)) %dopar% {
      initials = get_initials(Y, Theta_0, meds, m, n, k, m0, k0)
      Theta_it = initials$Theta
      s_it = initials$s
      P_it = initials$P
      mu_res = artdeconv_single_solve(Y, Theta_0, Theta_it, s_it, P_it, m, n, k, m0, k0, meds, ranges, alpha1, alpha2, beta, ...)
      mu_res
    }
  } else {
    ## a list to save all results
    rl = vector('list', n_start)
    
    ## start the iterates
    for (i in seq_len(n_start)) {
      initials = get_initials(Y, Theta_0, meds, m, n, k, m0, k0)
      Theta_it = initials$Theta
      s_it = initials$s
      P_it = initials$P
      rl[[i]] = artdeconv_single_solve(Y, Theta_0, Theta_it, s_it, P_it, m, n, k, m0, k0, meds, ranges, alpha1, alpha2, beta, ...)
    }
  }
  
  ## get the residuals of Y of each iteration
  itr_res = sapply(rl, function(x) norm(Y - x$Y_hat, type = 'F'))
  
  ## find the best iteration by argmin(Yr[i])
  best_ind = which(itr_res == min(itr_res))
  best_itr = rl[[best_ind[1]]]
  return(best_itr)
}
