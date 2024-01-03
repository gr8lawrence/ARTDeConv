#' @useDynLib ARTdeConv, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Adaptive Regularized Tri-factor NMF Deconvolution
#' 
#' This function performs deconvolution of the bulk expression data using the ARTdeConv algorithm. 
#' The user should input 
#' 
#' @param Y the bulk matrix.
#' @param Theta_0 \eqn{m_0 \times K_0} partial reference matrix that we have prior knowledge about.
#' @param m0 the number of signature genes for well characterized cell types in a tissue.
#' @param k0 the number of well characterized cell types in a tissue.
#' @param meds a vector of length \eqn{K} of pre-specified medians of cell proportions.
#' @param ranges a vector of length \eqn{K} of pre-specified ranges of cell-type proportions (1/`ranges` is used as weights for regularization parameters).
#' @param alpha1 the tuning parameter for regularizing the known part of Theta.
#' @param alpha2 the tuning parameter for regularizing the unknown part of Theta.
#' @param beta the tuning parameter for regularizing P.
#' @param n_start an integer indicating the number of restarts of the algorithm; default value is 10.
#' @param parallel a logical value indicating whether or not to run [artdeconv] in parallel; default is `TRUE`. 
#' @param s_fixed a logical value indicating whether the cell size matrix \eqn{S} is coerced to be the identity matrix; default is `FALSE`.
#' @param ... other parameters that can be passed to the function (see [artdeconv_single_solve_cpp]).
#' 
#' @seealso [artdeconv_single_solve_cpp]
#' 
#' @return A list with following items:
#' * Y_hat: the estimated Y from the product of estimated factors;
#' * Theta_hat: the estimated signature matrix from deconvolution;
#' * s_hat: the estimated "cell size" vector from deconvolution;
#' * P_hat the estimated Theta from deconvolution.
#' 
#' @importFrom foreach foreach
#' @importFrom foreach %dopar% 
#' @export

artdeconv <- function(Y, Theta_0, m0, k0, meds, ranges, alpha1, alpha2, beta, n_start = 10, parallel = TRUE, s_fixed = FALSE, ...) {
  ## parameters
  m <- nrow(Y)
  n <- ncol(Y)
  k <- ncol(Theta_0)
  
  ## the parallel version 
  if (parallel) {
    ## run all restarts in parallel
    rl <- foreach::foreach(i = seq_len(n_start)) %dopar% {
      initials <- get_initials(Y, Theta_0, meds, m, n, k, m0, k0)
      if (s_fixed) {
        mu_res <- artdeconv_single_solve_s_fixed_cpp(Y, Theta_0, initials$Theta, initials$P, m0, k0, meds, ranges, alpha1, alpha2, beta, ...)
      } else {
        mu_res <- artdeconv_single_solve_cpp(Y, Theta_0, initials$Theta, initials$s, initials$P, m0, k0, meds, ranges, alpha1, alpha2, beta, ...)
      }
      mu_res
    }
  } else {
    ## a list to save all results
    rl <- vector('list', n_start)
    
    ## start the iterates
    for (i in seq_len(n_start)) {
      initials <- get_initials(Y, Theta_0, meds, m, n, k, m0, k0)
      if (s_fixed) {
        rl[[i]] <- artdeconv_single_solve_s_fixed_cpp(Y, Theta_0, initials$Theta, initials$P, m0, k0, meds, ranges, alpha1, alpha2, beta, ...)
      } else {
        rl[[i]] <- artdeconv_single_solve_cpp(Y, Theta_0, initials$Theta, initials$s, initials$P, m0, k0, meds, ranges, alpha1, alpha2, beta, ...)
      }
    }
  }
  ## get the residuals of Y of each iteration
  itr_res <- sapply(rl, function(x) norm(Y - x$Y_hat, type = 'F'))
  
  ## find the best iteration by argmin(Yr[i])
  best_ind <- which(itr_res == min(itr_res))
  best_itr <- rl[[best_ind[1]]]
  return(best_itr)
}
