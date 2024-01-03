#' Function to Generate One Set of Initial Values For The ARTdeConv Algorithm
#'
#' Generate one set of \eqn{(\mathbf{\Theta}, \mathbf{s}, \mathbf{P})} based on the problem sizes and 
#' known part of \eqn{\mathbf{\Theta}}. Details can be found in Liu et al.
#'
#' @param Y_star The bulk expression matrix.
#' @param Theta_0 The \eqn{m \times K} reference signature matrix. The last \eqn{K_0} columns for cell types with unknown reference gene expression should be padded with 0s.
#' @param meds A vector of length \eqn{K} of the pre-specified medians of cell type proportions.
#' @param m The total number of genes in the deconvolution.
#' @param k The total number of cell types (CTs) in the deconvolution.
#' @param m0 The total number of genes of the known CTs in the deconvolution (\eqn{\leq m}, the default is \eqn{m}).
#' @param k0 The total number of known CTs in the deconvolution (\eqn{\leq K}).
#' 
#' @return A list of the initial values for \eqn{(\mathbf{\Theta}, \mathbf{s}, \mathbf{P})}:
#'  * Theta: \eqn{\mathbf{\Theta}^0};
#'  * s: \eqn{\mathbf{s}^0};
#'  * P: \eqn{\mathbf{P}^0}.
#' 
#' @export

get_initials <- function(Y_star, Theta_0, meds, m, n, k, m0, k0) {
  
  ## get Theta
  Theta_it <- matrix(0, m, k)
  for (ii in seq(m)) Theta_it[ii, ] <- abs(rnorm(k, mean(Y_star[ii, ]), sd(Y_star[ii, ])))
  Theta_it[seq_len(m0), seq_len(k0)] <- Theta_0[seq_len(m0), seq_len(k0)]
  
  ## get P
  P_it <- matrix(rep(meds, n), k, n) + matrix(rnorm(k * n, 0, 0.1), k, n) # need to add noises so columns of P_it are not collinear
  P_it[P_it < 0] <- 0.01 # so that no initial of Y gets stuck in 
  
  ## get s (solve a quadratic problem of Theta_it and P_it)
  u <- vector(mode = 'double', length = k)
  for (i in seq(k)) u[i] <- sum(outer(Theta_it[, i], P_it[i, ]) * Y_star)
  V <- t(Theta_it) %*% Theta_it * P_it %*% t(P_it)
  s_it <- as.vector(solve(V + diag(1e-16, ncol(V))) %*% u) # make sure V is positive definite
  s_it[s_it < 0] <- mean(s_it[s_it > 0])
  return(list(Theta = Theta_it, 
              s = s_it, 
              P = P_it))
}