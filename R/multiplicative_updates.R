#' function to update the proportion matrix through multiplicative update
#' 
#' @param Y the bulk matrix
#' @param Theta the signature matrix (assumed fixed here)
#' @param s the diagonal "cell size" vector (assumed fixed here)
#' @param P_old the proportion matrix before the update
#' @param meds the vector of known medians
#' @param beta the tuning parameter for regularizing P
#' @param wt the inverse weights of cell types; must be a a vector of length k if using non-uniform weights or a scalar if using uniform weights
#' 
#' @return the $K \times n$ proportion matrix after the current update

mu_P <- function(Y, Theta, s, meds, P_old, m, n, beta, wt = 1) {
  # calculating parts
  G = Theta %*% diag(s) 
  num_matrix = t(G) %*% Y + m * n * beta * (1/wt) * meds
  denom_matrix = t(G) %*% G %*% P_old + m * n * beta * (1/wt) * P_old # (1/wt) * P = diag(1/wt) %*% P 
  # updates (using the element-wise products)
  P_new = P_old * num_matrix/denom_matrix
  return(P_new)
}

#' function to update the signature matrix through multiplicative update
#' 
#' @param Y the bulk matrix
#' @param s the diagonal "cell size" vector (assumed fixed here)
#' @param P the proportions matrix (assumed fixed here)
#' @param Theta_0 signature matrix for the fixed part of Theta
#' @param Delta indicator matrix for non-missing signature matrix expression values
#' @param Delta_c complement of the indicator matrix
#' @param m0 number of signature genes for well characterized cell types in a tissue
#' @param k0 number of well characterized cell types in a tissue
#' @param Theta_old the signature matrix before the update
#' @param alpha1 the tuning parameter for regularizing the known part of Theta
#' @param alpha2 the tuning parameter for regularizing the unknown part of Theta
#' 
#' @return the $m \times K$ signature matrix after the update

mu_Theta <- function(Y, s, P, Theta_0, Delta, Delta_c, Theta_old, m, n, alpha1, alpha2) {
  # calculating parts
  C = diag(s) %*% P 
  num_matrix = Y %*% t(C) + m * n * alpha1 * Delta * Theta_0
  denom_matrix = Theta_old %*% C %*% t(C) + m * n * (alpha1 * Delta * Theta_old  + alpha2 * Delta_c * Theta_old)
  # updating
  Theta_new = Theta_old * num_matrix/denom_matrix
  return(Theta_new)
}

#' function to update the size vector through multiplicative update
#' 
#' @param Y the bulk matrix
#' @param Theta the signature matrix (assumed fixed here)
#' @param s_old the "cell size" vector before the update
#' @param P the proportions matrix (assumed fixed here)
#' 
#' @return the "cell size" vector of length \eqn{K} after the update

mu_s <- function(Y, Theta, s_old, P) {
  k <- ncol(Theta)
  u <- vector(mode = 'double', length = k)
  for (i in seq(k)) u[i] <- sum(outer(Theta[, i], P[i, ]) * Y)
  V <- t(Theta) %*% Theta * P %*% t(P)
  s_new <- c(u * s_old / (V %*% s_old)) # updating by minimizing the surrogate function
  return(s_new)
}