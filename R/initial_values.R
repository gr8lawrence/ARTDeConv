#' Function to Generate One Set of Initial Values For The ARTdeConv Algorithm
#'
#' Generate one set of \eqn{(\mathbf{\Theta}, \mathbf{s}, \mathbf{P})} based on the problem sizes and 
#' known part of \eqn{\mathbf{\Theta}}. Details can be found in Liu et al.
#'
#' @param Y_star bulk expression matrix
#' @param Theta_0 the known part of \eqn{\mathbf{\Theta}}
#' @param meds a vector of length \eqn{K} of the pre-specified medians of cell type proportions
#' @param m the total number of genes in the deconvolution
#' @param k the total number of cell types (CTs) in the deconvolution
#' @param m0 the total number of genes of the known CTs in the deconvolution (\eqn{\leq m})
#' @param k0 the total number of known CTs in the deconvolution (\eqn{\leq K})
#' 
#' @return A list of the initial values for \eqn{(\mathbf{\Theta}, \mathbf{s}, \mathbf{P})}.

get_initials <- function(Y_star, Theta_0, meds, m, n, k, m0, k0) {
  ## get Theta
  Theta_it = matrix(0, m, k)
  for (ii in seq(m)) Theta_it[ii, ] = abs(rnorm(k, mean(Y_star[ii, ]), sd(Y_star[ii, ])))
  Theta_it[seq_len(m0), seq_len(k0)] = Theta_0[seq_len(m0), seq_len(k0)]
  
  ## get P
  P_it = matrix(rep(meds, n), k, n) + matrix(rnorm(k * n, 0, 0.1), k, n) # need to add noises so columns of P_it are not collinear
  P_it[P_it < 0] = 0.01 # so that no initial of Y gets stuck in 
  
  ## get s (solve a quadratic problem of Theta_it and P_it)
  k = ncol(Theta_it)
  G_ls = lapply(seq_len(k), function(i) outer(Theta_it[, i], P_it[i, ])) # a list of rank-1 outer products Z_i
  u = sapply(G_ls, inner_product_TL, B = Y_star) # a vector of inner products between Y and Z_i
  G_ls_vec = sapply(G_ls, c) # vectorize Z_i and put the vectors into a matrix
  V = Reduce('+', lapply(seq_len(nrow(G_ls_vec)), function(i) outer(G_ls_vec[i, ], G_ls_vec[i, ]))) # sum over all outer products of the rows of the above matrix
  s_it = as.vector(solve(V + diag(1e-16, ncol(V))) %*% u) # make sure V is positive definite
  s_it[s_it < 0] = mean(s_it[s_it > 0])
  
  return(list(Theta = Theta_it, 
              s = s_it, 
              P = P_it))
}