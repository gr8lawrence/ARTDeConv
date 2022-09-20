source("a02_utility_fun.R")

# Update functions --------------------------------------------------------

## function to update the proportion matrix through minimizing the surrogate function
## IN: 
## Y - the bulk matrix
## Theta - the signature matrix
## s - the diagonal "size" vector
## P_old - the proportion matrix before the update
## meds - the vector of known medians
## beta - the tuning parameter for regularizing P
## wt - the inverse weights of cell types; must be a scalar if using uniform weights or a vector of length k
## OUT: 
## P_new - the proportion matrix after the current update
MU_update_P_adaptive_TL <- function(Y, Theta, s, meds, P_old, m, n, beta, wt = 1) {
  # calculating parts
  G = Theta %*% diag(s) 
  num_matrix = t(G) %*% Y + m * n * beta * (1/wt) * meds
  denom_matrix = t(G) %*% G %*% P_old + m * n * beta * (1/wt) * P_old # (1/wt) * P = diag(1/wt) %*% P 
  # updates (using the element-wise products)
  P_new = P_old * num_matrix/denom_matrix
  return(P_new)
}

## function to update the signature matrix through minimizing the surrogate function
## REQUIRE:
## get_Delta - function for getting the Delta-indicator matrix
## IN: 
## Y - the bulk matrix
## s - the diagonal "size" vector
## P - the proportions matrix
## Theta_0 - reference matrix for the fixed part of Theta
## Delta - delta indicator matrix
## Delta_c - complement of the delta indicator matrix
## m0 - number of signature genes for well characterized cell types in a tissue
## k0 - number of well characterized cell types in a tissue
## Theta_old - the signature matrix before the update
## alpha1 - the tuning parameter for regularizing the known part of Theta
## alpha2 - the tuning parameter for regularizing the unknown part of Theta
## OUT: 
## Theta_new - the signature matrix after the update
MU_update_Theta_TL <- function(Y, s, P, Theta_0, Delta, Delta_c, Theta_old, m, n, alpha1, alpha2) {
  # calculating parts
  C = diag(s) %*% P 
  num_matrix = Y %*% t(C) + m * n * alpha1 * Delta * Theta_0
  denom_matrix = Theta_old %*% C %*% t(C) + m * n * (alpha1 * Delta * Theta_old  + alpha2 * Delta_c * Theta_old)
  # updating
  Theta_new = Theta_old * num_matrix/denom_matrix
  return(Theta_new)
}

## function to update the size vector through minimizing the surrogate function for all K cell types
## IN: 
## Y - the bulk matrix
## Theta - the signature matrix
## s_old - the size vector before the update
## P - the proportions matrix
## OUT: 
## s_new - the size vector after the update
MU_update_s_TL <- function(Y, Theta, s_old, P) {
  k = ncol(Theta)
  Z_ls = lapply(seq_len(k), function(i) outer(Theta[, i], P[i, ])) # a list of rank-1 outer products Z_i
  u = sapply(Z_ls, inner_product_TL, B = Y) # a vector of inner products between Y and Z_i
  Z_ls_vec = sapply(Z_ls, c) # vectorize Z_i and put the vectors into a matrix
  V = Reduce('+', lapply(seq_len(nrow(Z_ls_vec)), function(i) outer(Z_ls_vec[i, ], Z_ls_vec[i, ]))) # sum over all outer products of the rows of the above matrix
  s_new = c(u * s_old / (V %*% s_old)) # updating by minimizing the surrogate function
  return(s_new)
}

# for testing
# A = matrix(1:6, 2, 3)
# B = matrix(7:12, 2, 3)
# s_old = c(1, 1)
# k = 2
# Z_ls = lapply(seq_len(k), function(i) outer(A[, i], B[i, ]))
# V2 = matrix(0, k, k)
# for (ii in seq_len(k)) {
#   for (jj in ii:k) {
#     V2[ii, jj] = V2[jj, ii] = sum(Z_ls[[ii]] * Z_ls[[jj]])
#   }
# }
# V2 should be the same as V

# New main algorithm ------------------------------------------------------

## function to update the signature matrix through minimizing the surrogate function
## REQUIRE:
## get_Delta - function for getting the Delta-indicator matrix
## update functions for Theta, s, and P
## IN: 
## Y - the bulk matrix
## Theta_0 - reference matrix for the fixed part of Theta
## Theta_it - the initial reference matrix
## P_it - the proportions matrix
## m0 - number of signature genes for well characterized cell types in a tissue
## k0 - number of well characterized cell types in a tissue
## meds - pre-specified medians of cell-type proportions
## ranges - pre-specified ranges of cell-type proportions (1/ranges as weights for penalization on prop)
## alpha1 - the tuning parameter for regularizing the known part of Theta
## alpha2 - the tuning parameter for regularizing the unknown part of Theta
## beta - the tuning parameter for regularizing P
## max_iter - maximal number of iterations of the algorithm
## tol - tolerance of difference for indicating convergence
## OUT: 
## Theta_new - the signature matrix after the update
MU_solve <- function(Y, Theta_0, Theta_it, s_it, P_it, m, n, k, m0, k0, meds, ranges, alpha1, alpha2, beta, max_iter = 1e5, tol = 1e-6, fixed_s = FALSE) {
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
    Theta_new = MU_update_Theta_TL(Y, s_old, P_old, Theta_0, Delta, Delta_c, Theta_old, m, n, alpha1, alpha2)
    ## update P 
    P_new = MU_update_P_adaptive_TL(Y, Theta_new, s_old, meds, P_old, m, n, beta, ranges)
    ## updating s depending on whether we want s to be fixed
    if (!fixed_s) {
      s_new = MU_update_s_TL(Y, Theta_new, s_old, P_new)
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

# The main function that includes restarts (does not include tuning)
main_solve <- function(Y, Theta_0, m0, k0, meds, ranges, alpha1, alpha2, beta, n_start = 1, parallel = TRUE, ...) {
  ## parameters
  m = nrow(Y)
  n = ncol(Y)
  k = ncol(Theta_0)
  
  ## the parallel version 
  if (parallel) {
    ## load required packages
    require(doParallel)
    require(foreach)
    
    ## run all restarts in parallel
    # registerDoParallel(4)
    rl = foreach(i = seq_len(n_start)) %dopar% {
      initials = get_initials(Y, Theta_0, meds, m, n, k, m0, k0)
      Theta_it = initials$Theta
      s_it = initials$s
      P_it = initials$P
      mu_res = MU_solve(Y, Theta_0, Theta_it, s_it, P_it, m, n, k, m0, k0, meds, ranges, alpha1, alpha2, beta, ...)
      mu_res
    }
    # registerDoSEQ()
  } else {
    ## a list to save all results
    rl = vector('list', n_start)
    
    ## start the iterates
    for (i in seq_len(n_start)) {
      initials = get_initials(Y, Theta_0, meds, m, n, k, m0, k0)
      Theta_it = initials$Theta
      s_it = initials$s
      P_it = initials$P
      rl[[i]] = MU_solve(Y, Theta_0, Theta_it, s_it, P_it, m, n, k, m0, k0, meds, ranges, alpha1, alpha2, beta, ...)
    }
  }
  
  ## get the residuals of Y of each iteration
  itr_res = sapply(rl, function(x) norm(Y - x$Y_hat, type = 'F'))
  # print(itr_res)
  
  ## find the best iteration by argmin(Yr[i])
  best_ind = which(itr_res == min(itr_res))
  best_itr = rl[[best_ind[1]]]
  return(best_itr)
}
