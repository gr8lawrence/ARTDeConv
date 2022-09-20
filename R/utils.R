library(extraDistr)
# Guide matrix (Delta) and related ----------------------------------------
### return a delta matrix according to input parameters
get_Delta <- function(k, k0, m, m0) {
  Delta = rbind(
    matrix(rep(c(rep(1, k0), rep(0, k - k0)), m0), m0, k, byrow=TRUE),
    matrix(0, m - m0, k)
  )
  return(Delta)
}


# Function to generate initial values -------------------------------------
get_initials <- function(Y_star, Theta_0, meds, m, n, k, m0, k0) {
  ## get Theta
  Theta_it = matrix(0, m, k)
  for (ii in seq(m)) Theta_it[ii, ] = abs(rnorm(k, mean(Y_star[ii, ]), sd(Y_star[ii, ])))
  Theta_it[seq_len(m0), seq_len(k0)] = Theta_0[seq_len(m0), seq_len(k0)]
  
  ## get s
  s_it = rep(1, k) # TODO: may need a better initializer
  
  ## get P
  P_it = matrix(rep(meds, n), k, n) + matrix(rnorm(k * n, 0, 0.1), k, n) # need to add noises so columns of P_it are not collinear
  P_it[P_it < 0] = 0.01 # so that no initial of Y gets stuck in 
  
  return(list(Theta = Theta_it, s = s_it, P = P_it))
}


# Function to calculate matrix inner products -----------------------------

## IN:
## A, B - two matrices of the same dimensions
## OUT:
## inner - the inner product of A and B
## Exit and returns an error message if:
## A and B are not of the same dimensions
inner_product_TL <- function(A, B) {
  if (prod(dim(A) == dim(B)) == 0) stop("Dimensions of two matrices must be the same")
  inner = sum(A * B)
  return(inner)
}
