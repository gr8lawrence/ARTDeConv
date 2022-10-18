#' Objective function for the ARTDeConv algorithm

artdeconv_obj_fun <- function(Y, Y_hat, Theta_hat, P_hat, m, n, k, m0, k0, Theta_0, meds, ranges, alpha_1, alpha_2, beta) {
  Delta <- get_Delta(k, k0, m, m0)
  return(1/(2 * m * n) * norm(Y - Y_hat, type = 'F')^2 + (1/2) * alpha_1 * norm(Delta * Theta_hat - Theta_0, type = 'F')^2 + 
           alpha_2 * norm((1 - Delta) * Theta_hat, type = 'F')^2 + (1/2) * beta * norm((1/ranges) *(P_hat - meds), type = 'F')^2)
}