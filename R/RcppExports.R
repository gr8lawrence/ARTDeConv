# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Multiplicative Update of P
mu_P_cpp <- function(Y, Theta, s, meds, P, beta, wt) {
    .Call(`_ARTdeConv_mu_P_cpp`, Y, Theta, s, meds, P, beta, wt)
}

#' Multiplicative Update of Theta
mu_Theta_cpp <- function(Y, s, P, Theta_0, Delta, Delta_c, Theta, alpha1, alpha2) {
    .Call(`_ARTdeConv_mu_Theta_cpp`, Y, s, P, Theta_0, Delta, Delta_c, Theta, alpha1, alpha2)
}

#' Multiplicative Update of S
mu_s_cpp <- function(Y, Theta, s, P) {
    .Call(`_ARTdeConv_mu_s_cpp`, Y, Theta, s, P)
}

#' Function for getting the Delta matrix
get_Delta_cpp <- function(k, k0, m, m0) {
    .Call(`_ARTdeConv_get_Delta_cpp`, k, k0, m, m0)
}

#' Function for getting the complement of the Delta matrix
get_Delta_c_cpp <- function(Delta) {
    .Call(`_ARTdeConv_get_Delta_c_cpp`, Delta)
}

#' ARTdeConv objective function
obj_fun_cpp <- function(Y, Y_hat, Theta_hat, P_hat, m0, k0, Theta_0, Delta, Delta_c, meds, ranges, alpha_1, alpha_2, beta) {
    .Call(`_ARTdeConv_obj_fun_cpp`, Y, Y_hat, Theta_hat, P_hat, m0, k0, Theta_0, Delta, Delta_c, meds, ranges, alpha_1, alpha_2, beta)
}

#' Core ARTdeConv Function For One Set of Initial Values (in C++)
#'
#' This is the core function that will run ARTdeConv once under a set of initial values. It is written in C++ 
#' and integrated in R through Rcpp and RcppArmadillo. Users running this function should provide their own 
#' initial values in the input. Otherwise, the required parameters are the same as the main ARTdeConv function
#' with restarts.
#'
#' @param Y The bulk matrix.
#' @param Theta_0 The \eqn{m \times K} reference signature matrix. The last \eqn{K_0} columns for cell types with unknown reference gene expression should be padded with 0s. 
#' @param Theta_it The initial value of Theta.
#' @param P_it The initial value of P.
#' @param m0 The number of gene features we have prior knowledge about in `Theta_0`. 
#' @param k0 The number of cell types whose cell-type level expression we have prior knowledge about in `Theta_0`.
#' @param meds The vector of median cell proportions in the population.
#' @param ranges The vector of the range interval lengths of proportions of different cell types.
#' @param alpha1 The tuning parameter for regularizing the part of Theta with prior knowledge.
#' @param alpha2 The tuning parameter for regularizing the part of Theta without prior knowledge.
#' @param beta The tuning parameter for regularizing P.
#' @param max_iter The maximal number of iterations this core function will run. The default is `1e5`.
#' @param tol The tolerance parameter for the convergence criterion of ARTdeConv. The default is `1e-5`.
#'
#' @return A list with the following items:
#' * Y_hat: the estimated bulk expression from deconvolution; 
#' * Theta_hat: the estimated signature matrix from deconvolution;
#' * s_hat: the vector of estimated mRNA amounts;
#' * P_hat: the estimated proportion (the results of primary interests);
#' * res_v: the vector of bulk expression residuals after each iteration;
#' * obj_v: the vector of objective function values after each iteration;
#' * fixed_s: the indicator of whether the mRNA amount vector is fixed to be 1 for all cell types;
#' * weights: the indicator of whether adaptive penalization is used;
#' * n_iter: the number of iterations the algorithm has run.
#' 
#' @export
artdeconv_single_solve_cpp <- function(Y, Theta_0, Theta_it, s_it, P_it, m0, k0, meds, ranges, alpha1, alpha2, beta, max_iter = 1e5L, tol = 1e-5) {
    .Call(`_ARTdeConv_artdeconv_single_solve_cpp`, Y, Theta_0, Theta_it, s_it, P_it, m0, k0, meds, ranges, alpha1, alpha2, beta, max_iter, tol)
}

#' Core ARTdeConv Function For One Set of Initial Values Assuming a Fixed S (in C++)
#'
#' This is the core function that will run ARTdeConv once under a set of initial values, but with the legacy bi-factor assumption (that the scale matrix S is fixed). 
#' It is written in C++ and integrated in R through Rcpp and RcppArmadillo. It can technically increase the speed since S will not be updated,
#' but should be only used when the user is confident that S can be represented by the identity matrix. Users running this function should provide their own 
#' initial values in the input. Otherwise, the required parameters are the same as the main ARTdeConv function
#' with restarts.
#'
#' @param Y The bulk matrix.
#' @param Theta_0 The \eqn{m \times K} reference signature matrix. The last \eqn{K_0} columns for cell types with unknown reference gene expression should be padded with 0s. 
#' @param Theta_it The initial value of Theta.
#' @param P_it The initial value of P.
#' @param m0 The number of gene features we have prior knowledge about in `Theta_0`. 
#' @param k0 The number of cell types whose cell-type level expression we have prior knowledge about in `Theta_0`.
#' @param meds The vector of median cell proportions in the population.
#' @param ranges The vector of the range interval lengths of proportions of different cell types.
#' @param alpha1 The tuning parameter for regularizing the part of Theta with prior knowledge.
#' @param alpha2 The tuning parameter for regularizing the part of Theta without prior knowledge.
#' @param beta The tuning parameter for regularizing P.
#' @param max_iter The maximal number of iterations this core function will run. The default is `1e5`.
#' @param tol The tolerance parameter for the convergence criterion of ARTdeConv. The default is `1e-5`.
#'
#' @return A list with the following items:
#' * Y_hat: the estimated bulk expression from deconvolution; 
#' * Theta_hat: the estimated signature matrix from deconvolution;
#' * s_hat: the vector of estimated mRNA amounts;
#' * P_hat: the estimated proportion (the results of primary interests);
#' * res_v: the vector of bulk expression residuals after each iteration;
#' * obj_v: the vector of objective function values after each iteration;
#' * fixed_s: the indicator of whether the mRNA amount vector is fixed to be 1 for all cell types;
#' * weights: the indicator of whether adaptive penalization is used;
#' * n_iter: the number of iterations the algorithm has run.
#' 
#' @export
artdeconv_single_solve_s_fixed_cpp <- function(Y, Theta_0, Theta_it, P_it, m0, k0, meds, ranges, alpha1, alpha2, beta, max_iter = 1e5L, tol = 1e-5) {
    .Call(`_ARTdeConv_artdeconv_single_solve_s_fixed_cpp`, Y, Theta_0, Theta_it, P_it, m0, k0, meds, ranges, alpha1, alpha2, beta, max_iter, tol)
}

