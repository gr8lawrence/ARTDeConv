#' Cross Validation Function For Finding Optimal Tuning Parameters For ARTdeConv
#' 
#' @param Y the bulk matrix
#' @param Theta_0 \eqn{m_0 \times K_0} reference matrix for the fixed part of the signature matrix
#' @param m0 the number of signature genes for well characterized cell types in a tissue
#' @param k0 the number of well characterized cell types in a tissue
#' @param meds a vector of length \eqn{K} of pre-specified medians of cell proportions 
#' @param ranges a vector of length \eqn{K} of pre-specified ranges of cell-type proportions (1/ranges as weights for regularization parameters)
#' @param alpha1_range a vector of tuning parameters for \eqn{\alpha_1} 
#' @param alpha2_range a vector of tuning parameters for \eqn{\alpha_2} 
#' @param beta_range a vector of tuning parameters for \eqn{\beta}
#' @param n_fold the number of folds of the cross validation; default is 5
#' @param n_start the number of restart for ARTdeConv in each cross validation fold; default is 10
#' @param parallel a logical value of whether to run ARTdeConv in parallel; default is `TRUE`
#' @param verbose a logical value: if `TRUE`, will print a message of all optimal tuning parameters after the cross validaiton is done; default is `TRUE`
#' @param ... other parameters for the [artdeconv] function
#' 
#' @importFrom nnls nnls
#' 
#' @seealso [artdeconv]
#' 
#' @export

cv_artdeconv <- function(Y, Theta_0, m0, k0, meds, ranges, alpha1_range, alpha2_range, beta_range, n_fold = 5, n_start = 10, parallel = TRUE, verbose = TRUE, ...) {
  ## get the dimensions from input matrices
  m = nrow(Y)
  n = ncol(Y)
  K = ncol(Theta_0)
  ## get all parameter combinations
  tune_grid = expand.grid(alpha1 = alpha1_range, alpha2 = alpha2_range, beta = beta_range)
  tune_grid$cv_error = NaN
  ## the loop for CV
  for (i in seq(nrow(tune_grid))) {
    # i = 1
    a1 = tune_grid$alpha1[i]
    a2 = tune_grid$alpha2[i]
    b = tune_grid$beta
    fold_id = sample(x = rep(seq(n_fold), ceiling(n/n_fold)), size = n) # assign fold id
    cv_error = 0 # initiate the CV error
    ## the loop to calculate total CV error for this param combination
    for (j in seq(n_fold)) {
      holdout_id = (fold_id == j) # locate the ids for the test set
      Y_test = Y[, holdout_id] # the test set
      Y_train = Y[, !holdout_id] # the training set
      train_res = artdeconv(Y = Y_train, Theta_0 = Theta_0, m0 = m0, k0 = k0, meds = meds, ranges = ranges, alpha1 = a1, alpha2 = a2, beta = b, n_start = n_start, parallel = parallel, ...)
      Theta_train = train_res$Theta_hat
      s_train = train_res$s_hat
      G_train = Theta_train %*% diag(s_train)
      P_test = apply(Y_test, 2, function(bb) {nnls::nnls(A = G_train, b = bb)$x}) # get P_test using NNLS since it and P_train share the same Theta and s
      cv_error = cv_error + norm(Y_test - G_train %*% P_test, type = 'F')^2
    }
    tune_grid$cv_error[i] = cv_error/n_fold # assign the averaged CV error
  }
  ## return the best parameter
  best_param_id = which(tune_grid$cv_error == min(tune_grid$cv_error))
  alpha1_best = tune_grid$alpha1[best_param_id]
  alpha2_best = tune_grid$alpha2[best_param_id]
  beta_best = tune_grid$beta[best_param_id]
  if (verbose) message(paste('The best tuning parameters are alpha1 = ', alpha1_best,
                             ', alpha2 = ', alpha2_best,
                             ', beta = ', beta_best,
                             '; they can be extracted from the returned list object.'))
  return(list(alpha1 = alpha1_best,
              alpha2 = alpha2_best,
              beta = beta_best))
}