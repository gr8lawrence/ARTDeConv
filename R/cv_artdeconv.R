#' Cross Validation Function For Finding Optimal Tuning Parameters For ARTdeConv
#' 
#' @param Y The bulk matrix.
#' @param Theta_0 The \eqn{m \times K} reference signature matrix. The last \eqn{K_0} columns for cell types with unknown reference gene expression should be padded with 0s.
#' @param m0 The number of signature genes for well characterized cell types in a tissue. Must satisfy \eqn{m_0 \leq m}. The default value is \eqn{m} (`nrow(Theta_0)`).
#' @param k0 The number of well characterized cell types in a tissue. Must satisfy \eqn{K_0 < K}.
#' @param meds A vector of length \eqn{K} of pre-specified medians of cell proportions.
#' @param ranges A vector of length \eqn{K} of pre-specified ranges of cell-type proportions (1/ranges as weights for regularization parameters).
#' @param alpha1_range A vector of tuning parameters for \eqn{\alpha_1}.
#' @param alpha2_range A vector of tuning parameters for \eqn{\alpha_2}.
#' @param beta_range A vector of tuning parameters for \eqn{\beta}.
#' @param n_fold The number of folds of the cross validation; default is 5.
#' @param n_start The number of restart for ARTdeConv in each cross validation fold; default is 10.
#' @param parallel A logical value of whether to run ARTdeConv in parallel; default is `TRUE`.
#' @param verbose A logical value: if `TRUE`, will print a message of all optimal tuning parameters after the cross validaiton is done; default is `TRUE`.
#' @param seed A numerical seed to ensure the reproducibility of results; default is 100.
#' @param seed_control A logical value dictating whether seed is controlled in the foreach loop during CV; default is `TRUE`.
#' @param ... Other parameters for the [artdeconv] function.
#' 
#' @return A list with following items:
#' * alpha_1: the chosen \eqn{\alpha_1} value from the CV;
#' * alpha_2: the chosen \eqn{\alpha_2} value from the CV;
#' * beta: the chosen \eqn{\beta} value from the CV;
#' * error: the CV error under the best tuning parameter selections. 
#' 
#' @importFrom nnls nnls
#' @importFrom progress progress_bar
#' 
#' @seealso [artdeconv]
#' 
#' @export

cv_artdeconv <- function(Y, Theta_0, m0 = nrow(Theta_0), k0, meds, ranges, alpha1_range, alpha2_range, beta_range, n_fold = 5, n_start = 10, parallel = TRUE, verbose = TRUE, seed = 100, seed_control = TRUE, ...) {
  set.seed(seed)
  
  ## get all parameter combinations
  tune_grid <- expand.grid(alpha1 = alpha1_range, alpha2 = alpha2_range, beta = beta_range)
  
  ## print the total number of grid value combinations if verbose = TRUE
  if (verbose) message(
    paste("In the tuning grid, there are", 
          length(alpha1_range), "value(s) for alpha1,", 
          length(alpha2_range), "for alpha2, and",
          length(beta_range), "for beta, for a total of",
          nrow(tune_grid), "combinations.")
  )
  
  ## determine the fold ids
  fold_id <- sample(x = rep(seq(n_fold), ceiling(ncol(Y)/n_fold)), size = ncol(Y)) # assign fold id
  
  ## the outer loop for each combination of tuning param
  message('Begin cross validation...')
  cv_errors_matrix <- foreach::foreach(i = 1:nrow(tune_grid), .combine = 'cbind') %:%  
    foreach::foreach(j = 1:n_fold, .combine = 'c') %dopar% {
      if (seed_control) set.seed(seed + j + (i - 1) * n_fold)
      ## the inner loop to calculate total CV error for this param combination
      holdout_id <- (fold_id == j) # locate the ids for the test set
      train_res <- artdeconv(Y = Y[, !holdout_id], Theta_0 = Theta_0, m0 = m0, k0 = k0, meds = meds, ranges = ranges, alpha1 = tune_grid$alpha1[i], alpha2 = tune_grid$alpha2[i], beta = tune_grid$beta[i], n_start = n_start, parallel = FALSE, ...)
      G_train <- train_res$Theta_hat %*% diag(c(train_res$s_hat))
      P_test <- apply(Y[, holdout_id], 2, function(bb) {nnls::nnls(A = G_train, b = bb)$x}) # get P_test using NNLS since it and P_train share the same Theta and s
      cv_error <- norm(Y[, holdout_id] - G_train %*% P_test, type = 'F')^2/(n_fold * nrow(Y) * ncol(Y))
      cv_error
    } 
  message('Finished cross validation...')
  
  ## return the best parameter
  cv_errors <- colSums(cv_errors_matrix) # the sum of each column corresponds to the total CV error for a combination of pars
  best_param_id <- which(cv_errors == min(cv_errors))
  alpha1_best <- tune_grid$alpha1[best_param_id]
  alpha2_best <- tune_grid$alpha2[best_param_id]
  beta_best <- tune_grid$beta[best_param_id]
  if (verbose) message(paste('The best tuning parameters are alpha1 = ', alpha1_best,
                             ', alpha2 = ', alpha2_best,
                             ', beta = ', beta_best,
                             '; they can be extracted from the returned list object.'))
  return(list(alpha1 = alpha1_best,
              alpha2 = alpha2_best,
              beta = beta_best,
              error = min(cv_errors)))
}