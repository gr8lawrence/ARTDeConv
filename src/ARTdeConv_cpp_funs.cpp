// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
using namespace Rcpp;
using namespace arma;

//' Multiplicative Update of P
//' @export 
// [[Rcpp::export]]
arma::mat mu_P_cpp(const arma::mat Y, const arma::mat Theta, const arma::vec s, const arma::vec meds, arma::mat P, const double beta, const arma::vec wt) {
  
  // dimensions
  int m = Y.n_rows;
  int n = Y.n_cols;
  int K = Theta.n_cols;
  
  // create diagonal S matrix from s vector
  arma::mat S = arma::diagmat(s);
  
  // create the inverse wt vector and put it into the matrix
  arma::mat W_inv = arma::diagmat(1 / wt);
  
  // create G = Theta * S
  arma::mat G = Theta * S; 
  
  // calculate parts in the multiplicative update (by columns of Y and P)
  for(int i = 0; i < n; ++i) {
    P.col(i) = P.col(i) % (G.t() * Y.col(i) + m * n * beta * W_inv * meds) / (G.t() * G * P.col(i) + m * n * beta * W_inv * P.col(i));
  }
  return P;
}

//' Multiplicative Update of Theta
//' @export 
// [[Rcpp::export]]
arma::mat mu_Theta_cpp(const arma::mat Y, const arma::vec s, const arma::mat P, const arma::mat Theta_0, const arma::mat Delta, const arma::mat Delta_c, arma::mat Theta, const double alpha1, const double alpha2) {
  
  // dimensions
  int m = Y.n_rows;
  int n = Y.n_cols;
  int K = P.n_rows;
  
  // create diagonal S matrix from s vector
  arma::mat S = arma::diagmat(s);
  
  // create C = S * P
  arma::mat C = S * P;
  
  // calculate parts in the multiplicative update (by columns of Y and P)
  Theta = Theta % (Y * C.t() + m * n * alpha1 * Delta % Theta_0) / (Theta * C * C.t() + m * n * (alpha1 * Delta % Theta + alpha2 * Delta_c % Theta));
  return Theta;
}

//' Multiplicative Update of S
//' @export 
// [[Rcpp::export]]
arma::vec mu_s_cpp(const arma::mat Y, const arma::mat Theta, arma::vec s, const arma::mat P) {
  
  // dimensions
  int m = Y.n_rows;
  int n = Y.n_cols;
  int K = Theta.n_cols;
  
  // make a new matrix
  arma::vec u(K);
  arma::mat Z_k(m, n);
  
  // calculate the update for s
  for(int k = 0; k < K; ++k) {
    
    // compute Z_k as the outer product of The k-th col of Theta and k-th row of P
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        Z_k(i, j) = Theta(i, k) * P(k, j);
      }
    }
    u(k) = arma::trace(Y.t() * Z_k);
  }
  arma::mat V = (Theta.t() * Theta) % (P * P.t());
  s = s % u / (V * s);
  return s;
}

//' Function for getting the Delta matrix
//' @export 
// [[Rcpp::export]]
arma::mat get_Delta_cpp (const int k, const int k0, const int m, const int m0) {
  arma::mat Delta(m, k, fill::zeros);
  for(int i = 0; i < m0; ++i) {
    for(int j = 0; j < k0; ++j) {
      Delta(i, j) = 1;
    }
  }
  return Delta;
}

//' Function for getting the complement of the Delta matrix
//' @export 
// [[Rcpp::export]]
arma::mat get_Delta_c_cpp (const arma::mat Delta) {
  int m = Delta.n_rows;
  int k = Delta.n_cols;
  arma::mat Delta_c(m, k, fill::ones);
  Delta_c -= Delta;
  return Delta_c;
}

//' ARTdeConv objective function
//' @export 
// [[Rcpp::export]]
double obj_fun_cpp(const arma::mat Y, const arma::mat Y_hat, const arma::mat Theta_hat, const arma::mat P_hat, const int m0, const int k0, const arma::mat Theta_0, const arma::mat Delta, const arma::mat Delta_c, const arma::vec meds, const arma::vec ranges, const double alpha_1, const double alpha_2, const double beta) {
  
  // dimensions
  int m = Y.n_rows;
  int n = Y.n_cols;
  int k = Theta_hat.n_cols;
  
  // objective function without the P-part (which needs to be summed by column)
  double obj = 1/(2 * (double) m * (double) n) * pow(norm(Y - Y_hat, "fro"), 2) +
    0.5 * alpha_1 * pow(norm(Delta % Theta_hat - Theta_0, "fro"), 2) +
    0.5 * alpha_2 * pow(norm(Delta_c % Theta_hat, "fro"), 2);
  
  // calculate the objective function with the P-part
  arma::mat W_inv = arma::diagmat(1 / ranges);
  for(int i = 0; i < n; ++i) {
    obj += 0.5 * beta * pow(norm(W_inv * (P_hat.col(i) - meds), 2), 2);
  }
  return obj;
}

//' Core ARTdeConv Function For One Set of Initial Values (in C++)
//'
//' This is the core function that will run ARTdeConv once under a set of initial values. It is written in C++ 
//' and integrated in R through Rcpp and RcppArmadillo. Users running this function should provide their own 
//' initial values in the input. Otherwise, the required parameters are the same as the main ARTdeConv function
//' with restarts.
//'
//' @param Y the bulk matrix.
//' @param Theta_0 the partial reference matrix that we have prior knowledge about. The unknown part of 
//' @param Theta_it the initial value of Theta.
//' @param s_it the initial value of the main diagonal of S.
//' @param P_it the initial value of P.
//' @param m0 the number of gene features we have prior knowledge about in `Theta_0`. 
//' @param k0 the number of cell types whose cell-type level expression we have prior knowledge about in `Theta_0`.
//' @param meds the vector of median cell proportions in the population.
//' @param ranges the vector of the range interval lengths of proportions of different cell types.
//' @param alpha1 the tuning parameter for regularizing the part of Theta with prior knowledge.
//' @param alpha2 the tuning parameter for regularizing the part of Theta without prior knowledge.
//' @param beta the tuning parameter for regularizing P.
//' @param max_iter the maximal number of iterations this core function will run. The default is `1e5`.
//' @param tol the tolerance parameter for the convergence criterion of ARTdeConv. The default is `1e-5`.
//' @export 
// [[Rcpp::export]]
Rcpp::List artdeconv_single_solve_cpp(const arma::mat Y, const arma::mat Theta_0, const arma::mat Theta_it, const arma::vec s_it, const arma::mat P_it, const int m0, const int k0, const arma::vec meds, const arma::vec ranges, const double alpha1, const double alpha2, const double beta, const int max_iter = 1e5, const double tol = 1e-5) {
  
  // initiate the clock
  arma::wall_clock timer;
  timer.tic();

  // dimensions
  int m = Y.n_rows;
  int n = Y.n_cols;
  int k = Theta_it.n_cols;

  // get Delta and Delta_c
  arma::mat Delta = get_Delta_cpp(k, k0, m, m0);
  arma::mat Delta_c = get_Delta_c_cpp(Delta);

  // pass the initiate values to "old" values
  arma::mat Theta_old = Theta_it;
  arma::vec s_old = s_it;
  arma::mat S_old = arma::diagmat(s_old);
  arma::mat P_old = P_it;
  arma::mat Y_old = Theta_old * S_old * P_old;
  double obj_old = obj_fun_cpp(Y, Y_old, Theta_old, P_old, m0, k0, Theta_0, Delta, Delta_c, meds, ranges, alpha1, alpha2, beta);

  // declare some new values 
  arma::mat Theta_new(m, k);
  arma::vec s_new(k);
  arma::mat S_new(k, k, fill::zeros);
  arma::mat P_new(k, n);
  arma::mat Y_new(m, n);
  double obj_new = std::numeric_limits<double>::infinity();

  // initiate a residual and objective function vector
  arma::vec res_v(max_iter);
  arma::vec obj_v(max_iter);
  double delt_obj = std::numeric_limits<double>::infinity(); // set the initial delt_obj to infinity

  // start the loop
  int tt = 0; // declare the iteration counter for outputs
  for (int t = 1; t < max_iter; ++t) {
    if (delt_obj <= tol) {
      
      // if delt_Y gets under the tolerance, break out of the loop
      break; 
    } else {
      
      // run the MU updates
      Theta_new = mu_Theta_cpp(Y, s_old, P_old, Theta_0, Delta, Delta_c, Theta_old, alpha1, alpha2);
      P_new = mu_P_cpp(Y, Theta_new, s_old, meds, P_old, beta, ranges);
      s_new = mu_s_cpp(Y, Theta_new, s_old, P_new);
      S_new.diag() = s_new;
      Y_new = Theta_new * S_new * P_new; // calculating the updated Y
      
      // calculate the updated objective function and its change
      obj_new = obj_fun_cpp(Y, Y_new, Theta_new, P_new, m0, k0, Theta_0, Delta, Delta_c, meds, ranges, alpha1, alpha2, beta);
      obj_v(t) = obj_new; 
      res_v(t) = 1/(2 * (double) m * (double) n) * pow(norm(Y_new - Y, "fro"), 2); // storing all residuals of Y
      delt_obj = abs(obj_new - obj_old)/obj_old; // calculate the change in objective function values
      
      // passing the updates to the next round
      Y_old = Y_new;
      Theta_old = Theta_new;
      s_old = s_new;
      P_old = P_new;
      obj_old = obj_new;
    }
    tt += 1; // update the iteration counter for outputs
  }

  // normalize the proportion matrix P so its row sums equal 1
  double col_sum = 0;
  for(int i = 0; i < n; ++i) {
    col_sum = sum(P_new.col(i));
    for(int j = 0; j < k; ++j) {
      P_new(j, i) = P_new(j, i)/col_sum;
    }
  }

  // record the time and print it
  double time_elapsed = timer.toc();
  
  // check if the weights are uniform
  arma::vec uni_v(k, fill::ones);
  std::string weights = "range-adaptive";
  double diff = sum(ranges - uni_v);
  if (diff == 0) {
    weights = "uniform";
  }

  // write the list for returned values;
  Rcpp::List ret = Rcpp::List::create(Named("Y_hat") = Y_new, 
                                      _["Theta_hat"] = Theta_new,
                                      _["s_hat"] = s_new,
                                      _["P_hat"] = P_new,
                                      _["res_v"] = res_v,
                                      _["obj_v"] = obj_v,
                                      _["fixed_s"] = "size-free",
                                      _["weights"] = weights,
                                      _["n_iter"] = tt);
  return ret;
}

//' Core ARTdeConv Function For One Set of Initial Values Assuming a Fixed S (in C++)
//'
//' This is the core function that will run ARTdeConv once under a set of initial values, but with the legacy bi-factor assumption (that the scale matrix S is fixed). 
//' It is written in C++ and integrated in R through Rcpp and RcppArmadillo. It can technically increase the speed since S will not be updated,
//' but should be only used when the user is confident that S can be represented by the identity matrix. Users running this function should provide their own 
//' initial values in the input. Otherwise, the required parameters are the same as the main ARTdeConv function
//' with restarts.
//'
//' @param Y the bulk matrix.
//' @param Theta_0 the partial reference matrix that we have prior knowledge about. The unknown part of 
//' @param Theta_it the initial value of Theta.
//' @param P_it the initial value of P.
//' @param m0 the number of gene features we have prior knowledge about in `Theta_0`. 
//' @param k0 the number of cell types whose cell-type level expression we have prior knowledge about in `Theta_0`.
//' @param meds the vector of median cell proportions in the population.
//' @param ranges the vector of the range interval lengths of proportions of different cell types.
//' @param alpha1 the tuning parameter for regularizing the part of Theta with prior knowledge.
//' @param alpha2 the tuning parameter for regularizing the part of Theta without prior knowledge.
//' @param beta the tuning parameter for regularizing P.
//' @param max_iter the maximal number of iterations this core function will run. The default is `1e5`.
//' @param tol the tolerance parameter for the convergence criterion of ARTdeConv. The default is `1e-5`.
//' @export 
// [[Rcpp::export]]
Rcpp::List artdeconv_single_solve_s_fixed_cpp(const arma::mat Y, const arma::mat Theta_0, const arma::mat Theta_it, const arma::mat P_it, const int m0, const int k0, const arma::vec meds, const arma::vec ranges, const double alpha1, const double alpha2, const double beta, const int max_iter = 1e5, const double tol = 1e-5) {
  
  // initiate the clock
  arma::wall_clock timer;
  timer.tic();
  // cout << "Finished timer initiation. \n";
  
  // dimensions
  int m = Y.n_rows;
  int n = Y.n_cols;
  int k = Theta_it.n_cols;
  // cout << "Finished dimnesions. \n";
  
  // get Delta and Delta_c
  arma::mat Delta = get_Delta_cpp(k, k0, m, m0);
  arma::mat Delta_c = get_Delta_c_cpp(Delta);
  // cout << "Finished getting Delta. \n";
  
  // pass the initiate values to "old" values
  arma::mat Theta_old = Theta_it;
  arma::vec s_old(k, fill::ones);
  // arma::mat S_old = arma::diagmat(s_old);
  arma::mat P_old = P_it;
  arma::mat Y_old = Theta_old * P_old;
  double obj_old = obj_fun_cpp(Y, Y_old, Theta_old, P_old, m0, k0, Theta_0, Delta, Delta_c, meds, ranges, alpha1, alpha2, beta);
  // cout << "Finished getting old values. \n";
  
  // declare some new values 
  arma::mat Theta_new(m, k);
  arma::mat P_new(k, n);
  arma::mat Y_new(m, n);
  double obj_new = std::numeric_limits<double>::infinity();
  // cout << "Finished declaring new values. \n";
  
  // initiate a residual and objective function vector
  arma::vec res_v(max_iter);
  arma::vec obj_v(max_iter);
  // int n_iter = 0; // initiate the counter
  double delt_obj = std::numeric_limits<double>::infinity(); // set the initial delt_obj to infinity
  // cout << "Finished declaring vectors for storing results. \n";
  
  // start the loop
  int tt = 0; // declare the iteration counter for outputs
  for (int t = 1; t < max_iter; ++t) {
    if (delt_obj <= tol) {
      
      // if delt_Y gets under the tolerance, break out of the loop
      break; 
    } else {
      
      // run the MU updates
      Theta_new = mu_Theta_cpp(Y, s_old, P_old, Theta_0, Delta, Delta_c, Theta_old, alpha1, alpha2);
      P_new = mu_P_cpp(Y, Theta_new, s_old, meds, P_old, beta, ranges);
      Y_new = Theta_new * P_new; // calculating the updated Y
      
      // calculate the updated objective function and its change
      obj_new = obj_fun_cpp(Y, Y_new, Theta_new, P_new, m0, k0, Theta_0, Delta, Delta_c, meds, ranges, alpha1, alpha2, beta);
      obj_v(t) = obj_new; 
      res_v(t) = 1/(2 * (double) m * (double) n) * pow(norm(Y_new - Y, "fro"), 2); // storing all residuals of Y
      delt_obj = abs(obj_new - obj_old)/obj_old; // calculate the change in objective function values
      
      // passing the updates to the next round
      Y_old = Y_new;
      Theta_old = Theta_new;
      P_old = P_new;
      obj_old = obj_new;
    }
    tt += 1; // update the iteration counter for outputs
    // cout << "Finished iteration " << t << ".\n";
  }
  // cout << "Finished all iterations .\n";
  
  // normalize the proportion matrix P so its row sums equal 1
  double col_sum = 0;
  for(int i = 0; i < n; ++i) {
    col_sum = sum(P_new.col(i));
    for(int j = 0; j < k; ++j) {
      P_new(j, i) = P_new(j, i)/col_sum;
    }
  }
  // cout << "Finished normalization. \n";
  
  // record the time and print it
  double time_elapsed = timer.toc();
  // cout << "number of seconds: " << time_elapsed << endl;
  
  // check if the weights are uniform
  arma::vec uni_v(k, fill::ones);
  std::string weights = "range-adaptive";
  double diff = sum(ranges - uni_v);
  if (diff == 0) {
    weights = "uniform";
  }
  // cout << "Finished weight check. \n";
  
  // write the list for returned values;
  Rcpp::List ret = Rcpp::List::create(Named("Y_hat") = Y_new, 
                                      _["Theta_hat"] = Theta_new,
                                      _["s_hat"] = s_old,
                                      _["P_hat"] = P_new,
                                      _["res_v"] = res_v,
                                      _["obj_v"] = obj_v,
                                      _["fixed_s"] = "size-fixed",
                                      _["weights"] = weights,
                                      _["n_iter"] = tt);
  return ret;
}
