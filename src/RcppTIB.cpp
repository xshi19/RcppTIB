// Rcpp functions of TIB

#include "RcppArmadillo.h"

// Rank 1 update for Cholesky
// en.wikipedia.org/wiki/Cholesky_decomposition#Rank-one_update
// [[Rcpp::export]]
arma::mat cholupdate(const arma::mat & L, const arma::vec & x) {
  arma::mat L_ = L;
  arma::vec x_ = x;
  int n = L_.n_rows;
  double r, c, s;
  for (int i = 0; i < n; i++) {
    r = sqrt(L_(i, i) * L_(i, i) + x_(i) * x_(i));
    c = r / L_(i, i);
    s = x_(i) / L_(i, i);
    if ((i + 1) < n) {
      for (int j = i; j < n; j++) {
        L_(j, i) = (L_(j, i) + s * x_(j)) / c;
        x_(j) = c * x_(j) - s * L_(j, i);
      }
    }
  }
  return L_;
}

// Forward subsititution of complex matrix
// [[Rcpp::export]]
arma::cx_mat forwardsolve(const arma::cx_mat & L, const arma::cx_mat & x) {
  arma::cx_mat L_ = L;
  arma::cx_mat x_ = x;
  arma::cx_mat y = arma::solve(arma::trimatl(L_), x_);
  return y;  
}

// Solve bi-diagonal matrix M of TIB
// The order of the vector M is the same as the Matlab sparse M
// [[Rcpp::export]]
arma::cx_mat solvetibM(const arma::cx_vec & M, const arma::cx_mat & x) {
  arma::cx_mat y = x;
  int k;
  y.row(0) = x.row(0) / M(0);
  for (int j = 1; j < x.n_rows; j++) {
    k = 2 * j;
    y.row(j) = (x.row(j) - M(k - 1) * y.row(j - 1)) / M(k);
  }
  return y;
}

// Multiply bi-diagonal matrix N of TIB
// The order of the vector N is the same as the Matlab sparse N
// [[Rcpp::export]]
arma::cx_mat timestibN(const arma::cx_vec & N, const arma::cx_mat & x) {
  arma :: cx_mat y = x;
  int k;
  y.row(0) = N(0) * x.row(0);
  for (int j = 1; j < x.n_rows; j++) {
    k = 2 * j;
    y.row(j) = N(k - 1) * x.row(j - 1) + N(k) * x.row(j);
  }
  return y;
}

// sptib in C++ 
// [[Rcpp::export]]
Rcpp::List sptib_cpp(const arma::cx_vec & lambda) {
  int n = lambda.n_elem;
  arma :: vec A(2 * n - 1);
  arma :: vec B(2 * n - 1);
  arma :: cx_vec M(A, B);
  arma :: cx_vec N(A, B);
  int k;

  M(0) = 1.0 + 0.0 * 1i;  
  N(0) = lambda(0);

  std :: complex<double> c1 = sqrt(1.0 - pow(abs(lambda(0)), 2));
  std :: complex<double> c2 = c1;

  for (int j = 1; j < n; j++) {
    k = 2 * j;
    c2 = sqrt(1.0 - pow(abs(lambda(j)), 2));
    M(k) = -M(k - 2);
    N(k) = M(k) * lambda(j);
    N(k - 1) = M(k - 2) * (c2 / c1);
    M(k - 1) = N(k - 1) * conj(lambda(j - 1));
    c1 = c2;
  }
  
  M = M / sqrt(1.0 - pow(abs(lambda(0)), 2));
  N = N / sqrt(1.0 - pow(abs(lambda(0)), 2));
  return Rcpp::List::create (Rcpp::Named("M") = M, 
                             Rcpp::Named("N") = N);
}


  
  
