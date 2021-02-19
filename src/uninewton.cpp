#include <RcppArmadillo.h>
#include <iostream>     // std::cout
#include <cmath>        // std::abs
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Uses the Newton method to find the root of a univariate function, denoted g.
//'
//' @param gp <func> First derivative of g.
//' @param gpp <func> Second derivative of g.
//' @param a <double> Initial guess for the root.
//' @param eps <double> Tolerance used to classify convergence.
//' @param maxIter <int> Maximum number of iterations to perform before stopping.
//' @return <list> Root of g and the number of iterations required for convergence.
//[[Rcpp::export]]
Rcpp::List uninewton(Rcpp::Function gp, Rcpp::Function gpp, Rcpp::NumericVector a, double eps, int maxIter){
  int i;
  double relConvCrit = 100;
  NumericVector x = a;
  NumericVector xOld (1);
  double gpx = 100;
  double gppx = 100;

  // update x based on updating equations
  for(i = 0; i < maxIter; i++){
    xOld[0] = x[0];

    // evaluate the derivative and double derivative at old value
    gpx = as<double>(gp(xOld));
    gppx = as<double>(gpp(xOld));
    // update rule to get new value
    x[0] = xOld[0] - gpx/gppx;

    // check for convergence and stop if convergence has been reached
    relConvCrit = std::abs(x[0] - xOld[0]) / std::abs(xOld[0]);
    if(relConvCrit < eps){
      break;
    }
  }
  // throw error if no convergence
  if(relConvCrit > eps){
    throw Rcpp::exception("convergence failed! Increase maxIter or eps.", false);
  }

  return List::create(
    Rcpp::Named("estimate") = x[0],
                               Rcpp::Named("iter") = i);
}
