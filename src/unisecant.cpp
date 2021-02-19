#include <RcppArmadillo.h>
#include <iostream>     // std::cout
#include <cmath>        // std::abs
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Uses the secant method to find the root of a univariate function, denoted g.
//'
//' @param gp <func> First derivative of g.
//' @param a <double> Initial guess for the left endpoint.
//' @param b <double> Initial guess for the right endpoint.
//' @param eps <double> Tolerance used to classify convergence.
//' @param maxIter <int> Maximum number of iterations to perform before stopping.
//' @return <list> Root of g and the number of iterations required for convergence.
//[[Rcpp::export]]
Rcpp::List unisecant(Rcpp::Function gp, Rcpp::NumericVector a, Rcpp::NumericVector b, double eps, int maxIter){
  int i;
  double relConvCrit = 100;
  NumericVector xNew (1);
  NumericVector xOld (1);
  double x = 1;
  double gpxNew = 1;
  double gpxOld = 1;

  // set midpoint based on initial values
  xOld[0] = as<double>(a);
  xNew[0] = as<double>(b);

  // update midpoint and interval based on updating equations
  for(i = 0; i < maxIter; i++){
    gpxOld = as<double>(gp(xOld));
    gpxNew = as<double>(gp(xNew));
    x = xNew[0];
    // update rule
    xNew[0] = xNew[0] - gpxNew*((xNew[0] - xOld[0])/(gpxNew - gpxOld));
    xOld[0] = x;
    // check for convergence and stop if convergence has been reached
    relConvCrit = std::abs(xNew[0] - xOld[0]) / std::abs(xOld[0]);
    if(relConvCrit < eps){
      break;
    }
  }
  // throw error if no convergence
  if(relConvCrit > eps){
    throw Rcpp::exception("convergence failed! Increase maxIter or eps.", false);
  }

  return List::create(
    Rcpp::Named("estimate") = xNew[0],
    Rcpp::Named("iter") = i);
}
