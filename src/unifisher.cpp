#include <RcppArmadillo.h>
#include <iostream>     // std::cout
#include <cmath>        // std::abs
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Uses the Fisher scoring method to find the root of a univariate likelihood function, denoted l.
//'
//' @param lp <func> First derivative of l, function of only the parameter.
//' @param fisherInfo <func> Fisher information, function of only the parameter.
//' @param n <int> Sample size.
//' @param a <double> Initial guess for parameter.
//' @param eps <double> Tolerance used to classify convergence.
//' @param maxIter <int> Maximum number of iterations to perform before stopping.
//' @return <list> Estimate for parameter and the number of iterations required for convergence.
//[[Rcpp::export]]
Rcpp::List unifisher(Rcpp::Function lp, Rcpp::Function fisherInfo, int n, Rcpp::NumericVector a, double eps, int maxIter){
  int i;
  double relConvCrit = 100;
  NumericVector x = a;
  NumericVector xOld (1);
  double lpx = 100;
  double fisherInfox = 100;

  // update x based on updating equations
  for(i = 0; i < maxIter; i++){
    xOld[0] = x[0];

    // evaluate derivative of likelihood at old value
    lpx = as<double>(lp(xOld));
    // evaluate fisher info at old value
    fisherInfox = as<double>(fisherInfo(xOld));
    // update rule to get new value
    x[0] = xOld[0] + (lpx/fisherInfox)/(double)n;

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
