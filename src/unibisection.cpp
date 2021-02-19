#include <RcppArmadillo.h>
#include <iostream>     // std::cout
#include <cmath>        // std::abs
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Uses the bisection method to find the root of a univariate function, denoted g.
//'
//' @param gp <func> First derivative of g.
//' @param a <double> Initial guess for the left endpoint.
//' @param b <double> Initial guess for the right endpoint.
//' @param eps <double> Tolerance used to classify convergence.
//' @param maxIter <int> Maximum number of iterations to perform before stopping.
//' @return <list> Root of g and the number of iterations required for convergence.
//[[Rcpp::export]]
Rcpp::List unibisection(Rcpp::Function gp, Rcpp::NumericVector a, Rcpp::NumericVector b, double eps, int maxIter){
  int i;
  double relConvCrit = 100;
  NumericVector xmid (1);
  NumericVector xmidOld (1);

  // check if bisection method is appropriate.
  if(as<double>(gp(a)) * as<double>(gp(b)) > 0){
    throw Rcpp::exception("cannot apply IVT. Pick new a and b.", false);
  }

  // set midpoint based on initial values
  xmid[0] = (as<double>(a) + as<double>(b))/2;

  // update midpoint and interval based on updating equations
  for(i = 0; i < maxIter; i++){
    // check the concavity and reassign a or b depending
    if(as<double>(gp(a)) * as<double>(gp(xmid)) <= 0){
      b[0] = xmid[0];
    } else{
      a[0] = xmid[0];
    }
    // set old midpoint
    xmidOld[0] = xmid[0];
    // recalculate midpoint with new a and b
    xmid[0] = (as<double>(a) + as<double>(b))/2;
    // check for convergence and stop if convergence has been reached
    relConvCrit = std::abs(xmid[0] - xmidOld[0]) / std::abs(xmidOld[0]);
    if(relConvCrit < eps){
      break;
    }
  }
  // throw error if no convergence
  if(relConvCrit > eps){
    throw Rcpp::exception("convergence failed! Increase maxIter or eps.", false);
  }

  return List::create(
    Rcpp::Named("estimate") = xmid[0],
    Rcpp::Named("iter") = i);
}
