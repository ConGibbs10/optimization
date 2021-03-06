% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{unibisection}
\alias{unibisection}
\title{Uses the bisection method to find the root of a univariate function, denoted g.}
\usage{
unibisection(gp, a, b, eps, maxIter)
}
\arguments{
\item{gp}{<func> First derivative of g.}

\item{a}{<double> Initial guess for the left endpoint.}

\item{b}{<double> Initial guess for the right endpoint.}

\item{eps}{<double> Tolerance used to classify convergence.}

\item{maxIter}{<int> Maximum number of iterations to perform before stopping.}
}
\value{
<list> Root of g and the number of iterations required for convergence.
}
\description{
Uses the bisection method to find the root of a univariate function, denoted g.
}
