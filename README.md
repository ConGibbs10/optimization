
<!-- README.md is generated from README.Rmd. Please edit that file -->

# optimization

<!-- badges: start -->
<!-- badges: end -->

The goal of optimization is to provided functions for finding roots.

## Installation

Install this package from GitHub. You will need to install `devtools` if
you do not already have it.

``` r
# install.packages('devtools')
devtools::install_github('ConGibbs10/optimization')
#> Skipping install of 'optimization' from a github remote, the SHA1 (2495cd32) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

## Example

This is a basic example which shows you how to find roots with
`optimization`. Suppose
*X*<sub>1</sub>, …, *X*<sub>25</sub> ∼ Cauchy(*θ*, 1). Then, one can
derive the log-likelihood and its first and second derivatives.
Furthermore, one can derive the Fisher information. These functions,
univariate in theta, are coded in `R` below.

``` r
library(optimization)

# log-likelihood
ll <- function(theta){
  obs <- c(-8.86, -6.82, -4.03, -2.84, 0.14, 0.19, 0.24, 0.27, 0.49, 0.62, 0.76, 1.09, 
         1.18, 1.32, 1.36, 1.58, 1.58, 1.78, 2.13, 2.15, 2.36, 4.05, 4.11, 4.12, 6.83)
  n <- length(obs)
  ll <- -n*log(pi) - sum(log(1 + (obs - theta)^2))
  return(ll)
}

# first derivative of log-likelihood
d1ll <- function(theta){
  obs <- c(-8.86, -6.82, -4.03, -2.84, 0.14, 0.19, 0.24, 0.27, 0.49, 0.62, 0.76, 1.09, 
         1.18, 1.32, 1.36, 1.58, 1.58, 1.78, 2.13, 2.15, 2.36, 4.05, 4.11, 4.12, 6.83)
  return(sum((2*(obs - theta))/(1 + (obs - theta)^2)))
}

# second derivative of log-likelihood
d2ll <- function(theta){
  obs <- c(-8.86, -6.82, -4.03, -2.84, 0.14, 0.19, 0.24, 0.27, 0.49, 0.62, 0.76, 1.09, 
         1.18, 1.32, 1.36, 1.58, 1.58, 1.78, 2.13, 2.15, 2.36, 4.05, 4.11, 4.12, 6.83)
  return(sum((-2*(1 - (obs - theta)^2)) / ((1 + (obs - theta)^2)^2)))
}

# fisher information
fi <- function(theta){1/2}
```

Let’s take a look at the derivative of the log-likelihood.

``` r
vd1ll <- Vectorize(d1ll)
root <- uniroot(vd1ll, lower = -5, upper = 5)$root
ggplot(data = data.frame(x = 0), aes(x = x)) +
  stat_function(fun = vd1ll) + 
  geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') +
  geom_vline(xintercept = root, color = 'red', linetype = 'dashed') +
  geom_text(data = data.frame(x = root, y = -10.75),
            aes(x = x, y = y, label = format(root, digits = 4)),
            inherit.aes = FALSE, angle = 90, vjust = -.25) +
  scale_x_continuous(bquote(theta), limits = c(-5, 5), breaks = seq(-5, 5, by = 1)) +
  scale_y_continuous(expression(paste("l'(", theta, ")"))) +
  theme_bw()
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

I’ll now demonstrate how to find the root with each of the included
numerical methods.

### Bisection

``` r
mresults <- list()
mresults$bisection <- 
  unibisection(gp = d1ll,          # derivative of log-likelihood
               a = 1,              # left endpoint
               b = 2,              # right endpoint
               eps = 1e-6,         # tolerance 
               maxIter = 10000)    # maximum number of iterations to perform
```

### Secant

``` r
mresults$secant <- 
  unisecant(gp = d1ll,        # derivative of log-likelihood
            a = 0,            # left endpoint
            b = 2,            # right endpoint
            eps = 1e-6,       # tolerance 
            maxIter = 10000)  # maximum number of iterations to perform
```

### Newton-Raphson

``` r
mresults$newton <- 
  uninewton(gp = d1ll,        # derivative of log-likelihood
            gpp = d2ll,       # second derivative of log-likelihood
            a = 1,            # initial guess
            eps = 1e-6,       # tolerance       
            maxIter = 10000)  # maximum number of iterations to perform
```

### Fisher Scoring

``` r
mresults$fisher <- 
  unifisher(lp = d1ll,        # derivative of log-likelihood
            fisherInfo = fi,  # fisher information
            n = 25,           # sample size
            a = 1,            # initial guess
            eps = 1e-6,       # tolerance       
            maxIter = 10000)  # maximum number of iterations to perform)
```

## Results

After approximating the root using each method, the results can be
summarized in the following table:

    #> 
    #> -----------------------------------
    #>   Method     Estimate   Iterations 
    #> ----------- ---------- ------------
    #>  Bisection    1.188         18     
    #> 
    #>   Secant      1.188         3      
    #> 
    #>   Newton      1.188         2      
    #> 
    #>   Fisher      1.188         3      
    #> -----------------------------------

Every method approximates the root to the same precision as `uniroot()`.
