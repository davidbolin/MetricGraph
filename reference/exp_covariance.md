# Exponential covariance function

Evaluates the exponential covariance function \$\$C(h) = \sigma^2
\exp\\-kappa h\\\$\$

## Usage

``` r
exp_covariance(h, theta)
```

## Arguments

- h:

  Distances to evaluate the covariance function at.

- theta:

  A vector `c(sigma, kappa)`, where `sigma` is the standard deviation
  and `kappa` is a range-like parameter.

## Value

A vector with the values of the covariance function.
