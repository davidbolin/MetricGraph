# Glance at a `graph_lme` object

Glance accepts a `graph_lme` object and returns a
[`tidyr::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with exactly one row of model summaries. The summaries are the square
root of the estimated variance of the measurement error, residual
degrees of freedom, AIC, BIC, log-likelihood, the type of latent model
used in the fit and the total number of observations.

## Usage

``` r
# S3 method for class 'graph_lme'
glance(x, ...)
```

## Arguments

- x:

  A `graph_lme` object.

- ...:

  Additional arguments. Currently not used.

## Value

A
[`tidyr::tibble()`](https://tibble.tidyverse.org/reference/tibble.html)
with exactly one row and columns:

- `nobs` Number of observations used.

- `sigma` the square root of the estimated residual variance

- `logLik` The log-likelihood of the model.

- `AIC` Akaike's Information Criterion for the model.

- `BIC` Bayesian Information Criterion for the model.

- `deviance` Deviance of the model.

- `df.residual` Residual degrees of freedom.

- `model.type` Type of latent model fitted.

## See also

[augment.graph_lme](https://davidbolin.github.io/MetricGraph/reference/augment.graph_lme.md)
