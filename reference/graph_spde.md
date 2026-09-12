# 'INLA' implementation of Whittle-Matérn fields for metric graphs

This function creates an 'INLA' object that can be used in 'INLA' or
'inlabru' to fit Whittle-Matérn fields on metric graphs.

## Usage

``` r
graph_spde(
  graph_object,
  alpha = 1,
  LGCP = FALSE,
  directional = FALSE,
  stationary_endpoints = "all",
  parameterization = c("matern", "spde"),
  start_range = NULL,
  prior_range = NULL,
  range_lower_bound = NULL,
  range_upper_bound = NULL,
  range_prec_inc = 1,
  start_kappa = NULL,
  prior_kappa = NULL,
  kappa_lower_bound = NULL,
  kappa_upper_bound = NULL,
  kappa_prec_inc = 1,
  start_sigma = NULL,
  prior_sigma = NULL,
  sigma_lower_bound = NULL,
  sigma_upper_bound = NULL,
  sigma_prec_inc = 1,
  start_tau = NULL,
  prior_tau = NULL,
  tau_lower_bound = NULL,
  tau_upper_bound = NULL,
  tau_prec_inc = 1,
  factor_start_range = 0.3,
  type_start_range_bbox = "diag",
  shared_lib = "detect",
  debug = FALSE,
  verbose = 0
)
```

## Arguments

- graph_object:

  A `metric_graph` object.

- alpha:

  The order of the SPDE.

- LGCP:

  Logical. If TRUE, the model will be used as a latent model for a
  Log-Gaussian Cox Process. Default is FALSE.

- directional:

  Should a directional model be used? Currently only implemented for
  `alpha=1`.

- stationary_endpoints:

  Which vertices of degree 1 should contain stationary boundary
  conditions? Set to "all" for all vertices of degree 1, "none" for none
  of the vertices of degree 1, or pass the indices of the vertices of
  degree 1 for which stationary conditions are desired.

- parameterization:

  Which parameterization to be used? The options are 'matern' (sigma and
  range) and 'spde' (sigma and kappa).

- start_range:

  Starting value for range parameter.

- prior_range:

  a `list` containing the elements `meanlog` and `sdlog`, that is, the
  mean and standard deviation of the range parameter on the log scale.
  Will not be used if prior.kappa is non-null. If bounds are specified,
  it can also contain `mean` and `prec` for the beta prior.

- range_lower_bound:

  Lower bound for the range parameter. Default is NULL (no lower bound,
  or 0).

- range_upper_bound:

  Upper bound for the range parameter. Default is NULL (no upper bound).
  If specified, a beta prior will be used.

- range_prec_inc:

  Amount to increase the precision in the beta prior distribution for
  the range parameter. Default is 1. Similar to nu.prec.inc in
  rspde.matern().

- start_kappa:

  Starting value for kappa.

- prior_kappa:

  a `list` containing the elements `meanlog` and `sdlog`, that is, the
  mean and standard deviation of kappa on the log scale. If bounds are
  specified, it can also contain `mean` and `prec` for the beta prior.

- kappa_lower_bound:

  Lower bound for the kappa parameter. Default is NULL (no lower bound,
  or 0).

- kappa_upper_bound:

  Upper bound for the kappa parameter. Default is NULL (no upper bound).
  If specified, a beta prior will be used.

- kappa_prec_inc:

  Amount to increase the precision in the beta prior distribution for
  the kappa parameter. Default is 1.

- start_sigma:

  Starting value for sigma.

- prior_sigma:

  a `list` containing the elements `meanlog` and `sdlog`, that is, the
  mean and standard deviation of sigma on the log scale. If bounds are
  specified, it can also contain `mean` and `prec` for the beta prior.

- sigma_lower_bound:

  Lower bound for the sigma parameter. Default is NULL (no lower bound,
  or 0).

- sigma_upper_bound:

  Upper bound for the sigma parameter. Default is NULL (no upper bound).
  If specified, a beta prior will be used.

- sigma_prec_inc:

  Amount to increase the precision in the beta prior distribution for
  the sigma parameter. Default is 1.

- start_tau:

  Starting value for tau.

- prior_tau:

  a `list` containing the elements `meanlog` and `sdlog`, that is, the
  mean and standard deviation of tau on the log scale. If bounds are
  specified, it can also contain `mean` and `prec` for the beta prior.

- tau_lower_bound:

  Lower bound for the tau parameter. Default is NULL (no lower bound, or
  0).

- tau_upper_bound:

  Upper bound for the tau parameter. Default is NULL (no upper bound).
  If specified, a beta prior will be used.

- tau_prec_inc:

  Amount to increase the precision in the beta prior distribution for
  the tau parameter. Default is 1.

- factor_start_range:

  Factor to multiply the max/min dimension of the bounding box to obtain
  a starting value for range. Default is 0.3.

- type_start_range_bbox:

  Which dimension from the bounding box should be used? The options are
  'diag', the default, 'max' and 'min'.

- shared_lib:

  Which shared lib to use for the cgeneric implementation? If "detect",
  the default, it will use the compiled rSPDE library if one is
  available locally, and otherwise the models built into the 'INLA'
  binary. If 'INLA', it will use the models built into the 'INLA'
  binary. If 'rSPDE', it will use the local installation of the rSPDE
  package, which requires rSPDE to have been installed from source with
  its cgeneric library compiled. Otherwise, you can directly supply the
  path of the .so (or .dll) file.

- debug:

  Should debug be displayed?

- verbose:

  Level of verbosity. 0 is silent, 1 prints basic information, 2 prints
  more.

## Value

An 'INLA' object.

## Details

This function is used to construct a Matern SPDE model on a metric
graph. The latent field \\u\\ is the solution of the SPDE
\$\$(\kappa^2 - \Delta)^\alpha u = \sigma W,\$\$ where \\W\\ is Gaussian
white noise on the metric graph. This model implements exactly the cases
in which \\\alpha = 1\\ or \\\alpha = 2\\. For a finite element
approximation for general \\\alpha\\ we refer the reader to the 'rSPDE'
package and to the Whittle–Matérn fields with general smoothness
vignette.

We also have the alternative parameterization \\\rho =
\frac{\sqrt{8(\alpha-0.5)}}{\kappa}\\, which can be interpreted as a
range parameter.

Let \\\kappa_0\\ and \\\sigma_0\\ be the starting values for \\\kappa\\
and \\\sigma\\, we write \\\sigma = \exp\\\theta_1\\\\ and \\\kappa =
\exp\\\theta_2\\\\. We assume priors on \\\theta_1\\ and \\\theta_2\\ to
be normally distributed with mean, respectively, \\\log(\sigma_0)\\ and
\\\log(\kappa_0)\\, and variance 10. Similarly, if we let \\\rho_0\\ be
the starting value for \\\rho\\, then we write \\\rho =
\exp\\\theta_2\\\\ and assume a normal prior for \\\theta_2\\, with mean
\\\log(\rho_0)\\ and variance 10.
