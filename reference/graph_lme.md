# Metric graph linear mixed effects models

Fitting linear mixed effects model in metric graphs. The random effects
can be Gaussian Whittle-Matern fields, discrete Gaussian Markov random
fields based on the graph Laplacian, as well as Gaussian random fields
with isotropic covariance functions.

## Usage

``` r
graph_lme(
  formula,
  graph,
  model = list(type = "linearModel"),
  which_repl = NULL,
  optim_method = "L-BFGS-B",
  possible_methods = "L-BFGS-B",
  model_options = list(),
  BC = 1,
  previous_fit = NULL,
  fix_coeff = FALSE,
  parallel = FALSE,
  n_cores = parallel::detectCores() - 1,
  optim_controls = list(),
  improve_hessian = FALSE,
  hessian_args = list(),
  check_euclidean = TRUE
)
```

## Arguments

- formula:

  Formula object describing the relation between the response variables
  and the fixed effects.

- graph:

  A `metric_graph` object.

- model:

  The random effects model that will be used (it also includes the
  option of not having any random effects). It can be either a
  character, whose options are 'lm', for linear models without random
  effects; 'WM1' and 'WM2' for Whittle-Matern models with \\\alpha\\=1
  and 2, with exact precision matrices, respectively; 'WM' for
  Whittle-Matern models where one also estimates the smoothness
  parameter via finite-element method; 'isoExp' for a model with
  isotropic exponential covariance; 'GL1' and 'GL2' for a SPDE model
  based on graph Laplacian with \\\alpha\\ = 1 and 2, respectively.
  'WMD1' is the directed Whittle-Matern with \\\alpha\\=1. There is also
  the option to provide it as a list containing the elements `type`,
  which can be `linearModel`, `WhittleMatern`, `graphLaplacian` or
  `isoCov`. `linearModel` corresponds to a linear model without random
  effects. For `WhittleMatern` models, that is, if the list contains
  `type = 'WhittleMatern'`, one can choose between a finite element
  approximation of the precision matrix by adding `fem = TRUE` to the
  list, or to use the exact precision matrix (by setting `fem = FALSE`).
  If `fem` is `FALSE`, there is also the parameter `alpha`, to determine
  the order of the SPDE, which is either 1 or 2. If `fem` is `FALSE` and
  `alpha` is not specified, then the default value of `alpha=1` will be
  used. If `fem` is `TRUE` and one does not specify `alpha`, it will be
  estimated from the data. However, if one wants to have `alpha` fixed
  to some value, the user can specify either `alpha` or `nu` in the
  list. See the vignettes for examples. Finally, for type
  'WhittleMatern', there is an optional argument, `rspde_order`, that
  chooses the order of the rational approximation. By default
  `rspde_order` is 2. Finally, if one wants to fit a nonstationary
  model, then `fem` necessarily needs to be `TRUE`, and one needs to
  also supply the matrices `B.tau` and `B.kappa` or `B.range` and
  `B.sigma`. For `graph-Laplacian` models, the list must also contain a
  parameter `alpha` (which is 1 by default). For `isoCov` models, the
  list must contain a parameter `cov_function`, containing the
  covariance function. The function accepts a string input for the
  following covariance functions: 'exp_covariance', 'WM1', 'WM2', 'GL1',
  'GL2'. For another covariance function, the function itself must be
  provided as the `cov_function` argument. The default is
  'exp_covariance', the exponential covariance. We also have
  covariance-based versions of the Whittle-Matern and graph Laplacian
  models, however they are much slower, they are the following (string)
  values for 'cov_function': 'alpha1' and 'alpha2' for Whittle-Matern
  fields, and 'GL1' and 'GL2' for graph Laplacian models. Finally, for
  `Whittle-Matern` models, there is an additional parameter `version`,
  which can be either 1 or 2, to tell which version of the likelihood
  should be used. Version is 1 by default.

- which_repl:

  Vector or list containing which replicates to consider in the model.
  If `NULL` all replicates will be considered.

- optim_method:

  The method to be used with `optim` function.

- possible_methods:

  Which methods to try in case the optimization fails or the hessian is
  not positive definite. The options are 'Nelder-Mead', 'L-BFGS-B',
  'BFGS', 'CG' and 'SANN'. By default only 'L-BFGS-B' is considered.

- model_options:

  A list containing additional options to be used in the model.
  Currently, it is possible to fix parameters during the estimation or
  change the starting values of the parameters. The general structure of
  the elements of the list is `fix_parname` and `start_parname`, where
  `parname` stands for the name of the parameter. If `fix_parname` is
  not `NULL`, then the model with be fitted with the `parname` being
  fixed at the value that was passed. If `start_parname` is not `NULL`,
  the model will be fitted using the value passed as starting value for
  `parname`. the For 'WM' models, the possible elements of the list are:
  `fix_sigma_e`, `start_sigma_e`, `fix_nu`, `start_nu`, `fix_sigma`,
  `start_sigma`, `fix_range`, `start_range`. Alternatively, one can use
  `fix_sigma_e`, `start_sigma_e`, `fix_nu`, `start_nu`, `fix_tau`,
  `start_tau`, `fix_kappa`, `start_kappa`. For 'WM1', 'WM2', 'isoExp',
  'GL1' and 'GL2' models, the possible elements of the list are
  `fix_sigma_e`, `start_sigma_e`, `fix_sigma`, `start_sigma`,
  `fix_range`, `start_range`. Alternatively, one can use `fix_sigma_e`,
  `start_sigma_e`, `fix_tau`, `start_tau`, `fix_kappa`, `start_kappa`.
  For 'isoCov' models, the possible values are `fix_sigma_e`,
  `start_sigma_e`, `fix_par_vec`, `start_par_vec`. Observe that contrary
  to the other models, for 'isoCov' models, both `fix_par_vec` and
  `start_par_vec` should be given as vectors of the size of the
  dimension of the vector for the input of the covariance function
  passed to the 'isoCov' model. Furthermore, for 'isoCov' models,
  `fix_par_vec` is a logical vector, indicating which parameters to be
  fixed, and the values will be kept fixed to the values given to
  `start_par_vec`, one can also use `fix_sigma_e` and `start_sigma_e`
  for controlling the std. deviation of the measurement error.

- BC:

  For `WhittleMatern` models, decides which boundary condition to use
  (0,1). Here, 0 is Neumann boundary conditions and 1 specifies
  stationary boundary conditions.

- previous_fit:

  An object of class `graph_lme`. Use the fitted coefficients as
  starting values.

- fix_coeff:

  If using a previous fit, should all coefficients be fixed at the
  starting values?

- parallel:

  logical. Indicating whether to use `optimParallel()` or not.

- n_cores:

  Number of cores to be used if parallel is true.

- optim_controls:

  Additional controls to be passed to
  [`optim()`](https://rdrr.io/r/stats/optim.html) or `optimParallel()`.

- improve_hessian:

  Should a more precise estimate of the hessian be obtained? Turning on
  might increase the overall time.

- hessian_args:

  List of controls to be used if `improve_hessian` is `TRUE`. The list
  can contain the arguments to be passed to the `method.args` argument
  in the `hessian` function. See the help of the `hessian` function in
  'numDeriv' package for details. Observet that it only accepts the
  "Richardson" method for now, the method "complex" is not supported.

- check_euclidean:

  Check if the graph used to compute the resistance distance has
  Euclidean edges? The graph used to compute the resistance distance has
  the observation locations as vertices.

## Value

A list containing the fitted model.
