# Whittle--Matérn fields with general smoothness

## Introduction

In this vignette we will introduce how to fit Whittle–Matérn fields with
general smoothness based on finite element and rational approximations.
The theory for this approach is provided in [Bolin, Kovács, et al.
(2023)](https://doi.org/10.1090/mcom/3929) and [Bolin, Simas, et al.
(2023)](https://doi.org/10.1080/10618600.2023.2231051). For the
implementation, we make use of the [`rSPDE`
package](https://davidbolin.github.io/rSPDE/) for the rational
approximations.

These models are thus implemented using finite element approximations.
Such approximations are not needed for integer smoothness parameters,
and for the details about the exact models we refer to the vignettes

- [Whittle–Matérn fields on metric
  graphs](https://davidbolin.github.io/MetricGraph/articles/random_fields.md)
- [INLA and inlabru
  interfaces](https://davidbolin.github.io/MetricGraph/articles/inla_interface.md)

For details on the construction of metric graphs, see [Working with
metric
graphs](https://davidbolin.github.io/MetricGraph/articles/metric_graphs.md)

For further details on data manipulation on metric graphs, see [Data
manipulation on metric
graphs](https://davidbolin.github.io/MetricGraph/articles/metric_graph_data.md)

## Constructing the graph and the mesh

We begin by loading the `rSPDE` and `MetricGraph` packages:

``` r

  library(rSPDE)
  library(MetricGraph)
```

As an example, we consider the following metric graph

``` r

  edge1 <- rbind(c(0,0),c(1,0))
  edge2 <- rbind(c(0,0),c(0,1))
  edge3 <- rbind(c(0,1),c(-1,1))
  theta <- seq(from=pi,to=3*pi/2,length.out = 20)
  edge4 <- cbind(sin(theta),1+ cos(theta))
  edges = list(edge1, edge2, edge3, edge4)
  graph <- metric_graph$new(edges = edges)
  graph$plot()
```

![](fem_models_files/figure-html/unnamed-chunk-2-1.png)

To construct a FEM approximation of a Whittle–Matérn field with general
smoothness, we must first construct a mesh on the graph.

``` r

  graph$build_mesh(h = 0.1)
  graph$plot(mesh=TRUE)
```

![](fem_models_files/figure-html/unnamed-chunk-3-1.png)

In the command `build_mesh`, the argument `h` decides the largest
spacing between nodes in the mesh. As can be seen in the plot, the mesh
is very coarse, so let’s reduce the value of `h` and rebuild the mesh:

``` r

graph$build_mesh(h = 0.01)
```

We are now ready to specify the model
``` math
(\kappa^2 - \Delta)^{\alpha/2} \tau u = \mathcal{W}
```
for the Whittle–Matérn field $`u`$. For this, we use the
`matern.operators` function from the `rSPDE` package:

``` r

  sigma <- 1.3
  range <- 0.15
  nu <- 0.8 

  rspde.order <- 2
  op <- matern.operators(nu = nu, range = range, sigma = sigma, 
                         parameterization = "matern",
                         m = rspde.order, graph = graph)                     
```

As can be seen in the code, we specify $`\kappa`$ via the practical
correlation range $`\sqrt{8\nu}/\kappa`$. Also, the model is not
parametrized by $`\tau, \alpha`$ but instead by $`\sigma, \nu`$. Here,
`sigma` denotes the standard deviation of the field and `nu` is the
smoothness parameter, which is related to $`\alpha`$ via the relation
$`\alpha = \nu + 1/2`$. The object `op` contains the matrices needed for
evaluating the distribution of the stochastic weights in the FEM
approximation.

Let us simulate the field $`u`$ at the mesh locations and plot the
result:

``` r

u <- simulate(op)
df_u <- data.frame(u = as.vector(u), edge_number = graph$mesh$VtE[,1],
                   distance_on_edge = graph$mesh$VtE[,2])
df_u <- graph$process_data(data = df_u, normalized = TRUE)
graph$plot_function(data = "u", newdata = df_u, type = "plotly")
```

If we want to evaluate $`u(s)`$ at some locations $`s_1,\ldots, s_n`$,
we need to multiply the weights with the FEM basis functions
$`\varphi_i(s)`$ evaluated at the locations. For this, we can construct
the observation matrix $`\boldsymbol{\mathrm{A}}`$, with elements
$`A_{ij} = \varphi_j(s_i)`$, which links the FEM basis functions to the
locations. This can be done by the function `fem_basis` in the metric
graph object. To illustrate this, let us simulate some observation
locations on the graph and construct the matrix:

``` r

obs.per.edge <- 100
obs.loc <- NULL
for(i in 1:graph$nE) {
  obs.loc <- rbind(obs.loc,
                   cbind(rep(i,obs.per.edge), runif(obs.per.edge)))
}
n.obs <- obs.per.edge*graph$nE
A <- graph$fem_basis(obs.loc)
```

In the code, we generate $`100`$ observation locations per edge in the
graph, drawn at random. It can be noted that we assume that the
observation locations are given in the format $`(e, d)`$ where $`e`$
denotes the edge of the observation and $`d`$ is the position on the
edge, i.e., the relative distance from the first vertex of the edge.

To compute the precision matrix from the covariance-based rational
approximation one can use the
[`precision()`](https://davidbolin.github.io/rSPDE/reference/precision.CBrSPDEobj.html)
method on object returned by the
[`matern.operators()`](https://davidbolin.github.io/rSPDE/reference/matern.operators.html)
function:

``` r

  Q <- precision(op)
```

As an illustration of the model, let us compute the covariance function
between the process at $`s=(2,0.1)`$, that is, the point at edge 2 and
distance on edge 0.1, and all the other mesh points. To this end, we can
use the helper function `cov_function_mesh` that is contained in the
`op` object:

``` r

  c_cov <- cov_function_mesh(op, matrix(c(2,0.1),1,2))
  df_c_cov <- data.frame(cov = as.vector(c_cov), edge_number = graph$mesh$VtE[,1],
                         distance_on_edge = graph$mesh$VtE[,2])
  df_c_cov <- graph$process_data(data = df_c_cov, normalized = TRUE)
  graph$plot_function(data = "cov", newdata = df_c_cov, type = "plotly")
```

## Using the model for inference

There is built-in support for computing log-likelihood functions and
performing kriging prediction in the `rSPDE` package which we can use
for the graph model. To illustrate this, we use the simulation to create
some noisy observations of the process. We generate the observations as
$`Y_i = 1 + 2x_{i1} - 3 x_{i2} + u(s_i) + \varepsilon_i`$, where
$`\varepsilon_i \sim N(0,\sigma_e^2)`$ is Gaussian measurement noise,
$`x_1`$ and $`x_2`$ are covariates generated from the relative positions
of the observations on the graph.

``` r

    sigma.e <- 0.1

    x1 <- obs.loc[,1]
    x2 <- obs.loc[,2]

    Y <- 1 + 2*x1 - 3*x2 + as.vector(A %*% u + sigma.e * rnorm(n.obs))
```

Let us now fit the model. To this end we will use the
[`graph_lme()`](https://davidbolin.github.io/MetricGraph/reference/graph_lme.md)
function (that, for the finite element models, acts as a wrapper for the
[`rspde_lme()`](https://davidbolin.github.io/rSPDE/reference/rspde_lme.html)
function from the `rSPDE` package). To this end, let us now assemble the
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) with the
observations, the observation locations and the covariates:

``` r

df_data <- data.frame(y = Y, edge_number = obs.loc[,1],
                        distance_on_edge = obs.loc[,2],
                        x1 = x1, x2 = x2)
```

Let us now add the data to the graph object and plot it:

``` r

graph$add_observations(data = df_data, normalized = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are normalized by the length of the edge.

``` r

graph$plot(data = "y")
```

![](fem_models_files/figure-html/unnamed-chunk-12-1.png)

We can now fit the model. To this end, we use the
[`graph_lme()`](https://davidbolin.github.io/MetricGraph/reference/graph_lme.md)
function and set the model to `'WM`’.

``` r

fit <- graph_lme(y ~ x1 + x2, graph = graph, model = "WM")
```

Let us obtain a summary of the model:

``` r

summary(fit)
```

    ## 
    ## Latent model - Whittle-Matern
    ## 
    ## Call:
    ## graph_lme(formula = y ~ x1 + x2, graph = graph, model = "WM")
    ## 
    ## Fixed effects:
    ##             Estimate Std.error z-value Pr(>|z|)    
    ## (Intercept)   0.2646    0.7101   0.373   0.7094    
    ## x1            2.1377    0.1989  10.747   <2e-16 ***
    ## x2           -2.1069    0.7074  -2.979   0.0029 ** 
    ## 
    ## Random effects:
    ##        Estimate Std.error z-value
    ## alpha  1.265733  0.019863  63.722
    ## tau    0.056240  0.005443  10.333
    ## kappa 15.712349  2.553094   6.154
    ## 
    ## Random effects (Matern parameterization):
    ##       Estimate Std.error z-value
    ## nu     0.76573   0.01986  38.550
    ## sigma  1.32354   0.13564   9.758
    ## range  0.15752   0.02441   6.453
    ## 
    ## Measurement error:
    ##          Estimate Std.error z-value
    ## std. dev 0.097018  0.006282   15.44
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
    ## 
    ## Log-Likelihood:  -126.561 
    ## Number of function calls by 'optim' = 65
    ## Optimization method used in 'optim' = L-BFGS-B
    ## 
    ## Time used to:     fit the model =  23.40653 secs

An improved estimate of the Hessian can be obtained by setting
`improve_hessian` to `TRUE`, which improves the precision of the
standard errors.

``` r

fit <- graph_lme(y ~ x1 + x2, graph = graph, model = "WM", improve_hessian = TRUE)
```

Let us obtain a summary of the model fitted with the improved Hessian:

``` r

summary(fit)
```

    ## 
    ## Latent model - Whittle-Matern
    ## 
    ## Call:
    ## graph_lme(formula = y ~ x1 + x2, graph = graph, model = "WM", 
    ##     improve_hessian = TRUE)
    ## 
    ## Fixed effects:
    ##             Estimate Std.error z-value Pr(>|z|)    
    ## (Intercept)   0.2646    0.7112   0.372   0.7098    
    ## x1            2.1377    0.1995  10.715   <2e-16 ***
    ## x2           -2.1069    0.7074  -2.978   0.0029 ** 
    ## 
    ## Random effects:
    ##       Estimate Std.error z-value
    ## alpha  1.26573   0.05265  24.041
    ## tau    0.05624   0.01402   4.012
    ## kappa 15.71235   3.03851   5.171
    ## 
    ## Random effects (Matern parameterization):
    ##       Estimate Std.error z-value
    ## nu     0.76573   0.05265  14.544
    ## sigma  1.32354   0.13564   9.758
    ## range  0.15752   0.02441   6.453
    ## 
    ## Measurement error:
    ##          Estimate Std.error z-value
    ## std. dev 0.097018  0.006288   15.43
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
    ## 
    ## Log-Likelihood:  -126.561 
    ## Number of function calls by 'optim' = 65
    ## Optimization method used in 'optim' = L-BFGS-B
    ## 
    ## Time used to:     fit the model =  21.04701 secs 
    ##   compute the Hessian = 3.23019 secs

We can also obtain additional information by using the function
[`glance()`](https://davidbolin.github.io/MetricGraph/reference/glance.graph_lme.md):

``` r

glance(fit)
```

    ## # A tibble: 1 × 9
    ##    nobs  sigma logLik   AIC   BIC deviance df.residual model         alpha
    ##   <int>  <dbl>  <dbl> <dbl> <dbl>    <dbl>       <dbl> <chr>         <dbl>
    ## 1   400 0.0970  -127.  267.  295.     253.         393 WhittleMatern  1.27

Let us compare the values of the parameters of the latent model with the
true ones:

``` r

print(data.frame(sigma = c(sigma, fit$alt_par_coeff$coeff["sigma"]), 
                   range = c(range, fit$alt_par_coeff$coeff["range"]),
                   nu = c(nu, fit$alt_par_coeff$coeff["nu"]),
                   row.names = c("Truth", "Estimates")))
```

    ##              sigma     range        nu
    ## Truth     1.300000 0.1500000 0.8000000
    ## Estimates 1.323538 0.1575225 0.7657331

### Kriging

Given that we have estimated the parameters, let us compute the kriging
predictor of the field given the observations at the mesh nodes.

We will perform kriging with the
[`predict()`](https://rdrr.io/r/stats/predict.html) method. To this end,
we need to provide a `data.frame` containing the prediction locations,
as well as the values of the covariates at the prediction locations.

``` r

  df_pred <- data.frame(edge_number = graph$mesh$VtE[,1],
                        distance_on_edge = graph$mesh$VtE[,2],
                        x1 = graph$mesh$VtE[,1],
                        x2 = graph$mesh$VtE[,2])

  u.krig <- predict(fit, newdata = df_pred, normalized = TRUE)
```

The estimate is shown in the following figure

``` r

  df_krig <- data.frame(mean = as.vector(u.krig$mean), edge_number = graph$mesh$VtE[,1],
                        distance_on_edge = graph$mesh$VtE[,2])
  df_krig <- graph$process_data(data = df_krig, normalized = TRUE)
  graph$plot_function(data = "mean", newdata = df_krig, type = "plotly")
```

We can also use the
[`augment()`](https://davidbolin.github.io/MetricGraph/reference/augment.graph_lme.md)
function to easily plot the predictions. Let us a build a 3d plot now
and add the observed values on top of the predictions:

``` r

df_pred <- graph$process_data(data = df_pred, normalized = TRUE)
p <- augment(fit, newdata = df_pred, normalized = TRUE) %>%
          graph$plot_function(data = ".fitted", type = "plotly")

graph$plot(data = "y", p = p, type = "plotly")
```

## Further details on `graph_lme`

The `graph_lme` function provides flexibility in model fitting by
allowing users to fix certain parameters at specific values or set
custom starting values for the optimization process. This can be useful
when you have prior knowledge about some parameters or when you want to
improve convergence by providing better starting points.

### Fixing parameters

Parameters can be fixed by using the `model_options` argument with
elements of the form `fix_parname = value`, where `parname` is the name
of the parameter you want to fix. The parameters that can be fixed in
the `model_options` list are:

- `fix_sigma_e`: Fix the standard deviation of the noise parameter
  $`\sigma_\varepsilon`$
- `fix_sigma`: Fix the standard deviation parameter $`\sigma`$
- `fix_range`: Fix the range parameter
- `fix_alpha`: Fix the fractional power $`\alpha`$

### Setting starting values

Similarly, you can set starting values for parameters using elements of
the form `start_parname = value` in the `model_options` list. This is
particularly useful when the default starting values might be far from
the optimal values, which could lead to slow convergence or convergence
to a local minimum.

- `start_sigma_e`: Starting value for the standard deviation of the
  noise parameter $`\sigma_\varepsilon`$
- `start_sigma`: Starting value for the standard deviation parameter
  $`\sigma`$
- `start_range`: Starting value for the range parameter
- `start_alpha`: Starting value for the fractional power $`\alpha`$

### Example: Fixing and setting starting values

Let us demonstrate how to use these parameters with a graph model
example. We will fix the standard deviation $`\sigma`$ to 1 and set the
starting value for the range parameter to 0.5, on the previous example.

``` r

# Fit model with fixed sigma and start_range
fit_fixed <- graph_lme(y ~ x1 + x2, graph = graph, model = "WM",
                      model_options = list(
                        fix_sigma = 1,      # Fix sigma to 1
                        start_range = 0.5   # Set starting value for range
                      ))

# Print summary of the model
summary(fit_fixed)
```

    ## 
    ## Latent model - Whittle-Matern
    ## 
    ## Call:
    ## graph_lme(formula = y ~ x1 + x2, graph = graph, model = "WM", 
    ##     model_options = list(fix_sigma = 1, start_range = 0.5))
    ## 
    ## Fixed effects:
    ##             Estimate Std.error z-value Pr(>|z|)    
    ## (Intercept)   0.2461    0.4741   0.519    0.604    
    ## x1            2.1407    0.1340  15.973  < 2e-16 ***
    ## x2           -2.1519    0.5063  -4.250 2.14e-05 ***
    ## 
    ## Random effects:
    ##               Estimate Std.error z-value
    ## nu             0.72823   0.01055   69.00
    ## sigma (fixed)  1.00000        NA      NA
    ## range          0.11440   0.00872   13.12
    ## 
    ## Random effects (SPDE parameterization):
    ##       Estimate Std.error z-value
    ## alpha  1.22823   0.01055   116.4
    ## tau    0.06768        NA      NA
    ## kappa 21.09912        NA      NA
    ## 
    ## Measurement error:
    ##          Estimate Std.error z-value
    ## std. dev 0.097303  0.006312   15.42
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
    ## 
    ## Log-Likelihood:  -132.2996 
    ## Number of function calls by 'optim' = 38
    ## Optimization method used in 'optim' = L-BFGS-B
    ## 
    ## Time used to:     fit the model =  15.97415 secs

``` r

# Compare log-likelihoods
cat("Log-likelihood with fixed sigma:", logLik(fit_fixed), "\n")
```

    ## Log-likelihood with fixed sigma: -132.2996

``` r

cat("Log-likelihood with all parameters estimated:", logLik(fit), "\n")
```

    ## Log-likelihood with all parameters estimated: -126.561

## Fitting a model with replicates

Let us now illustrate how to simulate a data set with replicates and
then fit a model to such data. To simulate a latent model with
replicates, all we do is set the `nsim` argument to the number of
replicates.

``` r

  n.rep <- 30
  u.rep <- simulate(op, nsim = n.rep)
```

Now, let us generate the observed values $`Y`$:

``` r

  sigma.e <- 0.3
  Y.rep <- A %*% u.rep + sigma.e * matrix(rnorm(n.obs * n.rep), ncol = n.rep)
```

Note that $`Y`$ is a matrix with 20 columns, each column containing one
replicate. We need to turn `y` into a vector and create an auxiliary
vector `repl` indexing the replicates of `y`:

``` r

y_vec <- as.vector(Y.rep)
repl <- rep(1:n.rep, each = n.obs)                       

df_data_repl <- data.frame(y = y_vec,
                              edge_number = rep(obs.loc[,1], n.rep),
                              distance_on_edge = rep(obs.loc[,2], n.rep), 
                              repl = repl)
```

Let us clear the previous observations and add the new data to the
graph:

``` r

graph$add_observations(data = df_data_repl, normalized = TRUE, 
                            group = "repl", clear_obs = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are normalized by the length of the edge.

We can now fit the model in the same way as before by using the
[`rspde_lme()`](https://davidbolin.github.io/rSPDE/reference/rspde_lme.html)
function. Note that we can optimize in parallel by setting `parallel` to
`TRUE`. If we do not specify which replicate to consider, in the
`which_repl` argument, all replicates will be considered.

``` r

fit_repl <- graph_lme(y ~ -1, graph = graph, model = "WM", parallel = TRUE)
```

Now, let us see a summary of the fit:

``` r

summary(fit_repl)
```

    ## 
    ## Latent model - Whittle-Matern
    ## 
    ## Call:
    ## graph_lme(formula = y ~ -1, graph = graph, model = "WM", parallel = TRUE)

    ## 
    ## No fixed effects.

    ## 
    ## Random effects:
    ##        Estimate Std.error z-value
    ## alpha  1.289682  0.006362  202.72
    ## tau    0.052035  0.001474   35.30
    ## kappa 15.901647  0.507588   31.33
    ## 
    ## Random effects (Matern parameterization):
    ##       Estimate Std.error z-value
    ## nu    0.789682  0.006362  124.13
    ## sigma 1.313555  0.024716   53.15
    ## range 0.158063  0.004862   32.51
    ## 
    ## Measurement error:
    ##          Estimate Std.error z-value
    ## std. dev   0.3014    0.0030   100.5
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
    ## 
    ## Log-Likelihood:  -9837.691 
    ## Number of function calls by 'optim' = 117
    ## Optimization method used in 'optim' = L-BFGS-B
    ## 
    ## Time used to:     fit the model =  1.15531 mins 
    ##   set up the parallelization = 2.71371 secs

Let us also take a glance of the fit:

``` r

glance(fit_repl)
```

    ## # A tibble: 1 × 9
    ##    nobs sigma logLik    AIC    BIC deviance df.residual model         alpha
    ##   <int> <dbl>  <dbl>  <dbl>  <dbl>    <dbl>       <dbl> <chr>         <dbl>
    ## 1 12000 0.301 -9838. 19683. 19713.   19675.       11996 WhittleMatern  1.29

Let us compare the values of the parameters of the latent model with the
true ones:

``` r

print(data.frame(sigma = c(sigma, fit_repl$alt_par_coeff$coeff["sigma"]), 
                   range = c(range, fit_repl$alt_par_coeff$coeff["range"]),
                   nu = c(nu, fit_repl$alt_par_coeff$coeff["nu"]),
                   row.names = c("Truth", "Estimates")))
```

    ##              sigma     range        nu
    ## Truth     1.300000 0.1500000 0.8000000
    ## Estimates 1.313555 0.1580625 0.7896819

Let us do kriging. We will use the same prediction locations as in the
previous example. Let us get prediction for replicate 10, then add the
original observations on top of them. Observe that even though the
group/replicate might be a number, it is stored as character in the data
frame, so we need to pass it as a character in the `which_repl`
argument.

``` r

p <- augment(fit_repl, which_repl = 10, newdata = df_pred, normalized = TRUE) %>% 
          graph$plot_function(data = ".fitted", type = "plotly")

graph$plot(data = "y", group = 10, type = "plotly", p = p)
```

## Using the R-INLA implementation

We also have an `R-INLA` implementation of the rational SPDE approach
for metric graphs.

We begin by defining the model by using the
[`rspde.metric_graph()`](https://davidbolin.github.io/rSPDE/reference/rspde.metric_graph.html)
function. This function contains the same arguments as the function
[`rspde.matern()`](https://davidbolin.github.io/rSPDE/reference/rspde.matern.html).
We refer the reader to the [R-INLA implementation of the rational SPDE
approach](https://davidbolin.github.io/rSPDE/articles/rspde_inla.html)
vignette for further details.

We begin by clearing the previous observations and adding the
observations (for the case without replicates) to the graph:

``` r

graph$clear_observations()
graph$add_observations(data = df_data, normalized = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are normalized by the length of the edge.

Let us create the model object:

``` r

  library(INLA)
  rspde_model <- rspde.metric_graph(graph)
```

By default, the order of the rational approximation is 2.

We can now create the auxiliary quantities that will be needed with the
[`graph_data_rspde()`](https://davidbolin.github.io/rSPDE/reference/graph_data_rspde.html)
function:

``` r

  data_rspde <- graph_data_rspde(rspde_model, name = "field")
```

The remaining is standard: we create the formula object, the stack
object, and then fit the model by using the
[`inla()`](https://rdrr.io/pkg/INLA/man/inla.html) function. So, first
we create the formula object:

``` r

  f.s <- y ~ -1 + Intercept + x1 + x2 + f(field, model = rspde_model)
```

Now we create the `inla.stack` object. To such an end, observe that
`data_rspde` contains the dataset as the `data` component, the index as
the `index` component and the so-called `A` matrix as the `basis`
component. We will now create the stack using these components:

``` r

  stk.dat <- inla.stack(
    data = data_rspde[["data"]]["y"], A = list(data_rspde[["basis"]],1), tag = "est",
    effects =
      list(c(
        data_rspde[["index"]],
        list(Intercept = 1)), list(x1 = data_rspde[["data"]]["x1"] ,
                                      x2 = data_rspde[["data"]]["x2"])
      )
    )
```

Finally, we can fit the model:

``` r

  rspde_fit <- inla(f.s, data = inla.stack.data(stk.dat),
    control.inla = list(int.strategy = "eb"),
    control.predictor = list(A = inla.stack.A(stk.dat), compute = TRUE),
    num.threads = "1:1"
  )
```

We can use the same functions as the `rspde` fitted models in `inla`.
For instance, we can see the results in the original scale by creating
the `result` object:

``` r

  result_fit <- rspde.result(rspde_fit, "field", rspde_model)
  summary(result_fit)
```

    ##             mean        sd 0.025quant 0.5quant 0.975quant     mode
    ## std.dev 1.397280 0.1600590   1.114000 1.386190   1.741860 1.361870
    ## range   0.186825 0.0514245   0.108250 0.179105   0.308829 0.164048
    ## nu      0.690081 0.0952478   0.509418 0.688366   0.882027 0.686564

Let us compare with the true values:

``` r

  result_df <- data.frame(
    parameter = c("std.dev", "range", "nu"),
    true = c(sigma, range, nu),
    mean = c(
      result_fit$summary.std.dev$mean,
      result_fit$summary.range$mean,
      result_fit$summary.nu$mean
    ),
    mode = c(
      result_fit$summary.std.dev$mode,
      result_fit$summary.range$mode,
      result_fit$summary.nu$mode
    )
  )
  print(result_df)
```

    ##   parameter true      mean      mode
    ## 1   std.dev 1.30 1.3972820 1.3618734
    ## 2     range 0.15 0.1868252 0.1640483
    ## 3        nu 0.80 0.6900807 0.6865640

We can also plot the posterior marginal densities with the help of the
[`gg_df()`](https://davidbolin.github.io/MetricGraph/reference/gg_df.metric_graph_spde_result.md)
function:

``` r

  posterior_df_fit <- gg_df(result_fit)

  library(ggplot2)

  ggplot(posterior_df_fit) + geom_line(aes(x = x, y = y)) + 
  facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](fem_models_files/figure-html/unnamed-chunk-40-1.png)

### Kriging with the `R-INLA` implementation

We will do kriging on the mesh locations:

``` r

  pred_loc <- graph$mesh$VtE
```

Let us now add the observations for prediction:

``` r

graph$add_observations(data = data.frame(y=rep(NA,nrow(pred_loc)), 
                                x1 = graph$mesh$VtE[,1],
                                x2 = graph$mesh$VtE[,2],
                                edge_number = pred_loc[,1], 
                                distance_on_edge = pred_loc[,2]), 
                                normalized = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are normalized by the length of the edge.

Let us now create a new model and, then, compute the auxiliary
components at the prediction locations. To this end, we set the argument
`only_pred` to `TRUE`, in which it will return the `data.frame`
containing the `NA` data.

``` r

  rspde_model_prd <- rspde.metric_graph(graph) 
  data_rspde_prd <- graph_data_rspde(rspde_model_prd, only_pred = TRUE)
```

Let us build the prediction stack using the components of
`data_rspde_prd` and gather it with the estimation stack.

``` r

  ef.prd <- 
    list(c(data_rspde_prd[["index"]], list(Intercept = 1)), 
          list(x1 = data_rspde_prd[["data"]][["x1"]],
                x2 = data_rspde_prd[["data"]][["x2"]]))
  stk.prd <- inla.stack(
    data = data.frame(y = data_rspde_prd[["data"]][["y"]]),
    A = list(data_rspde_prd[["basis"]],1), tag = "prd",
    effects = ef.prd
  )
  stk.all <- inla.stack(stk.dat, stk.prd)
```

Let us obtain the predictions:

``` r

rspde_fitprd <- inla(f.s,
  data = inla.stack.data(stk.all),
  control.predictor = list(
    A = inla.stack.A(stk.all),
    compute = TRUE, link = 1
  ),
  control.compute = list(
    return.marginals = FALSE,
    return.marginals.predictor = FALSE
  ),
  control.inla = list(int.strategy = "eb"),
  num.threads = "1:1"
)
```

Let us now extract the indices of the predicted nodes and store the
means:

``` r

id.prd <- inla.stack.index(stk.all, "prd")$data
m.prd <- rspde_fitprd$summary.fitted.values$mean[id.prd]
```

Finally, let us plot the predicted values. To this end we will use the
`plot_function()` graph method.

``` r

  df_prd <- data.frame(mean = m.prd, edge_number = pred_loc[,1],
                       distance_on_edge = pred_loc[,2])
  df_prd <- graph$process_data(data = df_prd, normalized = TRUE)
  graph$plot_function(data = "mean", newdata = df_prd, type = "plotly")
```

## Using `R-INLA` implementation to fit models with replicates

Let us begin by cloning the graph and clearing the observations on the
cloned graph:

``` r

graph_rep <- graph$clone()
graph_rep$clear_observations()
```

We will now add the data with replicates to the graph:

``` r

graph_rep$add_observations(data = data.frame(y=as.vector(Y.rep), 
                          edge_number = rep(obs.loc[,1], n.rep), 
                          distance_on_edge = rep(obs.loc[,2], n.rep),
                          repl = rep(1:n.rep, each = n.obs)), 
                          group = "repl",
                          normalized = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are normalized by the length of the edge.

Let us create a new `rspde` model object:

``` r

rspde_model_rep <- rspde.metric_graph(graph_rep)
```

To fit the model with replicates we need to create the auxiliary
quantities with the
[`graph_data_rspde()`](https://davidbolin.github.io/rSPDE/reference/graph_data_rspde.html)
function, where we set the `repl` argument in the function
`graph_data_spde` to `.all` since we want to use all replicates:

``` r

data_rspde_rep <- graph_data_rspde(rspde_model_rep, 
                      name = "field", repl = ".all",
                      repl_col = "repl")
```

Let us now create the corresponding `inla.stack` object:

``` r

st.dat.rep <- inla.stack(
  data = data_rspde_rep[["data"]],
  A = data_rspde_rep[["basis"]],
  effects = data_rspde_rep[["index"]]
)
```

Observe that we need the response variable `y` to be a vector. We can
now create the `formula` object, remembering that since we gave the name
argument `field`, when creating the index, we need to pass `field.repl`
to the `formula`:

``` r

f.rep <-
  y ~ -1 + f(field,
    model = rspde_model_rep,
    replicate = field.repl
  )
```

We can, finally, fit the model:

``` r

rspde_fit_rep <-
  inla(f.rep,
    data = inla.stack.data(st.dat.rep),
    family = "gaussian",
    control.predictor =
      list(A = inla.stack.A(st.dat.rep)),
    num.threads = "1:1"
  )
```

We can obtain the estimates in the original scale with the
[`rspde.result()`](https://davidbolin.github.io/rSPDE/reference/rspde.result.html)
function:

``` r

  result_fit_rep <- rspde.result(rspde_fit_rep, "field", rspde_model_rep)
  summary(result_fit_rep)
```

    ##             mean         sd 0.025quant 0.5quant 0.975quant     mode
    ## std.dev 1.325270 0.02459720   1.278100 1.324840   1.374690 1.323680
    ## range   0.169338 0.00778464   0.154011 0.169403   0.184548 0.169864
    ## nu      0.701944 0.02699590   0.652507 0.700527   0.758250 0.695990

Let us compare with the true values of the parameters:

``` r

  result_rep_df <- data.frame(
    parameter = c("std.dev", "range", "nu"),
    true = c(sigma, range, nu),
    mean = c(
      result_fit_rep$summary.std.dev$mean,
      result_fit_rep$summary.range$mean,
      result_fit_rep$summary.nu$mean
    ),
    mode = c(
      result_fit_rep$summary.std.dev$mode,
      result_fit_rep$summary.range$mode,
      result_fit_rep$summary.nu$mode
    )
  )
  print(result_rep_df)
```

    ##   parameter true      mean      mode
    ## 1   std.dev 1.30 1.3252690 1.3236826
    ## 2     range 0.15 0.1693380 0.1698636
    ## 3        nu 0.80 0.7019436 0.6959898

We can also plot the posterior marginal densities with the help of the
[`gg_df()`](https://davidbolin.github.io/MetricGraph/reference/gg_df.metric_graph_spde_result.md)
function:

``` r

  posterior_df_fit_rep <- gg_df(result_fit_rep)

  ggplot(posterior_df_fit_rep) + geom_line(aes(x = x, y = y)) + 
  facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](fem_models_files/figure-html/unnamed-chunk-57-1.png)

## Using `inlabru` implementation

The `inlabru` package allows us to fit models and do kriging in a
straighforward manner, without having to handle `A` matrices, indices
nor `inla.stack` objects. Therefore, we suggest the reader to use this
implementation when using our implementation to fit real data.

Let us clear the graph, since it contains `NA` observations we used for
prediction, add the observations again, and create a new `rSPDE` model
object:

``` r

graph$clear_observations()
graph$add_observations(data = df_data, 
                          normalized = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are normalized by the length of the edge.

``` r

rspde_model <- rspde.metric_graph(graph)
```

Let us now load the `inlabru` package and create the component (which is
`inlabru`’s formula-like object). Let us begin by building the auxiliary
data to be used with the
[`graph_data_rspde()`](https://davidbolin.github.io/rSPDE/reference/graph_data_rspde.html)
function, where we pass the name of the location variable in the above
formula as the `loc_name` argument, which in this case is `"loc"`:

``` r

data_rspde_bru <- graph_data_rspde(rspde_model, bru = TRUE)
```

Now, we create the component to be used in `inlabru`, in which we pass
the `index` element from the `data_rspde_bru` object as index locations:

``` r

    library(inlabru)
    cmp <-
    y ~ -1 + Intercept(1) + x1 + x2 + field(
                          cbind(.edge_number, .distance_on_edge), 
                          model = rspde_model
                          )                   
```

Now, we can directly fit the model, by using the `data` element of
`data_rspde_bru`:

``` r

  rspde_bru_fit <-
    bru(cmp,
        data=data_rspde_bru[["data"]],
        options = list(num.threads = "1:1")
    )
```

Let us now obtain the estimates of the parameters in the original scale
by using the
[`rspde.result()`](https://davidbolin.github.io/rSPDE/reference/rspde.result.html)
function:

``` r

  result_bru_fit <- rspde.result(rspde_bru_fit, "field", rspde_model)
  summary(result_bru_fit)
```

    ##             mean        sd 0.025quant 0.5quant 0.975quant     mode
    ## std.dev 1.398990 0.1706930   1.102360 1.384940   1.771370 1.352720
    ## range   0.187693 0.0523565   0.109026 0.179330   0.313103 0.163428
    ## nu      0.690243 0.0920886   0.514697 0.688932   0.874979 0.688221

Let us compare with the true values of the parameters:

``` r

  result_bru_df <- data.frame(
    parameter = c("std.dev", "range", "nu"),
    true = c(sigma, range, nu),
    mean = c(
      result_bru_fit$summary.std.dev$mean,
      result_bru_fit$summary.range$mean,
      result_bru_fit$summary.nu$mean
    ),
    mode = c(
      result_bru_fit$summary.std.dev$mode,
      result_bru_fit$summary.range$mode,
      result_bru_fit$summary.nu$mode
    )
  )
  print(result_bru_df)
```

    ##   parameter true      mean      mode
    ## 1   std.dev 1.30 1.3989924 1.3527157
    ## 2     range 0.15 0.1876933 0.1634278
    ## 3        nu 0.80 0.6902427 0.6882206

We can also plot the posterior marginal densities with the help of the
[`gg_df()`](https://davidbolin.github.io/MetricGraph/reference/gg_df.metric_graph_spde_result.md)
function:

``` r

  posterior_df_bru_fit <- gg_df(result_bru_fit)

  ggplot(posterior_df_bru_fit) + geom_line(aes(x = x, y = y)) + 
  facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](fem_models_files/figure-html/unnamed-chunk-64-1.png)

### Kriging with the `inlabru` implementation

It is very easy to do kriging with the `inlabru` implementation. We
simply need to provide the prediction locations to the
[`predict()`](https://rdrr.io/r/stats/predict.html) method.

In this example we will use the mesh locations. To this end we will use
the `get_mesh_locations()` method. We also set `bru=TRUE` to obtain a
data frame suitable to be used with `inlabru`. In this case, the mesh
locations will be returned as a `data.frame` with the location columns
`.edge_number` and `.distance_on_edge`. We will, then, add the
covariates `x1` and `x2` to the data frame:

``` r

  prd_loc <- graph$get_mesh_locations(bru = TRUE)
  prd_loc[["x1"]] <- prd_loc[,1]
  prd_loc[["x2"]] <- prd_loc[,2]  
```

Now, we can simply provide these locations to the `predict` method along
with the fitted object `rspde_bru_fit`:

``` r

  y_pred <- predict(rspde_bru_fit, newdata=prd_loc, 
                        ~Intercept + x1 + x2 + field)
```

Let us now prepare the predictions so we can plot them easily by using
the
[`process_rspde_predictions()`](https://davidbolin.github.io/MetricGraph/reference/process_rspde_predictions.md)
function:

``` r

y_pred <- process_rspde_predictions(y_pred, graph = graph, PtE = prd_loc)
```

Finally, let us plot the predicted values. To this end we will use the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method on
`y_pred`:

``` r

  plot(y_pred) 
```

![](fem_models_files/figure-html/unnamed-chunk-68-1.png)

We can also create the 3d plot, together with the true data:

``` r

p <- graph$plot(data = "y", type = "plotly")
plot(y_pred, type = "plotly", p = p)
```

## Using inlabru to fit models with replicates

We can also use our `inlabru` implementation to fit models with
replicates. We will consider the same data that was generated above,
where the number of replicates is 30.

For this implementation we will use the `rspde_model_rep` object.

We can now create the component, passing the vector with the indices of
the replicates as the `replicate` argument. To obtain the auxiliary
data, we will pass `repl` argument we use the function
[`graph_data_rspde()`](https://davidbolin.github.io/rSPDE/reference/graph_data_rspde.html),
where we set it to `.all`, since we want all replicates. Further, we
also set the argument `bru` to `TRUE`.

``` r

data_rspde_rep <- graph_data_rspde(rspde_model_rep, repl = ".all", 
                                    bru = TRUE, repl_col = "repl")
```

We can now define the `bru` component formula, passing the `repl` as the
`replicate` argument:

``` r

  cmp_rep <-
    y ~ -1 + field(cbind(.edge_number, .distance_on_edge), 
                              model = rspde_model_rep,
                              replicate = repl)
```

Now, we are ready to fit the model:

``` r

  rspde_bru_fit_rep <-
    bru(cmp_rep,
        data=data_rspde_rep[["data"]],
        options=list(
        family = "gaussian",
        num.threads = "1:1")
    )
```

We can obtain the estimates in the original scale with the
[`rspde.result()`](https://davidbolin.github.io/rSPDE/reference/rspde.result.html)
function:

``` r

  result_bru_fit_rep <- rspde.result(rspde_bru_fit_rep, "field", rspde_model_rep)
  summary(result_bru_fit_rep)
```

    ##             mean         sd 0.025quant 0.5quant 0.975quant     mode
    ## std.dev 1.325270 0.02459720   1.278100 1.324840   1.374690 1.323680
    ## range   0.169338 0.00778464   0.154011 0.169403   0.184548 0.169864
    ## nu      0.701944 0.02699590   0.652507 0.700527   0.758250 0.695990

Let us compare with the true values of the parameters:

``` r

  result_bru_rep_df <- data.frame(
    parameter = c("std.dev", "range", "nu"),
    true = c(sigma, range, nu),
    mean = c(
      result_bru_fit_rep$summary.std.dev$mean,
      result_bru_fit_rep$summary.range$mean,
      result_bru_fit_rep$summary.nu$mean
    ),
    mode = c(
      result_bru_fit_rep$summary.std.dev$mode,
      result_bru_fit_rep$summary.range$mode,
      result_bru_fit_rep$summary.nu$mode
    )
  )
  print(result_bru_rep_df)
```

    ##   parameter true      mean      mode
    ## 1   std.dev 1.30 1.3252690 1.3236826
    ## 2     range 0.15 0.1693380 0.1698636
    ## 3        nu 0.80 0.7019436 0.6959898

We can also plot the posterior marginal densities with the help of the
[`gg_df()`](https://davidbolin.github.io/MetricGraph/reference/gg_df.metric_graph_spde_result.md)
function:

``` r

  posterior_df_bru_fit_rep <- gg_df(result_bru_fit_rep)

  ggplot(posterior_df_bru_fit_rep) + geom_line(aes(x = x, y = y)) + 
  facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](fem_models_files/figure-html/unnamed-chunk-75-1.png)

Let us now do prediction for observations of replicate `10`. We start by
building the data list with the prediction locations:

``` r

  data_prd_repl <- graph$get_mesh_locations(bru = TRUE)
  data_prd_repl[["repl"]] <- rep(10, nrow(data_prd_repl))
```

Let us now obtain predictions for this replicate:

``` r

  y_pred <- predict(rspde_bru_fit_rep, 
                      newdata=data_prd_repl, 
                      ~field_eval(cbind(.edge_number, .distance_on_edge), 
                                    replicate = repl))
```

Let us now process the predictions:

``` r

  y_pred <- process_rspde_predictions(y_pred, graph = graph, PtE = data_prd_repl)
```

We can now plot the predictions along with the observed values for
replicate `10`:

``` r

p <- plot(y_pred, type = "plotly")
graph_rep$plot(data = "y", group = 10, type = "plotly", p = p)
```

## An example with a non-stationary model

Our goal now is to show how one can fit model with non-stationary
$`\sigma`$ (std. deviation) and non-stationary $`\rho`$ (a range
parameter). One can also use the parameterization in terms of
non-stationary SPDE parameters $`\kappa`$ and $`\tau`$.

We follow the same structure as `INLA`. However, `INLA` only allows one
to specify `B.tau` and `B.kappa` matrices, and, in `INLA`, if one wants
to parameterize in terms of range and standard deviation one needs to do
it manually. Here we provide the option to directly provide the matrices
`B.sigma` and `B.range`.

The usage of the matrices `B.tau` and `B.kappa` are identical to the
corresponding ones in
[`inla.spde2.matern()`](https://rdrr.io/pkg/INLA/man/inla.spde2.matern.html)
function. The matrices `B.sigma` and `B.range` work in the same way, but
they parameterize the stardard deviation and range, respectively.

The columns of the `B` matrices correspond to the same parameter. The
first column does not have any parameter to be estimated, it is a
constant column.

So, for instance, if one wants to share a parameter with both `sigma`
and `range` (or with both `tau` and `kappa`), one simply let the
corresponding column to be nonzero on both `B.sigma` and `B.range` (or
on `B.tau` and `B.kappa`).

### Creating the graph and adding data

For this example we will consider the `pems` data contained in the
`MetricGraph` package. The data consists of traffic speed observations
on highways in the city of San Jose, California. The variable `y`
contains the traffic speeds.

``` r

 pems_graph <- metric_graph$new(edges = pems$edges)
 pems_graph$add_observations(data = pems$data)
 pems_graph$prune_vertices( )
 pems_graph$build_mesh(h=0.1)
```

The summary of this graph:

``` r

summary(pems_graph)
```

    ## A metric graph object with:
    ## 
    ## Vertices:
    ##   Total: 347 
    ##   Degree 1: 11;  Degree 2: 16;  Degree 3: 315;  Degree 4: 5; 
    ##   With incompatible directions:  17 
    ## 
    ## Edges: 
    ##   Total: 504 
    ##   Lengths: 
    ##       Min: 0.01040037  ; Max: 7.663457  ; Total: 470.6126 
    ##   Weights: 
    ##       Columns: .weights 
    ##   That are circles:  0 
    ## 
    ## Graph units: 
    ##   Vertices unit:  degree  ; Lengths unit:  km 
    ## 
    ## Longitude and Latitude coordinates:  TRUE
    ##   Which spatial package:  sp 
    ##   CRS:  +proj=longlat +datum=WGS84 +no_defs
    ## 
    ## Some characteristics of the graph:
    ##   Connected: TRUE
    ##   Has loops: FALSE
    ##   Has multiple edges: TRUE
    ##   Is a tree: FALSE
    ##   Distance consistent: FALSE
    ##   Has Euclidean edges: FALSE
    ## 
    ## Computed quantities inside the graph: 
    ##   Laplacian:  FALSE  ; Geodesic distances:  TRUE 
    ##   Resistance distances:  FALSE  ; Finite element matrices:  FALSE 
    ## 
    ## Mesh: 
    ##   Max h_e:  0.0999534  ; Min n_e:  0 
    ## 
    ## Data: 
    ##   Columns:  y 
    ##   Groups:  .group 
    ## 
    ## Tolerances: 
    ##   vertex-vertex:  0.001 
    ##   vertex-edge:  0.001 
    ##   edge-edge:  0

Observe that it is a non-Euclidean graph.

We now define a non-stationary covariate based on the position along
each edge. This covariate is designed to indicate whether the traffic
speed observation was taken close to an intersection. More precisely, we
identify points whose relative position on the edge is either below 0.1
or above 0.9, corresponding to locations near the endpoints of the edge.
Thus, the covariate is equal to TRUE for points close to intersections
and FALSE otherwise.

``` r

cov_pos <- (pems_graph$mesh$VtE[,2] > 0.9) | (pems_graph$mesh$VtE[,2] < 0.1)
```

We will now build the non-stationary matrices to be used:

``` r

 B.sigma = cbind(0, 1, 0, cov_pos, 0)
 B.range = cbind(0, 0, 1,  0, cov_pos)
```

Let us also obtain the same covariate for the observations:

``` r

cov_obs <- pems$data[[".distance_on_edge"]]
cov_obs <- (cov_obs > 0.9) | (cov_obs < 0.1)
```

Let add this covariate to the data:

``` r

pems_graph$add_observations(data = pems_graph$mutate(cov_obs = cov_obs),
                            clear_obs = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are NOT normalized by the length of the edge.

    ## The unit for edge lengths is km

    ## The current tolerance for removing distant observations is (in km): 3.83172861400444

### Fitting the model with `graph_lme`

We are now in position to fit this model using the
[`graph_lme()`](https://davidbolin.github.io/MetricGraph/reference/graph_lme.md)
function. We will also add `cov_obs` as a covariate for the model.

``` r

fit <- graph_lme(y ~ 1, graph = pems_graph, model = list(type = "WhittleMatern", 
                    B.sigma = B.sigma, B.range = B.range, fem = TRUE))
```

Let us now obtain a summary of the fitted model:

``` r

summary(fit)
```

    ## 
    ## Latent model - Generalized Whittle-Matern
    ## 
    ## Call:
    ## graph_lme(formula = y ~ 1, graph = pems_graph, model = list(type = "WhittleMatern", 
    ##     B.sigma = B.sigma, B.range = B.range, fem = TRUE))
    ## 
    ## Fixed effects:
    ##             Estimate Std.error z-value Pr(>|z|)    
    ## (Intercept)   51.526     2.524   20.41   <2e-16 ***
    ## 
    ## Random effects:
    ##        Estimate Std.error z-value
    ## nu      2.27550   0.08526  26.688
    ## theta1  2.87528   0.32294   8.904
    ## theta2  1.85945   0.20965   8.870
    ## theta3  0.39400   1.88828   0.209
    ## theta4  0.12797   0.84706   0.151
    ## 
    ## Measurement error:
    ##          Estimate Std.error z-value
    ## std. dev   6.9554    0.3517   19.78
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
    ## 
    ## Log-Likelihood:  -1206.891 
    ## Number of function calls by 'optim' = 501
    ## Optimization method used in 'optim' = Nelder-Mead
    ## 
    ## Time used to:     fit the model =  55.68221 secs

Let us plot the range parameter along the mesh, so we can see how it is
varying:

``` r

est_range <- exp(B.range[,-1]%*%fit$coeff$random_effects[2:5])
df_range <- data.frame(range = est_range, edge_number = pems_graph$mesh$VtE[,1],
                       distance_on_edge = pems_graph$mesh$VtE[,2])
df_range <- pems_graph$process_data(data = df_range, normalized = TRUE)
pems_graph$plot_function(data = "range", newdata = df_range, vertex_size = 0,
                    type = "mapview", mapview_caption = "Range")
```

Similarly, we have for sigma:

``` r

est_sigma <- exp(B.sigma[,-1]%*%fit$coeff$random_effects[2:5])
df_sigma <- data.frame(sigma = est_sigma, edge_number = pems_graph$mesh$VtE[,1],
                       distance_on_edge = pems_graph$mesh$VtE[,2])
df_sigma <- pems_graph$process_data(data = df_sigma, normalized = TRUE)
pems_graph$plot_function(data = "sigma", newdata = df_sigma, vertex_size = 0,
                    type = "mapview", mapview_caption = "Sigma")
```

Our goal now is to plot the estimated marginal standard deviation of
this model. To this end, we start by creating the non-stationary Matérn
operator using the `rSPDE` package:

``` r

rspde_object_ns <- rSPDE::spde.matern.operators(graph = pems_graph,
                                                parameterization = "matern",
                                                B.sigma = B.sigma,
                                                B.range = B.range,
                                                theta = fit$coeff$random_effects[2:5],
                                                nu = fit$coeff$random_effects[1])
```

Now, we compute the estimated marginal standard deviation:

``` r

est_cov_matrix <- covariance_mesh(rspde_object_ns)
est_std_dev <- sqrt(Matrix::diag(est_cov_matrix))
```

We can now plot:

``` r

df_std <- data.frame(std = est_std_dev, edge_number = pems_graph$mesh$VtE[,1],
                     distance_on_edge = pems_graph$mesh$VtE[,2])
df_std <- pems_graph$process_data(data = df_std, normalized = TRUE)
pems_graph$plot_function(data = "std", newdata = df_std, vertex_size = 0,
          type = "mapview", mapview_caption = "Std. dev")
```

### Fixing parameters in non-stationary models

In non-stationary models, the parameters are labeled as theta1, theta2,
etc., corresponding to the coefficients in the B matrices. Similar to
the stationary case, we can fix individual parameters or set starting
values, but these options must be set in the `model_options` list
argument using `fix_theta1`, `fix_theta2`, etc., or `start_theta1`,
`start_theta2`, etc.

For example, if we want to fix the first coefficient theta1 to 0 in the
previous example:

``` r

# Fit model with fixed theta1 parameter
fit_ns_fixed_theta1 <- graph_lme(y ~ 1, 
                               graph = pems_graph, 
                               model = list(type = "WhittleMatern", 
                                          B.sigma = B.sigma, 
                                          B.range = B.range, 
                                          fem = TRUE),
                               model_options = list(fix_theta1 = 0.5))  # Fix theta1 to 0.5
```

    ## Warning in rSPDE::rspde_lme(formula = formula, loc =
    ## cbind(df_data[[".edge_number"]], : optim method L-BFGS-B failed to provide a
    ## positive-definite Hessian. Another optimization method was used.

``` r

summary(fit_ns_fixed_theta1)
```

    ## 
    ## Latent model - Generalized Whittle-Matern
    ## 
    ## Call:
    ## graph_lme(formula = y ~ 1, graph = pems_graph, model = list(type = "WhittleMatern", 
    ##     B.sigma = B.sigma, B.range = B.range, fem = TRUE), model_options = list(fix_theta1 = 0.5))
    ## 
    ## Fixed effects:
    ##             Estimate Std.error z-value Pr(>|z|)    
    ## (Intercept)  50.8495    0.8783    57.9   <2e-16 ***
    ## 
    ## Random effects:
    ##                Estimate Std.error z-value
    ## nu             0.266500  0.004959  53.743
    ## theta1 (fixed) 0.500000        NA      NA
    ## theta2         2.087741  0.511435   4.082
    ## theta3         2.481165  0.378906   6.548
    ## theta4         0.924256  0.885353   1.044
    ## 
    ## Measurement error:
    ##          Estimate Std.error z-value
    ## std. dev  13.7723    0.5875   23.44
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
    ## 
    ## Log-Likelihood:  -1328.889 
    ## Number of function calls by 'optim' = 443
    ## Optimization method used in 'optim' = Nelder-Mead
    ## 
    ## Time used to:     fit the model =  32.24226 secs

Similarly, we can provide starting values for the entire theta vector
with `start_theta`:

``` r

# Fit model with starting values for theta parameters
fit_ns_start <- graph_lme(y ~ 1, 
                        graph = pems_graph, 
                        model = list(type = "WhittleMatern", 
                                   B.sigma = B.sigma, 
                                   B.range = B.range, 
                                   fem = TRUE),
                        model_options = list(start_theta = c(0.4, 0.7, 1.0, 0.2)))  # Starting values for theta vector
```

    ## Warning in rSPDE::rspde_lme(formula = formula, loc =
    ## cbind(df_data[[".edge_number"]], : optim method L-BFGS-B failed to provide a
    ## positive-definite Hessian. Another optimization method was used.

``` r

summary(fit_ns_start)
```

    ## 
    ## Latent model - Generalized Whittle-Matern
    ## 
    ## Call:
    ## graph_lme(formula = y ~ 1, graph = pems_graph, model = list(type = "WhittleMatern", 
    ##     B.sigma = B.sigma, B.range = B.range, fem = TRUE), model_options = list(start_theta = c(0.4, 
    ##     0.7, 1, 0.2)))
    ## 
    ## Fixed effects:
    ##             Estimate Std.error z-value Pr(>|z|)    
    ## (Intercept)   51.190     2.789   18.36   <2e-16 ***
    ## 
    ## Random effects:
    ##        Estimate Std.error z-value
    ## nu       1.4767    0.6954   2.123
    ## theta1   3.0923    0.4429   6.982
    ## theta2   2.2318    0.5233   4.265
    ## theta3  -0.4515    1.0539  -0.428
    ## theta4  -0.3982    0.8347  -0.477
    ## 
    ## Measurement error:
    ##          Estimate Std.error z-value
    ## std. dev   7.1685    0.3959   18.11
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 
    ## 
    ## Log-Likelihood:  -1207.766 
    ## Number of function calls by 'optim' = 502
    ## Optimization method used in 'optim' = Nelder-Mead
    ## 
    ## Time used to:     fit the model =  47.0822 secs

### Fitting the inlabru rSPDE model

Let us then fit the same model using `inlabru` now. We start by defing
the `rSPDE` model with the
[`rspde.metric_graph()`](https://davidbolin.github.io/rSPDE/reference/rspde.metric_graph.html)
function:

``` r

rspde_model_nonstat <- rspde.metric_graph(pems_graph,
  B.sigma = B.sigma,
  B.range = B.range,
  parameterization = "matern") 
```

Let us now create the
[`data.frame()`](https://rdrr.io/r/base/data.frame.html) and the vector
with the replicates indexes:

``` r

 data_rspde_bru_ns <- graph_data_rspde(rspde_model_nonstat, bru = TRUE)
```

Let us create the component and fit.

``` r

cmp_nonstat <-
  y ~ -1 + Intercept(1) + field(
    cbind(.edge_number, .distance_on_edge),
    model = rspde_model_nonstat
  )


rspde_fit_nonstat <-
  bru(cmp_nonstat,
    data = data_rspde_bru_ns[["data"]],
    family = "gaussian",
    options = list(num.threads = "1:1")
  )
```

We can get the summary:

``` r

summary(rspde_fit_nonstat)
```

    ## inlabru version: 2.15.0 
    ## INLA version: 26.08.22 
    ## Latent components:
    ## Intercept: main = linear(1)
    ## field: main = cgeneric(cbind(.edge_number, .distance_on_edge))
    ## Observation models:
    ##   Model tag: <No tag>
    ##     Family: 'gaussian'
    ##     Data class: 'metric_graph_data', 'data.frame'
    ##     Response class: 'numeric'
    ##     Predictor: y ~ Intercept + field
    ##     Additive/Linear/Rowwise: TRUE/TRUE/TRUE
    ##     Used components: effect[Intercept, field], latent[] 
    ## Time used:
    ##     Pre = 0.161, Running = 60.3, Post = 0.332, Total = 60.8 
    ## Fixed effects:
    ##             mean    sd 0.025quant 0.5quant 0.975quant  mode kld
    ## Intercept 50.606 2.875     44.838   50.621     56.291 50.62   0
    ## 
    ## Random effects:
    ##   Name     Model
    ##     field CGeneric
    ## 
    ## Model hyperparameters:
    ##                                         mean    sd 0.025quant 0.5quant
    ## Precision for the Gaussian observations 0.02 0.002      0.016     0.02
    ## Theta1 for field                        3.04 0.166      2.725     3.03
    ## Theta2 for field                        2.02 0.179      1.677     2.01
    ## Theta3 for field                        2.67 0.189      2.289     2.68
    ## Theta4 for field                        1.66 0.154      1.342     1.67
    ## Theta5 for field                        1.23 0.130      1.027     1.22
    ##                                         0.975quant  mode
    ## Precision for the Gaussian observations      0.025 0.019
    ## Theta1 for field                             3.378 3.008
    ## Theta2 for field                             2.380 1.993
    ## Theta3 for field                             3.031 2.699
    ## Theta4 for field                             1.948 1.691
    ## Theta5 for field                             1.527 1.155
    ## 
    ## Marginal log-Likelihood:  -1242.54 
    ##  is computed 
    ## Posterior summaries for the linear predictor and the fitted values are computed
    ## (Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')

We can obtain outputs with respect to parameters in the original scale
by using the function
[`rspde.result()`](https://davidbolin.github.io/rSPDE/reference/rspde.result.html):

``` r

result_fit_nonstat <- rspde.result(rspde_fit_nonstat, "field", rspde_model_nonstat)
summary(result_fit_nonstat)
```

    ##                  mean        sd 0.025quant 0.5quant 0.975quant    mode
    ## Theta1.matern 3.03611 0.1661290    2.72466  3.03110    3.37801 3.00840
    ## Theta2.matern 2.01598 0.1787190    1.67706  2.01176    2.38046 1.99289
    ## Theta3.matern 2.67387 0.1885230    2.28897  2.67847    3.03091 2.69906
    ## Theta4.matern 1.66129 0.1541420    1.34210  1.66655    1.94813 1.69082
    ## nu            1.54736 0.0441706    1.47293  1.54312    1.64200 1.52303

We can also plot the posterior densities. To this end we will use the
[`gg_df()`](https://davidbolin.github.io/MetricGraph/reference/gg_df.metric_graph_spde_result.md)
function, which creates `ggplot2` user-friendly data frames:

``` r

posterior_df_fit <- gg_df(result_fit_nonstat)

ggplot(posterior_df_fit) + geom_line(aes(x = x, y = y)) + 
facet_wrap(~parameter, scales = "free") + labs(y = "Density")
```

![](fem_models_files/figure-html/plot_post_nonstat-1.png)

## References

Bolin, David, Mihály Kovács, Vivek Kumar, and Alexandre B. Simas. 2023.
“Regularity and Numerical Approximation of Fractional Elliptic
Differential Equations on Compact Metric Graphs.” *Mathematics of
Computation*.

Bolin, David, Alexandre B. Simas, and Zhen Xiong. 2023.
“Covariance-Based Rational Approximations of Fractional SPDEs for
Computationally Efficient Bayesian Inference.” *Journal of Computational
and Graphical Statistics*.
