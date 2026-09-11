# Profile likelihood and REML for Whittle--Matérn fields

## Introduction

When fitting a Whittle–Matérn field model
``` math
  Y = X\beta + u(s) + \varepsilon, \quad \varepsilon \sim N(0, \sigma_\varepsilon^2 I),
```
the log-likelihood is quadratic in the fixed-effect coefficients
$`\beta`$ for any fixed covariance parameters
$`\theta = (\log\sigma_\varepsilon, \log(1/\tau), \log\kappa)`$.
Profiling $`\beta`$ out analytically (GLS) reduces the numerical
optimization to $`\theta \in \mathbb{R}^3`$ regardless of how many fixed
effects $`p`$ there are:
``` math
  2\ell_p(\theta) = 2\ell_0(\theta) + h^\top H^{-1} h,
```
where $`\ell_0(\theta)`$ is the no-covariate likelihood and $`H`$, $`h`$
are the GLS information matrix and right-hand side. The profiled
maximizer is $`\hat\beta(\theta) = H^{-1}h`$.

The **restricted likelihood** (REML) further penalises by
$`-\frac{1}{2}\log|H_\theta|`$, removing the contribution of
$`\hat\beta`$ to the curvature.

This vignette uses the Middle Fork river network from the `SSN2` package
(**Peterson2010?**) (introduced in the [River graph
example](https://davidbolin.github.io/MetricGraph/articles/river_example.md))
with **simulated covariates and response** to demonstrate the
profile-likelihood API at realistic scales ($`n \in \{200, 400\}`$,
$`p \in \{5, 10\}`$) using the directional $`\alpha = 1`$ Whittle–Matérn
model.

## River graph

We load the Middle Fork network (163 river segments, 261 km) exactly as
in the river example vignette.

``` r

library(SSN2)
copy_lsn_to_temp()
path <- file.path(tempdir(), "MiddleFork04.ssn")
mf04p <- ssn_import(path = path, predpts = c("pred1km", "CapeHorn"),
                    overwrite = TRUE)
graph <- metric_graph$new(mf04p, check_connected = FALSE)
cat("nE:", graph$nE, " | total length:", round(sum(graph$edge_lengths), 1), "km\n")
```

    ## nE: 163  | total length: 260.9 km

## Parameterization helpers

The profile-likelihood functions use the SPDE parameterization
$`\theta = (\log\sigma_\varepsilon, \log(1/\tau), \log\kappa)`$. For the
directional $`\alpha = 1`$ model ($`\nu = 0.5`$), the conversion to the
user-facing $`(\sigma, \text{range})`$ is:

``` r

nu   <- 0.5
C2   <- sqrt(gamma(nu) / ((4 * pi)^0.5 * gamma(nu + 0.5)))   # ≈ 0.707

range_to_kappa <- function(r)         sqrt(8 * nu) / r
kappa_to_range <- function(k)         sqrt(8 * nu) / k
sigma_to_tau   <- function(sig, k)    C2 / (sig * k^nu)
tau_to_sigma   <- function(tau, k)    C2 / (tau * k^nu)

make_theta <- function(sigma_e, sig, r) {
  k   <- range_to_kappa(r)
  tau <- sigma_to_tau(sig, k)
  c(log(sigma_e), log(1 / tau), log(k))
}
```

## Simulation setup

We sample observation locations proportionally to edge length, ensuring
at least one location per edge to give adequate spatial coverage.

``` r

sample_locs <- function(n, seed = 139) {
  set.seed(seed)
  nE    <- graph$nE
  extra <- n - nE
  extra_edges <- if (extra > 0)
    sample(nE, extra, replace = TRUE,
           prob = graph$edge_lengths / sum(graph$edge_lengths))
  else integer(0)
  edges <- c(seq_len(nE), extra_edges)
  cbind(edge     = edges[order(edges)],
        distance = runif(length(edges)))
}
```

We simulate from the **directional $`\alpha = 1`$ model** with

``` math
  \beta_0 = 5,\quad \beta_j = 1\ (j = 1,\ldots,p-1),\quad
  \sigma = 1.5,\quad \text{range} = 1\text{ km},\quad \sigma_\varepsilon = 0.5.
```

``` r

sigma_e_true <- 0.5
sigma_true   <- 1.5
range_true   <- 1       # km
theta_true   <- make_theta(sigma_e_true, sigma_true, range_true)

simulate_data <- function(PtE, p, seed = 139) {
  n    <- nrow(PtE)
  set.seed(seed)
  Xcov <- matrix(rnorm(n * (p - 1)), n, p - 1,
                 dimnames = list(NULL, paste0("x", seq_len(p - 1))))
  Xfull <- cbind(1, Xcov)
  beta <- c(intercept = 5, setNames(rep(1, p - 1), paste0("x", seq_len(p - 1))))
  u    <- sample_spde(range = range_true, sigma = sigma_true, alpha = 1,
                      directional = TRUE, graph = graph, PtE = PtE)
  y    <- as.vector(Xfull %*% beta) + u + sigma_e_true * rnorm(n)
  list(y = y, X = Xfull, Xcov = Xcov, beta = beta)
}
```

## Profile likelihood — $`n = 200`$, $`p = 5`$

``` r

PtE_200 <- sample_locs(200)
dat_200_5 <- simulate_data(PtE_200, p = 5)

graph$clear_observations()
graph$add_observations(
  data = data.frame(y = dat_200_5$y,
                    edge_number      = PtE_200[, "edge"],
                    distance_on_edge = PtE_200[, "distance"]),
  normalized = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are normalized by the length of the edge.

    ## The unit for edge lengths is km

    ## The current tolerance for removing distant observations is (in km): 3.16483618132507

``` r

pcomp_200_5 <- MetricGraph:::precompute_alpha1_directional(
  graph, manual_y = dat_200_5$y, X_cov = dat_200_5$X, repl = NULL)
cat("n =", nrow(PtE_200), " | p =", ncol(dat_200_5$X), "\n")
```

    ## n = 200  | p = 5

### Verifying the core identity

``` r

lp <- MetricGraph:::likelihood_alpha1_directional_profile_precompute(
  theta_true, precomputed_data = pcomp_200_5,
  parameterization = "spde", reml = FALSE)

bhat_true <- MetricGraph:::profile_beta_estimate(
  theta_true, model = "alpha1_directional", graph = graph,
  precomputed_data = pcomp_200_5, X_cov = dat_200_5$X,
  parameterization = "spde")$beta

lj <- MetricGraph:::likelihood_alpha1_directional(
  c(theta_true, bhat_true), graph = graph,
  manual_y = dat_200_5$y, X_cov = dat_200_5$X,
  repl = NULL, parameterization = "spde")

cat("Profile log-lik:          ", round(lp, 6), "\n")
```

    ## Profile log-lik:           -347.2387

``` r

cat("Joint at beta_hat:        ", round(lj, 6), "\n")
```

    ## Joint at beta_hat:         -347.2387

``` r

cat("Difference (machine eps): ", format(abs(lp - lj), scientific = TRUE), "\n")
```

    ## Difference (machine eps):  2.728484e-12

### Optimization and parameter recovery

``` r

neg_lp_fn <- function(pcomp) {
  function(theta) {
    val <- tryCatch(
      -MetricGraph:::likelihood_alpha1_directional_profile_precompute(
        theta, precomputed_data = pcomp,
        parameterization = "spde", reml = FALSE),
      error = function(e) Inf)
    if (!is.finite(val)) Inf else val
  }
}

theta_init <- make_theta(0.4, 1.2, 0.8)

opt_200_5 <- optim(theta_init, neg_lp_fn(pcomp_200_5), method = "BFGS",
                   control = list(reltol = 1e-10, maxit = 500))

khat <- exp(opt_200_5$par[3])
that <- 1 / exp(opt_200_5$par[2])

res_200_5 <- MetricGraph:::profile_beta_estimate(
  opt_200_5$par, model = "alpha1_directional", graph = graph,
  precomputed_data = pcomp_200_5, X_cov = dat_200_5$X,
  parameterization = "spde")

cov_table_200_5 <- data.frame(
  sigma_e = c(sigma_e_true, exp(opt_200_5$par[1])),
  sigma   = c(sigma_true,   tau_to_sigma(that, khat)),
  range   = c(range_true,   kappa_to_range(khat)),
  row.names = c("Truth", "Estimate"))
print(round(cov_table_200_5, 3))
```

    ##          sigma_e sigma range
    ## Truth      0.500 1.500 1.000
    ## Estimate   0.602 1.328 0.929

``` r

beta_table_200_5 <- rbind(dat_200_5$beta, res_200_5$beta)
rownames(beta_table_200_5) <- c("Truth", "Estimate")
print(round(beta_table_200_5, 3))
```

    ##          intercept    x1    x2    x3    x4
    ## Truth        5.000 1.000 1.000 1.000 1.000
    ## Estimate     5.089 0.901 0.996 1.006 1.107

``` r

cat("Beta RMSE:", round(sqrt(mean((res_200_5$beta - dat_200_5$beta)^2)), 4), "\n")
```

    ## Beta RMSE: 0.0766

## Profile likelihood — $`n = 200`$, $`p = 10`$

``` r

dat_200_10 <- simulate_data(PtE_200, p = 10)

graph$clear_observations()
graph$add_observations(
  data = data.frame(y = dat_200_10$y,
                    edge_number      = PtE_200[, "edge"],
                    distance_on_edge = PtE_200[, "distance"]),
  normalized = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are normalized by the length of the edge.

    ## The unit for edge lengths is km

    ## The current tolerance for removing distant observations is (in km): 3.16483618132507

``` r

pcomp_200_10 <- MetricGraph:::precompute_alpha1_directional(
  graph, manual_y = dat_200_10$y, X_cov = dat_200_10$X, repl = NULL)

opt_200_10 <- optim(theta_init, neg_lp_fn(pcomp_200_10), method = "BFGS",
                    control = list(reltol = 1e-10, maxit = 500))

khat10 <- exp(opt_200_10$par[3]); that10 <- 1 / exp(opt_200_10$par[2])
res_200_10 <- MetricGraph:::profile_beta_estimate(
  opt_200_10$par, model = "alpha1_directional", graph = graph,
  precomputed_data = pcomp_200_10, X_cov = dat_200_10$X,
  parameterization = "spde")

cov_table_200_10 <- data.frame(
  sigma_e = c(sigma_e_true, exp(opt_200_10$par[1])),
  sigma   = c(sigma_true,   tau_to_sigma(that10, khat10)),
  range   = c(range_true,   kappa_to_range(khat10)),
  row.names = c("Truth", "Estimate"))
print(round(cov_table_200_10, 3))
```

    ##          sigma_e sigma range
    ## Truth      0.500 1.500 1.000
    ## Estimate   0.588 1.486 0.819

``` r

cat("Beta RMSE (p=10):", round(sqrt(mean((res_200_10$beta - dat_200_10$beta)^2)), 4), "\n")
```

    ## Beta RMSE (p=10): 0.0991

## Profile likelihood — $`n = 400`$, $`p = 5`$ and $`p = 10`$

``` r

PtE_400 <- sample_locs(400)
cat("n =", nrow(PtE_400), "| edges covered:", length(unique(PtE_400[, "edge"])),
    "of", graph$nE, "\n")
```

    ## n = 400 | edges covered: 163 of 163

``` r

dat_400_5 <- simulate_data(PtE_400, p = 5)

graph$clear_observations()
graph$add_observations(
  data = data.frame(y = dat_400_5$y,
                    edge_number      = PtE_400[, "edge"],
                    distance_on_edge = PtE_400[, "distance"]),
  normalized = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are normalized by the length of the edge.

    ## The unit for edge lengths is km

    ## The current tolerance for removing distant observations is (in km): 3.16483618132507

``` r

pcomp_400_5 <- MetricGraph:::precompute_alpha1_directional(
  graph, manual_y = dat_400_5$y, X_cov = dat_400_5$X, repl = NULL)

opt_400_5 <- optim(theta_init, neg_lp_fn(pcomp_400_5), method = "BFGS",
                   control = list(reltol = 1e-10, maxit = 500))

k5 <- exp(opt_400_5$par[3]); t5 <- 1 / exp(opt_400_5$par[2])
res_400_5 <- MetricGraph:::profile_beta_estimate(
  opt_400_5$par, model = "alpha1_directional", graph = graph,
  precomputed_data = pcomp_400_5, X_cov = dat_400_5$X,
  parameterization = "spde")
```

``` r

dat_400_10 <- simulate_data(PtE_400, p = 10)

graph$clear_observations()
graph$add_observations(
  data = data.frame(y = dat_400_10$y,
                    edge_number      = PtE_400[, "edge"],
                    distance_on_edge = PtE_400[, "distance"]),
  normalized = TRUE)
```

    ## Adding observations...

    ## Assuming the observations are normalized by the length of the edge.

    ## The unit for edge lengths is km

    ## The current tolerance for removing distant observations is (in km): 3.16483618132507

``` r

pcomp_400_10 <- MetricGraph:::precompute_alpha1_directional(
  graph, manual_y = dat_400_10$y, X_cov = dat_400_10$X, repl = NULL)

opt_400_10 <- optim(theta_init, neg_lp_fn(pcomp_400_10), method = "BFGS",
                    control = list(reltol = 1e-10, maxit = 500))

k10 <- exp(opt_400_10$par[3]); t10 <- 1 / exp(opt_400_10$par[2])
res_400_10 <- MetricGraph:::profile_beta_estimate(
  opt_400_10$par, model = "alpha1_directional", graph = graph,
  precomputed_data = pcomp_400_10, X_cov = dat_400_10$X,
  parameterization = "spde")
```

## Summary across scenarios

``` r

scenarios <- list(
  list(lab="n=200 p=5",  opt=opt_200_5,  res=res_200_5,  beta=dat_200_5$beta),
  list(lab="n=200 p=10", opt=opt_200_10, res=res_200_10, beta=dat_200_10$beta),
  list(lab="n=400 p=5",  opt=opt_400_5,  res=res_400_5,  beta=dat_400_5$beta),
  list(lab="n=400 p=10", opt=opt_400_10, res=res_400_10, beta=dat_400_10$beta))

summary_df <- do.call(rbind, lapply(scenarios, function(s) {
  k   <- exp(s$opt$par[3]); t <- 1 / exp(s$opt$par[2])
  data.frame(
    scenario  = s$lab,
    sigma_e   = round(exp(s$opt$par[1]), 3),
    sigma     = round(tau_to_sigma(t, k), 3),
    range_km  = round(kappa_to_range(k), 3),
    beta_rmse = round(sqrt(mean((s$res$beta - s$beta)^2)), 4),
    conv      = s$opt$convergence)
}))
print(summary_df, row.names = FALSE)
```

    ##    scenario sigma_e sigma range_km beta_rmse conv
    ##   n=200 p=5   0.602 1.328    0.929    0.0766    0
    ##  n=200 p=10   0.588 1.486    0.819    0.0991    0
    ##   n=400 p=5   1.109 1.048    0.880    0.0599    0
    ##  n=400 p=10   1.036 1.132    1.037    0.0624    0

True values: $`\sigma_\varepsilon = 0.5`$, $`\sigma = 1.5`$, range
$`= 1`$ km. The covariance parameters show moderate variability across
scenarios, which is expected given the sparse observation design on a
261 km network — see the [River graph
example](https://davidbolin.github.io/MetricGraph/articles/river_example.md)
for a discussion of model identifiability. The fixed-effect vector
$`\beta`$ is consistently recovered with RMSE below 0.10 across all four
scenarios.

## Comparison with `graph_lme` ($`n = 400`$, $`p = 5`$)

The high-level `graph_lme(..., model = "wmd1")` jointly optimizes over
$`(\theta, \beta)`$. To compare on a common scale we evaluate
`likelihood_alpha1_directional` (the non-profiled joint likelihood) at
both sets of parameter estimates.

``` r

# Restore n=400, p=5 observations, now including covariates in graph data
graph$clear_observations()
graph$add_observations(
  data = data.frame(y                = dat_400_5$y,
                    as.data.frame(dat_400_5$Xcov),
                    edge_number      = PtE_400[, "edge"],
                    distance_on_edge = PtE_400[, "distance"]),
  normalized = TRUE)

lme_formula <- as.formula(paste("y ~", paste(colnames(dat_400_5$Xcov), collapse = " + ")))
fit_lme <- graph_lme(lme_formula, graph = graph, model = "wmd1")

lme_se    <- fit_lme$coeff$measurement_error
lme_sigma <- fit_lme$matern_coeff$random_effects["sigma"]
lme_range <- fit_lme$matern_coeff$random_effects["range"]
lme_beta  <- fit_lme$coeff$fixed_effects
lme_theta <- make_theta(lme_se, lme_sigma, lme_range)

lme_vs_profile <- data.frame(
  sigma_e   = c(sigma_e_true, exp(opt_400_5$par[1]), lme_se),
  sigma     = c(sigma_true,   tau_to_sigma(t5, k5),  lme_sigma),
  range     = c(range_true,   kappa_to_range(k5),    lme_range),
  row.names = c("Truth", "Profile", "graph_lme"))
print(round(lme_vs_profile, 4))
```

    ##           sigma_e  sigma  range
    ## Truth      0.5000 1.5000 1.0000
    ## Profile    1.1085 1.0480 0.8800
    ## graph_lme  0.4827 1.4559 0.6756

``` r

# Evaluate the joint log-likelihood at both sets of estimates using the same
# function so the values are on an identical scale.
lj_profile <- MetricGraph:::likelihood_alpha1_directional(
  c(opt_400_5$par, res_400_5$beta),
  graph = graph, manual_y = dat_400_5$y,
  X_cov = dat_400_5$X, repl = NULL, parameterization = "spde")

lj_lme <- MetricGraph:::likelihood_alpha1_directional(
  c(lme_theta, lme_beta),
  graph = graph, manual_y = dat_400_5$y,
  X_cov = dat_400_5$X, repl = NULL, parameterization = "spde")

cat(sprintf("Joint loglik at profile (theta_hat, beta_hat): %10.4f\n", lj_profile))
```

    ## Joint loglik at profile (theta_hat, beta_hat):  -723.5477

``` r

cat(sprintf("Joint loglik at graph_lme (theta_hat, beta):   %10.4f\n", lj_lme))
```

    ## Joint loglik at graph_lme (theta_hat, beta):    -743.3861

``` r

cat(sprintf("Difference (profile - graph_lme):              %+10.4f\n", lj_profile - lj_lme))
```

    ## Difference (profile - graph_lme):                +19.8384

## REML ($`n = 400`$, $`p = 10`$)

``` r

lp_val <- MetricGraph:::likelihood_alpha1_directional_profile_precompute(
  opt_400_10$par, precomputed_data = pcomp_400_10,
  parameterization = "spde", reml = FALSE)
reml_val <- MetricGraph:::likelihood_alpha1_directional_profile_precompute(
  opt_400_10$par, precomputed_data = pcomp_400_10,
  parameterization = "spde", reml = TRUE)
penalty <- -0.5 * as.numeric(determinant(res_400_10$H)$modulus)

cat("Profile log-lik:         ", round(lp_val,  4), "\n")
```

    ## Profile log-lik:          -719.5523

``` r

cat("REML log-lik:            ", round(reml_val, 4), "\n")
```

    ## REML log-lik:             -745.5261

``` r

cat("Penalty -½ log|H|:       ", round(penalty, 4), "\n")
```

    ## Penalty -½ log|H|:        -25.9738

``` r

cat("REML = Profile + Penalty:", isTRUE(all.equal(reml_val, lp_val + penalty)), "\n")
```

    ## REML = Profile + Penalty: TRUE

Optimizing the REML criterion for $`n = 400`$, $`p = 10`$:

``` r

neg_reml_fn <- function(pcomp) {
  function(theta) {
    val <- tryCatch(
      -MetricGraph:::likelihood_alpha1_directional_profile_precompute(
        theta, precomputed_data = pcomp,
        parameterization = "spde", reml = TRUE),
      error = function(e) Inf)
    if (!is.finite(val)) Inf else val
  }
}

opt_reml <- optim(theta_init, neg_reml_fn(pcomp_400_10), method = "BFGS",
                  control = list(reltol = 1e-10, maxit = 500))

kr <- exp(opt_reml$par[3]); tr <- 1 / exp(opt_reml$par[2])
reml_cov <- data.frame(
  sigma_e = c(sigma_e_true, exp(opt_reml$par[1])),
  sigma   = c(sigma_true,   tau_to_sigma(tr, kr)),
  range   = c(range_true,   kappa_to_range(kr)),
  row.names = c("Truth", "REML"))
print(round(reml_cov, 3))
```

    ##       sigma_e sigma range
    ## Truth   0.500 1.500 1.000
    ## REML    1.058 1.137 1.058

## References
