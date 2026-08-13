# likelihood_directional_ou_covariance.R
#
# Gaussian log-likelihood (with beta profiled out) for the "proper global OU"
# directional model, computed via the dense closed-form covariance of
# R/covariance_directional_ou.R (directional_ou_covariance()) rather than the
# sparse edge/precision-matrix route used by likelihood_alpha1_directional*
# (R/graph_likelihoods.R, R/graph_likelihoods_v2.R). This is a slow-but-simple
# reference implementation -- O(n^2) to build Sigma and O(n^3) to factor it --
# useful for cross-checking the sparse likelihood and for small graphs; it is
# not a replacement for the sparse path on anything but tiny fixtures (see
# DIRECTIONAL_OU_MAX_POINTS in R/covariance_directional_ou.R).
#
# theta convention (pinned, matches check_theta_profile() /
# profile_lik_core_alpha1_directional() in R/graph_likelihoods_v2.R):
#   theta = (log sigma_e, log reciprocal_tau, log kappa)
# with reciprocal_tau <- exp(theta[2]) and tau <- 1 / reciprocal_tau -- tau is
# NOT exp(theta[2]) directly. This tau is what gets passed to
# directional_ou_covariance()'s `tau` argument (its pinned convention:
# stationary variance = 1/(2*kappa*tau^2), same as MetricGraph:::r_1()).
#
# Uses the existing profile core/finish contract from R/graph_likelihoods_v2.R
# (check_theta_profile(), profile_lik_finish()): any *_core_* function need
# only return list(base, H, h, n_cov) with base = log-lik at beta=0 on the raw
# y, H = X'Sigma^-1 X, h = X'Sigma^-1 y; profile_lik_finish() completes the
# square over beta.


#' Core quantities for the dense directional-OU-covariance profiled
#' log-likelihood
#'
#' Builds the dense covariance `Sigma = directional_ou_covariance(...) +
#' sigma_e^2 I` at the observation locations and evaluates the no-covariate
#' log-likelihood with the raw response `y` (`base`), plus the profiled
#' fixed-effect accumulators `H = X' Sigma^-1 X` and `h = X' Sigma^-1 y`, via
#' a single Cholesky factor of `Sigma`. Output is consumable by
#' [profile_lik_finish()] exactly like `profile_lik_core_alpha1()` /
#' `profile_lik_core_alpha2()` / `profile_lik_core_alpha1_directional()`.
#'
#' @param theta `(log sigma_e, log reciprocal_tau, log kappa)`; length 3,
#'  beta is profiled out and must not be part of `theta`. Note
#'  `reciprocal_tau <- exp(theta[2])`, `tau <- 1 / reciprocal_tau` -- this
#'  `tau` (not `exp(theta[2])`) is what is passed to
#'  [directional_ou_covariance()].
#' @param graph A `metric_graph` object with directional edge weights and
#'  vertex weight functions set, as required by
#'  [directional_ou_covariance()].
#' @param data_name Name of the response variable attached to `graph`
#'  (ignored if `manual_y` is supplied).
#' @param manual_y Manual response vector, in the same row order as
#'  `graph$get_PtE()` (required, instead of `data_name`, whenever `PtE` is
#'  explicitly supplied -- see `PtE` below).
#' @param X_cov Optional fixed-effect design matrix, `n x p`, in the same row
#'  order as `y`. `NULL` (default) means zero covariates (`n_cov = 0`).
#' @param PtE Evaluation points passed to [directional_ou_covariance()].
#'  `NULL` (default) uses `graph$get_PtE()`, which is also the row order
#'  `data_name`/`X_cov` are assumed to follow. If `PtE` is explicitly
#'  supplied, `manual_y` must be supplied too (there is no safe way to
#'  subset a `graph`-attached data column to match an arbitrary
#'  caller-supplied `PtE`).
#' @param sigma_source Anchoring variance at source vertices, forwarded to
#'  [directional_ou_covariance()]. **Must be left at its default `NULL`** to
#'  match `likelihood_alpha1_directional_profile_precompute()` (which
#'  internally always uses `stationary_points = "all"`); a silently
#'  different `sigma_source` here would make any comparison between this
#'  likelihood and the sparse one meaningless without erroring, so this is a
#'  modeling requirement, not a mere default.
#' @param normalized Passed to [directional_ou_covariance()]; `TRUE`
#'  (default) means `PtE`/observation `distance_on_edge` are in `[0, 1]`.
#' @param force_dense Forwarded to [directional_ou_gls_core_from_sigma()];
#'  `TRUE` forces the dense `base::chol()` route regardless of Sigma's
#'  sparsity. A benchmarking/comparison knob, not a modeling choice --
#'  see that function's docs.
#' @return `list(base, H, h, n_cov)` -- see R/graph_likelihoods_v2.R's
#'  profile core/finish contract.
#' @noRd
directional_ou_covariance_loglik_core <- function(theta, graph, data_name = NULL,
                                                   manual_y = NULL, X_cov = NULL,
                                                   PtE = NULL, sigma_source = NULL,
                                                   normalized = TRUE,
                                                   force_dense = FALSE,
                                                   cpp = TRUE) {
  check_theta_profile(theta)

  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  tau <- 1 / reciprocal_tau
  kappa <- exp(theta[3])

  if (!is.null(PtE) && is.null(manual_y)) {
    stop("PtE was supplied explicitly but manual_y was not: a graph-attached ",
         "data column (data_name) cannot be safely subset to match an ",
         "arbitrary caller-supplied PtE. Supply manual_y (in the same row ",
         "order as PtE) whenever PtE is not left at its default NULL.")
  }

  if (is.null(manual_y)) {
    if (is.null(data_name)) {
      stop("Either data_name or manual_y must be not NULL")
    }
    y <- graph$.__enclos_env__$private$data[[data_name]]
  } else {
    y <- manual_y
  }

  n_cov <- if (is.null(X_cov)) 0L else ncol(X_cov)

  Sigma <- directional_ou_covariance(graph, kappa, tau, PtE = PtE,
                                     sigma_source = sigma_source,
                                     normalized = normalized,
                                     cpp = cpp)
  diag(Sigma) <- diag(Sigma) + sigma_e^2

  directional_ou_gls_core_from_sigma(Sigma, y, X_cov, n_cov, force_dense = force_dense)
}

# Below this fraction of nonzero entries, directional_ou_gls_core_from_sigma()
# routes through a sparse Cholesky (Matrix::Cholesky()) instead of
# base::chol(). The dendritic directional-OU covariance is structurally (not
# just numerically) sparse: directional_ou_covariance_dendritic_cpp() sets
# Sigma(i,j) <- 0 exactly, for any pair of points with no last common
# ancestor in the flow order (src/directional_ou_covariance.cpp, the
# `Sigma(i, j) = 0.0;` branch) -- a graph-topology fact independent of
# kappa/tau/sigma_e, not a numerical coincidence at one theta. On the real
# Mid-Columbia network this is ~97-98% exact zeros regardless of n (measured
# at n = 1000..5000), so a dense O(n^3) factorization is doing a large
# amount of avoidable work on structural zeros. 0.3 is a deliberately
# conservative cutoff -- comfortably above the ~0.02-0.03 nonzero fraction
# actually observed, so genuinely dense inputs (small/non-dendritic graphs,
# arbitrary Sigma passed to directional_ou_covariance_loglik_core()) are
# unaffected and keep taking the original dense path.
DIRECTIONAL_OU_GLS_SPARSE_NNZ_THRESHOLD <- 0.3

# Shared Cholesky/GLS tail of the dense directional-OU-covariance profiled
# log-likelihood: given a covariance matrix `Sigma` (nugget already added on
# the diagonal), the response `y`, and the design matrix `X_cov` (or `NULL`),
# factor once and return the profile core/finish contract `list(base, H, h,
# n_cov)`. Shared by [directional_ou_covariance_loglik_core()] (Sigma via
# [directional_ou_covariance()]) and
# [directional_ou_covariance_loglik_precompute()] (Sigma via
# [directional_ou_covariance_from_setup()]) -- the two differ only in how
# Sigma is assembled, not in what happens to it afterwards.
#
# Sigma's nonzero fraction decides the linear-algebra route (see
# DIRECTIONAL_OU_GLS_SPARSE_NNZ_THRESHOLD above); both routes compute
# exactly the same mathematical quantities (base = Gaussian log-density at
# beta = 0, H = X'Sigma^-1 X, h = X'Sigma^-1 y), verified to agree with the
# dense path to ~1e-12 (floating-point noise) on real Mid-Columbia Sigma
# matrices up to n = 5000 -- this is purely an implementation/performance
# choice, not a change to the statistical model.
#
# force_dense = TRUE skips the nnz_frac check and always takes the
# base::chol() route below, regardless of how sparse Sigma actually is --
# a benchmarking/comparison knob, not something a real caller should need.
# Note this
# does not avoid the DIRECTIONAL_OU_MAX_POINTS guard in
# directional_ou_covariance(): Sigma is always materialized as a dense
# matrix before this function ever sees it, whichever route runs.
#' @noRd
directional_ou_gls_core_from_sigma <- function(Sigma, y, X_cov, n_cov, force_dense = FALSE) {
  n <- length(y)
  nnz_frac <- sum(Sigma != 0) / (as.double(n) * ncol(Sigma))

  if (!force_dense && nnz_frac <= DIRECTIONAL_OU_GLS_SPARSE_NNZ_THRESHOLD) {
    Sigma_sp <- Matrix::forceSymmetric(Matrix::Matrix(Sigma, sparse = TRUE))
    ch <- Matrix::Cholesky(Sigma_sp, LDL = FALSE)

    logdet <- as.numeric(2 * Matrix::determinant(ch, sqrt = TRUE, logarithm = TRUE)$modulus)
    SigmaInv_y <- as.numeric(Matrix::solve(ch, y, system = "A"))
    base <- -0.5 * (n * log(2 * pi) + logdet + sum(y * SigmaInv_y))

    if (n_cov > 0) {
      SigmaInv_X <- as.matrix(Matrix::solve(ch, X_cov, system = "A"))
      H <- crossprod(X_cov, SigmaInv_X)
      h <- as.vector(crossprod(X_cov, SigmaInv_y))
    } else {
      H <- NULL
      h <- NULL
    }

    return(list(base = base, H = H, h = h, n_cov = n_cov))
  }

  Rchol <- base::chol(Sigma)  # upper triangular, t(Rchol) %*% Rchol == Sigma

  yt <- backsolve(Rchol, y, transpose = TRUE)

  base <- -0.5 * (n * log(2 * pi) + 2 * sum(log(diag(Rchol))) + sum(yt^2))

  if (n_cov > 0) {
    Xt <- backsolve(Rchol, X_cov, transpose = TRUE)
    H <- crossprod(Xt)
    h <- as.vector(crossprod(Xt, yt))
  } else {
    H <- NULL
    h <- NULL
  }

  list(base = base, H = H, h = h, n_cov = n_cov)
}


#' Profiled log-likelihood for the dense directional-OU-covariance model
#'
#' The dense-covariance counterpart of
#' `likelihood_alpha1_directional_profile()` /
#' `likelihood_alpha1_directional_profile_precompute()`
#' (R/graph_likelihoods_v2.R): same profiled-beta model (`theta = (log
#' sigma_e, log reciprocal_tau, log kappa)`, beta profiled out analytically),
#' but the covariance is assembled directly via
#' [directional_ou_covariance()] instead of the sparse edge/precision-matrix
#' route. Intended as a slow, simple reference/cross-check on small graphs
#' (see [directional_ou_covariance_loglik_core()] and
#' [directional_ou_covariance()]'s `DIRECTIONAL_OU_MAX_POINTS` guard), not a
#' performance replacement for the sparse likelihood.
#'
#' Out of scope, deliberately (not oversights):
#'   * No `repl` argument -- unlike the sparse likelihood, this v1 does not
#'     support multiple replicates.
#'   * No `profile_beta_estimate()`-style helper -- this v1 does not provide
#'     a way to recover `beta_hat` after optimizing `theta`; use
#'     [directional_ou_covariance_loglik_core()] directly and
#'     `profile_beta_solve()` if you need it.
#'
#' @inheritParams directional_ou_covariance_loglik_core
#' @param reml If `TRUE`, return the restricted (REML) log-likelihood
#'  instead of the profile log-likelihood.
#' @return A single finite numeric: the profile (or restricted)
#'  log-likelihood.
#' @noRd
directional_ou_covariance_loglik <- function(theta, graph, data_name = NULL,
                                             manual_y = NULL, X_cov = NULL,
                                             PtE = NULL, sigma_source = NULL,
                                             normalized = TRUE, reml = FALSE,
                                             cpp = TRUE) {
  core <- directional_ou_covariance_loglik_core(theta, graph, data_name = data_name,
                                                manual_y = manual_y, X_cov = X_cov,
                                                PtE = PtE, sigma_source = sigma_source,
                                                normalized = normalized,
                                                cpp = cpp)
  profile_lik_finish(core, reml)
}


#' Precompute the graph-structural part of the dense directional-OU-covariance
#' log-likelihood
#'
#' Caches everything [directional_ou_covariance_loglik_precompute()] needs
#' that does not depend on `theta`: the graph-structural (kappa/tau/
#' sigma_source-independent) setup piece
#' [directional_ou_setup_structure()], the resolved response `y`, the design
#' matrix `X_cov`, the evaluation points `PtE`, and `n_cov`. Split out so the
#' structural setup (topological edge order, beta_v inputs, Euler-tour
#' ancestor labels) is computed once and reused across repeated `theta`
#' evaluations (e.g. during likelihood optimization), mirroring
#' [directional_ou_setup_structure()]'s own rationale.
#'
#' @param graph A `metric_graph` object with directional edge weights and
#'  vertex weight functions set, as required by
#'  [directional_ou_covariance()].
#' @param data_name Name of the response variable attached to `graph`
#'  (ignored if `manual_y` is supplied).
#' @param manual_y Manual response vector, in the same row order as
#'  `graph$get_PtE()` (required, instead of `data_name`, whenever `PtE` is
#'  explicitly supplied -- see `PtE` below).
#' @param X_cov Optional fixed-effect design matrix, `n x p`, in the same row
#'  order as `y`. `NULL` (default) means zero covariates (`n_cov = 0`).
#' @param PtE Evaluation points, forwarded unevaluated to
#'  [directional_ou_covariance_from_setup()] at each
#'  [directional_ou_covariance_loglik_precompute()] call. `NULL` (default)
#'  means "use `graph$get_PtE()` at evaluation time" -- the same convention
#'  as [directional_ou_covariance()]'s own `PtE = NULL` default -- rather
#'  than freezing a snapshot of `graph$get_PtE()` here. If `PtE` is
#'  explicitly supplied, `manual_y` must be supplied too (there is no safe
#'  way to subset a `graph`-attached data column to match an arbitrary
#'  caller-supplied `PtE`), exactly as in
#'  [directional_ou_covariance_loglik_core()].
#' @return A list with components `structure` (from
#'  [directional_ou_setup_structure()]), `y`, `X_cov`, `PtE` (as passed in,
#'  possibly `NULL`), and `n_cov`.
#' @noRd
precompute_directional_ou_covariance <- function(graph, data_name = NULL,
                                                  manual_y = NULL,
                                                  X_cov = NULL, PtE = NULL) {
  if (!is.null(PtE) && is.null(manual_y)) {
    stop("PtE was supplied explicitly but manual_y was not: a graph-attached ",
         "data column (data_name) cannot be safely subset to match an ",
         "arbitrary caller-supplied PtE. Supply manual_y (in the same row ",
         "order as PtE) whenever PtE is not left at its default NULL.")
  }

  if (is.null(manual_y)) {
    if (is.null(data_name)) {
      stop("Either data_name or manual_y must be not NULL")
    }
    y <- graph$.__enclos_env__$private$data[[data_name]]
  } else {
    y <- manual_y
  }

  n_cov <- if (is.null(X_cov)) 0L else ncol(X_cov)

  structure <- directional_ou_setup_structure(graph)

  list(structure = structure, y = y, X_cov = X_cov, PtE = PtE, n_cov = n_cov)
}


#' Per-`theta` evaluator for the precomputed dense directional-OU-covariance
#' profiled log-likelihood
#'
#' The precomputed counterpart of
#' [directional_ou_covariance_loglik_core()]: given the output of
#' [precompute_directional_ou_covariance()], recomputes only the
#' `kappa`/`tau`/`sigma_source`-dependent piece of the setup
#' ([directional_ou_setup_numeric()]), builds `Sigma` via
#' [directional_ou_covariance_from_setup()] instead of rebuilding the
#' graph-structural setup from scratch, then does the same Cholesky-based
#' profiled GLS as [directional_ou_covariance_loglik_core()] (factored into
#' the shared [directional_ou_gls_core_from_sigma()] helper). Output is
#' consumable by [profile_lik_finish()], same as
#' [directional_ou_covariance_loglik_core()].
#'
#' @param theta `(log sigma_e, log reciprocal_tau, log kappa)`; length 3,
#'  beta is profiled out and must not be part of `theta`. Note
#'  `reciprocal_tau <- exp(theta[2])`, `tau <- 1 / reciprocal_tau` -- this
#'  `tau` (not `exp(theta[2])`) is what is passed to
#'  [directional_ou_setup_numeric()].
#' @param precomputed Output of [precompute_directional_ou_covariance()].
#' @param sigma_source Anchoring variance at source vertices, forwarded to
#'  [directional_ou_setup_numeric()]. **Must be left at its default `NULL`**
#'  for the same reason as in [directional_ou_covariance_loglik_core()]: a
#'  silently different `sigma_source` here would make any comparison with
#'  the sparse likelihood meaningless without erroring.
#' @param normalized Passed to [directional_ou_covariance_from_setup()];
#'  `TRUE` (default) means `precomputed$PtE`/observation
#'  `distance_on_edge` are in `[0, 1]`.
#' @param reml If `TRUE`, return the restricted (REML) log-likelihood
#'  instead of the profile log-likelihood.
#' @param cpp If `TRUE`, use [directional_ou_covariance_from_setup_cpp()]
#'  (the C++ in-tree covariance fill, with the mirrored fast R fill for
#'  out-trees and the generic R fallback for irregular trees) instead of
#'  [directional_ou_covariance_from_setup()] to build `Sigma`; `FALSE` always
#'  uses the pure-R path. Same
#'  spirit as `likelihood_alpha1_directional_profile_precompute()`'s own
#'  `cpp` argument (R/graph_likelihoods_v2.R), though that one is tri-state
#'  (`NULL`/not-`NULL`) while this one is a plain logical.
#' @param force_dense Forwarded to [directional_ou_gls_core_from_sigma()];
#'  `TRUE` forces the dense `base::chol()` route regardless of Sigma's
#'  sparsity. A benchmarking/comparison knob, not a modeling choice --
#'  see that function's docs. Does not affect how `Sigma` itself is built
#'  (still always a dense matrix, capped by
#'  [directional_ou_covariance()]'s `DIRECTIONAL_OU_MAX_POINTS` guard),
#'  only which linear-algebra route factors it.
#' @return A single finite numeric: the profile (or restricted)
#'  log-likelihood.
#' @noRd
directional_ou_covariance_loglik_precompute <- function(theta, precomputed,
                                                         sigma_source = NULL,
                                                         normalized = TRUE,
                                                         reml = FALSE,
                                                         cpp = TRUE,
                                                         force_dense = FALSE) {
  check_theta_profile(theta)

  sigma_e <- exp(theta[1])
  reciprocal_tau <- exp(theta[2])
  tau <- 1 / reciprocal_tau
  kappa <- exp(theta[3])

  structure <- precomputed$structure
  numeric_part <- directional_ou_setup_numeric_dispatch(
    structure, kappa, tau, sigma_source, cpp = cpp
  )
  setup <- directional_ou_compose_setup(structure, numeric_part)

  build_Sigma <- if (cpp) directional_ou_covariance_from_setup_cpp else directional_ou_covariance_from_setup
  Sigma <- build_Sigma(setup, PtE = precomputed$PtE, normalized = normalized)
  diag(Sigma) <- diag(Sigma) + sigma_e^2

  core <- directional_ou_gls_core_from_sigma(Sigma, precomputed$y,
                                             precomputed$X_cov,
                                             precomputed$n_cov,
                                             force_dense = force_dense)
  profile_lik_finish(core, reml)
}
