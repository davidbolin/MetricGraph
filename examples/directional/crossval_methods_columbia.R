# Full-graph model fitting and plug-in LOO comparison on the packaged
# Mid-Columbia component. Covariance parameters are fitted once per method;
# they are not refitted for each omitted observation.
#
# Inputs : packaged data("columbia_main_component"); edit the Settings block
#          below to change the run.
# Outputs: parameter/score/timing tables printed to the console; a checkpoint
#          file if out_file is set below.
#
# Run from the package root with:
#
#   Rscript examples/directional/crossval_methods_columbia.R


## Settings ------------------------------------------------------------------

# Set to a file path (e.g. "out/crossval.rds") to save a checkpoint; "" runs
# the comparison and prints the summary tables without writing anything.
out_file <- ""
save_result <- nzchar(out_file)

response_name <- "STREAM_AUG"
weights_name <- "h2oAreaKm2"
covariates <- c("ELEV", "SLOPE", "PRECIP")

# Order fixed by the profile likelihood's theta argument (see
# likelihood_alpha1_directional_profile_precompute /
# likelihood_alpha1_profile_precompute), not by this file; initial_thetas,
# lower/upper, and the optimization diagnostics below all assume it.
theta_names <- c("log_sigma_e", "log_reciprocal_tau", "log_kappa")

# These starts are fixed before reading the response and are identical for
# every method. Requiring independent starts to reproduce the same optimum
# guards against a favourable hand-picked or method-specific starting value.
initial_thetas <- rbind(
  baseline = c(log(0.5), log(1), log(0.5)),
  unit = c(log(1), log(1), log(1)),
  long_range = c(log(1), log(1), log(0.05))
)
# Log-scale parameter bounds; wide enough to not bind at any of the starts
# above while still keeping nlminb() away from numerically degenerate
# sigma_e/tau/kappa values.
lower <- rep(-10, length(theta_names))
upper <- rep(10, length(theta_names))
# nlminb() iteration budget per start.
max_iterations <- 300L
# At least this many of the independent starts in initial_thetas must
# reproduce the objective/theta within tolerance before an optimum is
# accepted.
minimum_reproducing_starts <- 2L
# nlminb() function-evaluation budget per start, as a multiple of
# max_iterations.
eval_max_per_iteration <- 5L
# Negative-log-likelihood units; two runs are treated as the same optimum if
# their final objectives are within this gap.
objective_tolerance <- 1e-4
# Log-scale units; reproducing starts must also agree on theta within this.
theta_tolerance <- 1e-2
# Max |numerical gradient| at the accepted optimum (see the gradient check
# below); nlminb's convergence == 0 alone does not guarantee a stationary
# point for this profile likelihood.
gradient_tolerance <- 5e-2
# Log-scale distance from lower/upper treated as "at boundary".
bound_tolerance <- 1e-7
# Relative step for the independent finite-difference gradient check.
finite_difference_relative_step <- 1e-4
# Max abs difference allowed between the profile and LOO fixed-effect
# estimates, which are computed by two independent code paths.
beta_consistency_tolerance <- 1e-6
# optimization_objective() returns this for a rejected/non-finite evaluation;
# valid_objective_ceiling is the cutoff used to recognize such a value later,
# kept well below the penalty so the two can never be confused.
rejected_evaluation_penalty <- 1e100
valid_objective_ceiling <- 1e50

methods <- list(
  K1 = list(
    label = "K1 directional",
    model = "directional",
    reversed = FALSE,
    directional_rule = "K1",
    BC = NA_integer_
  ),
  K2 = list(
    label = "K2 directional",
    model = "directional",
    reversed = FALSE,
    directional_rule = "K2",
    BC = NA_integer_
  ),
  WM1_BC1 = list(
    label = "Symmetric WM1, stationary endpoints (BC=1)",
    model = "wm1",
    reversed = FALSE,
    directional_rule = NA_character_,
    BC = 1L
  ),
  K1_REVERSED_CONTINUITY = list(
    label = "Reversed-direction K1 continuity",
    model = "directional",
    reversed = TRUE,
    directional_rule = "K1",
    BC = NA_integer_
  )
  # WM2/BC1 has no entry here: its full-component optimization is
  # substantially slower than the other comparisons above.
)
# Set id from each entry's own list name, rather than typing it twice, so
# the two can never drift apart.
methods <- Map(function(id, spec) c(list(id = id), spec), names(methods), methods)

## Data and main component ---------------------------------------------------

pkgload::load_all(".", reset = TRUE, export_all = TRUE)
source("examples/directional/columbia_full_graph_helpers.R")

cat(
  "Using local MetricGraph ",
  as.character(utils::packageVersion("MetricGraph")),
  " from ",
  getNamespaceInfo(asNamespace("MetricGraph"), "path"),
  "\n",
  sep = ""
)

component_timing <- system.time({
  columbia_component <- MetricGraph::columbia_main_component
})

graph_timing <- system.time({
  graph_original <- columbia_make_graph(
    component = columbia_component,
    reversed = FALSE,
    weights_name = weights_name
  )
  graph_reversed <- columbia_make_graph(
    component = columbia_component,
    reversed = TRUE,
    weights_name = weights_name
  )
})
graph_data <- graph_original$get_data(format = "tibble", drop_na = FALSE)
y_full <- graph_data[[response_name]]
observation_ids <- graph_data$columbia_obs_id

design_timing <- system.time({
  covariate_formula <- stats::reformulate(covariates)
  covariate_data <- stats::model.frame(
    covariate_formula,
    data = graph_data,
    na.action = stats::na.pass
  )
  X <- stats::model.matrix(covariate_formula, data = covariate_data)
})

fitting_rows <- !is.na(y_full)

reversed_data <- graph_reversed$get_data(format = "tibble", drop_na = FALSE)
reversed_rows <- match(reversed_data$columbia_obs_id, observation_ids)
X_reversed <- X[reversed_rows, , drop = FALSE]
y_reversed <- reversed_data[[response_name]]
observation_ids_reversed <- reversed_data$columbia_obs_id


## Fit and cross-validate each method ----------------------------------------

# nlminb()$evaluations may be absent or missing a key; NA marks that case
# for the diagnostics table instead of failing the whole run.
nlminb_evaluation_count <- function(run, key) {
  if (!is.null(run$evaluations) && key %in% names(run$evaluations)) {
    unname(run$evaluations[[key]])
  } else {
    NA_integer_
  }
}

# Build the model-specific precomputed likelihood inputs and the profile
# likelihood (+ its fixed arguments) that optimize_covariance_parameters()
# will evaluate repeatedly. Also reads response_name from the Settings block.
precompute_likelihood_inputs <- function(method, graph, X_method) {
  precompute_timing <- system.time({
    if (method$model == "directional") {
      precomputed <- MetricGraph:::precompute_alpha1_directional(
        graph = graph,
        data_name = response_name,
        X_cov = X_method,
        repl = NULL
      )
    } else {
      precomputed <- MetricGraph:::precompute_alpha1(
        graph = graph,
        data_name = response_name,
        X_cov = X_method,
        repl = NULL
      )
    }
  })
  if (method$model == "directional") {
    profile_likelihood <-
      MetricGraph:::likelihood_alpha1_directional_profile_precompute
    profile_arguments <- list(
      precomputed_data = precomputed,
      parameterization = "spde"
    )
  } else {
    profile_likelihood <-
      MetricGraph:::likelihood_alpha1_profile_precompute
    profile_arguments <- list(
      graph = graph,
      precomputeddata = precomputed,
      BC = method$BC,
      parameterization = "spde"
    )
  }
  list(
    precomputed = precomputed,
    profile_likelihood = profile_likelihood,
    profile_arguments = profile_arguments,
    timing = precompute_timing
  )
}

# Multi-start nlminb() optimization of the profile likelihood, followed by
# the acceptance policy: take the best converged, interior run, but only if
# no rejected/boundary run is materially better, and at least
# minimum_reproducing_starts runs agree with it within objective_tolerance
# and theta_tolerance. Returns the accepted run plus everything
# verify_local_optimum() and the result assembly need afterwards.
# Also reads initial_thetas, lower, upper, max_iterations,
# eval_max_per_iteration, bound_tolerance, valid_objective_ceiling,
# rejected_evaluation_penalty, minimum_reproducing_starts,
# objective_tolerance, theta_tolerance, and theta_names from the
# Settings block.
optimize_covariance_parameters <- function(
    method, profile_likelihood, profile_arguments) {
  rejected_profile_evaluations <- 0L
  last_profile_error <- NA_character_
  raw_profile <- function(theta) {
    do.call(
      profile_likelihood,
      c(list(theta = theta), profile_arguments)
    )
  }
  optimization_objective <- function(theta) {
    profile_error <- NULL
    value <- tryCatch(
      raw_profile(theta),
      error = function(error) {
        profile_error <<- conditionMessage(error)
        NA_real_
      }
    )
    if (length(value) != 1L || !is.finite(value)) {
      rejected_profile_evaluations <<-
        rejected_profile_evaluations + 1L
      last_profile_error <<- if (is.null(profile_error)) {
        "Profile likelihood was not one finite scalar."
      } else {
        profile_error
      }
      return(rejected_evaluation_penalty)
    }
    -as.numeric(value)
  }

  optimization_runs <- vector("list", nrow(initial_thetas))
  names(optimization_runs) <- rownames(initial_thetas)
  optimization_diagnostics <- data.frame()
  optimization_timing <- system.time({
    for (start_index in seq_len(nrow(initial_thetas))) {
      start_id <- rownames(initial_thetas)[start_index]
      start_theta <- initial_thetas[start_index, ]
      failures_before <- rejected_profile_evaluations
      start_objective <- optimization_objective(start_theta)
      run <- stats::nlminb(
        start = start_theta,
        objective = optimization_objective,
        lower = lower,
        upper = upper,
        control = list(
          iter.max = as.integer(max_iterations),
          eval.max = as.integer(eval_max_per_iteration * max_iterations)
        )
      )
      # Recompute at run$par rather than trusting nlminb's own run$objective:
      # this forces the final point through the same tryCatch/finite-value
      # check as every other evaluation, so final_objective below can never
      # disagree with rejected_evaluations about whether this point is valid.
      verified_objective <- optimization_objective(run$par)
      at_lower_bound <- abs(run$par - lower) <= bound_tolerance
      at_upper_bound <- abs(run$par - upper) <= bound_tolerance
      run$objective <- verified_objective
      run$at_lower_bound <- at_lower_bound
      run$at_upper_bound <- at_upper_bound
      optimization_runs[[start_id]] <- run

      function_evaluations <- nlminb_evaluation_count(run, "function")
      gradient_evaluations <- nlminb_evaluation_count(run, "gradient")
      # Column names for the per-parameter start/final values are derived
      # from theta_names so they can't drift from it (see reproducing_thetas
      # below, which selects columns by the same derived names).
      start_theta_cols <- stats::setNames(
        as.list(start_theta), paste0("start_", theta_names)
      )
      final_theta_cols <- stats::setNames(
        as.list(run$par), paste0("final_", theta_names)
      )
      optimization_diagnostics <- rbind(
        optimization_diagnostics,
        do.call(data.frame, c(
          list(start_id = start_id),
          start_theta_cols,
          list(start_objective = start_objective),
          final_theta_cols,
          list(
            final_objective = verified_objective,
            convergence = run$convergence,
            message = if (is.null(run$message)) {
              NA_character_
            } else {
              run$message
            },
            iterations = run$iterations,
            function_evaluations = function_evaluations,
            gradient_evaluations = gradient_evaluations,
            rejected_evaluations =
              rejected_profile_evaluations - failures_before,
            at_bound = any(at_lower_bound | at_upper_bound),
            row.names = NULL
          )
        ))
      )
    }
  })

  print(
    optimization_diagnostics[
      , c(
        "start_id", "final_objective", "convergence", "iterations",
        "rejected_evaluations", "at_bound"
      )
    ],
    digits = 8,
    row.names = FALSE
  )
  valid_objective <- is.finite(optimization_diagnostics$final_objective) &
    optimization_diagnostics$final_objective < valid_objective_ceiling
  eligible <- valid_objective &
    optimization_diagnostics$convergence == 0L &
    !optimization_diagnostics$at_bound
  if (!any(eligible)) {
    stop(method$id, " did not produce a converged, interior optimum.")
  }
  eligible_indices <- which(eligible)
  best_index <- eligible_indices[which.min(
    optimization_diagnostics$final_objective[eligible_indices]
  )]
  best_objective <- optimization_diagnostics$final_objective[best_index]
  materially_better_unconverged <- valid_objective & !eligible &
    optimization_diagnostics$final_objective <
      best_objective - objective_tolerance
  if (any(materially_better_unconverged)) {
    stop(
      method$id,
      " has an unconverged or boundary run better than the accepted optimum."
    )
  }
  reproducing_indices <- which(
    eligible &
      abs(optimization_diagnostics$final_objective - best_objective) <=
        objective_tolerance
  )
  if (length(reproducing_indices) < minimum_reproducing_starts) {
    stop(
      method$id,
      " optimum was not reproduced by at least ",
      minimum_reproducing_starts,
      " independent starts."
    )
  }
  reproducing_thetas <- as.matrix(optimization_diagnostics[
    reproducing_indices,
    paste0("final_", theta_names)
  ])
  theta_spread <- apply(
    reproducing_thetas,
    2,
    function(value) max(value) - min(value)
  )
  if (any(theta_spread > theta_tolerance)) {
    stop(
      method$id,
      " reproduced the objective but not stable covariance parameters."
    )
  }

  list(
    optimization = optimization_runs[[best_index]],
    optimization_diagnostics = optimization_diagnostics,
    optimization_objective = optimization_objective,
    best_index = best_index,
    reproducing_indices = reproducing_indices,
    theta_spread = theta_spread,
    rejected_profile_evaluations = rejected_profile_evaluations,
    last_profile_error = last_profile_error,
    timing = optimization_timing
  )
}

# nlminb's convergence == 0 does not guarantee a stationary point for this
# profile likelihood; verify independently with a numerical gradient, then
# confirm the Hessian is positive definite (a genuine local maximum). Also
# reads finite_difference_relative_step, gradient_tolerance, and theta_names
# from the Settings block.
verify_local_optimum <- function(
    method, optimization, optimization_objective) {
  gradient_step <-
    finite_difference_relative_step * pmax(1, abs(optimization$par))
  numerical_gradient <- vapply(seq_along(optimization$par), function(index) {
    theta_plus <- optimization$par
    theta_minus <- optimization$par
    theta_plus[index] <- theta_plus[index] + gradient_step[index]
    theta_minus[index] <- theta_minus[index] - gradient_step[index]
    (
      optimization_objective(theta_plus) -
        optimization_objective(theta_minus)
    ) / (2 * gradient_step[index])
  }, numeric(1))
  names(numerical_gradient) <- theta_names
  max_abs_gradient <- max(abs(numerical_gradient))
  if (!is.finite(max_abs_gradient) || max_abs_gradient > gradient_tolerance) {
    stop(
      method$id,
      " failed the numerical gradient check (max |gradient| = ",
      signif(max_abs_gradient, 5),
      ")."
    )
  }
  objective_hessian <- stats::optimHess(
    optimization$par,
    optimization_objective
  )
  hessian_eigenvalues <- eigen(
    (objective_hessian + t(objective_hessian)) / 2,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (any(!is.finite(hessian_eigenvalues)) ||
      min(hessian_eigenvalues) <= 0) {
    stop(method$id, " did not finish at a local likelihood maximum.")
  }
  list(
    numerical_gradient = numerical_gradient,
    max_abs_gradient = max_abs_gradient,
    objective_hessian = objective_hessian,
    hessian_eigenvalues = hessian_eigenvalues
  )
}

# Profile out the fixed effects at the accepted covariance parameters.
recover_fixed_effects <- function(method, graph, precomputed, optimization) {
  beta_arguments <- list(
    theta = optimization$par,
    model = if (method$model == "directional") {
      "alpha1_directional"
    } else {
      "alpha1"
    },
    graph = graph,
    precomputed_data = precomputed,
    parameterization = "spde"
  )
  if (method$model != "directional") {
    beta_arguments$BC <- method$BC
  }
  fixed_effect_timing <- system.time({
    beta_estimate <- do.call(
      MetricGraph:::profile_beta_estimate,
      beta_arguments
    )
  })
  list(beta_estimate = beta_estimate, timing = fixed_effect_timing)
}

# Plug-in LOO predictions at the accepted parameters, cross-checked against
# recover_fixed_effects()'s independent fixed-effect estimate. Also reads
# beta_consistency_tolerance from the Settings block.
compute_plugin_loo <- function(
    method, graph, X_method, y_method, optimization, precomputed,
    beta_estimate) {
  loo_arguments <- list(
    theta = optimization$par,
    precomputed_data = precomputed,
    parameterization = "spde"
  )
  if (method$model == "directional") {
    loo_function <- MetricGraph:::cv_core_alpha1_directional
    loo_arguments$method <- "selinv"
  } else {
    loo_function <- MetricGraph:::cv_core_alpha1
    loo_arguments$graph <- graph
    loo_arguments$y_resp_full <- y_method
    loo_arguments$BC <- method$BC
  }
  loo_timing <- system.time({
    loo <- do.call(loo_function, loo_arguments)
  })

  n_predictions <- length(loo$idx)
  if (n_predictions == 0L ||
      length(loo$mu) != n_predictions ||
      length(loo$var) != n_predictions ||
      anyNA(loo$idx) || anyDuplicated(loo$idx) ||
      any(loo$idx < 1L) || any(loo$idx > length(y_method))) {
    stop(method$id, " returned inconsistent LOO indices or dimensions.")
  }
  if (any(!is.finite(y_method[loo$idx])) ||
      any(!is.finite(loo$mu)) ||
      any(!is.finite(loo$var)) || any(loo$var <= 0)) {
    stop(method$id, " returned invalid LOO observations or predictions.")
  }
  if (length(beta_estimate$beta) != ncol(X_method) ||
      length(loo$beta_hat) != ncol(X_method)) {
    stop(method$id, " returned the wrong number of fixed effects.")
  }
  beta_difference <- max(abs(beta_estimate$beta - loo$beta_hat))
  if (!is.finite(beta_difference) ||
      beta_difference > beta_consistency_tolerance) {
    stop(method$id, " profile and LOO fixed-effect estimates disagree.")
  }

  list(loo = loo, timing = loo_timing)
}

results <- list()

for (method_id in names(methods)) {
  method <- methods[[method_id]]
  # precompute_likelihood_inputs(), recover_fixed_effects(), and
  # compute_plugin_loo() all branch on method$model == "directional" and
  # otherwise assume "wm1"; validate the full known domain here, once, so a
  # future model (e.g. the deferred WM2 variant noted above) can't silently
  # fall through the "wm1" branch instead of erroring.
  if (!method$model %in% c("directional", "wm1")) {
    stop(method$id, " has an unrecognized model: ", method$model)
  }
  method_started <- Sys.time()
  cat("\n", method$label, "\n", sep = "")

  if (method$reversed) {
    graph <- graph_reversed
    X_method <- X_reversed
    y_method <- y_reversed
    observation_ids_method <- observation_ids_reversed
  } else {
    graph <- graph_original
    X_method <- X
    y_method <- y_full
    observation_ids_method <- observation_ids
  }
  if (method$model == "directional") {
    # K1/K2 meaning: see the directional_rule comment above the methods list.
    weight_function <- switch(
      method$directional_rule,
      K1 = NULL,
      K2 = .columbia_k2_weights,
      stop(
        method$id, " has an unrecognized directional_rule: ",
        method$directional_rule
      )
    )
    if (is.null(weight_function)) {
      graph$setDirectionalWeightFunction()
    } else {
      graph$setDirectionalWeightFunction(f_in = weight_function)
    }
  }

  cat("  precomputing likelihood inputs...\n")
  likelihood_inputs <- precompute_likelihood_inputs(method, graph, X_method)
  precomputed <- likelihood_inputs$precomputed
  precompute_timing <- likelihood_inputs$timing

  cat("  optimizing three covariance parameters from independent starts...\n")
  fit <- optimize_covariance_parameters(
    method, likelihood_inputs$profile_likelihood,
    likelihood_inputs$profile_arguments
  )
  optimization <- fit$optimization
  optimization_diagnostics <- fit$optimization_diagnostics
  optimization_timing <- fit$timing

  verification <- verify_local_optimum(
    method, optimization, fit$optimization_objective
  )

  fixed_effects_fit <- recover_fixed_effects(
    method, graph, precomputed, optimization
  )
  beta_estimate <- fixed_effects_fit$beta_estimate
  fixed_effect_timing <- fixed_effects_fit$timing

  cat("  computing plug-in LOO predictions...\n")
  loo_fit <- compute_plugin_loo(
    method, graph, X_method, y_method, optimization, precomputed,
    beta_estimate
  )
  loo <- loo_fit$loo
  loo_timing <- loo_fit$timing

  loo_result <- data.frame(
    observation_id = observation_ids_method[loo$idx],
    observation_row = loo$idx,
    observed = y_method[loo$idx],
    mean = loo$mu,
    variance = loo$var,
    row.names = NULL
  )

  fit_total <- sum(c(
    unname(precompute_timing[["elapsed"]]),
    unname(optimization_timing[["elapsed"]]),
    unname(fixed_effect_timing[["elapsed"]])
  ))
  timings <- c(
    precompute = unname(precompute_timing[["elapsed"]]),
    parameter_optimization = unname(optimization_timing[["elapsed"]]),
    fixed_effect_recovery = unname(fixed_effect_timing[["elapsed"]]),
    fit_total = fit_total,
    loo = unname(loo_timing[["elapsed"]]),
    overall = as.numeric(difftime(
      Sys.time(),
      method_started,
      units = "secs"
    ))
  )

  theta <- stats::setNames(optimization$par, theta_names)
  parameters <- c(
    sigma_e = exp(theta[["log_sigma_e"]]),
    tau = exp(-theta[["log_reciprocal_tau"]]),
    kappa = exp(theta[["log_kappa"]]),
    # Practical range for alpha=1 (nu=0.5) Whittle-Matern fields:
    # sqrt(8*nu)/kappa, the same convention used elsewhere in the package
    # (e.g. kappa <- sqrt(8*0.5)/range in R/util.R).
    range = sqrt(8 * 0.5) / exp(theta[["log_kappa"]])
  )
  log_likelihood <- -optimization$objective
  fixed_effects <- stats::setNames(
    as.vector(beta_estimate$beta),
    colnames(X_method)
  )
  optimization_record <- list(
    selected_start = optimization_diagnostics$start_id[fit$best_index],
    theta = theta,
    objective = optimization$objective,
    convergence = optimization$convergence,
    message = optimization$message,
    iterations = optimization$iterations,
    evaluations = optimization$evaluations,
    at_lower_bound = stats::setNames(
      optimization$at_lower_bound,
      names(theta)
    ),
    at_upper_bound = stats::setNames(
      optimization$at_upper_bound,
      names(theta)
    ),
    numerical_gradient = verification$numerical_gradient,
    max_abs_gradient = verification$max_abs_gradient,
    hessian = verification$objective_hessian,
    hessian_eigenvalues = verification$hessian_eigenvalues,
    reproduced_by = optimization_diagnostics$start_id[fit$reproducing_indices],
    reproducing_theta_spread = stats::setNames(
      fit$theta_spread, names(theta)
    ),
    rejected_profile_evaluations = fit$rejected_profile_evaluations,
    last_profile_error = fit$last_profile_error,
    starts = optimization_diagnostics
  )

  results[[method_id]] <- list(
    theta = theta,
    parameters = parameters,
    fixed_effects = fixed_effects,
    log_likelihood = log_likelihood,
    optimization = optimization_record,
    loo = loo_result,
    timings = timings
  )

  cat(sprintf(
    paste0(
      "  fit %.2fs, LOO %.2fs, log-likelihood %.6f, ",
      "convergence %d, max |gradient| %.3g\n"
    ),
    timings[["fit_total"]],
    timings[["loo"]],
    log_likelihood,
    optimization$convergence,
    verification$max_abs_gradient
  ))

  # Free this method's graph/precomputed/Hessian objects before fitting the
  # next one; peak memory on the full Columbia component is otherwise high
  # enough that four methods' worth can accumulate across iterations.
  rm(
    graph, X_method, y_method, observation_ids_method,
    likelihood_inputs, precomputed, fit, optimization,
    optimization_diagnostics, verification, fixed_effects_fit,
    beta_estimate, loo_fit, loo, loo_result, parameters, theta,
    fixed_effects, log_likelihood, optimization_record
  )
  invisible(gc(verbose = FALSE))
}


## Summary and saving --------------------------------------------------------

# All Settings-block tunables that affect fitting/acceptance, captured by
# name (not hand-copied) so this can't drift from the constants above it.
fitting_settings <- mget(c(
  "initial_thetas", "lower", "upper", "max_iterations",
  "eval_max_per_iteration", "minimum_reproducing_starts",
  "objective_tolerance", "theta_tolerance", "gradient_tolerance",
  "bound_tolerance", "finite_difference_relative_step",
  "beta_consistency_tolerance", "rejected_evaluation_penalty",
  "valid_objective_ceiling"
))
attr(results, "audit") <- c(
  list(
    created_at = Sys.time(),
    package = "MetricGraph",
    package_version = as.character(utils::packageVersion("MetricGraph")),
    package_namespace_path = getNamespaceInfo(asNamespace("MetricGraph"), "path"),
    component_fingerprint = columbia_component$fingerprint,
    component_summary = columbia_component$summary,
    component_source = "https://doi.org/10.6084/m9.figshare.24132840.v1",
    component_license = "CC BY 4.0",
    response = response_name,
    covariates = covariates
  ),
  fitting_settings
)

if (save_result) {
  columbia_save_checkpoint(results, out_file)
}

parameter_table <- data.frame()
score_table <- data.frame()
timing_table <- data.frame()

for (method_id in names(results)) {
  method_result <- results[[method_id]]
  loo_result <- method_result$loo
  standard_deviation <- sqrt(loo_result$variance)
  residual <- loo_result$observed - loo_result$mean

  parameter_table <- rbind(
    parameter_table,
    data.frame(
      method = method_id,
      sigma_e = method_result$parameters[["sigma_e"]],
      tau = method_result$parameters[["tau"]],
      kappa = method_result$parameters[["kappa"]],
      range = method_result$parameters[["range"]],
      log_likelihood = method_result$log_likelihood,
      convergence = method_result$optimization$convergence,
      reproduced_starts = length(method_result$optimization$reproduced_by),
      max_abs_gradient = method_result$optimization$max_abs_gradient,
      row.names = NULL
    )
  )
  score_table <- rbind(
    score_table,
    data.frame(
      method = method_id,
      logscore = -mean(MetricGraph:::LS(
        loo_result$observed,
        loo_result$mean,
        standard_deviation
      )),
      crps = -mean(MetricGraph:::CRPS(
        loo_result$observed,
        loo_result$mean,
        standard_deviation
      )),
      scrps = -mean(MetricGraph:::SCRPS(
        loo_result$observed,
        loo_result$mean,
        standard_deviation
      )),
      mae = mean(abs(residual)),
      rmse = sqrt(mean(residual^2)),
      row.names = NULL
    )
  )
  timing_table <- rbind(
    timing_table,
    data.frame(
      method = method_id,
      precompute = method_result$timings[["precompute"]],
      parameter_optimization =
        method_result$timings[["parameter_optimization"]],
      fixed_effect_recovery =
        method_result$timings[["fixed_effect_recovery"]],
      fit_total = method_result$timings[["fit_total"]],
      loo = method_result$timings[["loo"]],
      overall = method_result$timings[["overall"]],
      row.names = NULL
    )
  )
}

cat("\nFitted parameters and log-likelihoods (higher is better)\n")
print(parameter_table, digits = 6, row.names = FALSE)
cat("\nPlug-in LOO scores (lower is better)\n")
print(score_table, digits = 6, row.names = FALSE)
cat("\nTimings in seconds\n")
print(timing_table, digits = 6, row.names = FALSE)
cat(
  "\nShared component, graph, and design matrix timings:",
  unname(component_timing[["elapsed"]]),
  unname(graph_timing[["elapsed"]]),
  unname(design_timing[["elapsed"]]),
  "\n"
)
cat("\nCompleted:", paste(names(results), collapse = ", "), "\n")

if (save_result) {
  cat("Saved:", out_file, "\n")
} else {
  cat(
    "Result was not saved; set out_file in the Settings block to persist it.\n"
  )
}
