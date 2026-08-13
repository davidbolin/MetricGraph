# Time one precompute and likelihood evaluation as observation count grows on
# the packaged Mid-Columbia component.
#
# Run from the package root: Rscript examples/directional/sim_speed_columbia.R
# Inputs: none; edit the Settings block below to change the run.
# Outputs: a timing table on stdout, and if output_file is set below, an .rds
# of the results plus "<stem>_precompute.png"/"<stem>_evaluate.png".


## Settings ------------------------------------------------------------------

# Arbitrary round-number defaults spanning two orders of magnitude; no
# benchmarking backs this specific grid.
n_obs_values <- c(100L, 500L, 1000L, 5000L)
# Balances timing stability (more repetitions) against wall time.
n_repetitions <- 5L
available_cores <- max(1L, parallel::detectCores(logical = FALSE))
# No point using more cores than repetitions; leave one core free for the OS.
n_cores <- min(n_repetitions, max(1L, available_cores - 1L))
# mclapply forks; forking is unavailable on Windows, so run serially.
if (.Platform$OS.type == "windows") {
  n_cores <- 1L
}

# Set to a file path (e.g. "out/sim_speed.rds") to save results and plots;
# "" runs the sweep and prints the table without writing anything.
output_file <- "examples/directional/results/sim_speed.rds"
# (log sigma_e, log reciprocal_tau, log kappa); mid-range values used only to
# fix a representative evaluation point for timing, not fitted.
theta_test <- c(log(0.5), log(1), log(0.5))
# Matches DIRECTIONAL_OU_MAX_POINTS in R/covariance_directional_ou.R -- keep
# in sync; dense n x n covariance is not viable beyond this.
maximum_covariance_n <- 10000L
seed <- 1L
# Keep observations off edge endpoints (0/1 coincide with vertices).
edge_margin <- 0.01


## Load the code under test --------------------------------------------------

pkgload::load_all(".", reset = TRUE, export_all = TRUE)
source("examples/directional/columbia_full_graph_helpers.R")


## Data and graphs -----------------------------------------------------------

component <- MetricGraph::columbia_main_component
graph_original <- columbia_make_graph(component, reversed = FALSE)
graph_reversed <- columbia_make_graph(component, reversed = TRUE)


## Methods under test --------------------------------------------------------

# K1/K2: profile-likelihood evaluation with the K1/K2 directional weights on
# the original (downstream-flowing) graph. continuity: K1 weights on the
# reversed graph, i.e. timing the reversed-direction continuity condition.
# K2_covariance: direct dense-covariance evaluation instead of the profile
# likelihood, capped at maximum_covariance_n observations.
methods <- list(
  K1 = list(
    graph = "original", weights = "K1",
    use_dense_covariance = FALSE, maximum_n = Inf
  ),
  K2 = list(
    graph = "original", weights = "K2",
    use_dense_covariance = FALSE, maximum_n = Inf
  ),
  continuity = list(
    graph = "reversed", weights = "K1",
    use_dense_covariance = FALSE, maximum_n = Inf
  ),
  K2_covariance = list(
    graph = "original", weights = "K2",
    use_dense_covariance = TRUE, maximum_n = maximum_covariance_n
  )
)
# Plot color per method, derived from `methods` so a renamed or added method
# can't silently drop off the plots in plot_stage() below.
method_palette <- c("black", "steelblue", "darkgreen", "orange")
if (length(method_palette) < length(methods)) {
  stop(
    "method_palette needs at least as many colors as methods; got ",
    length(method_palette), " colors for ", length(methods), " methods"
  )
}
method_colors <- stats::setNames(
  method_palette[seq_along(methods)], names(methods)
)


## Timing and plotting helpers -----------------------------------------------

run_repetitions <- function(indices, task, cores) {
  if (cores > 1L) {
    parallel::mclapply(indices, task, mc.cores = cores, mc.preschedule = FALSE)
  } else {
    lapply(indices, task)
  }
}

# Also reads n_cores and theta_test from the Settings block.
time_method <- function(graph, method, inputs) {
  # C is rebuilt from the new weight function on the next constraint call;
  # the "constraint matrix deleted" warning is expected every iteration here
  # and is not a sign of a stale graph.
  suppressWarnings(
    graph$setDirectionalWeightFunction(
      f_in = if (identical(method$weights, "K2")) {
        .columbia_k2_weights
      } else {
        NULL
      }
    )
  )

  runs <- run_repetitions(seq_along(inputs), function(repetition) {
    input <- inputs[[repetition]]
    graph$clear_observations()
    graph$add_observations(
      data = data.frame(
        y = input$y,
        edge_number = input$edge_number,
        distance_on_edge = input$distance_on_edge
      ),
      normalized = TRUE,
      verbose = 0
    )

    if (method$use_dense_covariance) {
      precompute_seconds <- system.time({
        precomputed <- MetricGraph:::precompute_directional_ou_covariance(
          graph, data_name = "y"
        )
      })[["elapsed"]]
      evaluate_seconds <- system.time({
        loglik_value <- MetricGraph:::directional_ou_covariance_loglik_precompute(
          theta_test, precomputed, cpp = TRUE, force_dense = FALSE
        )
      })[["elapsed"]]
    } else {
      precompute_seconds <- system.time({
        precomputed <- MetricGraph:::precompute_alpha1_directional(
          graph, data_name = "y"
        )
      })[["elapsed"]]
      evaluate_seconds <- system.time({
        loglik_value <-
          MetricGraph:::likelihood_alpha1_directional_profile_precompute(
            theta_test, precomputed, parameterization = "spde", cpp = TRUE
          )
      })[["elapsed"]]
    }

    c(
      precompute = precompute_seconds,
      evaluate = evaluate_seconds,
      loglik = loglik_value
    )
  }, cores = n_cores)

  # Only mclapply() repetitions (cores > 1) can produce "try-error" objects;
  # a failing lapply() repetition (cores == 1, e.g. always on Windows) raises
  # immediately instead and never reaches this check.
  failures <- vapply(runs, inherits, logical(1), what = "try-error")
  if (any(failures)) {
    stop("Timing repetition failed: ", runs[[which(failures)[1L]]])
  }
  loglik_values <- vapply(runs, `[[`, numeric(1), "loglik")
  if (!all(is.finite(loglik_values))) {
    stop(
      "Timing repetitions returned non-finite likelihoods: ",
      paste(loglik_values, collapse = ", ")
    )
  }

  c(
    precompute = median(vapply(runs, `[[`, numeric(1), "precompute")),
    evaluate = median(vapply(runs, `[[`, numeric(1), "evaluate"))
  )
}

# Also reads method_colors (built in the Methods-under-test section above)
# and timing_results (built later, in the Timing sweep section below) --
# both are in scope by the time this function is actually called, in Output.
plot_stage <- function(stage, file) {
  # A log-scale axis can't show a zero elapsed time.
  log_scale_floor <- .Machine$double.eps
  stage_results <- timing_results[timing_results$stage == stage, ]
  stage_results$elapsed_seconds <-
    pmax(stage_results$elapsed_seconds, log_scale_floor)

  plot <- ggplot2::ggplot(
    stage_results,
    ggplot2::aes(n_obs, elapsed_seconds, color = method)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::labs(x = "observations", y = "median elapsed seconds")

  ggplot2::ggsave(file, plot, width = 800, height = 550, units = "px", dpi = 96)
}


## Timing sweep --------------------------------------------------------------

set.seed(seed)
results <- list()
result_index <- 0L

message(sprintf(
  "n_obs = %s, repetitions = %d, cores = %d, seed = %d, MetricGraph %s",
  paste(n_obs_values, collapse = ","), n_repetitions, n_cores, seed,
  utils::packageVersion("MetricGraph")
))

for (n_obs in n_obs_values) {
  inputs <- lapply(seq_len(n_repetitions), function(repetition) {
    list(
      edge_number = sample.int(graph_original$nE, n_obs,
                               replace = TRUE),
      distance_on_edge = runif(
        n_obs, edge_margin, 1 - edge_margin
      ),
      y = rnorm(n_obs)
    )
  })

  message("n_obs = ", n_obs)
  for (method_name in names(methods)) {
    method <- methods[[method_name]]
    if (n_obs > method$maximum_n) {
      message("  ", method_name, ": skipped above ", method$maximum_n)
      next
    }
    graph <- if (identical(method$graph, "original")) {
      graph_original
    } else {
      graph_reversed
    }
    timing <- time_method(graph, method, inputs)
    result_index <- result_index + 1L
    results[[result_index]] <- data.frame(
      n_obs = n_obs,
      method = method_name,
      stage = names(timing),
      elapsed_seconds = as.numeric(timing),
      row.names = NULL
    )
    message(
      sprintf("  %-13s precompute %.4fs, evaluate %.4fs",
              method_name, timing[["precompute"]], timing[["evaluate"]])
    )
  }
}

timing_results <- do.call(rbind, results)
print(timing_results, row.names = FALSE)


## Output --------------------------------------------------------------------

if (nzchar(output_file)) {

  attr(timing_results, "run_config") <- list(
    n_obs_values = n_obs_values, n_repetitions = n_repetitions,
    n_cores = n_cores, seed = seed, theta_test = theta_test
  )
  attr(timing_results, "package_version") <-
    as.character(utils::packageVersion("MetricGraph"))
  saveRDS(timing_results, output_file)
  figure_stem <- tools::file_path_sans_ext(output_file)
  plot_stage("precompute", paste0(figure_stem, "_precompute.png"))
  plot_stage("evaluate", paste0(figure_stem, "_evaluate.png"))
}
