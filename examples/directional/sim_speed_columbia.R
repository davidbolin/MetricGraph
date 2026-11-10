# Time one precompute and likelihood evaluation as observation count grows on
# the packaged Mid-Columbia component.
#
# Run from the package root: Rscript examples/directional/sim_speed_columbia.R
# Inputs: none; edit the Settings block below to change the run.
# Outputs: a timing table on stdout, and if output_file is set below, an .rds
# of the results plus "<stem>_precompute.pdf"/"<stem>_evaluate.pdf".

## Settings ------------------------------------------------------------------

# Arbitrary round-number defaults spanning two orders of magnitude; no
# benchmarking backs this specific grid.
n_obs_values <- c(
  100L, 500L, 1000L, 2000L, 4000L,
  8000L, 10000L, 12000L, 15000L, 20000L
)
# Balances timing stability (more repetitions) against wall time.
n_repetitions <- 5L
available_cores <- max(1L, min(n_repetitions,parallel::detectCores(logical = FALSE)))
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
# Per-method dense-covariance caps. Continuity uses a lower benchmark cap
# because its covariance evaluation requires more working memory.
maximum_covariance_n <- c(K1 = 12000L, K2 = 12000L, continuity = 8000L)
# Raise the package's default covariance guard to the largest configured cap;
# forked workers inherit the option. Per-method caps below still skip
# continuity above its lower limit.
options(DIRECTIONAL_OU_MAX_POINTS = max(maximum_covariance_n))
seed <- 1L
# Keep observations off edge endpoints (0/1 coincide with vertices).
edge_margin <- 0.01
# Show the crowded low-observation timings in an inset whenever the main plot
# also contains observation counts above this threshold.
inset_max_n <- 1000L


## Load the code under test --------------------------------------------------

pkgload::load_all(".", reset = TRUE, export_all = TRUE)
source("examples/directional/columbia_full_graph_helpers.R")


## Data and graphs -----------------------------------------------------------

component <- MetricGraph::columbia_main_component
graph_original <- columbia_make_graph(component, reversed = FALSE)
graph_reversed <- columbia_make_graph(component, reversed = TRUE)

# Remove compatible degree-2 vertices before any benchmark timings. The
# default check_weights = TRUE preserves directional-weight boundaries.
graph_original$prune_vertices(verbose = 1)
graph_reversed$prune_vertices(verbose = 1)


## Methods under test --------------------------------------------------------

# K1/K2: profile-likelihood evaluation with the K1/K2 directional weights on
# the original (downstream-flowing) graph. continuity: K1 weights on the
# reversed graph, i.e. timing the reversed-direction continuity condition.
# Each family also gets a _covariance counterpart built from the same
# (graph, weights) pair, timing direct dense-covariance evaluation instead of
# the profile likelihood, with its observation cap taken from
# maximum_covariance_n.
base_configs <- list(
  K1         = list(graph = "original", weights = "K1"),
  K2         = list(graph = "original", weights = "K2"),
  continuity = list(graph = "reversed", weights = "K1")
)
methods <- c(
  lapply(base_configs, function(config) {
    c(config, use_dense_covariance = FALSE, maximum_n = Inf)
  }),
  stats::setNames(
    lapply(names(base_configs), function(method_name) {
      config <- base_configs[[method_name]]
      c(
        config,
        use_dense_covariance = TRUE,
        maximum_n = maximum_covariance_n[[method_name]]
      )
    }),
    paste0(names(base_configs), "_covariance")
  )
)
# Plotmath labels give K1/K2 proper subscripts and shorten the dense-
# covariance method suffix to "cov" in the legend.
method_labels <- c(
  K1 = expression(K[1]),
  K1_covariance = expression(K[1] ~ plain(cov)),
  K2 = expression(K[2]),
  K2_covariance = expression(K[2] ~ plain(cov)),
  continuity = expression(plain(continuity)),
  continuity_covariance = expression(plain(continuity) ~ plain(cov))
)
if (!setequal(names(method_labels), names(methods))) {
  stop("method_labels must contain exactly one label for every method")
}

# Use one color for each method family. Dashed lines and open points distinguish
# the dense-covariance variants from their corresponding profile methods.
method_palette <- c(
  K1 = "black",
  K2 = "#0072B2",
  continuity = "#009E73"
)
method_families <- sub("_covariance$", "", names(method_labels))
unknown_families <- setdiff(unique(method_families), names(method_palette))
if (length(unknown_families) > 0L) {
  stop("method_palette is missing: ", paste(unknown_families, collapse = ", "))
}
method_colors <- stats::setNames(
  method_palette[method_families], names(method_labels)
)
covariance_methods <- grepl("_covariance$", names(method_labels))
method_linetypes <- stats::setNames(
  ifelse(covariance_methods, "dashed", "solid"), names(method_labels)
)
method_shapes <- stats::setNames(
  ifelse(covariance_methods, 1, 16), names(method_labels)
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

# Also reads the method plotting controls above and timing_results (built later,
# in the Timing sweep section below). All are in scope when this function is
# called in Output.
plot_stage <- function(stage, file) {
  # A log-scale axis can't show a zero elapsed time.
  log_scale_floor <- .Machine$double.eps
  stage_results <- timing_results[timing_results$stage == stage, ]
  stage_results$elapsed_seconds <-
    pmax(stage_results$elapsed_seconds, log_scale_floor)

  plot <- ggplot2::ggplot(
    stage_results,
    ggplot2::aes(
      n_obs, elapsed_seconds,
      color = method, linetype = method, shape = method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 2.2, stroke = 0.8) +
    ggplot2::scale_x_continuous(
      labels = function(x) format(x, big.mark = ",", scientific = FALSE,
                                  trim = TRUE)
    ) +
    ggplot2::scale_y_log10(
      labels = function(x) format(x, scientific = FALSE, trim = TRUE)
    ) +
    ggplot2::scale_color_manual(
      values = method_colors,
      breaks = names(method_labels),
      labels = method_labels
    ) +
    ggplot2::scale_linetype_manual(
      values = method_linetypes,
      breaks = names(method_labels),
      labels = method_labels
    ) +
    ggplot2::scale_shape_manual(
      values = method_shapes,
      breaks = names(method_labels),
      labels = method_labels
    ) +
    ggplot2::labs(
      x = "# observations",
      y = "Median elapsed time (seconds)",
      color = "Method",
      linetype = "Method",
      shape = "Method"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  low_n <- stage_results$n_obs <= inset_max_n
  add_inset <- any(stage_results$n_obs > inset_max_n) &&
    length(unique(stage_results$n_obs[low_n])) >= 2L
  if (add_inset) {
    inset_results <- stage_results[low_n, ]
    inset_plot <- ggplot2::ggplot(
      inset_results,
      ggplot2::aes(
        n_obs, elapsed_seconds,
        color = method, linetype = method, shape = method
      )
    ) +
      ggplot2::geom_line(linewidth = 0.55) +
      ggplot2::geom_point(size = 1.7, stroke = 0.7) +
      ggplot2::scale_x_continuous(
        limits = c(0, inset_max_n),
        breaks = c(0, inset_max_n / 2, inset_max_n),
        labels = function(x) format(x, big.mark = ",", scientific = FALSE,
                                    trim = TRUE)
      ) +
      ggplot2::scale_y_log10(
        labels = function(x) format(x, scientific = FALSE, trim = TRUE)
      ) +
      ggplot2::scale_color_manual(values = method_colors) +
      ggplot2::scale_linetype_manual(values = method_linetypes) +
      ggplot2::scale_shape_manual(values = method_shapes) +
      ggplot2::labs(
        title = paste0(
          "# observations: 0-",
          format(inset_max_n, big.mark = ",", scientific = FALSE)
        ),
        x = NULL,
        y = NULL
      ) +
      ggplot2::theme_bw(base_size = 8) +
      ggplot2::theme(
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 8, hjust = 0.5),
        plot.margin = ggplot2::margin(2, 3, 2, 2)
      )

    x_limits <- range(stage_results$n_obs, finite = TRUE)
    log_y_limits <- range(log10(stage_results$elapsed_seconds), finite = TRUE)
    x_span <- diff(x_limits)
    log_y_span <- diff(log_y_limits)
    if (x_span > 0 && log_y_span > 0) {
      plot <- plot + ggplot2::annotation_custom(
        grob = ggplot2::ggplotGrob(inset_plot),
        xmin = x_limits[1L] + 0.55 * x_span,
        xmax = x_limits[1L] + 0.98 * x_span,
        ymin = 10^(log_y_limits[1L] + 0.05 * log_y_span),
        ymax = 10^(log_y_limits[1L] + 0.43 * log_y_span)
      )
    }
  }

  ggplot2::ggsave(
    file, plot,
    width = 8, height = 5.5, units = "in",
    device = grDevices::cairo_pdf
  )
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
  # Derived from the actual names so the log columns stay aligned regardless
  # of which methods are configured above.
  method_name_width <- max(nchar(names(methods)))
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
      sprintf("  %-*s precompute %.4fs, evaluate %.4fs",
              method_name_width, method_name,
              timing[["precompute"]], timing[["evaluate"]])
    )
  }
}

timing_results <- do.call(rbind, results)
print(timing_results, row.names = FALSE)


## Output --------------------------------------------------------------------

if (nzchar(output_file)) {

  attr(timing_results, "run_config") <- list(
    n_obs_values = n_obs_values, n_repetitions = n_repetitions,
    n_cores = n_cores, seed = seed, theta_test = theta_test,
    DIRECTIONAL_OU_MAX_POINTS = getOption("DIRECTIONAL_OU_MAX_POINTS")
  )
  attr(timing_results, "package_version") <-
    as.character(utils::packageVersion("MetricGraph"))
  saveRDS(timing_results, output_file)
  figure_stem <- tools::file_path_sans_ext(output_file)
  plot_stage("precompute", paste0(figure_stem, "_precompute.pdf"))
  plot_stage("evaluate", paste0(figure_stem, "_evaluate.pdf"))
}
