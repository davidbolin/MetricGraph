# Visual guide to the small graphs used by the directional covariance tests.
#
# Interactive use from the package root (one graph per plot-history page):
#   source("inst/dev/plot_directional_covariance_test_graphs.R")
#   result <- plot_directional_covariance_test_graphs()
#
# Command-line use from the package root writes a labelled, multi-page PDF:
#   Rscript inst/dev/plot_directional_covariance_test_graphs.R
#
# This script does not load or depend on testthat.

plot_directional_covariance_test_graphs <- function(
    package_root = ".", case_names = NULL, show_observations = TRUE,
    direction = TRUE, output_file = NULL) {
  package_root <- normalizePath(package_root, mustWork = TRUE)
  helper_path <- file.path(
    package_root, "tests", "testthat", "helper-subgraph-equivalence.R"
  )
  if (!file.exists(helper_path)) {
    stop(
      "Could not find the directional covariance fixture helper. ",
      "Run from the MetricGraph package root or set package_root explicitly."
    )
  }
  if (!requireNamespace("MetricGraph", quietly = TRUE)) {
    stop(
      "MetricGraph must be installed before running this script; ",
      "from the package root, run `R CMD INSTALL .` first."
    )
  }

  fixture_environment <- new.env(parent = asNamespace("MetricGraph"))
  sys.source(helper_path, envir = fixture_environment)
  cases <- fixture_environment$subgraph_equivalence_cases()
  available_case_names <- names(cases)
  if (is.null(case_names)) {
    case_names <- available_case_names
  } else {
    unknown_case_names <- setdiff(case_names, available_case_names)
    if (length(unknown_case_names) > 0L) {
      stop(
        "Unknown case_names: ", paste(unknown_case_names, collapse = ", "),
        ". Available cases: ", paste(available_case_names, collapse = ", "),
        "."
      )
    }
    cases <- cases[case_names]
  }
  plots <- Map(function(graph, case_name) {
    graph$plot(
      data = if (show_observations) "temp" else NULL,
      direction = direction
    ) +
      ggplot2::ggtitle(case_name) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }, cases, case_names)
  names(plots) <- case_names

  if (!is.null(output_file)) {
    output_file <- normalizePath(output_file, mustWork = FALSE)
    grDevices::pdf(output_file, width = 8, height = 6, onefile = TRUE)
    output_device <- grDevices::dev.cur()
    on.exit({
      if (!is.null(output_device) &&
          output_device %in% grDevices::dev.list()) {
        grDevices::dev.off(output_device)
      }
    }, add = TRUE)
  } else {
    output_device <- NULL
  }

  for (plot in plots) {
    print(plot)
  }
  if (!is.null(output_device)) {
    grDevices::dev.off(output_device)
    output_device <- NULL
  }

  invisible(list(cases = cases, plots = plots, output_file = output_file))
}

script_arguments <- grep(
  "^--file=", commandArgs(trailingOnly = FALSE), value = TRUE
)
running_directly <- length(script_arguments) == 1L &&
  identical(
    basename(sub("^--file=", "", script_arguments)),
    "plot_directional_covariance_test_graphs.R"
  )
if (running_directly) {
  output_file <- file.path(
    getwd(), "directional_covariance_test_graphs.pdf"
  )
  plot_directional_covariance_test_graphs(output_file = output_file)
  message("Wrote ", output_file)
}
rm(script_arguments, running_directly)
