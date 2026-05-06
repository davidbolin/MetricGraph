## benchmark_helpers.R
## Shared helpers for the MetricGraph benchmark suite.

if (!requireNamespace("bench", quietly = TRUE)) {
  stop("The 'bench' package is required. Install with: install.packages('bench')")
}

library(MetricGraph)
library(bench)
library(sf)

# ---------------------------------------------------------------------------
# Graph generators
# ---------------------------------------------------------------------------

#' Build a planar n x n square grid as a list of edge matrices (2x2).
#' Produces n^2 vertices and 2*n*(n-1) edges.
make_grid_edges <- function(n) {
  xs <- ys <- seq_len(n)
  out <- vector("list", 2L * n * (n - 1L))
  k <- 1L
  for (j in seq_len(n)) {
    for (i in seq_len(n - 1L)) {
      out[[k]] <- rbind(c(xs[i], ys[j]), c(xs[i + 1L], ys[j]))
      k <- k + 1L
    }
  }
  for (i in seq_len(n)) {
    for (j in seq_len(n - 1L)) {
      out[[k]] <- rbind(c(xs[i], ys[j]), c(xs[i], ys[j + 1L]))
      k <- k + 1L
    }
  }
  out
}

#' Build a lon/lat grid as an sf object with CRS 4326,
#' mimicking graphs imported from osmdata.
#' Uses slightly perturbed coordinates for realism.
make_longlat_grid <- function(n, lon0 = 2.34, lon1 = 2.36,
                              lat0 = 48.85, lat1 = 48.87) {
  set.seed(42)
  lons <- seq(lon0, lon1, length.out = n)
  lats <- seq(lat0, lat1, length.out = n)
  # Add small perturbation to interior vertices to mimic real OSM data
  perturb <- (lon1 - lon0) / (n - 1) * 0.05
  edges <- list()
  k <- 1L
  for (j in seq_len(n)) {
    for (i in seq_len(n - 1L)) {
      p1 <- c(lons[i], lats[j])
      p2 <- c(lons[i + 1L], lats[j])
      if (i > 1 && i < n && j > 1 && j < n) {
        p1 <- p1 + runif(2, -perturb, perturb)
      }
      edges[[k]] <- rbind(p1, p2)
      k <- k + 1L
    }
  }
  for (i in seq_len(n)) {
    for (j in seq_len(n - 1L)) {
      p1 <- c(lons[i], lats[j])
      p2 <- c(lons[i], lats[j + 1L])
      edges[[k]] <- rbind(p1, p2)
      k <- k + 1L
    }
  }
  ls_list <- lapply(edges, function(e) st_linestring(unname(e)))
  st_sf(geometry = st_sfc(ls_list, crs = 4326))
}

#' Generate random observations on a graph.
#' Returns a data.frame with columns: y, edge_number, distance_on_edge.
generate_observations <- function(graph, n_obs_per_edge, sigma = 1.3,
                                  range = 0.2, sigma_e = 0.1, alpha = 1) {
  PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
    cbind(rep(i, n_obs_per_edge), runif(n_obs_per_edge))))
  u <- sample_spde(range = range, sigma = sigma, alpha = alpha,
                   graph = graph, PtE = PtE)
  y <- u + sigma_e * rnorm(n_obs_per_edge * graph$nE)
  data.frame(y = y, edge_number = PtE[, 1], distance_on_edge = PtE[, 2])
}

# ---------------------------------------------------------------------------
# Timing utilities
# ---------------------------------------------------------------------------

#' Run a benchmark by timing an expression n_rep times.
#' The expression is evaluated in the caller's environment.
#' Returns a list with median, min, max, n_itr, mem_alloc, timings.
run_bench <- function(expr, n_rep = 3, env = parent.frame()) {
  expr_q <- substitute(expr)

  # Memory measurement using bench for the first run
  mem <- tryCatch(
    bench::bench_process_memory(),
    error = function(e) NULL
  )

  timings <- numeric(n_rep)
  for (i in seq_len(n_rep)) {
    gc(FALSE)
    t0 <- proc.time()
    eval(expr_q, envir = env)
    timings[i] <- (proc.time() - t0)[["elapsed"]]
  }

  mem_after <- tryCatch(
    bench::bench_process_memory(),
    error = function(e) NULL
  )

  mem_diff <- if (!is.null(mem) && !is.null(mem_after)) {
    bench::as_bench_bytes(max(0, as.numeric(mem_after["current"]) -
                                  as.numeric(mem["current"])))
  } else {
    bench::as_bench_bytes(0)
  }

  list(
    median = median(timings),
    min = min(timings),
    max = max(timings),
    n_itr = n_rep,
    mem_alloc = mem_diff,
    timings = timings
  )
}

# ---------------------------------------------------------------------------
# Output formatting
# ---------------------------------------------------------------------------

#' Format a bench_time or numeric seconds value as a human-readable string.
#' Shows ms for sub-second, s for seconds, min for minutes.
format_time <- function(x) {
  secs <- as.numeric(x)
  if (secs < 0.0001) {
    sprintf("%.1f us", secs * 1e6)
  } else if (secs < 1) {
    sprintf("%.1f ms", secs * 1000)
  } else if (secs < 60) {
    sprintf("%.2f s", secs)
  } else {
    sprintf("%.1f min", secs / 60)
  }
}

#' Format bench results as a readable table.
#' Takes a list of named bench results and combines them.
format_bench_results <- function(results_list) {
  df <- do.call(rbind, lapply(names(results_list), function(nm) {
    r <- results_list[[nm]]
    data.frame(
      scenario = nm,
      median = format_time(r$median),
      min = format_time(r$min),
      mem_alloc = format(r$mem_alloc),
      n_itr = as.integer(r$n_itr),
      stringsAsFactors = FALSE
    )
  }))
  rownames(df) <- NULL
  df
}

#' Print a section header for benchmark output.
print_section <- function(title) {
  width <- 70
  cat("\n", strrep("=", width), "\n", sep = "")
  cat(" ", title, "\n", sep = "")
  cat(strrep("=", width), "\n\n", sep = "")
}

#' Print a formatted results data.frame.
print_results <- function(df) {
  # Use format for alignment
  out <- capture.output(print(df, right = FALSE, row.names = FALSE))
  cat(paste(out, collapse = "\n"), "\n\n")
}

#' Print system info header.
print_system_info <- function() {
  cat(strrep("=", 70), "\n")
  cat(" MetricGraph Benchmark Suite\n")
  cat(strrep("=", 70), "\n")
  cat(sprintf(" R version: %s\n", R.version.string))
  cat(sprintf(" Platform:  %s\n", R.version$platform))
  cat(sprintf(" MetricGraph: %s\n", packageVersion("MetricGraph")))
  cat(sprintf(" Date:      %s\n", Sys.Date()))
  cat(strrep("=", 70), "\n\n")
}

# ---------------------------------------------------------------------------
# Size configurations
# ---------------------------------------------------------------------------

# Grid sizes for non-fitting benchmarks (graph construction, mesh, data, plotting)
GRID_SIZES_FULL <- c(5, 10, 20, 30, 50, 100)
GRID_SIZES_QUICK <- c(5, 10, 20)

# Grid sizes for model fitting (smaller due to optimization cost)
FIT_SIZES_FULL <- c(5, 10, 20)
FIT_SIZES_QUICK <- c(5, 10)

# Mesh resolutions to benchmark
MESH_H_FULL <- c(0.5, 0.1, 0.05)
MESH_H_QUICK <- c(0.5, 0.1)

# Observation densities (per edge)
OBS_DENSITIES_FULL <- c(10, 50, 100, 200)
OBS_DENSITIES_QUICK <- c(10, 50)

#' Get the appropriate sizes based on mode.
get_grid_sizes <- function(quick = FALSE) {
  if (quick) GRID_SIZES_QUICK else GRID_SIZES_FULL
}

get_fit_sizes <- function(quick = FALSE) {
  if (quick) FIT_SIZES_QUICK else FIT_SIZES_FULL
}

get_mesh_h <- function(quick = FALSE) {
  if (quick) MESH_H_QUICK else MESH_H_FULL
}

get_obs_densities <- function(quick = FALSE) {
  if (quick) OBS_DENSITIES_QUICK else OBS_DENSITIES_FULL
}
