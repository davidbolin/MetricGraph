## bench_data_operations.R
## Benchmark data operations: add_observations, observation_to_vertex,
## compute_geodist, compute_resdist.

bench_data_operations <- function(quick = FALSE, include_slow = FALSE) {
  print_section("Data Operations")

  sizes <- get_grid_sizes(quick)
  obs_densities <- get_obs_densities(quick)

  results_add <- list()
  results_o2v <- list()
  results_geodist <- list()
  results_resdist <- list()

  # --- add_observations at varying densities (plain graphs) ---
  cat("  add_observations (plain graphs):\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    nE <- 2 * n * (n - 1)

    for (nobs in obs_densities) {
      total_obs <- nobs * nE
      # Skip very large combos that would be excessive
      if (total_obs > 500000) next

      label <- sprintf("add_obs_plain_%dx%d_%dpe", n, n, nobs)
      cat(sprintf("    %dx%d, %d obs/edge (total=%d) ... ", n, n, nobs, total_obs))

      graph <- metric_graph$new(edges = edges, verbose = 0)
      graph$build_mesh(h = 0.5)

      PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
        cbind(rep(i, nobs), runif(nobs))))
      y <- rnorm(total_obs)
      df_obs <- data.frame(y = y, edge_number = PtE[, 1],
                           distance_on_edge = PtE[, 2])

      bm <- run_bench({
        g <- graph$clone()
        g$add_observations(data = df_obs, normalized = TRUE, verbose = 0)
      })
      results_add[[label]] <- bm
      results_add[[label]]$nV <- n * n
      results_add[[label]]$nE <- nE
      results_add[[label]]$n_obs <- total_obs
      results_add[[label]]$type <- "plain"
      cat(format_time(bm$median), "\n")
    }
  }

  # --- add_observations (CRS graphs) ---
  cat("\n  add_observations (CRS graphs):\n")
  for (n in sizes) {
    nE <- 2 * n * (n - 1)

    for (nobs in obs_densities) {
      total_obs <- nobs * nE
      if (total_obs > 500000) next

      label <- sprintf("add_obs_crs_%dx%d_%dpe", n, n, nobs)
      cat(sprintf("    %dx%d, %d obs/edge (total=%d) ... ", n, n, nobs, total_obs))

      sf_edges <- make_longlat_grid(n)
      graph <- metric_graph$new(edges = sf_edges, verbose = 0)
      graph$build_mesh(h = 0.5)

      PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
        cbind(rep(i, nobs), runif(nobs))))
      y <- rnorm(total_obs)
      df_obs <- data.frame(y = y, edge_number = PtE[, 1],
                           distance_on_edge = PtE[, 2])

      bm <- run_bench({
        g <- graph$clone()
        g$add_observations(data = df_obs, normalized = TRUE, verbose = 0)
      })
      results_add[[label]] <- bm
      results_add[[label]]$nV <- n * n
      results_add[[label]]$nE <- nE
      results_add[[label]]$n_obs <- total_obs
      results_add[[label]]$type <- "crs"
      cat(format_time(bm$median), "\n")
    }
  }

  # --- observation_to_vertex ---
  cat("\n  observation_to_vertex (plain graphs):\n")
  n_o2v <- if (quick) c(10, 20) else c(10, 20, 30, 50)
  nobs_o2v <- 20
  for (n in n_o2v) {
    edges <- make_grid_edges(n)
    nE <- 2 * n * (n - 1)
    total_obs <- nobs_o2v * nE
    label <- sprintf("o2v_plain_%dx%d_%dpe", n, n, nobs_o2v)
    cat(sprintf("    %dx%d, %d obs/edge (total=%d) ... ", n, n, nobs_o2v, total_obs))

    graph <- metric_graph$new(edges = edges, verbose = 0)
    PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
      cbind(rep(i, nobs_o2v), runif(nobs_o2v))))
    y <- rnorm(total_obs)
    graph$add_observations(
      data = data.frame(y = y, edge_number = PtE[, 1],
                        distance_on_edge = PtE[, 2]),
      normalized = TRUE, verbose = 0
    )

    bm <- run_bench({
      g <- graph$clone()
      g$observation_to_vertex()
    })
    results_o2v[[label]] <- bm
    results_o2v[[label]]$nV <- n * n
    results_o2v[[label]]$n_obs <- total_obs
    cat(format_time(bm$median), "\n")
  }

  # --- compute_geodist ---
  cat("\n  compute_geodist:\n")
  geodist_sizes <- if (quick) c(5, 10, 20) else c(5, 10, 20, 30, 50)
  for (n in geodist_sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)
    label <- sprintf("geodist_%dx%d", n, n)
    cat(sprintf("    %dx%d (V=%d) ... ", n, n, n * n))
    bm <- run_bench(graph$compute_geodist())
    results_geodist[[label]] <- bm
    results_geodist[[label]]$nV <- n * n
    cat(format_time(bm$median), "\n")
  }

  # --- compute_resdist (slow — only with --slow) ---
  if (include_slow) {
    cat("\n  compute_resdist:\n")
    resdist_sizes <- if (quick) c(5, 10, 20) else c(5, 10, 20, 30, 50)
    for (n in resdist_sizes) {
      edges <- make_grid_edges(n)
      graph <- metric_graph$new(edges = edges, verbose = 0)
      label <- sprintf("resdist_%dx%d", n, n)
      cat(sprintf("    %dx%d (V=%d) ... ", n, n, n * n))
      bm <- run_bench(graph$compute_resdist())
      results_resdist[[label]] <- bm
      results_resdist[[label]]$nV <- n * n
      cat(format_time(bm$median), "\n")
    }
  } else {
    cat("\n  (Skipping compute_resdist — use --slow to include)\n")
  }

  # --- Summary tables ---
  df_add <- do.call(rbind, lapply(names(results_add), function(nm) {
    r <- results_add[[nm]]
    data.frame(
      scenario = nm,
      type = r$type,
      nV = r$nV,
      nE = r$nE,
      n_obs = r$n_obs,
      median = format_time(r$median),
      mem_alloc = format(r$mem_alloc),
      stringsAsFactors = FALSE
    )
  }))
  rownames(df_add) <- NULL

  df_o2v <- do.call(rbind, lapply(names(results_o2v), function(nm) {
    r <- results_o2v[[nm]]
    data.frame(
      scenario = nm,
      nV = r$nV,
      n_obs = r$n_obs,
      median = format_time(r$median),
      mem_alloc = format(r$mem_alloc),
      stringsAsFactors = FALSE
    )
  }))
  rownames(df_o2v) <- NULL

  cat("\n  add_observations summary:\n")
  print_results(df_add)

  cat("  observation_to_vertex summary:\n")
  print_results(df_o2v)

  invisible(list(add_obs = df_add, o2v = df_o2v))
}
