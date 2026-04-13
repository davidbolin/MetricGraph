## bench_graph_spde.R
## Benchmark the INLA interface functions from graph_spde.R:
## graph_spde(), graph_data_spde().
## Requires INLA to be installed.

bench_graph_spde <- function(quick = FALSE) {
  print_section("graph_spde Interface (INLA)")

  if (!requireNamespace("INLA", quietly = TRUE)) {
    cat("  INLA not installed — skipping graph_spde benchmarks.\n")
    return(invisible(NULL))
  }

  sizes <- get_grid_sizes(quick)
  n_obs_per_edge <- 20
  results_spde <- list()
  results_data <- list()

  # -------------------------------------------------------------------------
  # graph_spde(): model setup (alpha=1, non-directional)
  # -------------------------------------------------------------------------
  cat("  graph_spde (alpha=1):\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    nE <- 2 * n * (n - 1)
    total_obs <- n_obs_per_edge * nE

    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()
    obs <- generate_observations(graph, n_obs_per_edge)
    graph$add_observations(data = obs, normalized = TRUE, verbose = 0)

    label <- sprintf("spde_a1_%dx%d", n, n)
    cat(sprintf("    %dx%d (nV=%d, nE=%d, nObs=%d) ... ", n, n, n*n, nE, total_obs))
    t0 <- proc.time()
    tryCatch({
      spde_obj <- graph_spde(graph, alpha = 1)
      elapsed <- (proc.time() - t0)[["elapsed"]]
      results_spde[[label]] <- list(median = elapsed, alpha = 1,
                                    nV = n*n, nE = nE, n_obs = total_obs)
      cat(format_time(elapsed), "\n")
    }, error = function(e) {
      elapsed <- (proc.time() - t0)[["elapsed"]]
      results_spde[[label]] <<- list(median = elapsed, alpha = 1,
                                     nV = n*n, nE = nE, n_obs = total_obs,
                                     error = conditionMessage(e))
      cat(sprintf("FAILED (%s): %s\n", format_time(elapsed), conditionMessage(e)))
    })
  }

  cat("\n  graph_spde (alpha=2):\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    nE <- 2 * n * (n - 1)
    total_obs <- n_obs_per_edge * nE

    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()
    # Use alpha=1 for data generation (alpha=2 sampling can fail on large graphs)
    obs <- generate_observations(graph, n_obs_per_edge, alpha = 1)
    graph$add_observations(data = obs, normalized = TRUE, verbose = 0)

    label <- sprintf("spde_a2_%dx%d", n, n)
    cat(sprintf("    %dx%d (nV=%d, nE=%d, nObs=%d) ... ", n, n, n*n, nE, total_obs))
    t0 <- proc.time()
    tryCatch({
      spde_obj <- graph_spde(graph, alpha = 2)
      elapsed <- (proc.time() - t0)[["elapsed"]]
      results_spde[[label]] <- list(median = elapsed, alpha = 2,
                                    nV = n*n, nE = nE, n_obs = total_obs)
      cat(format_time(elapsed), "\n")
    }, error = function(e) {
      elapsed <- (proc.time() - t0)[["elapsed"]]
      results_spde[[label]] <<- list(median = elapsed, alpha = 2,
                                     nV = n*n, nE = nE, n_obs = total_obs,
                                     error = conditionMessage(e))
      cat(sprintf("FAILED (%s): %s\n", format_time(elapsed), conditionMessage(e)))
    })
  }

  # -------------------------------------------------------------------------
  # graph_spde scaling with observation density (fixed graph size)
  # -------------------------------------------------------------------------
  cat("\n  graph_spde (alpha=1) varying obs density on 20x20 grid:\n")
  n <- 20
  edges <- make_grid_edges(n)
  nE <- 2 * n * (n - 1)
  obs_densities <- if (quick) c(10, 50) else c(10, 50, 100, 200)

  for (nobs in obs_densities) {
    total_obs <- nobs * nE
    if (total_obs > 500000) next

    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()
    obs <- generate_observations(graph, nobs)
    graph$add_observations(data = obs, normalized = TRUE, verbose = 0)

    label <- sprintf("spde_a1_20x20_%dpe", nobs)
    cat(sprintf("    %d obs/edge (total=%d) ... ", nobs, total_obs))
    t0 <- proc.time()
    tryCatch({
      spde_obj <- graph_spde(graph, alpha = 1)
      elapsed <- (proc.time() - t0)[["elapsed"]]
      results_spde[[label]] <- list(median = elapsed, alpha = 1,
                                    nV = n*n, nE = nE, n_obs = total_obs)
      cat(format_time(elapsed), "\n")
    }, error = function(e) {
      elapsed <- (proc.time() - t0)[["elapsed"]]
      results_spde[[label]] <<- list(median = elapsed, alpha = 1,
                                     nV = n*n, nE = nE, n_obs = total_obs,
                                     error = conditionMessage(e))
      cat(sprintf("FAILED (%s): %s\n", format_time(elapsed), conditionMessage(e)))
    })
  }

  # -------------------------------------------------------------------------
  # graph_data_spde(): data preparation / basis construction
  # -------------------------------------------------------------------------
  cat("\n  graph_data_spde (alpha=1, varying graph size):\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    nE <- 2 * n * (n - 1)
    total_obs <- n_obs_per_edge * nE

    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()
    obs <- generate_observations(graph, n_obs_per_edge)
    graph$add_observations(data = obs, normalized = TRUE, verbose = 0)

    label_data <- sprintf("data_spde_a1_%dx%d", n, n)
    cat(sprintf("    %dx%d (nObs=%d) ... ", n, n, total_obs))
    t0_spde <- proc.time()
    spde_obj <- NULL
    tryCatch({
      spde_obj <- graph_spde(graph, alpha = 1)
    }, error = function(e) NULL)

    if (!is.null(spde_obj)) {
      t0 <- proc.time()
      tryCatch({
        data_obj <- graph_data_spde(spde_obj, name = "field")
        elapsed <- (proc.time() - t0)[["elapsed"]]
        results_data[[label_data]] <- list(median = elapsed,
                                           nV = n*n, nE = nE, n_obs = total_obs)
        cat(format_time(elapsed), "\n")
      }, error = function(e) {
        elapsed <- (proc.time() - t0)[["elapsed"]]
        results_data[[label_data]] <<- list(median = elapsed,
                                            nV = n*n, nE = nE, n_obs = total_obs,
                                            error = conditionMessage(e))
        cat(sprintf("FAILED (%s): %s\n", format_time(elapsed), conditionMessage(e)))
      })
    } else {
      cat("SKIPPED (graph_spde failed)\n")
    }
  }

  # graph_data_spde with replicates
  cat("\n  graph_data_spde (alpha=1, 10x10, varying replicates):\n")
  n <- 10
  edges <- make_grid_edges(n)
  nE <- 2 * n * (n - 1)
  n_obs_repl <- 15
  repl_values <- if (quick) c(1, 5) else c(1, 5, 10, 30)

  for (n_repl in repl_values) {
    total_obs <- n_obs_repl * nE * n_repl
    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()

    PtE <- do.call(rbind, lapply(seq_len(nE), function(i)
      cbind(rep(i, n_obs_repl), runif(n_obs_repl))))
    u <- sample_spde(range = 0.2, sigma = 1.3, alpha = 1,
                     graph = graph, PtE = PtE, nsim = n_repl)
    y <- u + 0.1 * matrix(rnorm(nrow(PtE) * n_repl), ncol = n_repl)

    df_repl <- data.frame(y, edge_number = PtE[, 1],
                          distance_on_edge = PtE[, 2])
    if (n_repl == 1) {
      colnames(df_repl)[1] <- "y"
      graph$add_observations(data = df_repl, normalized = TRUE, verbose = 0)
    } else {
      y_cols <- paste0("y.", seq_len(n_repl))
      colnames(df_repl)[1:n_repl] <- y_cols
      df_long <- tidyr::pivot_longer(df_repl, cols = all_of(y_cols),
                                     names_to = "repl", values_to = "y")
      graph$add_observations(data = df_long, normalized = TRUE,
                             group = "repl", verbose = 0)
    }

    label <- sprintf("data_spde_repl%d_10x10", n_repl)
    cat(sprintf("    n_repl=%d (total obs=%d) ... ", n_repl, total_obs))
    spde_obj <- NULL
    tryCatch(spde_obj <- graph_spde(graph, alpha = 1), error = function(e) NULL)

    if (!is.null(spde_obj)) {
      t0 <- proc.time()
      tryCatch({
        data_obj <- graph_data_spde(spde_obj, name = "field")
        elapsed <- (proc.time() - t0)[["elapsed"]]
        results_data[[label]] <- list(median = elapsed, n_repl = n_repl,
                                      nV = n*n, nE = nE, n_obs = total_obs)
        cat(format_time(elapsed), "\n")
      }, error = function(e) {
        elapsed <- (proc.time() - t0)[["elapsed"]]
        results_data[[label]] <<- list(median = elapsed, n_repl = n_repl,
                                       nV = n*n, nE = nE, n_obs = total_obs,
                                       error = conditionMessage(e))
        cat(sprintf("FAILED (%s): %s\n", format_time(elapsed), conditionMessage(e)))
      })
    } else {
      cat("SKIPPED (graph_spde failed)\n")
    }
  }

  # -------------------------------------------------------------------------
  # Summary tables
  # -------------------------------------------------------------------------
  df_spde <- do.call(rbind, lapply(names(results_spde), function(nm) {
    r <- results_spde[[nm]]
    data.frame(scenario = nm, alpha = r$alpha,
               nV = r$nV, nE = r$nE, n_obs = r$n_obs,
               time = format_time(r$median),
               status = if (is.null(r$error)) "OK" else "FAIL",
               stringsAsFactors = FALSE)
  }))
  rownames(df_spde) <- NULL

  df_data <- do.call(rbind, lapply(names(results_data), function(nm) {
    r <- results_data[[nm]]
    data.frame(scenario = nm,
               nV = r$nV, nE = r$nE, n_obs = r$n_obs,
               n_repl = if (!is.null(r$n_repl)) r$n_repl else 1L,
               time = format_time(r$median),
               status = if (is.null(r$error)) "OK" else "FAIL",
               stringsAsFactors = FALSE)
  }))
  rownames(df_data) <- NULL

  cat("\n  graph_spde summary:\n")
  print_results(df_spde)

  cat("  graph_data_spde summary:\n")
  print_results(df_data)

  invisible(list(graph_spde = df_spde, graph_data_spde = df_data))
}
