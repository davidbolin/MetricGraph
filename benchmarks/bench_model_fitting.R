## bench_model_fitting.R
## Benchmark graph_lme model fitting across graph sizes and model types.

bench_model_fitting <- function(quick = FALSE, include_slow = !quick) {
  print_section("Model Fitting (graph_lme)")

  sizes <- get_fit_sizes(quick)
  # isoexp and GL1 are much slower; only include when include_slow = TRUE
  models <- if (include_slow) {
    c("WM1", "WM2", "GL1", "isoexp")
  } else {
    c("WM1", "WM2")
  }
  if (!include_slow) {
    cat("  (Skipping GL1 and isoexp — use --slow to include)\n\n")
  }
  n_obs_per_edge <- 20
  results <- list()

  # --- Model fitting across sizes and model types ---
  cat("  Fitting models (", n_obs_per_edge, " obs/edge):\n", sep = "")

  for (n in sizes) {
    edges <- make_grid_edges(n)
    nE <- 2 * n * (n - 1)
    total_obs <- n_obs_per_edge * nE

    for (model in models) {
      label <- sprintf("%dx%d_%s", n, n, model)
      cat(sprintf("    %dx%d, model=%s (nObs=%d) ... ", n, n, model, total_obs))

      # Build a fresh graph with data for each benchmark
      graph <- metric_graph$new(edges = edges, verbose = 0)
      graph$build_mesh(h = 0.5)

      # For isoexp, need to check Euclidean
      if (model == "isoexp") {
        tryCatch(graph$check_euclidean(), error = function(e) NULL)
      }

      # Generate observations
      alpha <- if (model == "WM2") 2 else 1
      obs_data <- generate_observations(graph, n_obs_per_edge, alpha = alpha)
      graph$add_observations(data = obs_data, normalized = TRUE, verbose = 0)

      # For GL models, observations must be snapped to vertices first
      if (grepl("^GL", model)) {
        graph$observation_to_vertex()
        graph$compute_laplacian()
      }

      # Time the model fit (single iteration since fitting can be slow)
      t0 <- proc.time()
      tryCatch({
        res <- graph_lme(y ~ -1, graph = graph, model = model)
        elapsed <- (proc.time() - t0)[["elapsed"]]
        results[[label]] <- list(
          median = elapsed,
          model = model,
          nV = n * n,
          nE = nE,
          n_obs = total_obs,
          converged = TRUE
        )
        cat(format_time(elapsed), "\n")
      }, error = function(e) {
        elapsed <- (proc.time() - t0)[["elapsed"]]
        results[[label]] <<- list(
          median = elapsed,
          model = model,
          nV = n * n,
          nE = nE,
          n_obs = total_obs,
          converged = FALSE
        )
        cat(sprintf("FAILED (%s): %s\n", format_time(elapsed),
                    conditionMessage(e)))
      })
    }
  }

  # --- Replicate models on medium graph ---
  cat("\n  Replicate models (WM1, 10x10 grid):\n")
  n_repl_values <- if (quick) c(1, 5) else c(1, 5, 10)
  n <- 10
  edges <- make_grid_edges(n)
  nE <- 2 * n * (n - 1)
  n_obs_repl <- 15

  for (n_repl in n_repl_values) {
    total_obs <- n_obs_repl * nE * n_repl
    label <- sprintf("repl_%d_%dx%d_WM1", n_repl, n, n)
    cat(sprintf("    n_repl=%d (total obs=%d) ... ", n_repl, total_obs))

    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)

    # Generate replicate observations
    PtE <- do.call(rbind, lapply(seq_len(graph$nE), function(i)
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

    t0 <- proc.time()
    tryCatch({
      res <- graph_lme(y ~ -1, graph = graph, model = "WM1")
      elapsed <- (proc.time() - t0)[["elapsed"]]
      results[[label]] <- list(
        median = elapsed,
        model = "WM1",
        nV = n * n,
        nE = nE,
        n_obs = total_obs,
        n_repl = n_repl,
        converged = TRUE
      )
      cat(format_time(elapsed), "\n")
    }, error = function(e) {
      elapsed <- (proc.time() - t0)[["elapsed"]]
      results[[label]] <<- list(
        median = elapsed,
        model = "WM1",
        nV = n * n,
        nE = nE,
        n_obs = total_obs,
        n_repl = n_repl,
        converged = FALSE
      )
      cat(sprintf("FAILED (%s): %s\n", format_time(elapsed),
                  conditionMessage(e)))
    })
  }

  # --- Summary table ---
  df <- do.call(rbind, lapply(names(results), function(nm) {
    r <- results[[nm]]
    data.frame(
      scenario = nm,
      model = r$model,
      nV = r$nV,
      nE = r$nE,
      n_obs = r$n_obs,
      n_repl = if (!is.null(r$n_repl)) r$n_repl else 1L,
      time = format_time(r$median),
      converged = r$converged,
      stringsAsFactors = FALSE
    )
  }))
  rownames(df) <- NULL

  cat("\n  Model fitting summary:\n")
  print_results(df)

  invisible(df)
}
