## bench_spde.R
## Benchmark SPDE functions: spde_precision, sample_spde, spde_covariance.

bench_spde <- function(quick = FALSE) {
  print_section("SPDE Functions")

  sizes <- get_grid_sizes(quick)
  results_precision <- list()
  results_sample <- list()
  results_covariance <- list()

  kappa <- 20
  tau <- 1
  sigma <- 1.3
  range_val <- 0.2

  # -------------------------------------------------------------------------
  # spde_precision: build Q matrix for alpha = 1 and alpha = 2
  # -------------------------------------------------------------------------
  cat("  spde_precision (alpha=1):\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()
    n_mesh <- nrow(graph$mesh$V)

    label <- sprintf("Q_a1_%dx%d", n, n)
    cat(sprintf("    %dx%d (mesh=%d) ... ", n, n, n_mesh))
    bm <- run_bench(spde_precision(kappa = kappa, tau = tau,
                                   alpha = 1, graph = graph))
    bm$n_mesh <- n_mesh
    results_precision[[label]] <- bm
    cat(format_time(bm$median), "\n")
  }

  cat("\n  spde_precision (alpha=2):\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()
    n_mesh <- nrow(graph$mesh$V)

    label <- sprintf("Q_a2_%dx%d", n, n)
    cat(sprintf("    %dx%d (mesh=%d) ... ", n, n, n_mesh))
    bm <- run_bench(spde_precision(kappa = kappa, tau = tau,
                                   alpha = 2, graph = graph))
    bm$n_mesh <- n_mesh
    results_precision[[label]] <- bm
    cat(format_time(bm$median), "\n")
  }

  # -------------------------------------------------------------------------
  # sample_spde: draw samples at mesh locations
  # -------------------------------------------------------------------------
  cat("\n  sample_spde (alpha=1, type='mesh', method='Q'):\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()
    n_mesh <- nrow(graph$mesh$V)

    label <- sprintf("sample_a1_%dx%d", n, n)
    cat(sprintf("    %dx%d (mesh=%d) ... ", n, n, n_mesh))
    bm <- run_bench(sample_spde(range = range_val, sigma = sigma,
                                alpha = 1, graph = graph,
                                type = "mesh", method = "Q"))
    bm$n_mesh <- n_mesh
    results_sample[[label]] <- bm
    cat(format_time(bm$median), "\n")
  }

  cat("\n  sample_spde (alpha=2, type='mesh', method='Q'):\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()
    n_mesh <- nrow(graph$mesh$V)

    label <- sprintf("sample_a2_%dx%d", n, n)
    cat(sprintf("    %dx%d (mesh=%d) ... ", n, n, n_mesh))
    bm <- run_bench(sample_spde(range = range_val, sigma = sigma,
                                alpha = 2, graph = graph,
                                type = "mesh", method = "Q"))
    bm$n_mesh <- n_mesh
    results_sample[[label]] <- bm
    cat(format_time(bm$median), "\n")
  }

  # sample at observation locations (PtE)
  cat("\n  sample_spde (alpha=1, at PtE locations, varying n_obs):\n")
  n <- if (quick) 10 else 20
  edges <- make_grid_edges(n)
  nE <- 2 * n * (n - 1)
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$build_mesh(h = 0.5)
  graph$compute_fem()

  obs_counts <- if (quick) c(10, 50) else c(10, 50, 100, 200)
  for (nobs in obs_counts) {
    total <- nobs * nE
    PtE <- do.call(rbind, lapply(seq_len(nE), function(i)
      cbind(rep(i, nobs), runif(nobs))))
    label <- sprintf("sample_PtE_%dx%d_%dpe", n, n, nobs)
    cat(sprintf("    %dx%d, %d obs/edge (total=%d) ... ", n, n, nobs, total))
    bm <- run_bench(sample_spde(range = range_val, sigma = sigma,
                                alpha = 1, graph = graph, PtE = PtE))
    bm$n_obs <- total
    results_sample[[label]] <- bm
    cat(format_time(bm$median), "\n")
  }

  # sample with multiple replicates
  cat("\n  sample_spde (alpha=1, mesh, nsim varying):\n")
  n <- if (quick) 10 else 20
  edges <- make_grid_edges(n)
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$build_mesh(h = 0.5)
  graph$compute_fem()
  n_mesh <- nrow(graph$mesh$V)

  nsim_values <- if (quick) c(1, 10) else c(1, 10, 50, 100)
  for (nsim in nsim_values) {
    label <- sprintf("sample_nsim%d_%dx%d", nsim, n, n)
    cat(sprintf("    %dx%d (mesh=%d), nsim=%d ... ", n, n, n_mesh, nsim))
    bm <- run_bench(sample_spde(range = range_val, sigma = sigma,
                                alpha = 1, graph = graph,
                                type = "mesh", method = "Q",
                                nsim = nsim))
    bm$nsim <- nsim
    bm$n_mesh <- n_mesh
    results_sample[[label]] <- bm
    cat(format_time(bm$median), "\n")
  }

  # -------------------------------------------------------------------------
  # spde_covariance: covariance from a single point to all mesh nodes
  # -------------------------------------------------------------------------
  cat("\n  spde_covariance (alpha=1):\n")
  cov_sizes <- if (quick) c(5, 10, 20) else c(5, 10, 20, 30, 50)
  for (n in cov_sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()
    n_mesh <- nrow(graph$mesh$V)

    # Pick a point near the center of the graph
    P <- c(1, 0.5)

    label <- sprintf("cov_a1_%dx%d", n, n)
    cat(sprintf("    %dx%d (mesh=%d) ... ", n, n, n_mesh))
    bm <- run_bench(spde_covariance(P, kappa = kappa, tau = tau,
                                    alpha = 1, graph = graph))
    bm$n_mesh <- n_mesh
    results_covariance[[label]] <- bm
    cat(format_time(bm$median), "\n")
  }

  cat("\n  spde_covariance (alpha=2):\n")
  for (n in cov_sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.5)
    graph$compute_fem()
    n_mesh <- nrow(graph$mesh$V)
    P <- c(1, 0.5)

    label <- sprintf("cov_a2_%dx%d", n, n)
    cat(sprintf("    %dx%d (mesh=%d) ... ", n, n, n_mesh))
    bm <- run_bench(spde_covariance(P, kappa = kappa, tau = tau,
                                    alpha = 2, graph = graph))
    bm$n_mesh <- n_mesh
    results_covariance[[label]] <- bm
    cat(format_time(bm$median), "\n")
  }

  # -------------------------------------------------------------------------
  # Summary tables
  # -------------------------------------------------------------------------
  df_prec <- do.call(rbind, lapply(names(results_precision), function(nm) {
    r <- results_precision[[nm]]
    data.frame(scenario = nm, n_mesh = r$n_mesh,
               median = format_time(r$median),
               mem_alloc = format(r$mem_alloc),
               stringsAsFactors = FALSE)
  }))
  rownames(df_prec) <- NULL

  df_sample <- do.call(rbind, lapply(names(results_sample), function(nm) {
    r <- results_sample[[nm]]
    data.frame(scenario = nm,
               n_mesh = if (!is.null(r$n_mesh)) r$n_mesh else NA_integer_,
               n_obs = if (!is.null(r$n_obs)) r$n_obs else NA_integer_,
               nsim = if (!is.null(r$nsim)) r$nsim else 1L,
               median = format_time(r$median),
               mem_alloc = format(r$mem_alloc),
               stringsAsFactors = FALSE)
  }))
  rownames(df_sample) <- NULL

  df_cov <- do.call(rbind, lapply(names(results_covariance), function(nm) {
    r <- results_covariance[[nm]]
    data.frame(scenario = nm, n_mesh = r$n_mesh,
               median = format_time(r$median),
               mem_alloc = format(r$mem_alloc),
               stringsAsFactors = FALSE)
  }))
  rownames(df_cov) <- NULL

  cat("\n  spde_precision summary:\n")
  print_results(df_prec)

  cat("  sample_spde summary:\n")
  print_results(df_sample)

  cat("  spde_covariance summary:\n")
  print_results(df_cov)

  invisible(list(precision = df_prec, sample = df_sample, covariance = df_cov))
}
