## bench_plotting.R
## Benchmark plotting functions at varying graph sizes.

bench_plotting <- function(quick = FALSE) {
  print_section("Plotting")

  sizes <- get_grid_sizes(quick)
  results_plot <- list()
  results_plotfn <- list()

  # --- graph$plot() for plain graphs ---
  cat("  graph$plot() (plain graphs):\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)

    label <- sprintf("plot_plain_%dx%d", n, n)
    cat(sprintf("    %dx%d (V=%d, E=%d) ... ", n, n, graph$nV, graph$nE))
    bm <- run_bench(graph$plot())
    results_plot[[label]] <- bm
    results_plot[[label]]$nV <- graph$nV
    results_plot[[label]]$nE <- graph$nE
    results_plot[[label]]$type <- "plain"
    cat(format_time(bm$median), "\n")
  }

  # --- graph$plot() for CRS graphs ---
  cat("\n  graph$plot() (CRS graphs):\n")
  for (n in sizes) {
    sf_edges <- make_longlat_grid(n)
    graph <- metric_graph$new(edges = sf_edges, verbose = 0)

    label <- sprintf("plot_crs_%dx%d", n, n)
    cat(sprintf("    %dx%d (V=%d, E=%d) ... ", n, n, graph$nV, graph$nE))
    bm <- run_bench(graph$plot())
    results_plot[[label]] <- bm
    results_plot[[label]]$nV <- graph$nV
    results_plot[[label]]$nE <- graph$nE
    results_plot[[label]]$type <- "crs"
    cat(format_time(bm$median), "\n")
  }

  # --- plot_function with sampled field ---
  # Use smaller sizes since plot_function needs mesh + sampled data
  plotfn_sizes <- if (quick) c(5, 10) else c(5, 10, 20, 30, 50)
  cat("\n  graph$plot_function() (plain graphs, sampled field on mesh):\n")
  for (n in plotfn_sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.1)

    # Sample a field on the mesh
    X <- sample_spde(range = 0.2, sigma = 1.3, alpha = 1,
                     graph = graph, type = "mesh")
    n_mesh <- nrow(graph$mesh$V)

    label <- sprintf("plotfn_plain_%dx%d", n, n)
    cat(sprintf("    %dx%d (mesh=%d) ... ", n, n, n_mesh))
    bm <- run_bench(suppressWarnings(graph$plot_function(X = X)))
    results_plotfn[[label]] <- bm
    results_plotfn[[label]]$nV <- graph$nV
    results_plotfn[[label]]$n_mesh <- n_mesh
    results_plotfn[[label]]$type <- "plain"
    cat(format_time(bm$median), "\n")
  }

  # --- plot_function for CRS graphs ---
  cat("\n  graph$plot_function() (CRS graphs, sampled field on mesh):\n")
  for (n in plotfn_sizes) {
    sf_edges <- make_longlat_grid(n)
    graph <- metric_graph$new(edges = sf_edges, verbose = 0)
    graph$build_mesh(h = 0.1)

    X <- sample_spde(range = 0.2, sigma = 1.3, alpha = 1,
                     graph = graph, type = "mesh")
    n_mesh <- nrow(graph$mesh$V)

    label <- sprintf("plotfn_crs_%dx%d", n, n)
    cat(sprintf("    %dx%d (mesh=%d) ... ", n, n, n_mesh))
    bm <- run_bench(suppressWarnings(graph$plot_function(X = X)))
    results_plotfn[[label]] <- bm
    results_plotfn[[label]]$nV <- graph$nV
    results_plotfn[[label]]$n_mesh <- n_mesh
    results_plotfn[[label]]$type <- "crs"
    cat(format_time(bm$median), "\n")
  }

  # --- Summary tables ---
  df_plot <- do.call(rbind, lapply(names(results_plot), function(nm) {
    r <- results_plot[[nm]]
    data.frame(
      scenario = nm,
      type = r$type,
      nV = r$nV,
      nE = r$nE,
      median = format_time(r$median),
      mem_alloc = format(r$mem_alloc),
      stringsAsFactors = FALSE
    )
  }))
  rownames(df_plot) <- NULL

  df_plotfn <- do.call(rbind, lapply(names(results_plotfn), function(nm) {
    r <- results_plotfn[[nm]]
    data.frame(
      scenario = nm,
      type = r$type,
      nV = r$nV,
      n_mesh = r$n_mesh,
      median = format_time(r$median),
      mem_alloc = format(r$mem_alloc),
      stringsAsFactors = FALSE
    )
  }))
  rownames(df_plotfn) <- NULL

  cat("\n  graph$plot() summary:\n")
  print_results(df_plot)

  cat("  graph$plot_function() summary:\n")
  print_results(df_plotfn)

  invisible(list(plot = df_plot, plot_function = df_plotfn))
}
