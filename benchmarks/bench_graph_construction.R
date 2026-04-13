## bench_graph_construction.R
## Benchmark graph construction at varying sizes, with and without CRS.

bench_graph_construction <- function(quick = FALSE) {
  print_section("Graph Construction")

  sizes <- get_grid_sizes(quick)
  results <- list()

  # --- Plain edge-list construction ---
  cat("  Plain edge-list graphs:\n")
  for (n in sizes) {
    label <- sprintf("plain_%dx%d", n, n)
    edges <- make_grid_edges(n)
    nV <- n * n
    nE <- 2 * n * (n - 1)
    cat(sprintf("    %dx%d grid (V=%d, E=%d) ... ", n, n, nV, nE))
    bm <- run_bench(metric_graph$new(edges = edges, verbose = 0))
    results[[label]] <- bm
    results[[label]]$nV <- nV
    results[[label]]$nE <- nE
    results[[label]]$type <- "plain"
    cat(format_time(bm$median), "\n")
  }

  # --- CRS (longlat sf) construction ---
  cat("\n  CRS 4326 (sf) graphs:\n")
  for (n in sizes) {
    label <- sprintf("crs_%dx%d", n, n)
    sf_edges <- make_longlat_grid(n)
    nV <- n * n
    nE <- 2 * n * (n - 1)
    cat(sprintf("    %dx%d grid (V=%d, E=%d) ... ", n, n, nV, nE))
    bm <- run_bench(metric_graph$new(edges = sf_edges, verbose = 0))
    results[[label]] <- bm
    results[[label]]$nV <- nV
    results[[label]]$nE <- nE
    results[[label]]$type <- "crs"
    cat(format_time(bm$median), "\n")
  }

  # --- Summary table ---
  df <- do.call(rbind, lapply(names(results), function(nm) {
    r <- results[[nm]]
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
  rownames(df) <- NULL

  cat("\n  Summary:\n")
  print_results(df)

  invisible(df)
}
