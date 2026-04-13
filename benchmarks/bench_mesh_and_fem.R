## bench_mesh_and_fem.R
## Benchmark mesh construction, FEM assembly, and Laplacian computation.

bench_mesh_and_fem <- function(quick = FALSE) {
  print_section("Mesh Construction & FEM Assembly")

  sizes <- get_grid_sizes(quick)
  h_values <- get_mesh_h(quick)

  results_mesh <- list()
  results_fem <- list()
  results_laplacian <- list()

  # --- Plain graphs: mesh + FEM ---
  cat("  Plain edge-list graphs:\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)

    for (h in h_values) {
      label <- sprintf("plain_%dx%d_h%.2f", n, n, h)

      # Mesh construction
      cat(sprintf("    %dx%d, h=%.2f: build_mesh ... ", n, n, h))
      g_clone <- graph$clone()
      bm_mesh <- run_bench(g_clone$build_mesh(h = h))
      n_mesh <- nrow(g_clone$mesh$V)
      results_mesh[[label]] <- bm_mesh
      results_mesh[[label]]$nV <- n * n
      results_mesh[[label]]$h <- h
      results_mesh[[label]]$n_mesh <- n_mesh
      results_mesh[[label]]$type <- "plain"
      cat(format_time(bm_mesh$median), sprintf(" (%d mesh nodes)\n", n_mesh))

      # FEM assembly (requires mesh)
      cat(sprintf("    %dx%d, h=%.2f: compute_fem ... ", n, n, h))
      bm_fem <- run_bench(g_clone$compute_fem())
      results_fem[[label]] <- bm_fem
      results_fem[[label]]$nV <- n * n
      results_fem[[label]]$h <- h
      results_fem[[label]]$n_mesh <- n_mesh
      results_fem[[label]]$type <- "plain"
      cat(format_time(bm_fem$median), "\n")
    }
  }

  # --- CRS graphs: mesh + FEM ---
  cat("\n  CRS 4326 (sf) graphs:\n")
  for (n in sizes) {
    sf_edges <- make_longlat_grid(n)
    graph <- metric_graph$new(edges = sf_edges, verbose = 0)

    for (h in h_values) {
      label <- sprintf("crs_%dx%d_h%.2f", n, n, h)

      cat(sprintf("    %dx%d, h=%.2f: build_mesh ... ", n, n, h))
      g_clone <- graph$clone()
      bm_mesh <- run_bench(g_clone$build_mesh(h = h))
      n_mesh <- nrow(g_clone$mesh$V)
      results_mesh[[label]] <- bm_mesh
      results_mesh[[label]]$nV <- n * n
      results_mesh[[label]]$h <- h
      results_mesh[[label]]$n_mesh <- n_mesh
      results_mesh[[label]]$type <- "crs"
      cat(format_time(bm_mesh$median), sprintf(" (%d mesh nodes)\n", n_mesh))

      cat(sprintf("    %dx%d, h=%.2f: compute_fem ... ", n, n, h))
      bm_fem <- run_bench(g_clone$compute_fem())
      results_fem[[label]] <- bm_fem
      results_fem[[label]]$nV <- n * n
      results_fem[[label]]$h <- h
      results_fem[[label]]$n_mesh <- n_mesh
      results_fem[[label]]$type <- "crs"
      cat(format_time(bm_fem$median), "\n")
    }
  }

  # --- Laplacian computation ---
  cat("\n  Laplacian computation:\n")
  for (n in sizes) {
    edges <- make_grid_edges(n)
    graph <- metric_graph$new(edges = edges, verbose = 0)
    graph$build_mesh(h = 0.1)

    label <- sprintf("laplacian_%dx%d", n, n)
    cat(sprintf("    %dx%d: compute_laplacian ... ", n, n))
    bm <- run_bench(graph$compute_laplacian())
    results_laplacian[[label]] <- bm
    results_laplacian[[label]]$nV <- n * n
    cat(format_time(bm$median), "\n")
  }

  # --- Summary tables ---
  df_mesh <- do.call(rbind, lapply(names(results_mesh), function(nm) {
    r <- results_mesh[[nm]]
    data.frame(
      scenario = nm,
      type = r$type,
      nV = r$nV,
      h = r$h,
      n_mesh = r$n_mesh,
      median = format_time(r$median),
      mem_alloc = format(r$mem_alloc),
      stringsAsFactors = FALSE
    )
  }))
  rownames(df_mesh) <- NULL

  df_fem <- do.call(rbind, lapply(names(results_fem), function(nm) {
    r <- results_fem[[nm]]
    data.frame(
      scenario = nm,
      type = r$type,
      n_mesh = r$n_mesh,
      median = format_time(r$median),
      mem_alloc = format(r$mem_alloc),
      stringsAsFactors = FALSE
    )
  }))
  rownames(df_fem) <- NULL

  cat("\n  Mesh construction summary:\n")
  print_results(df_mesh)

  cat("  FEM assembly summary:\n")
  print_results(df_fem)

  invisible(list(mesh = df_mesh, fem = df_fem))
}
