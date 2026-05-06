## bench_prediction.R
## Benchmark prediction and cross-validation.

bench_prediction <- function(quick = FALSE, include_slow = FALSE) {
  print_section("Prediction & Cross-Validation")

  results_pred <- list()
  results_cv <- list()

  # --- Setup: fit a model on 10x10 grid ---
  cat("  Setting up fitted model (WM1, 10x10 grid) ...\n")
  n <- 10
  edges <- make_grid_edges(n)
  graph <- metric_graph$new(edges = edges, verbose = 0)
  graph$build_mesh(h = 0.1)

  obs_data <- generate_observations(graph, n_obs_per_edge = 20)
  graph$add_observations(data = obs_data, normalized = TRUE, verbose = 0)
  fit <- graph_lme(y ~ -1, graph = graph, model = "WM1")
  cat("  Model fitted.\n\n")

  # --- Prediction at varying number of locations ---
  n_pred_values <- if (quick) c(100, 500) else c(100, 500, 1000, 5000)

  cat("  predict() at varying number of locations:\n")
  for (n_pred in n_pred_values) {
    label <- sprintf("predict_n%d", n_pred)
    cat(sprintf("    n_pred=%d ... ", n_pred))

    # Generate prediction locations spread across edges
    pred_edge <- sample(seq_len(graph$nE), n_pred, replace = TRUE)
    pred_dist <- runif(n_pred)
    pred_df <- data.frame(edge_number = pred_edge,
                          distance_on_edge = pred_dist)

    bm <- run_bench(predict(fit, newdata = pred_df, normalized = TRUE))
    results_pred[[label]] <- bm
    results_pred[[label]]$n_pred <- n_pred
    cat(format_time(bm$median), "\n")
  }

  # --- Predict at mesh locations ---
  cat(sprintf("    mesh locations (n=%d) ... ", nrow(graph$mesh$VtE)))
  mesh_df <- data.frame(edge_number = graph$mesh$VtE[, 1],
                        distance_on_edge = graph$mesh$VtE[, 2])
  bm <- run_bench(predict(fit, newdata = mesh_df, normalized = TRUE))
  results_pred[["predict_mesh"]] <- bm
  results_pred[["predict_mesh"]]$n_pred <- nrow(graph$mesh$VtE)
  cat(format_time(bm$median), "\n")

  # --- Cross-validation (slow — only with --slow) ---
  if (include_slow) {
    cat("\n  posterior_crossvalidation:\n")
    cat("    LOO cross-validation ... ")
    t0 <- proc.time()
    tryCatch({
      cv <- posterior_crossvalidation(list("WM1" = fit), mode = "loo",
                                      factor = 1000)
      elapsed <- (proc.time() - t0)[["elapsed"]]
      results_cv[["cv_loo"]] <- list(median = elapsed)
      cat(format_time(elapsed), "\n")
    }, error = function(e) {
      elapsed <- (proc.time() - t0)[["elapsed"]]
      results_cv[["cv_loo"]] <<- list(median = elapsed, error = conditionMessage(e))
      cat(sprintf("FAILED (%s): %s\n", format_time(elapsed), conditionMessage(e)))
    })
  } else {
    cat("\n  (Skipping cross-validation — use --slow to include)\n")
  }

  # --- Summary ---
  df_pred <- do.call(rbind, lapply(names(results_pred), function(nm) {
    r <- results_pred[[nm]]
    data.frame(
      scenario = nm,
      n_pred = r$n_pred,
      median = format_time(r$median),
      mem_alloc = format(r$mem_alloc),
      stringsAsFactors = FALSE
    )
  }))
  rownames(df_pred) <- NULL

  cat("\n  Prediction summary:\n")
  print_results(df_pred)

  invisible(list(pred = df_pred))
}
