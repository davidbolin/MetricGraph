#!/usr/bin/env Rscript
## run_benchmarks.R
## Main entry point for the MetricGraph benchmark suite.
##
## Usage:
##   Rscript benchmarks/run_benchmarks.R            # Default: all tests at full sizes,
##                                                   #   excluding slow models (isoexp, GL1)
##                                                   #   and compute_resdist
##   Rscript benchmarks/run_benchmarks.R --quick     # Quick mode (small sizes, fast models only)
##   Rscript benchmarks/run_benchmarks.R --slow      # Include slow tests (isoexp, GL1, resdist)
##   Rscript benchmarks/run_benchmarks.R --save      # Save results to RDS

args <- commandArgs(trailingOnly = TRUE)
quick <- "--quick" %in% args
save_results <- "--save" %in% args
include_slow <- "--slow" %in% args

# Source all helpers and benchmark scripts
# Determine script directory robustly across source() and Rscript
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    # Fallback for Rscript: use --file arg or default to "benchmarks"
    cmd_args <- commandArgs(FALSE)
    file_arg <- grep("^--file=", cmd_args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      "benchmarks"
    }
  }
)
source(file.path(script_dir, "benchmark_helpers.R"))
source(file.path(script_dir, "bench_graph_construction.R"))
source(file.path(script_dir, "bench_mesh_and_fem.R"))
source(file.path(script_dir, "bench_data_operations.R"))
source(file.path(script_dir, "bench_model_fitting.R"))
source(file.path(script_dir, "bench_prediction.R"))
source(file.path(script_dir, "bench_plotting.R"))
source(file.path(script_dir, "bench_spde.R"))
source(file.path(script_dir, "bench_graph_spde.R"))

# Print header
print_system_info()
if (quick) {
  cat("  Mode: QUICK (reduced sizes, fast models only)\n\n")
} else if (include_slow) {
  cat("  Mode: FULL + SLOW (all sizes, all models incl. isoexp/GL1/resdist)\n\n")
} else {
  cat("  Mode: DEFAULT (all sizes, excluding isoexp/GL1/resdist)\n\n")
}

set.seed(12345)
all_results <- list()
total_t0 <- proc.time()

# Run each benchmark section
all_results$construction <- bench_graph_construction(quick = quick)
all_results$mesh_fem     <- bench_mesh_and_fem(quick = quick)
all_results$data_ops     <- bench_data_operations(quick = quick,
                                                    include_slow = include_slow)
all_results$model_fit    <- bench_model_fitting(quick = quick,
                                                   include_slow = include_slow)
all_results$prediction   <- bench_prediction(quick = quick,
                                                include_slow = include_slow)
all_results$plotting     <- bench_plotting(quick = quick)
all_results$spde         <- bench_spde(quick = quick)
all_results$graph_spde   <- bench_graph_spde(quick = quick)

total_elapsed <- (proc.time() - total_t0)[["elapsed"]]

cat(strrep("=", 70), "\n")
cat(sprintf(" Total benchmark time: %.1f seconds\n", total_elapsed))
cat(strrep("=", 70), "\n")

# Optionally save results
if (save_results) {
  out_file <- file.path(script_dir,
    sprintf("benchmark_results_%s.rds", format(Sys.Date(), "%Y%m%d")))
  saveRDS(all_results, out_file)
  cat(sprintf("\nResults saved to: %s\n", out_file))
}
