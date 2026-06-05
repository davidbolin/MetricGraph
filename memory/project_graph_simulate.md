---
name: project-graph-simulate
description: Status of the bridge simulator C++ acceleration phase (simulation.md phase 3)
metadata: 
  node_type: memory
  type: project
  originSessionId: ea33ddff-e850-48ca-9932-3ca088415e18
---

Phase complete as of 2026-06-04. C++ acceleration, extended method, and benchmark complete.

**Why:** Implement the C++ acceleration plan from `jonas_local/simulation.md` (phase 3).

**What was done (this session):**
- `src/edge_simulation.cpp` created: `draw_edge_direct_cpp` and `draw_edge_kriging_cpp` (Rcpp/Eigen, all 4 cases: direct/kriging × alpha=1/2). Uses `Eigen::LLT`, `R::norm_rand()`, unique-lag cache for kriging alpha=2. Matches R output to machine precision (< 1e-14 max diff).
- `Rcpp::compileAttributes()` run → `src/RcppExports.cpp` and `R/RcppExports.R` updated.
- `R/graph_simulate.R` updated: added `impl = c("cpp","R")` parameter (cpp is default), added `method="extended"` (`.simulate_extended_wm()` helper: inserts PtE as graph vertices, single sparse Cholesky draw for alpha=1; CoB-transformed precision + vertex value extraction for alpha=2).
- `tests/testthat/test-graph-simulate.R`: 33 tests, all pass (R vs C++ equivalence, extended method shape + distributional tests, impl validation).
- `examples/fast_simulation/run_study.R` updated: warm-up iteration, n_rep=5 medians, n_pts ∈ {8,...,2048}, all 3 methods.
- `examples/fast_simulation/benchmark.R` created: per-edge R vs C++ speedup (m ∈ {8,...,2048}), whole-field comparison (n_pts ∈ {8,...,512}).
- `NEWS.md` updated, docs regenerated via `devtools::document()`.

**Pre-existing test failures (unrelated):** `test_cross_validation.R` and `test_fem_predict_manual.R` have 29 failures with snapshot files in `tests/testthat/_problems/`. These predate the bridge simulation work and are NOT introduced by it.

**How to apply:** The package is in a working state. To run the study: `Rscript examples/fast_simulation/run_study.R`. To run the benchmark: `Rscript examples/fast_simulation/benchmark.R`.

[[feedback-no-confirmation]]
