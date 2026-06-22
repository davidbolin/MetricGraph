# examples/fast_simulation — third-party spectral simulator

This folder contains an **example/benchmark study only**. None of the code here
is part of the `MetricGraph` package build — nothing in this directory is sourced
from `R/` or `src/`, and it is not exported, documented, or shipped on CRAN. It
exists solely to run the timing/realization study in `run_study*.R`.

## Provenance of the ported code

The spectral simulator added to `run_study_temp.R` (Algorithm 1, exponential
covariance) is **ported from a third-party reference implementation**, not written
from scratch:

- **Source repository:** FastSimNetworks — <https://github.com/alfredoalegria/FastSimNetworks>
- **Paper:** Alfredo Alegría, Xavier Emery, Tobia Filosi & Emilio Porcu (2026),
  *Computationally Efficient Algorithms for Simulating Isotropic Gaussian Random
  Fields on Graphs with Euclidean Edges*, **Journal of Computational and Graphical
  Statistics** 35(2), 937–950. DOI:
  [10.1080/10618600.2025.2574535](https://doi.org/10.1080/10618600.2025.2574535).
  Open Access (CC BY 4.0).
- **Algorithm:** Algorithm 1, the spectral method (paper Sections 3.1 & 4.2).

### Files taken from FastSimNetworks

| Local file (here) | Upstream file | Notes |
|-------------------|---------------|-------|
| `BB.c` | `basic_functions/c/BB.c` | Brownian-bridge sampler, copied verbatim; compiled locally with `R CMD SHLIB BB.c`. |
| (inlined in `run_study_temp.R`) | `basic_functions/r/BB.R` | `.C("BB", …)` interface wrapper. |
| (inlined in `run_study_temp.R`) | `basic_functions/functions.R` | only `genSites()` and `laplacian()`. |
| (inlined in `run_study_temp.R`) | `sim_algorithms/auxiliary_process.R` | `SimAux()` — auxiliary field `Z_G`. |
| (inlined in `run_study_temp.R`) | `sim_algorithms/SpectralSim.R` | `SimSpec()` — Algorithm 1. |

Upstream Algorithms 2/3 (`Dilution1Sim.R`, `Dilution2Sim.R`) and the
variogram/resistance-metric diagnostics are **not** ported — only the
exponential-covariance spectral method is needed here.

The upstream code is released under the terms of that repository; this port is
included here for reproducing the paper's study within `MetricGraph`'s examples.
Please cite the paper above when using these results.

## What is original to this repo

`run_study.R` and the MetricGraph-based methods (`direct`, `kriging`, `extended`)
are part of this project. `run_study_temp.R` is `run_study.R` plus the ported
spectral row.

## Files

| File | Description |
|------|-------------|
| `run_study.R` | Existing study: `direct`, `kriging`, `extended` (Whittle–Matérn). |
| `run_study_temp.R` | *(to be added)* copy of `run_study.R` + ported spectral method. |
| `BB.c` | *(to be added)* upstream Brownian-bridge C routine. |
| `README.md` | This file. |

See `../../jonas_local/update_simulate.md` for the full porting guide and the
mathematical background.
