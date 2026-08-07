# run_study_temp.R  —  Chicago street-network simulation study + spectral method
#
# Extends run_study.R with Algorithm 1 (spectral / SimSpec) from:
#   Alegría, Emery, Filosi & Porcu (2026), JCGS 35:2, 937-950.
#   https://github.com/alfredoalegria/FastSimNetworks
#
# Run from the package root:
#   Rscript examples/fast_simulation/run_study_temp.R

library(MetricGraph)
library(Matrix)
library(spatstat.data)
library(spatstat.geom)

# ── Script directory (for output paths and BB.c compilation) ───────────────
.sd <- local({
  args <- commandArgs(trailingOnly = FALSE)
  m    <- regmatches(args, regexpr("(?<=--file=).*", args, perl = TRUE))
  if (length(m)) dirname(normalizePath(m)) else "."
})

# ── Compile and load BB.c (Brownian-bridge C routine) ──────────────────────
.bb_so <- file.path(.sd, paste0("BB", .Platform$dynlib.ext))
if (!file.exists(.bb_so)) {
  message("Compiling BB.c ...")
  system2("R", c("CMD", "SHLIB", "-o", .bb_so, file.path(.sd, "BB.c")))
}
dyn.load(.bb_so)

# ── Inlined helpers from FastSimNetworks ────────────────────────────────────

BB <- function(nEdges, lenEdges, nPointsE, vecNorm) {
  storage.mode(lenEdges) <- "double"
  storage.mode(nPointsE) <- "integer"
  storage.mode(vecNorm)  <- "double"
  num <- sum(nPointsE)
  o <- .C("BB", nEdges = as.integer(nEdges), lenEdges = lenEdges,
                 nPointsE = nPointsE, vecNorm = vecNorm, vec = double(num))
  o$vec
}

genSites <- function(nPointsE) {
  e <- rep(seq_along(nPointsE), nPointsE)
  t <- (sequence(nPointsE) - 1) / rep(nPointsE - 1, nPointsE)
  cbind(t, e)
}

laplacian <- function(adjM, vDist) {
  c <- adjM / vDist
  diag(c) <- 0
  c0 <- rowSums(c)
  c0[1] <- c0[1] + 1
  -c + diag(c0)
}

SimAux <- function(lenEdges, nEdges, nPointsE, nsites, nVert,
                   cholILM, from, to, sites) {
  bridges  <- BB(nEdges, lenEdges, nPointsE, rnorm(sum(nPointsE)))
  zVert    <- t(cholILM) %*% rnorm(nVert)
  zInterp  <- (1 - sites[, 1]) * zVert[from[sites[, 2]]] +
                   sites[, 1]  * zVert[to[sites[, 2]]]
  zAux     <- bridges + zInterp
  list(bridges, zVert, zInterp, zAux)
}

SimSpec <- function(L, nPointsE, nCopies, scaleParam) {
  cholILM <- chol(solve(laplacian(L$m, L$dpath)))
  sites   <- genSites(nPointsE)
  nsites  <- nrow(sites)
  nEdges  <- nsegments(L$lines)
  zSim    <- rep(0, nsites)
  for (m in seq_len(nCopies)) {
    V    <- scaleParam
    U    <- 2 * pi * runif(1)
    W    <- runif(1)
    zAux <- SimAux(lengths_psp(L$lines), nEdges, nPointsE,
                   nsites, nrow(cholILM), cholILM, L$from, L$to, sites)[[4]]
    zSim <- zSim + sqrt(-log(W)) * cos(V * zAux + U)
  }
  data.frame(e = sites[, 2], t = sites[, 1],
             values = zSim * sqrt(2 / nCopies))
}

# ============================================================
# 1. Build the Chicago street network
# ============================================================
data(chicago)
g   <- linnet.to.graph(chicago$domain, crs = sf::st_crs(NA))
Net <- chicago$domain  # already has dense $m and $dpath

stopifnot(nsegments(Net$lines) == g$nE)
cat(sprintf("Chicago graph: %d vertices, %d edges\n", g$nV, g$nE))
cat(sprintf("Edge lengths: min=%.1f  mean=%.1f  max=%.1f  total=%.0f\n",
            min(g$edge_lengths), mean(g$edge_lengths),
            max(g$edge_lengths), sum(g$edge_lengths)))

# ── Model parameters: practical range = 10 % of network length, sigma = 1 ──
range_val <- 0.1 * sum(g$edge_lengths)
sigma_val <- 1

.kt <- function(alpha, range, sigma) {
  nu    <- alpha - 0.5
  kappa <- sqrt(8 * nu) / range
  tau   <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
                               (4 * pi)^0.5 * gamma(nu + 0.5)))
  list(kappa = kappa, tau = tau)
}
p1 <- .kt(1L, range_val, sigma_val); kappa1 <- p1$kappa; tau1 <- p1$tau
p2 <- .kt(2L, range_val, sigma_val); kappa2 <- p2$kappa; tau2 <- p2$tau
cat(sprintf("Parameters: range=%.0f  sigma=%.1f\n", range_val, sigma_val))
cat(sprintf("  alpha=1: kappa=%.6f  tau=%.4f\n", kappa1, tau1))
cat(sprintf("  alpha=2: kappa=%.6f  tau=%.4f\n", kappa2, tau2))

# ============================================================
# 2. Timing study  (Table 4)
# ============================================================
n_pts_vals  <- c(8L, 16L, 32L, 64L, 128L, 256L, 512L, 1024L)
n_rep       <- 5L
master_seed <- 42L
BC          <- 1L
methods     <- c("direct", "kriging", "extended", "spectral")
alphas      <- c(1L, 2L)
spectral_scale_param <- 0.2
isotropic_exact_pts  <- head(n_pts_vals, 3L)
cat(sprintf("\nn_rep: %d\n", n_rep))
cat(sprintf("Edges: %d,  n_pts sweep: %s\n", g$nE, paste(n_pts_vals, collapse = " ")))

make_PtE <- function(n_pts) {
  t_norm <- seq(0, 1, length.out = n_pts + 2L)[-c(1L, n_pts + 2L)]
  do.call(rbind, lapply(seq_len(g$nE), function(e) cbind(e, t_norm)))
}

simulate_isotropic_exact <- function(PtE, scale_param = spectral_scale_param) {
  D <- g$compute_resdist_PtE(PtE, normalized = TRUE)
  Sigma <- Matrix::forceSymmetric(exp(-0.5 * scale_param^2 * D))
  R <- tryCatch(
    chol(Sigma),
    error = function(e) {
      stop(
        "Exact isotropic covariance Cholesky failed: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  as.numeric(t(R) %*% rnorm(nrow(PtE)))
}

results <- list()
cond_i  <- 0L
not_pd  <- list()

for (alpha in alphas) {
  kappa <- if (alpha == 1L) kappa1 else kappa2
  tau   <- if (alpha == 1L) tau1   else tau2

  for (meth in methods) {
    # spectral only supports exponential covariance (alpha = 1)
    if (meth == "spectral" && alpha == 2L) next

    key <- paste(meth, alpha, sep = "_")

    for (n_pts in n_pts_vals) {
      cond_i <- cond_i + 1L

      cat(sprintf("  method=%-8s alpha=%d n_pts=%4d", meth, alpha, n_pts))
      flush.console()

      if (meth == "spectral") {
        nPointsE_s <- rep(n_pts, nsegments(Net$lines))

        # warm-up
        SimSpec(Net, nPointsE_s, nCopies = 1000L,
                scaleParam = spectral_scale_param)

        times_ms <- numeric(n_rep)
        for (r in seq_len(n_rep)) {
          t0 <- proc.time()[3]
          SimSpec(Net, nPointsE_s, nCopies = 1000L,
                  scaleParam = spectral_scale_param)
          times_ms[r] <- (proc.time()[3] - t0) * 1000
        }
        med_ms <- median(times_ms)
        iqr_ms <- IQR(times_ms)

        cat(sprintf("  median=%.1f ms  IQR=%.1f ms  n_loc=%d\n",
                    med_ms, iqr_ms, g$nE * n_pts))

        results[[cond_i]] <- data.frame(
          method    = meth,
          alpha     = alpha,
          n_pts     = n_pts,
          n_loc     = g$nE * n_pts,
          median_ms = med_ms,
          iqr_ms    = iqr_ms,
          stringsAsFactors = FALSE
        )
        next
      }

      # ── existing methods ─────────────────────────────────────────────────
      PtE <- make_PtE(n_pts)

      if (isTRUE(not_pd[[key]])) {
        cat(sprintf("  median=NA (not PD, skipped)  n_loc=%d\n", g$nE * n_pts))
        results[[cond_i]] <- data.frame(
          method = meth, alpha = alpha, n_pts = n_pts,
          n_loc = g$nE * n_pts, median_ms = NA_real_, iqr_ms = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }

      simulate(g, alpha = alpha, method = meth,
               kappa = kappa, tau = tau, PtE = PtE, BC = BC,
               seed = master_seed)

      if (meth == "extended") {
        t0  <- proc.time()[3]
        res <- simulate(g, nsim = n_rep, alpha = alpha, method = meth,
                        kappa = kappa, tau = tau, PtE = PtE, BC = BC)
        if (anyNA(res)) {
          not_pd[[key]] <- TRUE
          med_ms <- NA_real_
        } else {
          med_ms <- (proc.time()[3] - t0) * 1000 / n_rep
        }
        iqr_ms <- 0
      } else {
        times_ms <- numeric(n_rep)
        for (r in seq_len(n_rep)) {
          t0 <- proc.time()[3]
          simulate(g, alpha = alpha, method = meth,
                   kappa = kappa, tau = tau, PtE = PtE, BC = BC,
                   seed = master_seed + r)
          times_ms[r] <- (proc.time()[3] - t0) * 1000
        }
        med_ms <- median(times_ms)
        iqr_ms <- IQR(times_ms)
      }

      cat(sprintf("  %s  IQR=%.1f ms  n_loc=%d\n",
                  if (is.na(med_ms)) "median=NA (not PD)" else sprintf("median=%.1f ms", med_ms),
                  iqr_ms, g$nE * n_pts))

      results[[cond_i]] <- data.frame(
        method    = meth,
        alpha     = alpha,
        n_pts     = n_pts,
        n_loc     = g$nE * n_pts,
        median_ms = med_ms,
        iqr_ms    = iqr_ms,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Exact isotropic exponential covariance benchmark. This is kept separate from
# the alpha loop because it is not a Whittle--Matern model.
for (n_pts in isotropic_exact_pts) {
  cond_i <- cond_i + 1L
  PtE <- make_PtE(n_pts)

  cat(sprintf("  method=%-15s n_pts=%4d", "isotropic_exact", n_pts))
  flush.console()

  set.seed(master_seed)
  warmup <- simulate_isotropic_exact(PtE)
  stopifnot(length(warmup) == nrow(PtE), all(is.finite(warmup)))
  rm(warmup)

  times_ms <- numeric(n_rep)
  for (r in seq_len(n_rep)) {
    gc()
    set.seed(master_seed + r)
    t0 <- proc.time()[3]
    u_iso <- simulate_isotropic_exact(PtE)
    times_ms[r] <- (proc.time()[3] - t0) * 1000
    stopifnot(length(u_iso) == nrow(PtE), all(is.finite(u_iso)))
    rm(u_iso)
  }

  med_ms <- median(times_ms)
  iqr_ms <- IQR(times_ms)
  cat(sprintf("  median=%.1f ms  IQR=%.1f ms  n_loc=%d\n",
              med_ms, iqr_ms, nrow(PtE)))

  results[[cond_i]] <- data.frame(
    method    = "isotropic_exact",
    alpha     = NA_integer_,
    n_pts     = n_pts,
    n_loc     = nrow(PtE),
    median_ms = med_ms,
    iqr_ms    = iqr_ms,
    stringsAsFactors = FALSE
  )
}

timing <- do.call(rbind, results)
write.csv(timing, file.path(.sd, "study_timing.csv"), row.names = FALSE)
saveRDS(timing,   file.path(.sd, "study_timing.rds"))
cat("Saved study_timing.csv and study_timing.rds\n")

# ── Format and print Table 4 ─────────────────────────────────────────────────
make_row_label <- function(method, alpha)
  sprintf("%-8s  alpha=%d",
          ifelse(method == "direct",   "Method A",
          ifelse(method == "kriging",  "Method B",
          ifelse(method == "extended", "Extended", "Spectral"))),
          alpha)

art_pts    <- c(32L, 64L, 128L, 256L, 512L, 1024L)
art_timing <- timing[timing$n_pts %in% art_pts, ]
n_locs_art <- sort(unique(art_timing$n_loc))
col_names  <- formatC(n_locs_art, format = "d", big.mark = ",")

row_labels <- c(
  make_row_label(rep(c("direct", "kriging", "extended"), each = 2), rep(c(1L, 2L), 3)),
  make_row_label("spectral", 1L)
)
art_timing$row_label <- make_row_label(art_timing$method, art_timing$alpha)

mat <- matrix(NA_real_, nrow = length(row_labels), ncol = length(n_locs_art),
              dimnames = list(row_labels, col_names))
for (i in seq_len(nrow(art_timing))) {
  r  <- art_timing$row_label[i]
  cc <- formatC(art_timing$n_loc[i], format = "d", big.mark = ",")
  if (r %in% row_labels && cc %in% col_names && !is.na(art_timing$median_ms[i]))
    mat[r, cc] <- art_timing$median_ms[i] / 1000
}

cat("\n=== Timing table: median wall-clock time (seconds), n_rep =", n_rep, "===\n")
cat(sprintf("Graph: Chicago (%d edges)\n\n", g$nE))
cat(sprintf("%-22s", ""))
cat(paste(sprintf("%9s", col_names), collapse = ""), "\n")
cat(strrep("-", 22 + 9 * length(n_locs_art)), "\n")
for (r in row_labels) {
  cat(sprintf("%-22s", r))
  cat(paste(sprintf("%9s",
    ifelse(is.na(mat[r, ]), "  —  ",
           ifelse(mat[r, ] < 1, sprintf("%.3f", mat[r, ]),
                               sprintf("%.2f",  mat[r, ])))),
    collapse = ""), "\n")
}

iso_timing <- timing[timing$method == "isotropic_exact", ]
iso_col_names <- formatC(iso_timing$n_loc, format = "d", big.mark = ",")
cat("\n=== Exact isotropic exponential: end-to-end median time (seconds) ===\n\n")
cat(sprintf("%-22s", ""))
cat(paste(sprintf("%9s", iso_col_names), collapse = ""), "\n")
cat(strrep("-", 22 + 9 * nrow(iso_timing)), "\n")
cat(sprintf("%-22s", "Isotropic exact"))
cat(paste(sprintf("%9s",
                  ifelse(iso_timing$median_ms < 1000,
                         sprintf("%.3f", iso_timing$median_ms / 1000),
                         sprintf("%.2f", iso_timing$median_ms / 1000))),
          collapse = ""), "\n")

wide <- as.data.frame(mat)
wide <- cbind(method_alpha = row_labels, wide)
write.csv(wide, file.path(.sd, "study_table.csv"), row.names = FALSE)
cat("Saved study_table.csv\n")

# ============================================================
# 3. Realizations figure  (Figure 1)
# ============================================================
n_pts_fig <- 200L
t_norm    <- seq(0, 1, length.out = n_pts_fig + 2L)[-c(1L, n_pts_fig + 2L)]
PtE_fig   <- do.call(rbind,
                     lapply(seq_len(g$nE), function(e) cbind(e, t_norm)))

for (alpha in c(1L, 2L)) {
  kappa <- if (alpha == 1L) kappa1 else kappa2
  tau   <- if (alpha == 1L) tau1   else tau2

  u <- simulate(g, alpha = alpha, method = "kriging",
                kappa = kappa, tau = tau, PtE = PtE_fig, seed = alpha * 100L)
  cat(sprintf("alpha=%d: n_loc=%d  range=[%.2f,%.2f]\n",
              alpha, length(u), min(u), max(u)))

  df <- data.frame(edge_number      = PtE_fig[, 1],
                   distance_on_edge = PtE_fig[, 2], u = u)
  sc <- ggplot2::scale_color_viridis_c(option = "D", limits = range(u))
  gr <- g$clone()
  gr$add_observations(data = df, normalized = TRUE, suppress_warnings = TRUE)
  p <- gr$plot_function(data = "u", scale_color = sc, vertex_size = 0) +
         ggplot2::theme_void(base_size = 10) +
         ggplot2::guides(color = ggplot2::guide_colorbar(title = NULL))

  out_file <- file.path(.sd, sprintf("chicago_field_alpha%d.png", alpha))
  png(out_file, width = 600, height = 550, res = 120)
  print(p)
  dev.off()
  cat(sprintf("Saved chicago_field_alpha%d.png\n", alpha))
}

# ── Spectral realization panel ──────────────────────────────────────────────
# n_pts_fig + 2 so that after dropping the t=0/1 endpoints we keep n_pts_fig
# interior points per edge (matching the kriging panels above).
nPointsE_fig_s <- rep(n_pts_fig + 2L, nsegments(Net$lines))
spec_df <- SimSpec(Net, nPointsE_fig_s, nCopies = 1000L,
                   scaleParam = spectral_scale_param)
spec_df <- spec_df[spec_df$t > 0 & spec_df$t < 1, ]
cat(sprintf("spectral: n_loc=%d  range=[%.2f,%.2f]\n",
            nrow(spec_df), min(spec_df$values), max(spec_df$values)))

df_spec <- data.frame(edge_number      = spec_df$e,
                      distance_on_edge = spec_df$t,
                      u                = spec_df$values)
sc_spec <- ggplot2::scale_color_viridis_c(option = "D", limits = range(df_spec$u))
gr_spec <- g$clone()
gr_spec$add_observations(data = df_spec, normalized = TRUE, suppress_warnings = TRUE)
p_spec <- gr_spec$plot_function(data = "u", scale_color = sc_spec, vertex_size = 0) +
            ggplot2::theme_void(base_size = 10) +
            ggplot2::guides(color = ggplot2::guide_colorbar(title = NULL))

out_file_spec <- file.path(.sd, "chicago_field_spectral.png")
png(out_file_spec, width = 600, height = 550, res = 120)
print(p_spec)
dev.off()
cat(sprintf("Saved chicago_field_spectral.png\n"))
