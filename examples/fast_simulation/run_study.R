# run_study.R  —  Chicago street-network simulation study (package-only driver)
#
# Reproduces Table 4 and Figure 1 from the paper using the MetricGraph package.
# Dependencies: MetricGraph (>= dev), spatstat.data, ggplot2, Matrix.
#
# Run from the package root:
#   Rscript examples/fast_simulation/run_study.R

library(MetricGraph)
library(Matrix)
library(spatstat.data)

# ── Script directory (for output paths) ────────────────────────────────────
.sd <- local({
  args <- commandArgs(trailingOnly = FALSE)
  m    <- regmatches(args, regexpr("(?<=--file=).*", args, perl = TRUE))
  if (length(m)) dirname(normalizePath(m)) else "."
})

# ============================================================
# 1. Build the Chicago street network
# ============================================================
data(chicago)
g <- linnet.to.graph(chicago$domain, crs = sf::st_crs(NA))
cat(sprintf("Chicago graph: %d vertices, %d edges\n", g$nV, g$nE))
cat(sprintf("Edge lengths: min=%.1f  mean=%.1f  max=%.1f  total=%.0f\n",
            min(g$edge_lengths), mean(g$edge_lengths),
            max(g$edge_lengths), sum(g$edge_lengths)))

# ── Model parameters: practical range = 10 % of network length, sigma = 1 ──
.rng       <- 0.1 * sum(g$edge_lengths)
range_val  <- .rng
sigma_val  <- 1

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
n_rep       <- 5L         # timed repetitions per condition (warm-up excluded)
master_seed <- 42L
BC          <- 1L
methods     <- c("direct", "kriging", "extended")
alphas      <- c(1L, 2L)
cat(sprintf("\nn_rep: %d\n", n_rep))
cat(sprintf("Edges: %d,  n_pts sweep: %s\n", g$nE, paste(n_pts_vals, collapse = " ")))

# Helper: build PtE for a given n_pts
make_PtE <- function(n_pts) {
  t_norm <- seq(0, 1, length.out = n_pts + 2L)[-c(1L, n_pts + 2L)]
  do.call(rbind, lapply(seq_len(g$nE), function(e) cbind(e, t_norm)))
}

results <- list()
cond_i  <- 0L

not_pd <- list()   # tracks first n_pts where each method/alpha went non-PD

for (alpha in alphas) {
  kappa <- if (alpha == 1L) kappa1 else kappa2
  tau   <- if (alpha == 1L) tau1   else tau2

  for (meth in methods) {
    key <- paste(meth, alpha, sep = "_")
    for (n_pts in n_pts_vals) {
      cond_i <- cond_i + 1L
      PtE    <- make_PtE(n_pts)

      cat(sprintf("  method=%-8s alpha=%d n_pts=%4d", meth, alpha, n_pts))
      flush.console()

      # If a smaller n_pts already failed, skip and record NA
      if (isTRUE(not_pd[[key]])) {
        cat(sprintf("  median=NA (not PD, skipped)  n_loc=%d\n", g$nE * n_pts))
        results[[cond_i]] <- data.frame(
          method = meth, alpha = alpha, n_pts = n_pts,
          n_loc = g$nE * n_pts, median_ms = NA_real_, iqr_ms = NA_real_,
          stringsAsFactors = FALSE
        )
        next
      }

      # Warm-up: one untimed run
      simulate(g, alpha = alpha, method = meth,
               kappa = kappa, tau = tau, PtE = PtE, BC = BC,
               seed = master_seed)

      # Timed repetitions (all serial)
      # Extended: one nsim=n_rep batch divided by n_rep (amortised graph setup)
      # Direct/kriging: n_rep individual calls, report median
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

timing <- do.call(rbind, results)
write.csv(timing, file.path(.sd, "study_timing.csv"), row.names = FALSE)
saveRDS(timing,   file.path(.sd, "study_timing.rds"))
cat("Saved study_timing.csv and study_timing.rds\n")

# ── Format and print Table 4 (article n_pts subset) ──────────────────────────
timing$label <- sprintf("%-8s  alpha=%d",
                        ifelse(timing$method == "direct",   "Method A",
                        ifelse(timing$method == "kriging",  "Method B", "Extended")),
                        timing$alpha)
# Article columns: n_pts in {32,64,128,256,512,1024}
art_pts    <- c(32L, 64L, 128L, 256L, 512L, 1024L)
art_timing <- timing[timing$n_pts %in% art_pts, ]
n_locs_art <- sort(unique(art_timing$n_loc))
col_names  <- formatC(n_locs_art, format = "d", big.mark = ",")

make_row_label <- function(method, alpha)
  sprintf("%-8s  alpha=%d",
          ifelse(method == "direct", "Method A",
          ifelse(method == "kriging", "Method B", "Extended")),
          alpha)

row_labels <- make_row_label(
  rep(c("direct", "kriging", "extended"), each = 2),
  rep(c(1L, 2L), 3)
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
