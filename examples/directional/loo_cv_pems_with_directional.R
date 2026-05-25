## --------------------------------------------------------------------------
## Leave-one-out cross-validation comparison on the pems data, extended
## with several variants of the directional Whittle--Matern model (WMD1).
##
## The non-directional baselines reproduce the model comparison from the
## `comparison` vignette (vignettes/comparison.rmd):
##   - WM alpha=1 (BC=0 and BC=1)
##   - WM alpha=2 (BC=0 and BC=1)
##   - graph Laplacian (alpha=1, alpha=2)
##   - isotropic exponential covariance
##
## On top of those we fit seven directional WMD1 variants that vary in two
## axes -- which graph the model sees, and which vertex condition it uses:
##
##   1. WMD1                                raw graph,        f_in = w/sum(w)
##   2. WMD1 var-stationary                 raw graph,        f_in = sqrt(w/sum(w))
##   3. WMD1 dir-corrected                  corrected graph,  f_in = w/sum(w)
##   4. WMD1 dir-corrected2                 corrected2 graph, f_in = w/sum(w)
##   5. WMD1 dir-corrected2 len-weighted    corrected2 graph, weights = edge length
##   6. WMD1 dir-corrected2 type-weighted   corrected2 graph, weights = 5 (long) / 1 (short)
##   7. WMD1 dir-corrected2 osm-weighted    corrected2 graph, weights from OSM highway class
##
## The graph variants are:
##   * `pems_graph`        : raw pems edges, mostly directionally consistent
##                           but with 16 pass-through vertices having two
##                           inbound or two outbound edges (artefacts of
##                           the source linestring coordinate order).
##   * `pems_graph_corr`   : pass-through inversions resolved by BFS
##                           propagation along degree-2 chains (22 edges
##                           flipped).
##   * `pems_graph_corr2`  : additional junction-edge flips chosen by
##                           per-edge log-likelihood evaluation at three
##                           topologically suspect deg>=3 vertices.
##
## The vertex condition is determined by the directional weight function
## f_in (controlled by graph$setDirectionalWeightFunction) and by the
## per-edge directional weights (controlled by graph$set_edge_weights).
## With unit weights the default f_in = w/sum(w) is the equal-weight
## average of incoming edge values; changing the weights changes the
## linear combination at each merge without touching f_in.
##
## The OSM-weighted variant requires the per-edge classification stored
## in `pems_osm_highway.rds` (created by `pems_fetch_osm_highway.R`); if
## that file isn't present, the script skips that fit.
## --------------------------------------------------------------------------

rm(list = ls())
library(MetricGraph)
set.seed(1)

## ---- 1. Build the pems graph and attach observations ---------------------

pems_graph <- metric_graph$new(edges = pems$edges, verbose = 0)
pems_graph$add_observations(data = pems$data, normalized = TRUE,
                            verbose = 0)
cat(sprintf("pems graph: %d edges, %d vertices, %d observations.\n",
            pems_graph$nE, pems_graph$nV,
            length(pems_graph$get_data()[[".group"]])))

## ---- 2. Fit the Whittle--Matern models (symmetric) -----------------------

cat("Fitting WM alpha=1 (BC=0)...\n")
fit_alpha1 <- graph_lme(y ~ 1, graph = pems_graph, BC = 0,
                        model = list(type = "WhittleMatern", alpha = 1),
                        optim_method = "Nelder-Mead")

cat("Fitting WM alpha=2 (BC=0)...\n")
fit_alpha2 <- graph_lme(y ~ 1, graph = pems_graph, BC = 0,
                        model = list(type = "WhittleMatern", alpha = 2),
                        optim_method = "Nelder-Mead")

cat("Fitting WM alpha=1 (BC=1)...\n")
fit_alpha1_bc <- graph_lme(y ~ 1, graph = pems_graph, BC = 1,
                           model = list(type = "WhittleMatern", alpha = 1),
                           optim_method = "Nelder-Mead")

cat("Fitting WM alpha=2 (BC=1)...\n")
fit_alpha2_bc <- graph_lme(y ~ 1, graph = pems_graph, BC = 1,
                           model = list(type = "WhittleMatern", alpha = 2),
                           optim_method = "Nelder-Mead")

## ---- 3. Fit the directional Whittle--Matern models (WMD1) ----------------
## We fit WMD1 with two different vertex conditions on the raw pems
## graph. The vertex condition is controlled by the graph's
## DirectionalWeightFunction_in (rebuilt from
## graph$setDirectionalWeightFunction). To keep the two fits cleanly
## separated, we clone the graph for the variance-stationary variant.

cat("Fitting WMD alpha=1 (default vertex condition)...\n")
fit_wmd1 <- graph_lme(y ~ 1, graph = pems_graph, BC = 0,
                      model = "WMD1",
                      optim_method = "Nelder-Mead")

cat("Fitting WMD alpha=1 (variance-stationary vertex condition)...\n")
pems_graph_vs <- pems_graph$clone()
pems_graph_vs$setDirectionalWeightFunction(
  f_in = function(x) sqrt(x / sum(x)))
fit_wmd1_vs <- graph_lme(y ~ 1, graph = pems_graph_vs, BC = 0,
                         model = "WMD1",
                         optim_method = "Nelder-Mead")

## ---- 3b. Direction-corrected pems graph ----------------------------------
## The raw pems edges are mostly directionally consistent (most pass-
## through vertices have one inbound and one outbound edge), but a small
## number have local inversions: two segments meeting at a pass-through
## vertex with both pointing in (or both pointing out). For a directed
## model, these inversions misrepresent the actual flow of traffic. We
## fix them by BFS-propagating orientation along each degree-2 chain,
## flipping the coordinate order of any edge that contradicts the
## propagated direction. At junctions (degree >= 3) we do nothing -- the
## directional model allows arbitrary in/out splits there.

propagate_orientation_flips <- function(g) {
  E   <- g$E
  nE  <- g$nE
  nV  <- g$nV
  adj <- vector("list", nV)
  for (e in seq_len(nE)) {
    adj[[E[e, 1]]] <- c(adj[[E[e, 1]]], e)
    adj[[E[e, 2]]] <- c(adj[[E[e, 2]]], e)
  }
  deg <- lengths(adj)

  flip    <- logical(nE)
  visited <- logical(nE)
  orient  <- function(e) if (flip[e]) c(E[e, 2], E[e, 1]) else E[e, ]

  for (seed in seq_len(nE)) {
    if (visited[seed]) next
    q <- seed
    while (length(q) > 0L) {
      cur <- q[1L]; q <- q[-1L]
      if (visited[cur]) next
      visited[cur] <- TRUE
      o <- orient(cur)
      for (v in o) {
        if (deg[v] != 2L) next
        other <- setdiff(adj[[v]], cur)
        if (length(other) != 1L || visited[other]) next
        cur_in_at_v   <- (v == o[2L])               # cur enters v?
        o_other       <- orient(other)
        other_in_at_v <- (v == o_other[2L])         # other enters v?
        ## Consistent at v means exactly one in and one out: opposite
        ## values. Same values => both in or both out => flip other.
        if (cur_in_at_v == other_in_at_v) {
          flip[other] <- !flip[other]
        }
        q <- c(q, other)
      }
    }
  }
  flip
}

direction_diagnostic <- function(g, label) {
  in_deg  <- g$get_degrees("indegree")
  out_deg <- g$get_degrees("outdegree")
  deg <- in_deg + out_deg
  consistent <- sum(deg == 2 & in_deg == 1 & out_deg == 1)
  inverted   <- sum(deg == 2 & ((in_deg == 2 & out_deg == 0) |
                                (in_deg == 0 & out_deg == 2)))
  cat(sprintf(
    "  %s: %d pass-through vertices (%d consistent, %d inverted).\n",
    label, consistent + inverted, consistent, inverted))
}

cat("Direction consistency on the raw pems graph:\n")
direction_diagnostic(pems_graph, "raw")

flip_vec <- propagate_orientation_flips(pems_graph)
cat(sprintf("  Edges to flip: %d / %d (%.1f%%).\n",
            sum(flip_vec), pems_graph$nE,
            100 * sum(flip_vec) / pems_graph$nE))

## Reverse the coordinate matrix of every flipped edge and rebuild a
## fresh metric_graph from the corrected edge list.
edge_list_corr <- lapply(seq_along(pems_graph$edges), function(i) {
  m <- pems_graph$edges[[i]]
  class(m) <- NULL
  m <- unname(m)
  if (flip_vec[i]) m[nrow(m):1, , drop = FALSE] else m
})
class(edge_list_corr) <- NULL

pems_graph_corr <- metric_graph$new(edges = edge_list_corr,
                                    longlat = TRUE, verbose = 0)
cat("Direction consistency on the corrected pems graph:\n")
direction_diagnostic(pems_graph_corr, "corrected")

pems_graph_corr$add_observations(data = pems$data, normalized = TRUE,
                                 verbose = 0)

cat("Fitting WMD alpha=1 (direction-corrected graph)...\n")
fit_wmd1_corr <- graph_lme(y ~ 1, graph = pems_graph_corr, BC = 0,
                           model = "WMD1",
                           optim_method = "Nelder-Mead")

## ---- 3c. Second correction: junction directions --------------------------
## The BFS correction in 3b enforces consistency at degree-2 (pass-
## through) vertices but leaves junctions alone. Some junctions in the
## corrected graph still look topologically odd: degree >= 3 vertices
## where all edges point in (pure sink) or all out (pure source) make
## little sense for a road network unless they're at the boundary, and
## a 4-way intersection with 3 edges in and 1 out is similarly unusual.
## For each such junction we test which incident edge, if flipped,
## maximises the WMD1 log-likelihood at the *fitted* parameters (a
## cheap evaluation since coefficients are fixed). We then greedily
## accumulate the flips that improve likelihood, building a second
## corrected graph for a fresh WMD1 fit.

flip_edges_and_rebuild <- function(base_edge_list, flip_idx) {
  el <- base_edge_list
  for (i in flip_idx) {
    el[[i]] <- el[[i]][nrow(el[[i]]):1, , drop = FALSE]
  }
  g_tmp <- metric_graph$new(edges = el, longlat = TRUE, verbose = 0)
  prop <- propagate_orientation_flips(g_tmp)
  if (any(prop)) {
    for (i in which(prop)) {
      el[[i]] <- el[[i]][nrow(el[[i]]):1, , drop = FALSE]
    }
  }
  metric_graph$new(edges = el, longlat = TRUE, verbose = 0)
}

eval_loglik_flipped <- function(base_edge_list, flip_idx,
                                obs_data, baseline_fit) {
  g <- flip_edges_and_rebuild(base_edge_list, flip_idx)
  g$add_observations(data = obs_data, normalized = TRUE,
                     suppress_warnings = TRUE, verbose = 0)
  fit <- graph_lme(y ~ 1, graph = g, BC = 0, model = "WMD1",
                   previous_fit = baseline_fit, fix_coeff = TRUE,
                   optim_method = "Nelder-Mead")
  as.numeric(logLik(fit))
}

in_deg_c  <- pems_graph_corr$get_degrees("indegree")
out_deg_c <- pems_graph_corr$get_degrees("outdegree")
deg_c     <- in_deg_c + out_deg_c
susp_v <- which(deg_c >= 3 &
                (in_deg_c == 0 | out_deg_c == 0 |
                 abs(in_deg_c - out_deg_c) >= 2))
cat(sprintf("Suspicious junctions (deg>=3, pure source/sink or |in-out|>=2): %d\n",
            length(susp_v)))
for (v in susp_v) {
  cat(sprintf("  vertex %d: in=%d out=%d  (%.5f, %.5f)\n",
              v, in_deg_c[v], out_deg_c[v],
              pems_graph_corr$V[v, 1], pems_graph_corr$V[v, 2]))
}

adj_c <- vector("list", pems_graph_corr$nV)
for (e in seq_len(pems_graph_corr$nE)) {
  adj_c[[pems_graph_corr$E[e, 1]]] <- c(adj_c[[pems_graph_corr$E[e, 1]]], e)
  adj_c[[pems_graph_corr$E[e, 2]]] <- c(adj_c[[pems_graph_corr$E[e, 2]]], e)
}

cat("Investigating per-junction flips...\n")
baseline_loglik <- as.numeric(logLik(fit_wmd1_corr))
accepted_flips  <- integer(0)
running_loglik  <- baseline_loglik
for (v in susp_v) {
  incident <- adj_c[[v]]
  liks <- vapply(incident, function(e) {
    eval_loglik_flipped(edge_list_corr,
                        c(accepted_flips, e),
                        pems$data, fit_wmd1_corr)
  }, numeric(1))
  best_idx <- which.max(liks)
  best_e   <- incident[best_idx]
  best_ll  <- liks[best_idx]
  improve  <- best_ll - running_loglik
  cat(sprintf(
    "  vertex %d (in=%d/out=%d): best flip = edge %d, loglik %.3f (delta %+.3f)\n",
    v, in_deg_c[v], out_deg_c[v], best_e, best_ll, improve))
  if (improve > 0) {
    accepted_flips <- c(accepted_flips, best_e)
    running_loglik <- best_ll
  }
}
cat(sprintf("Accepted %d junction flips: %s\n",
            length(accepted_flips),
            paste(accepted_flips, collapse = ", ")))

cat("Building second-corrected graph...\n")
pems_graph_corr2 <- flip_edges_and_rebuild(edge_list_corr, accepted_flips)
cat("Direction consistency on the second-corrected pems graph:\n")
direction_diagnostic(pems_graph_corr2, "corrected-2")

pems_graph_corr2$add_observations(data = pems$data, normalized = TRUE,
                                  verbose = 0)

cat("Fitting WMD alpha=1 (second-corrected graph)...\n")
fit_wmd1_corr2 <- graph_lme(y ~ 1, graph = pems_graph_corr2, BC = 0,
                            model = "WMD1",
                            optim_method = "Nelder-Mead")

## ---- 3d. Vertex conditions on the second-corrected graph -----------------
## All WMD1 fits so far on `pems_graph_corr` and `pems_graph_corr2` use
## the *default* vertex condition,
##     f_in(w) = w / sum(w),
## with unit per-edge directional weights, so at every internal vertex
## the outgoing edge values equal the equal-weight average of incoming
## edge values. The variance-stationary alternative was tried on the
## raw graph and was uniformly worse on the corrected2 graph.
##
## For traffic-speed data a more physically motivated choice is a
## *length-weighted* average: at a merge of an on-ramp and a main
## highway segment, the outgoing speed should reflect the longer
## segment (the main road) more than the shorter one (the ramp). This
## is implemented by keeping f_in = w / sum(w) but setting the per-edge
## directional weights equal to the edge lengths. (For unit weights the
## default f_in collapses to 1/n; with length weights it becomes a
## length-weighted average without any change to f_in itself.)

cat("Fitting WMD alpha=1 (dir-corrected2, length-weighted)...\n")
pems_graph_corr2_len <- pems_graph_corr2$clone()
pems_graph_corr2_len$set_edge_weights(
  weights = as.numeric(pems_graph_corr2_len$edge_lengths), verbose = 0)
fit_wmd1_corr2_len <- graph_lme(y ~ 1, graph = pems_graph_corr2_len, BC = 0,
                                model = "WMD1",
                                optim_method = "Nelder-Mead")

## A third option treats road type as a discrete variable. The pems sf
## table has no road-class column, so we use a length threshold as a
## proxy: edges longer than 0.5 km are taken as "main freeway"
## segments, shorter ones as ramps / short connectors. Main edges get
## directional weight 5 and ramp edges weight 1. Since f_in normalises
## by sum(w), what matters at each merge is the *ratio* of weights, so
## the absolute scale is arbitrary; the 5:1 ratio gives the main
## incoming road about 5x the influence of an incoming ramp.

cat("Fitting WMD alpha=1 (dir-corrected2, road-type-weighted)...\n")
pems_graph_corr2_type <- pems_graph_corr2$clone()
edge_len_corr2 <- as.numeric(pems_graph_corr2_type$edge_lengths)
is_highway     <- edge_len_corr2 > 0.5
type_weight    <- ifelse(is_highway, 5, 1)
cat(sprintf("  Classified %d highway / %d other edges (threshold 0.5 km).\n",
            sum(is_highway), sum(!is_highway)))
pems_graph_corr2_type$set_edge_weights(weights = type_weight, verbose = 0)
fit_wmd1_corr2_type <- graph_lme(y ~ 1, graph = pems_graph_corr2_type,
                                 BC = 0, model = "WMD1",
                                 optim_method = "Nelder-Mead")

## A fourth option uses actual OSM road-class tags. The bundled `pems`
## sf table doesn't carry road-class metadata, so we use a pre-computed
## per-edge classification stored in pems_osm_highway.rds, produced by
## the companion script `pems_fetch_osm_highway.R`. Edges are matched
## to OSM ways by nearest-feature distance; for the current `pems`
## graph every edge matches an OSM way within ~3 m. The resulting tag
## breakdown is dominated by motorway / motorway_link (the freeway
## mainlines and ramps) with smaller numbers of primary / secondary /
## trunk. We map the OSM road-class hierarchy to numeric weights and
## use them as directional weights with the default f_in.

if (requireNamespace("here", quietly = TRUE)) {
  osm_path <- here::here("examples", "directional", "pems_osm_highway.rds")
} else {
  ## Fall back to a path relative to the working directory.
  osm_path <- file.path("examples", "directional", "pems_osm_highway.rds")
}
if (file.exists(osm_path)) {
  osm_tags <- readRDS(osm_path)
  if (nrow(osm_tags) != pems_graph_corr2$nE) {
    warning("OSM tag file row count does not match graph edges; skipping OSM-weighted fit.")
    fit_wmd1_corr2_osm <- NULL
  } else {
    cat("Fitting WMD alpha=1 (dir-corrected2, OSM-class-weighted)...\n")
    cat("  Matched OSM highway class:\n")
    print(table(osm_tags$highway))
    osm_weight_map <- c(motorway       = 8,
                        trunk          = 6,
                        primary        = 4,
                        secondary      = 2,
                        motorway_link  = 1,
                        trunk_link     = 1,
                        primary_link   = 1,
                        secondary_link = 1,
                        unmatched      = 1)
    osm_w <- osm_weight_map[osm_tags$highway]
    osm_w[is.na(osm_w)] <- 1
    pems_graph_corr2_osm <- pems_graph_corr2$clone()
    pems_graph_corr2_osm$set_edge_weights(weights = as.numeric(osm_w),
                                          verbose = 0)
    fit_wmd1_corr2_osm <- graph_lme(y ~ 1, graph = pems_graph_corr2_osm,
                                    BC = 0, model = "WMD1",
                                    optim_method = "Nelder-Mead")
  }
} else {
  cat("OSM tags file not found at ", osm_path,
      "; skipping OSM-weighted fit.\n")
  cat("Run examples/directional/pems_fetch_osm_highway.R to create it.\n")
  fit_wmd1_corr2_osm <- NULL
}

## ---- 4. Fit the graph Laplacian models -----------------------------------
## Two-step fit (no fixed effects -> with fixed effects) matches the
## vignette; the no-fixed-effect fit is used as starting values to
## stabilise optimisation.

cat("Fitting GL alpha=1...\n")
fit_GL1 <- graph_lme(y ~ -1, graph = pems_graph,
                     model = list(type = "graphLaplacian", alpha = 1),
                     optim_method = "Nelder-Mead")
fit_GL1 <- graph_lme(y ~ 1, graph = pems_graph,
                     model = list(type = "graphLaplacian", alpha = 1),
                     previous_fit = fit_GL1,
                     optim_method = "Nelder-Mead")

cat("Fitting GL alpha=2...\n")
fit_GL2 <- graph_lme(y ~ 1, graph = pems_graph,
                     model = list(type = "graphLaplacian", alpha = 2),
                     previous_fit = fit_GL1,
                     optim_method = "Nelder-Mead")

## ---- 5. Fit the isotropic exponential covariance model -------------------

cat("Fitting isoExp...\n")
fit_isoexp <- graph_lme(y ~ 1, graph = pems_graph,
                        model = list(type = "isoCov"),
                        optim_method = "Nelder-Mead")

## ---- 6. Collect models and run LOO CV ------------------------------------

fitted_models_list <- list(
  "isoExp"      = fit_isoexp,
  "GL1"         = fit_GL1,
  "GL2"         = fit_GL2,
  "alpha=1"     = fit_alpha1,
  "alpha=1 bc"  = fit_alpha1_bc,
  "alpha=2"     = fit_alpha2,
  "alpha=2 bc"  = fit_alpha2_bc,
  "WMD1"                    = fit_wmd1,
  "WMD1 var-stationary"     = fit_wmd1_vs,
  "WMD1 dir-corrected"      = fit_wmd1_corr,
  "WMD1 dir-corrected2"     = fit_wmd1_corr2,
  "WMD1 dir-corrected2 len-weighted"  = fit_wmd1_corr2_len,
  "WMD1 dir-corrected2 type-weighted" = fit_wmd1_corr2_type)
if (!is.null(fit_wmd1_corr2_osm)) {
  fitted_models_list[["WMD1 dir-corrected2 osm-weighted"]] <- fit_wmd1_corr2_osm
}

cat("\n--- Negative log-likelihoods ---\n")
print(-sapply(fitted_models_list, logLik))

cat("\nRunning posterior_crossvalidation_loo...\n")
t0 <- Sys.time()
res <- posterior_crossvalidation_loo(fitted_models_list, tibble = FALSE)
cat(sprintf("LOO CV done in %.1fs.\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

cat("\n--- LOO CV scores (lower is better) ---\n")
print(round(res$scores, 3))

invisible(list(fits = fitted_models_list, cv = res))
