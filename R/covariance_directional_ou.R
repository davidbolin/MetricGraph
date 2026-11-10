# covariance_directional_ou.R
#
# Closed-form covariance of the "proper global OU" directional model
# (alpha = 1, directional weights set via set_edge_weights()/
# setDirectionalWeightFunction()) on a finite acyclic directed metric tree.
# Implements Theorem thm:cov-lca of jonas_local/tex/main_rssb.tex (see
# jonas_local/cov_lca_plan.md sec 1 for a self-contained restatement).
#
# Inputs: a metric_graph with directional weights + DirectionalWeightFunction
# set, rate kappa > 0, white-noise scale tau, evaluation points (PtE, default
# the observation locations), and an optional per-source anchoring variance
# sigma_source.
# Outputs: directional_ou_covariance() -- dense covariance matrix;
# directional_ou_variance() -- marginal variance vector.
#
# Math summary (see the plan for the derivation):
#   - At interior vertex v, u(v-out) = sum_in beta_v(out,in) u(v-in) a.s.
#     beta_v(e_out, e_in) = -in_weight_by_edge[e_in]/out_weight_by_edge[e_out]
#     (directional_weight_vectors() in R/util.R -- shared with the constraint
#     matrix builder, so this is consistent with the fitted precision matrix
#     by construction).
#   - Transfer factor for a forward path x -> y: A(x,y) = exp(-kappa*d(x,y)) *
#     product of beta_v crossed strictly between x and y.
#   - Marginal variance recursion along an edge of length h, tail variance v0:
#     r(head) = exp(-2*kappa*h)*v0 + (1-exp(-2*kappa*h))/(2*kappa*tau^2).
#   - Covariance: r(s,t) = r(a,a)*A(a,s)*A(a,t), a = last common ancestor of
#     s and t (unique maximal element of the upstream sets); r(s,t) = 0 if no
#     common ancestor exists.
#
# tau convention (pinned, see plan sec 3.5): stationary variance is
# 1/(2*kappa*tau^2) -- this matches R/spde_covariance.R's r_1(), NOT the
# "reciprocal_tau" convention used elsewhere in the package
# (R/graph_likelihoods.R etc., where tau = 1/reciprocal_tau).


## ---------------------------------------------------------------------------
## Public API
## ---------------------------------------------------------------------------

#' Closed-form covariance of the directional OU model
#'
#' Evaluates the covariance function of the "proper global OU" directional
#' model (Whittle-Matern alpha = 1 with directional vertex conditions) in
#' closed form, via the last-common-ancestor (LCA) representation, without
#' assembling or inverting the precision matrix.
#'
#' @details
#' The covariance between two points `s` and `t` is
#' \deqn{r(s,t) = r(a,a) \, A(a,s) \, A(a,t)}
#' where `a` is the last common ancestor (LCA) of `s` and `t` (or covariance
#' `0` if no common ancestor exists), `r(a,a)` is the marginal variance at
#' `a`, and
#' \deqn{A(x,y) = \exp(-\kappa \, d(x,y)) \prod_v \beta_v}
#' is the transfer factor for the (unique, forward) path from `x` to `y`,
#' the product running over the `beta_v` vertex weights of every vertex
#' crossed strictly between `x` and `y`.
#'
#' In plain terms: the last common ancestor of two points is the
#' furthest-downstream point that is upstream of both of them (it does not
#' exist if the two points lie on disjoint headwater branches, in which case
#' their covariance is exactly zero). The transfer factor `A(x,y)` is how
#' much of the value at `x` propagates forward to `y`: exponential decay
#' over the distance `d(x,y)`, attenuated by the directional weight
#' `beta_v` at every vertex the path passes through.
#'
#' This closed form requires `graph` to represent a directed *tree*: the
#' undirected skeleton must have no cycles, and the directed graph must
#' itself be acyclic. If `graph` has a directed cycle, is not a tree, or has
#' no directional weight function set, this errors (via the internal
#' `directional_ou_setup()` validation) -- see
#' `graph$setDirectionalWeightFunction()` to set one.
#'
#' @param graph A `metric_graph` object with directional edge weights set
#'  (`set_edge_weights(directional_weights = ...)`) and the vertex weight
#'  functions set (`setDirectionalWeightFunction()`).
#' @param kappa Rate parameter of the OU process (`kappa > 0`).
#' @param tau White-noise scale. The stationary marginal variance is
#'  `1/(2*kappa*tau^2)` -- the same convention as `MetricGraph:::r_1()`. Note
#'  this is NOT `reciprocal_tau` (used elsewhere in the package as
#'  `tau = 1/reciprocal_tau`); pass `tau` directly here.
#' @param PtE Evaluation points, a matrix with columns `(edge_number,
#'  distance_on_edge)`. Defaults to `graph$get_PtE()` (the observation
#'  locations).
#' @param PtE2 Second set of evaluation points for the cross-covariance
#'  `Cov(u(PtE), u(PtE2))`. Defaults to `PtE`, giving the symmetric
#'  `n x n` covariance matrix at the observation locations.
#' @param sigma_source Anchoring variance at the source vertices
#'  (indegree 0). `NULL` (default) anchors every source at the stationary
#'  variance `1/(2*kappa*tau^2)`, matching
#'  `Qalpha1_edges_custom(..., stationary_points = "all")`. Otherwise a
#'  numeric vector giving `Var(u)` at each source vertex, either named by
#'  vertex index (as character) or positional in the order
#'  `which(graph$get_degrees("indegree") == 0)`.
#' @param normalized If `TRUE` (default), `distance_on_edge` in `PtE`/`PtE2`
#'  is in `[0, 1]`; if `FALSE`, distances are absolute.
#' @param cpp If `TRUE` (default), use the C++ numeric setup and the shared C++
#'  covariance fill for oriented in-trees and out-trees. Irregular trees use
#'  the generic R fallback. Set to `FALSE` to use the pure-R implementation
#'  throughout.
#'
#' @return A dense `matrix` with `nrow(PtE)` rows and `nrow(PtE2)` columns.
#' @seealso [directional_ou_variance()] for the diagonal only.
#' @examples
#' edge1 <- rbind(c(0, 0), c(1, 0))
#' edge2 <- rbind(c(1, 0), c(2, 0))
#' edge3 <- rbind(c(1, 0), c(1, 1))
#' graph <- metric_graph$new(edges = list(edge1, edge2, edge3))
#' graph$set_edge_weights(weights = data.frame(w = c(1, 1, 1)),
#'                        directional_weights = "w")
#' graph$setDirectionalWeightFunction(f_in = function(x) sqrt(x / sum(x)))
#' Sigma <- directional_ou_covariance(graph, kappa = 1, tau = 1,
#'                                    PtE = rbind(c(1, 0.5), c(2, 0.5)))
#' Sigma
#' @export
directional_ou_covariance <- function(graph, kappa, tau,
                                      PtE = NULL, PtE2 = NULL,
                                      sigma_source = NULL,
                                      normalized = TRUE,
                                      cpp = TRUE) {
  setup <- directional_ou_setup(
    graph, kappa, tau, sigma_source, cpp = cpp
  )
  covariance_builder <- if (cpp) {
    directional_ou_covariance_from_setup_cpp
  } else {
    directional_ou_covariance_from_setup
  }
  covariance_builder(
    setup, PtE = PtE, PtE2 = PtE2, normalized = normalized
  )
}

# Resolve PtE/PtE2 defaults, convert to absolute distance, and apply the
# DIRECTIONAL_OU_MAX_POINTS dense-matrix size guard -- the boundary logic
# shared by directional_ou_covariance_from_setup() (R pairwise fill) and
# directional_ou_covariance_from_setup_cpp() (C++ oriented-tree fill dispatch).
# Extracted so the size guard (correctness-relevant: a dense n1 x n2 matrix
# beyond this is not viable) can't silently drift between the two paths.
#'
#' @param DIRECTIONAL_OU_MAX_POINTS Maximum number of points allowed in either
#'   `PtE` or `PtE2`. Defaults to the `DIRECTIONAL_OU_MAX_POINTS` option, or
#'   `10000L` when that option is unset. Larger inputs are rejected because the
#'   covariance calculation constructs a dense matrix.
#' @noRd
directional_ou_resolve_PtE_pair <- function(setup, PtE, PtE2, normalized,
                                            DIRECTIONAL_OU_MAX_POINTS =
                                              getOption(
                                                "DIRECTIONAL_OU_MAX_POINTS",
                                                10000L
                                              )) {
  graph <- setup$graph
  if (is.null(PtE)) {
    PtE <- graph$get_PtE()
  }
  if (is.null(PtE2)) {
    PtE2 <- PtE
  }

  # normalized=TRUE distances are in [0,1]; convert to absolute distance once,
  # here, at the boundary -- everything below this point works in absolute
  # edge-length units.
  edge_lengths <- graph$edge_lengths
  PtE_abs <- to_absolute_PtE(PtE, edge_lengths, normalized)
  PtE2_abs <- to_absolute_PtE(PtE2, edge_lengths, normalized)

  n1 <- nrow(PtE_abs)
  n2 <- nrow(PtE2_abs)
  if (length(DIRECTIONAL_OU_MAX_POINTS) != 1L ||
      !is.numeric(DIRECTIONAL_OU_MAX_POINTS) ||
      !is.finite(DIRECTIONAL_OU_MAX_POINTS) ||
      DIRECTIONAL_OU_MAX_POINTS < 1) {
    stop("option DIRECTIONAL_OU_MAX_POINTS must be one positive finite number")
  }
  if (n1 > DIRECTIONAL_OU_MAX_POINTS || n2 > DIRECTIONAL_OU_MAX_POINTS) {
    stop(sprintf(
      "directional_ou_covariance() builds a dense n1 x n2 matrix; got n1=%d, n2=%d, both must be <= %d. Supply a smaller PtE/PtE2.",
      n1, n2, DIRECTIONAL_OU_MAX_POINTS
    ))
  }

  list(PtE_abs = PtE_abs, PtE2_abs = PtE2_abs)
}

# Pairwise-fill covariance from an already-built setup (from
# directional_ou_setup(), or a cached directional_ou_setup_structure() plus
# a freshly-computed directional_ou_setup_numeric() merged into the same
# shape -- see directional_ou_setup()'s composer). Split out of
# directional_ou_covariance() so a caller that already has a setup (e.g. a
# precomputed likelihood evaluator reusing the graph-structural part across
# repeated kappa/tau) can skip rebuilding it.
#' @noRd
directional_ou_covariance_from_setup <- function(setup, PtE, PtE2 = NULL,
                                                  normalized = TRUE) {
  resolved <- directional_ou_resolve_PtE_pair(setup, PtE, PtE2, normalized)
  PtE_abs <- resolved$PtE_abs
  PtE2_abs <- resolved$PtE2_abs

  n1 <- nrow(PtE_abs)
  n2 <- nrow(PtE2_abs)
  Sigma <- matrix(0, nrow = n1, ncol = n2)
  same_point_set <- identical(PtE_abs, PtE2_abs)
  for (i in seq_len(n1)) {
    j_start <- if (same_point_set) i else 1L
    for (j in seq_len(n2)) {
      if (same_point_set && j < i) next  # fill by symmetry below
      Sigma[i, j] <- directional_ou_pair_covariance(
        setup, PtE_abs[i, ], PtE2_abs[j, ]
      )
    }
  }
  if (same_point_set) {
    Sigma[lower.tri(Sigma)] <- t(Sigma)[lower.tri(Sigma)]
  }
  Sigma
}

# Dispatch wrapper for the pairwise-fill covariance: routes both supported
# oriented-tree shapes through one O(n^2) C++ kernel -- converging in-trees
# (including the K1/K2 river graphs) use Euler ancestor containment, while
# diverging out-trees (including reversed continuity) use their cached
# Euler/RMQ LCA index. Only irregularly oriented trees retain the generic R
# graph-walk fallback.
#' @noRd
directional_ou_covariance_from_setup_cpp <- function(setup, PtE, PtE2 = NULL,
                                                       normalized = TRUE) {
  if (!(setup$tree_orientation %in% c("in", "out"))) {
    return(directional_ou_covariance_from_setup(setup, PtE, PtE2, normalized))
  }

  resolved <- directional_ou_resolve_PtE_pair(setup, PtE, PtE2, normalized)
  PtE_abs <- resolved$PtE_abs
  PtE2_abs <- resolved$PtE2_abs

  directional_ou_covariance_oriented_tree_cpp(
    E = setup$graph$E,
    edge_lengths = setup$graph$edge_lengths,
    tree_orientation = setup$tree_orientation,
    kappa = setup$kappa,
    sigma_stationary = setup$sigma_stationary,
    var_tail = setup$var_tail,
    logG_tail = setup$logG_tail,
    signG_tail = setup$signG_tail,
    enter = setup$enter,
    exit = setup$exit,
    depth = setup$depth,
    out_tree_lca_index = setup$out_tree_lca_index,
    PtE_abs = PtE_abs,
    PtE2_abs = PtE2_abs,
    same_point_set = identical(PtE_abs, PtE2_abs)
  )
}

#' Closed-form marginal variance of the directional OU model
#'
#' `r(x,x)` at each evaluation point -- the diagonal of
#' [directional_ou_covariance()], computed directly via the marginal
#' variance recursion along each edge, without any last-common-ancestor or
#' transfer-factor work.
#'
#' @inheritParams directional_ou_covariance
#' @return A numeric vector of length `nrow(PtE)`.
#' @seealso [directional_ou_covariance()] for the full closed-form
#'  covariance (including the `\deqn{r(s,t) = r(a,a) A(a,s) A(a,t)}` formula
#'  and the tree/acyclicity requirements).
#' @examples
#' edge1 <- rbind(c(0, 0), c(1, 0))
#' edge2 <- rbind(c(1, 0), c(2, 0))
#' edge3 <- rbind(c(1, 0), c(1, 1))
#' graph <- metric_graph$new(edges = list(edge1, edge2, edge3))
#' graph$set_edge_weights(weights = data.frame(w = c(1, 1, 1)),
#'                        directional_weights = "w")
#' graph$setDirectionalWeightFunction(f_in = function(x) sqrt(x / sum(x)))
#' directional_ou_variance(graph, kappa = 1, tau = 1,
#'                         PtE = rbind(c(1, 0.5), c(2, 0.5)))
#' @export
directional_ou_variance <- function(graph, kappa, tau, PtE = NULL,
                                    sigma_source = NULL, normalized = TRUE,
                                    cpp = TRUE) {
  if (is.null(PtE)) {
    PtE <- graph$get_PtE()
  }
  edge_lengths <- graph$edge_lengths
  PtE_abs <- to_absolute_PtE(PtE, edge_lengths, normalized)

  setup <- directional_ou_setup(
    graph, kappa, tau, sigma_source, cpp = cpp
  )
  directional_point_variance(setup, PtE_abs)
}


## ---------------------------------------------------------------------------
## Boundary helpers (public-function plumbing only, no graph math)
## ---------------------------------------------------------------------------

# Convert a PtE matrix to absolute distance-on-edge, once, at the boundary of
# the public functions.
#' @noRd
to_absolute_PtE <- function(PtE, edge_lengths, normalized) {
  PtE <- as.matrix(PtE)
  if (normalized) {
    PtE[, 2] <- PtE[, 2] * edge_lengths[PtE[, 1]]
  }
  PtE
}

# r(s,t) for a single pair of (already absolute-distance) points.
#' @noRd
directional_ou_pair_covariance <- function(setup, point_s, point_t) {
  lca <- directional_lca(setup, point_s, point_t)
  if (is.null(lca)) {
    return(0)
  }
  r_aa <- directional_point_variance(setup, matrix(lca$point, nrow = 1))
  log_A_s <- directional_log_transfer(setup, lca$point, point_s)
  log_A_t <- directional_log_transfer(setup, lca$point, point_t)
  r_aa * log_A_s$sign * log_A_t$sign * exp(log_A_s$log_abs + log_A_t$log_abs)
}


## ---------------------------------------------------------------------------
## Internal layer -- per-graph setup
## ---------------------------------------------------------------------------
##
## directional_ou_setup() precomputes everything the pairwise/pointwise
## evaluators above need, once per (graph, kappa, tau, sigma_source):
##   - a topological order of the edges (source -> outlet), which also
##     doubles as the acyclicity check the theorem requires;
##   - beta_v(e_out, e_in) inputs (out_by_edge/in_by_edge), shared with the
##     constraint-matrix builder via directional_weight_vectors();
##   - the vertex variance recursion var_tail/var_head (marginal variance of
##     u just before/after crossing each edge);
##   - root-normalised log-transfer accumulators logG_*/signG_* -- the
##     log|A(edge-head or edge-tail, outlet)| and its sign, needed by the
##     (not-yet-written) fast dendritic covariance path so it can form
##     A(x,y) = A(x,outlet)/A(y,outlet) without repeated products of beta's
##     along a path.

# Topological order of the edges: edge e depends on every edge f whose head
# is e's tail (E[f,2] == E[e,1]) -- f must be "crossed" before e can be. This
# is Kahn's algorithm on that edge-precedence graph. An edge's initial
# in-edge-graph-degree is exactly the indegree of its tail vertex
# (V_indegree[E[e,1]]), since every in-edge at that vertex is a predecessor.
#' @noRd
directional_edge_topo_order <- function(E, nE, nV, V_indegree) {
  out_edges_by_vertex <- split(seq_len(nE), factor(E[, 1], levels = seq_len(nV)))
  remaining_indeg <- V_indegree  # decremented as each in-edge at a vertex is processed

  queue <- which(V_indegree[E[, 1]] == 0)  # edges whose tail vertex is already a source
  topo_order <- integer(nE)
  n_done <- 0L
  head_ptr <- 1L
  while (head_ptr <= length(queue)) {
    e <- queue[head_ptr]
    head_ptr <- head_ptr + 1L
    n_done <- n_done + 1L
    topo_order[n_done] <- e

    v <- E[e, 2]
    remaining_indeg[v] <- remaining_indeg[v] - 1L
    if (remaining_indeg[v] == 0L) {
      queue <- c(queue, out_edges_by_vertex[[v]])
    }
  }

  if (n_done != nE) {
    stop("the directed graph has a cycle: not all edges could be topologically ordered. directional_ou_setup() requires an acyclic directed graph (Theorem thm:cov-lca).")
  }
  topo_order
}

# Number of connected components of the undirected skeleton (E[,1]-E[,2]
# edges, direction ignored), via union-find. Used for the tree-structure
# guard: a directed graph gives every point a unique last-common-ancestor
# only if it is a tree once direction is forgotten.
#' @noRd
count_undirected_components <- function(E, nE, nV) {
  parent <- seq_len(nV)
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]  # path halving
      x <- parent[x]
    }
    x
  }
  for (e in seq_len(nE)) {
    r1 <- find_root(E[e, 1])
    r2 <- find_root(E[e, 2])
    if (r1 != r2) parent[r1] <- r2
  }
  roots <- vapply(seq_len(nV), find_root, numeric(1))
  length(unique(roots))
}

#' Precompute the graph-structural (kappa/tau/sigma_source-independent) part
#' of the directional OU setup
#'
#' Validates that `graph` is a valid input for Theorem thm:cov-lca
#' (directional weight function set, acyclic, a tree), then computes a
#' topological edge order and the beta_v inputs (`out_by_edge`/`in_by_edge`).
#' This is the part of `directional_ou_setup()` that depends only on `graph`,
#' not on `kappa`/`tau`/`sigma_source` -- split out so it can be computed once
#' and reused across multiple `(kappa, tau, sigma_source)` calls (e.g. during
#' likelihood optimization).
#'
#' @param graph A `metric_graph` object with directional edge weights and
#'  vertex weight functions set, as in [directional_ou_covariance()].
#' @return A list with components `graph`, `E`, `nE`, `nV`, `V_indegree`,
#'  `V_outdegree`, `topo_order`, `out_by_edge`, `in_by_edge`, `is_dendritic`,
#'  `tree_orientation` (`"in"`/`"out"`/`"irregular"`),
#'  `in_edges_by_vertex`, `out_edges_by_vertex`, `enter`, `exit`, and `depth`
#'  (the last three `NULL` for irregular trees), plus `parent_edge` and
#'  `out_tree_lca_index` (`NULL` unless `tree_orientation == "out"`).
#' @noRd
directional_ou_setup_structure <- function(graph) {
  if (is.null(graph$DirectionalWeightFunction_in)) {
    stop("graph has no directional weight function set; call graph$setDirectionalWeightFunction() before using the directional OU covariance.")
  }

  E <- graph$E
  nE <- graph$nE
  nV <- graph$nV
  V_indegree <- graph$get_degrees("indegree")
  V_outdegree <- graph$get_degrees("outdegree")

  topo_order <- directional_edge_topo_order(E, nE, nV, V_indegree)

  n_components <- count_undirected_components(E, nE, nV)
  if (nE != nV - n_components) {
    stop(paste0(
      "graph is not a tree (nE = ", nE, ", nV = ", nV, ", undirected components = ",
      n_components, "); directional OU covariance requires a directed tree. ",
      "Vertices with indegree > 1 (candidates for a broken unique-path structure): ",
      paste(which(V_indegree > 1), collapse = ", ")
    ))
  }

  # Same accessor pattern as buildDirectionalConstraints() (R/metric_graph.R).
  dw <- graph$.__enclos_env__$private$directional_weights
  weight <- as.vector(graph$get_edge_weights()[[dw]])
  weight_vectors <- directional_weight_vectors(
    E = E, nE = nE, weight = weight,
    DirectionalWeightFunction_out = graph$DirectionalWeightFunction_out,
    DirectionalWeightFunction_in = graph$DirectionalWeightFunction_in
  )
  out_by_edge <- weight_vectors$out_by_edge
  in_by_edge <- weight_vectors$in_by_edge

  is_dendritic <- all(V_outdegree <= 1)
  is_out_tree <- !is_dendritic && all(V_indegree <= 1)
  tree_orientation <- if (is_dendritic) {
    "in"
  } else if (is_out_tree) {
    "out"
  } else {
    "irregular"
  }

  # Vertex -> incident-edge lookups, replacing which(E[,2]==v)/which(E[,1]==v)
  # scans (each O(nE)) with an O(1) list index. split() keys its result by
  # as.character() of the grouping values, so index by as.character(v), not
  # v directly; a vertex with no such edges simply has no key (index with
  # `[[...]]` returns NULL, not the empty-which()'s integer(0) -- callers
  # must guard, e.g. `%||%` a fallback of integer(0)).
  in_edges_by_vertex <- split(seq_len(nE), E[, 2])
  out_edges_by_vertex <- split(seq_len(nE), E[, 1])

  ancestor_labels <- directional_ou_ancestor_labels(
    E, V_indegree, V_outdegree, in_edges_by_vertex,
    out_edges_by_vertex, tree_orientation
  )
  out_tree_lca_index <- if (identical(tree_orientation, "out")) {
    directional_ou_out_tree_lca_index_cpp(
      ancestor_labels$parent_edge, ancestor_labels$depth
    )
  } else {
    NULL
  }

  list(
    graph = graph,
    E = E,
    nE = nE,
    nV = nV,
    V_indegree = V_indegree,
    V_outdegree = V_outdegree,
    topo_order = topo_order,
    out_by_edge = out_by_edge,
    in_by_edge = in_by_edge,
    is_dendritic = is_dendritic,
    tree_orientation = tree_orientation,
    in_edges_by_vertex = in_edges_by_vertex,
    out_edges_by_vertex = out_edges_by_vertex,
    enter = ancestor_labels$enter,
    exit = ancestor_labels$exit,
    depth = ancestor_labels$depth,
    parent_edge = ancestor_labels$parent_edge,
    out_tree_lca_index = out_tree_lca_index
  )
}

#' Euler-tour interval labels for O(1) ancestor tests on oriented trees
#'
#' In-trees are rooted at outlets and walked upstream; out-trees are rooted
#' at sources and walked downstream. The resulting interval-containment test
#' is used by `directional_edge_downstream_of()`. Irregularly oriented trees
#' return `NULL` labels and retain the generic graph-walk fallback.
#'
#' @param E The graph's edge matrix (columns `(tail, head)`).
#' @param V_indegree Indegree of each vertex.
#' @param V_outdegree Outdegree of each vertex.
#' @param in_edges_by_vertex List (as produced by `split()`, keyed by
#'  `as.character(vertex)`) of in-edges at each vertex.
#' @param out_edges_by_vertex Corresponding list of out-edges.
#' @param tree_orientation One of `"in"`, `"out"`, or `"irregular"`.
#' @return `list(enter, exit, depth, parent_edge)`. `parent_edge` is only
#'  populated for out-trees.
#' @noRd
directional_ou_ancestor_labels <- function(E, V_indegree, V_outdegree,
                                           in_edges_by_vertex,
                                           out_edges_by_vertex,
                                           tree_orientation) {
  if (!(tree_orientation %in% c("in", "out"))) {
    return(list(
      enter = NULL, exit = NULL, depth = NULL, parent_edge = NULL
    ))
  }

  nE <- nrow(E)
  enter <- integer(nE)
  exit <- integer(nE)
  depth <- integer(nE)
  counter <- 0L

  if (tree_orientation == "in") {
    roots <- which(V_outdegree[E[, 2]] == 0)
    children_of <- function(edge) {
      in_edges_by_vertex[[as.character(E[edge, 1])]]
    }
  } else {
    roots <- which(V_indegree[E[, 1]] == 0)
    children_of <- function(edge) {
      out_edges_by_vertex[[as.character(E[edge, 2])]]
    }
  }

  # Iterative, explicit-stack DFS. A recursive visit() (one R call per edge
  # along the deepest upstream chain) can exhaust the native C stack --an
  # uncatchable crash, not an R-level error -- well within plausible depth
  # for a real river network's finer headwater chains. Preallocated to
  # 2*nE frames (every edge is pushed exactly once as an "enter" frame and
  # once as an "exit" frame across the whole traversal) to avoid repeated
  # vector-growth copies. Pushing the exit frame before a node's (reversed)
  # children reproduces the recursive visit()'s pre/post-order numbering
  # exactly: children are fully unwound (enter *and* exit) before their
  # next sibling is popped, just as nested recursive calls would.
  stack_edge <- integer(2L * nE)
  stack_is_enter <- logical(2L * nE)
  stack_depth <- integer(2L * nE)
  sp <- 0L
  push <- function(edge, is_enter, edge_depth) {
    sp <<- sp + 1L
    stack_edge[sp] <<- edge
    stack_is_enter[sp] <<- is_enter
    stack_depth[sp] <<- edge_depth
  }

  for (root in roots) {
    push(root, TRUE, 0L)
    while (sp > 0L) {
      edge <- stack_edge[sp]
      is_enter <- stack_is_enter[sp]
      edge_depth <- stack_depth[sp]
      sp <- sp - 1L
      if (is_enter) {
        counter <- counter + 1L
        enter[edge] <- counter
        depth[edge] <- edge_depth
        push(edge, FALSE, edge_depth)
        for (child in rev(children_of(edge))) {
          push(child, TRUE, edge_depth + 1L)
        }
      } else {
        counter <- counter + 1L
        exit[edge] <- counter
      }
    }
  }

  parent_edge <- NULL
  if (tree_orientation == "out") {
    # In an out-tree every vertex has at most one incoming edge. Map each
    # vertex to that edge once, then index by edge tail to obtain all parent
    # edges without per-edge character conversion and named-list lookup.
    incoming_edge_by_vertex <- rep(NA_integer_, length(V_indegree))
    incoming_edge_by_vertex[E[, 2]] <- seq_len(nE)
    parent_edge <- incoming_edge_by_vertex[E[, 1]]
  }

  list(
    enter = enter, exit = exit, depth = depth, parent_edge = parent_edge
  )
}

#' Resolve `sigma_source` into a flat, per-vertex anchoring-variance vector
#'
#' Implements the `sigma_source` name-matching/positional-matching semantics
#' documented on [directional_ou_covariance()]: `NULL` anchors every source
#' vertex (indegree 0) at `sigma_stationary`; otherwise `sigma_source` is
#' matched either by name (character vertex index) or positionally, in the
#' order `which(V_indegree == 0)`. Extracted out of
#' `directional_ou_setup_numeric()` so the C++ port
#' (`directional_ou_setup_numeric_cpp()`) has a single R-side place to
#' resolve this from -- R closures/named-vector semantics don't cross the
#' C++ boundary, so this is evaluated once, up front, in R, exactly like
#' `directional_weight_vectors()` is for the constraint-matrix builder.
#'
#' @param structure Output of `directional_ou_setup_structure()`.
#' @param sigma_source `NULL` or a numeric vector, as documented on
#'  [directional_ou_covariance()].
#' @param sigma_stationary The stationary variance `1/(2*kappa*tau^2)`.
#' @return A length-`nV` numeric vector, indexed by vertex: entry `v` is the
#'  anchoring variance for vertex `v`. Only meaningful where
#'  `V_indegree[v] == 0`; non-source entries default to `sigma_stationary`
#'  but are otherwise unused -- callers must still guard on `V_indegree`.
#' @noRd
directional_ou_resolve_source_var <- function(structure, sigma_source, sigma_stationary) {
  V_indegree <- structure$V_indegree
  nV <- structure$nV

  source_vertices <- which(V_indegree == 0)
  named_source_var <- stats::setNames(rep(sigma_stationary, length(source_vertices)),
                                      as.character(source_vertices))
  if (!is.null(sigma_source)) {
    if (!is.null(names(sigma_source))) {
      bad_names <- setdiff(names(sigma_source), as.character(source_vertices))
      if (length(bad_names) > 0) {
        stop(sprintf(
          "directional_ou_resolve_source_var(): sigma_source name(s) %s do not match any source vertex (indegree 0); source vertices are %s",
          paste(bad_names, collapse = ", "), paste(source_vertices, collapse = ", ")
        ))
      }
      named_source_var[names(sigma_source)] <- sigma_source
    } else {
      if (length(sigma_source) != length(source_vertices)) {
        stop(sprintf(
          "directional_ou_resolve_source_var(): unnamed sigma_source has length %d, but there are %d source vertices (indegree 0); supply a length-%d vector (in the order of `which(V_indegree == 0)`) or name sigma_source by vertex index",
          length(sigma_source), length(source_vertices), length(source_vertices)
        ))
      }
      named_source_var[as.character(source_vertices)] <- sigma_source
    }
  }

  flat <- rep(sigma_stationary, nV)
  flat[source_vertices] <- named_source_var[as.character(source_vertices)]
  flat
}

#' Precompute the kappa/tau/sigma_source-dependent part of the directional OU
#' setup
#'
#' Given the graph-structural part from `directional_ou_setup_structure()`,
#' computes the vertex variance recursion (`var_tail`/`var_head`) and the
#' root-normalised log-transfer accumulators (`logG_*`/`signG_*`) -- the
#' quantities that change whenever `kappa`, `tau`, or `sigma_source` change,
#' even if `graph` doesn't.
#'
#' @param structure Output of `directional_ou_setup_structure()`.
#' @inheritParams directional_ou_covariance
#' @return A list with components `kappa`, `tau`, `sigma_stationary`,
#'  `var_tail`, `var_head`, `logG_head`, `logG_tail`, `signG_head`,
#'  `signG_tail`.
#' @noRd
directional_ou_setup_numeric <- function(structure, kappa, tau, sigma_source) {
  E <- structure$E
  nE <- structure$nE
  topo_order <- structure$topo_order
  V_indegree <- structure$V_indegree
  V_outdegree <- structure$V_outdegree
  out_by_edge <- structure$out_by_edge
  in_by_edge <- structure$in_by_edge
  edge_lengths <- structure$graph$edge_lengths
  in_edges_by_vertex <- structure$in_edges_by_vertex
  out_edges_by_vertex <- structure$out_edges_by_vertex

  sigma_stationary <- 1 / (2 * kappa * tau^2)

  # Per-source-vertex anchoring variance, broadcast across all out-edges of
  # that vertex (a source can have outdegree > 1). See
  # directional_ou_resolve_source_var() for the sigma_source name-matching
  # logic (extracted out so the C++ port has a single R-side place to get a
  # flat, resolved anchoring-variance vector from).
  source_var <- directional_ou_resolve_source_var(structure, sigma_source, sigma_stationary)

  var_tail <- numeric(nE)
  var_head <- numeric(nE)
  for (e in topo_order) {
    v <- E[e, 1]
    if (V_indegree[v] == 0) {
      var_tail[e] <- source_var[v]
    } else {
      E_in <- in_edges_by_vertex[[as.character(v)]]  # replaces which(E[,2]==v); V_indegree[v]>0 guarantees a non-NULL key here
      beta_sq <- (in_by_edge[E_in] / out_by_edge[e])^2  # beta_v(e,f) = -in_by_edge[f]/out_by_edge[e]
      var_tail[e] <- sum(beta_sq * var_head[E_in])
    }
    c2 <- exp(-2 * kappa * edge_lengths[e])
    var_head[e] <- c2 * var_tail[e] + (1 - c2) * sigma_stationary
  }

  # Root-normalised log-transfer accumulators. In-trees are normalised toward
  # their outlet and processed in reverse topological order. Out-trees are
  # normalised from their source and processed in forward topological order.
  # The latter is the orientation produced by reversing a dendritic river
  # graph for the continuity stand-in.
  logG_head <- numeric(nE)
  logG_tail <- numeric(nE)
  signG_head <- numeric(nE)
  signG_tail <- numeric(nE)
  if (identical(structure$tree_orientation, "out")) {
    for (e in topo_order) {
      v <- E[e, 1]
      if (V_indegree[v] == 0) {
        logG_tail[e] <- 0
        signG_tail[e] <- 1
      } else {
        f <- in_edges_by_vertex[[as.character(v)]][1]
        beta_ef <- -in_by_edge[f] / out_by_edge[e]
        logG_tail[e] <- logG_head[f] + log(abs(beta_ef))
        signG_tail[e] <- signG_head[f] * sign(beta_ef)
      }
      logG_head[e] <- logG_tail[e] - kappa * edge_lengths[e]
      signG_head[e] <- signG_tail[e]
    }
  } else {
    for (e in rev(topo_order)) {
      v <- E[e, 2]
      if (V_outdegree[v] == 0) {
        logG_head[e] <- 0
        signG_head[e] <- 1
      } else {
        # Off the in-tree path the first out-edge is arbitrary, so irregular
        # trees never consume these accumulators.
        g <- out_edges_by_vertex[[as.character(v)]][1]
        beta_ge <- -in_by_edge[e] / out_by_edge[g]
        logG_head[e] <- logG_tail[g] + log(abs(beta_ge))
        signG_head[e] <- signG_tail[g] * sign(beta_ge)
      }
      logG_tail[e] <- logG_head[e] - kappa * edge_lengths[e]
      signG_tail[e] <- signG_head[e]
    }
  }

  list(
    kappa = kappa,
    tau = tau,
    sigma_stationary = sigma_stationary,
    var_tail = var_tail,
    var_head = var_head,
    logG_head = logG_head,
    logG_tail = logG_tail,
    signG_head = signG_head,
    signG_tail = signG_tail
  )
}

#' Precompute per-edge quantities for the directional OU covariance
#'
#' Internal setup shared by [directional_ou_covariance()] and
#' [directional_ou_variance()]: validates that `graph` is a valid input for
#' Theorem thm:cov-lca (directional weight function set, acyclic, a tree),
#' then computes a topological edge order, the beta_v inputs, the vertex
#' variance recursion, and root-normalised log-transfer accumulators.
#'
#' A thin composer of `directional_ou_setup_structure()` (the
#' `graph`-only part) and `directional_ou_setup_numeric()` (the
#' `kappa`/`tau`/`sigma_source`-dependent part) -- kept split so the
#' structural part can be reused across repeated calls with different
#' `kappa`/`tau`/`sigma_source` (e.g. likelihood optimization).
#'
#' @inheritParams directional_ou_covariance
#' @return A list with components `graph`, `kappa`, `tau`, `sigma_stationary`,
#'  `V_indegree`, `V_outdegree`, `topo_order`, `out_by_edge`, `in_by_edge`,
#'  `var_tail`, `var_head`, `logG_head`, `logG_tail`, `signG_head`,
#'  `signG_tail`, `is_dendritic`, `tree_orientation`, `enter`, `exit`,
#'  `depth`, and `parent_edge`.
#' @noRd
# Compose the structural and numeric setup pieces in one place. Public
# covariance calls and cached likelihood calls use the same object shape.
#' @noRd
directional_ou_compose_setup <- function(structure, numeric_part) {
  list(
    graph = structure$graph,
    kappa = numeric_part$kappa,
    tau = numeric_part$tau,
    sigma_stationary = numeric_part$sigma_stationary,
    V_indegree = structure$V_indegree,
    V_outdegree = structure$V_outdegree,
    topo_order = structure$topo_order,
    out_by_edge = structure$out_by_edge,
    in_by_edge = structure$in_by_edge,
    var_tail = numeric_part$var_tail,
    var_head = numeric_part$var_head,
    logG_head = numeric_part$logG_head,
    logG_tail = numeric_part$logG_tail,
    signG_head = numeric_part$signG_head,
    signG_tail = numeric_part$signG_tail,
    is_dendritic = structure$is_dendritic,
    tree_orientation = structure$tree_orientation,
    enter = structure$enter,
    exit = structure$exit,
    depth = structure$depth,
    parent_edge = structure$parent_edge,
    out_tree_lca_index = structure$out_tree_lca_index,
    in_edges_by_vertex = structure$in_edges_by_vertex,
    out_edges_by_vertex = structure$out_edges_by_vertex
  )
}

# Graph-skeleton view of directional_ou_setup_structure()'s output, shaped
# for an eventual R->C++ boundary (matching the skeleton <- list(nE=,
# edge_lengths=, E=, get_degrees=) convention in
# R/graph_likelihoods_v2.R's profile_Q_alpha1_directional(), extended with
# the directional-specific fields a C++ core will also need). Consumed by
# directional_ou_covariance_loglik_precompute() (R/likelihood_directional_ou_covariance.R).
#' @noRd
directional_ou_skeleton <- function(structure) {
  list(
    nE = structure$nE,
    edge_lengths = structure$graph$edge_lengths,
    E = structure$E,
    get_degrees = structure$V_indegree == 0,
    V_indegree = structure$V_indegree,
    V_outdegree = structure$V_outdegree,
    out_by_edge = structure$out_by_edge,
    in_by_edge = structure$in_by_edge,
    topo_order = structure$topo_order,
    tree_orientation = structure$tree_orientation
  )
}

# Select the numeric setup implementation without duplicating source-variance
# resolution or setup composition in each caller.
#' @noRd
directional_ou_setup_numeric_dispatch <- function(structure, kappa, tau,
                                                  sigma_source, cpp = TRUE) {
  if (!cpp) {
    return(directional_ou_setup_numeric(
      structure, kappa, tau, sigma_source
    ))
  }

  sigma_stationary <- 1 / (2 * kappa * tau^2)
  source_var <- directional_ou_resolve_source_var(
    structure, sigma_source, sigma_stationary
  )
  directional_ou_setup_numeric_cpp(
    directional_ou_skeleton(structure), kappa, tau, source_var
  )
}

#' @noRd
directional_ou_setup <- function(graph, kappa, tau, sigma_source,
                                 cpp = TRUE) {
  structure <- directional_ou_setup_structure(graph)
  numeric_part <- directional_ou_setup_numeric_dispatch(
    structure, kappa, tau, sigma_source, cpp = cpp
  )
  directional_ou_compose_setup(structure, numeric_part)
}


## ---------------------------------------------------------------------------
## Internal layer -- point evaluators
## ---------------------------------------------------------------------------
##
## The pointwise/pairwise evaluators used by directional_ou_pair_covariance()
## and directional_ou_variance(): the last-common-ancestor (LCA) search, the
## marginal-variance recursion evaluated at an arbitrary point on an edge
## (not just its head), and the transfer factor A(from, to) between two
## points with from upstream of (or equal to) to -- via the O(1) dendritic
## fast path when available, else an explicit forward path walk.

#' Last common ancestor (LCA) of two points
#'
#' Finds the unique maximal element of `Lambda-up(point_s) \cap
#' Lambda-up(point_t)` (the ancestor sets, i.e. points reachable by walking
#' upstream against edge direction), needed by
#' `directional_ou_pair_covariance()` to anchor `r(s,t) = r(a,a) * A(a,s) *
#' A(a,t)`. Four cases, tried in order:
#'   1. Same edge: the LCA is the upstream one of the two (smaller distance).
#'   2. One point's edge is (weakly) downstream of the other's: the LCA is
#'      the upstream point itself.
#'   3. Neither contains the other: only possible when the graph is
#'      non-dendritic (some vertex has outdegree > 1 -- a true diffluence
#'      upstream of both points). In a dendritic tree, two points whose
#'      branches only merge downstream of both share no common ancestor at
#'      all (confluences, indegree > 1, don't create upstream ambiguity) --
#'      that's case 4, not case 3. Best-effort / untested here: no
#'      non-dendritic fixture exists yet.
#'   4. No common ancestor: returns `NULL` (the caller,
#'      `directional_ou_pair_covariance()`, treats this as covariance 0).
#'
#' @param setup A setup list from `directional_ou_setup()`.
#' @param point_s,point_t Length-2 `c(edge_number, absolute_distance)`
#'  vectors.
#' @return `list(point = c(edge_number, absolute_distance))` giving the LCA,
#'  or `NULL` if `point_s` and `point_t` have no common ancestor.
#' @noRd
directional_lca <- function(setup, point_s, point_t) {
  E <- setup$graph$E

  # Case 1: same edge -- the LCA is the upstream (smaller-distance) point.
  if (point_s[1] == point_t[1]) {
    return(list(point = c(point_s[1], min(point_s[2], point_t[2]))))
  }

  # Case 2: one point's edge is (weakly) downstream of the other's.
  if (directional_edge_downstream_of(setup, point_s[1], point_t[1])) {
    return(list(point = point_s))
  }
  if (directional_edge_downstream_of(setup, point_t[1], point_s[1])) {
    return(list(point = point_t))
  }

  # Case 3: neither contains the other -- non-dendritic diffluence case.
  if (!setup$is_dendritic) {
    lca <- if (identical(setup$tree_orientation, "out")) {
      directional_lca_edge_out_tree(setup, point_s[1], point_t[1])
    } else {
      directional_lca_edge_nondendritic(setup, point_s[1], point_t[1])
    }
    if (!is.null(lca)) {
      if (!is.null(lca$root)) {
        # The two branches share a common source vertex (indegree 0) but no
        # common upstream *edge* -- which, by construction, only happens
        # when that source has outdegree > 1 (a single out-edge would
        # itself be a shared ancestor edge, so directional_lca_edge_nondendritic()
        # would have returned an `edge`, not a `root`, above). There is no
        # upstream edge to anchor a single lca$point to, and
        # c(lca_edge, 0) only reaches the lca_edge branch itself via the
        # forward-BFS transfer walk, not its sibling branch(es). Narrow,
        # known gap -- not silently mishandled.
        stop(sprintf(
          "directional_lca(): LCA at a branching source vertex (indegree 0, outdegree > 1) is not yet supported -- vertex %d",
          lca$root
        ))
      }
      lca_edge <- lca$edge
      return(list(point = c(lca_edge, setup$graph$edge_lengths[lca_edge])))
    }
  }

  # Case 4: no common ancestor.
  NULL
}

# Is to_edge reachable by walking forward (via out-edges) from from_edge's
# head vertex? Same BFS-over-out-edges shape as directional_transfer_walk(),
# but only a boolean reachability check -- no transfer-factor bookkeeping.
# Kept as the non-dendritic fallback for directional_edge_downstream_of();
# not itself touched by the O(1) dendritic fast path below. Uses the
# out_edges_by_vertex bucket (built once in directional_ou_setup_structure())
# instead of which(E[,1]==v) -- the latter is an O(nE) scan repeated at
# every BFS step, which is invisible on small fixtures but makes this
# function O(nE) *per step* on a real, heavily-branched network (measured:
# single-digit seconds per call on the 18,668-edge Mid-Columbia graph before
# this fix).
#' @noRd
directional_edge_reachable_forward <- function(E, out_edges_by_vertex, from_edge, to_edge) {
  visited <- out_edges_by_vertex[[as.character(E[from_edge, 2])]]
  frontier <- visited
  while (!(to_edge %in% visited) && length(frontier) > 0) {
    next_frontier <- integer(0)
    for (g in frontier) {
      children <- setdiff(out_edges_by_vertex[[as.character(E[g, 2])]], visited)
      next_frontier <- c(next_frontier, children)
    }
    visited <- c(visited, next_frontier)
    frontier <- next_frontier
  }
  to_edge %in% visited
}

# Is to_edge reachable forward from from_edge (i.e. from_edge is upstream of
# or equal to to_edge)? O(1) via orientation-aware Euler interval containment
# for in/out trees; otherwise use the explicit BFS fallback.
#' @noRd
directional_edge_downstream_of <- function(setup, from_edge, to_edge) {
  if (!is.null(setup$enter)) {
    if (identical(setup$tree_orientation, "out")) {
      setup$enter[from_edge] <= setup$enter[to_edge] &&
        setup$exit[to_edge] <= setup$exit[from_edge]
    } else {
      setup$enter[to_edge] <= setup$enter[from_edge] &&
        setup$exit[from_edge] <= setup$exit[to_edge]
    }
  } else {
    directional_edge_reachable_forward(setup$graph$E, setup$out_edges_by_vertex, from_edge, to_edge)
  }
}

# Deepest common ancestor of edge_s and edge_t, for the non-dendritic case 3
# above (case 2 already ruled out one containing the other). Walks upstream
# (via in-edges) from each edge's tail vertex to build its ancestor edge-set
# -- upstream can itself branch at confluences (indegree > 1), so the BFS
# may fan out over multiple in-edges per step. Each walk also records the
# source vertex(es) (indegree 0) it dead-ends at, since a source has no key
# in in_edges_by_vertex to expand from.
#
# Returns list(edge = <edge or NULL>, root = <vertex or NULL>):
#   - If the two edge-sets share an edge, the deepest (latest in
#     setup$topo_order) is a genuine, edge-anchorable LCA: `edge` is set.
#   - Else, if the two walks share a dead-end source vertex, that source is
#     the true (vertex-level) LCA but has no upstream edge to anchor to:
#     `root` is set. (That source is guaranteed to have outdegree > 1: a
#     single out-edge would itself be a shared ancestor edge, which the
#     first branch above would already have caught.)
#   - Else NULL: the ancestor sets are genuinely disjoint (case 4, no
#     common ancestor at all).
#
# Uses the in_edges_by_vertex bucket instead of which(E[,2]==v) -- see
# directional_edge_reachable_forward()'s comment for why the which()-in-
# a-BFS-loop pattern is a real, not cosmetic, bottleneck on this graph.
#' @noRd
directional_lca_edge_nondendritic <- function(setup, edge_s, edge_t) {
  E <- setup$graph$E
  in_edges_by_vertex <- setup$in_edges_by_vertex
  V_indegree <- setup$V_indegree
  ancestors_of <- function(edge) {
    start_tail <- E[edge, 1]
    roots <- if (V_indegree[start_tail] == 0) as.character(start_tail) else character(0)
    visited <- in_edges_by_vertex[[as.character(start_tail)]]
    frontier <- visited
    while (length(frontier) > 0) {
      next_frontier <- integer(0)
      for (f in frontier) {
        tail_f <- E[f, 1]
        if (V_indegree[tail_f] == 0) {
          roots <- c(roots, as.character(tail_f))
        } else {
          next_frontier <- c(next_frontier, setdiff(in_edges_by_vertex[[as.character(tail_f)]], visited))
        }
      }
      visited <- c(visited, next_frontier)
      frontier <- next_frontier
    }
    list(edges = visited, roots = unique(roots))
  }

  anc_s <- ancestors_of(edge_s)
  anc_t <- ancestors_of(edge_t)

  common_edges <- intersect(anc_s$edges, anc_t$edges)
  if (length(common_edges) > 0) {
    lca_edge <- common_edges[which.max(match(common_edges, setup$topo_order))]
    return(list(edge = lca_edge, root = NULL))
  }

  common_roots <- intersect(anc_s$roots, anc_t$roots)
  if (length(common_roots) > 0) {
    return(list(edge = NULL, root = as.integer(common_roots[1])))
  }

  NULL
}

#' Last common ancestor edge on an out-tree
#'
#' Uses precomputed parent pointers and depths to align two edges, then climbs
#' them in lockstep until they meet. The caller has already ruled out either
#' edge being an ancestor of the other, so a match is the third edge above a
#' genuine diffluence.
#'
#' @param setup A setup list from `directional_ou_setup()` with
#'  `tree_orientation == "out"`.
#' @param edge_s,edge_t Edge numbers.
#' @return `list(edge, root)` for an edge-anchorable LCA, or `NULL` for
#'  different components. A branching-source-only LCA remains unsupported.
#' @noRd
directional_lca_edge_out_tree <- function(setup, edge_s, edge_t) {
  depth <- setup$depth
  parent_edge <- setup$parent_edge
  E <- setup$graph$E

  a <- edge_s
  b <- edge_t
  while (depth[a] > depth[b]) a <- parent_edge[a]
  while (depth[b] > depth[a]) b <- parent_edge[b]

  while (a != b) {
    if (depth[a] == 0L) {
      if (E[a, 1] == E[b, 1]) {
        stop(sprintf(
          paste0(
            "directional_lca(): LCA at a branching source vertex ",
            "(indegree 0, outdegree > 1) is not yet supported -- vertex %d"
          ),
          E[a, 1]
        ))
      }
      return(NULL)
    }
    a <- parent_edge[a]
    b <- parent_edge[b]
  }

  list(edge = a, root = NULL)
}

#' Marginal variance at arbitrary points (vectorized)
#'
#' `r(x,x)` at each row of `PtE_abs`, via the same intra-edge variance
#' recursion used to build `var_head` in `directional_ou_setup()`
#' (`r(head) = exp(-2*kappa*h)*v0 + (1-exp(-2*kappa*h))*sigma_stationary`),
#' evaluated at the point's own distance along the edge rather than only at
#' the edge head.
#'
#' @param setup A setup list from `directional_ou_setup()`.
#' @param PtE_abs A matrix with columns `(edge_number, absolute_distance)`,
#'  one row per point (may be a single row).
#' @return A numeric vector of length `nrow(PtE_abs)`.
#' @noRd
directional_point_variance <- function(setup, PtE_abs) {
  PtE_abs <- as.matrix(PtE_abs)
  e <- PtE_abs[, 1]
  t <- PtE_abs[, 2]
  decay <- exp(-2 * setup$kappa * t)
  as.numeric(decay * setup$var_tail[e] + (1 - decay) * setup$sigma_stationary)
}

#' Transfer factor A(from, to) by explicit forward path walk
#'
#' Fallback for the non-dendritic case (and a cross-check for the fast
#' dendritic path in `directional_log_transfer()`): finds the unique forward
#' edge path from `from`'s edge to `to`'s edge by BFS through the directed
#' tree, then accumulates the `-kappa * length` legs and the `beta_v`
#' crossings edge by edge. `from` must be upstream of (or equal to) `to` on
#' the unique directed tree path (caller's responsibility, not re-verified
#' here). O(depth) per call -- this is an explicit fallback, not the fast
#' path, and is not meant to be optimized further.
#'
#' @param setup A setup list from `directional_ou_setup()`.
#' @param from,to Length-2 `c(edge_number, absolute_distance)` vectors.
#' @return `list(sign, log_abs)` with `sign * exp(log_abs) == A(from, to)`.
#' @noRd
directional_transfer_walk <- function(setup, from, to) {
  kappa <- setup$kappa
  E <- setup$graph$E
  edge_lengths <- setup$graph$edge_lengths
  out_edges_by_vertex <- setup$out_edges_by_vertex

  if (from[1] == to[1]) {
    return(list(sign = 1, log_abs = -kappa * (to[2] - from[2])))
  }

  # BFS forward from from's edge (starting at its head vertex) until to's
  # edge is discovered, recording a parent-edge map for path reconstruction.
  # Uses the out_edges_by_vertex bucket instead of which(E[,1]==v) -- see
  # directional_edge_reachable_forward()'s comment for why the which()-in-
  # a-BFS-loop pattern is a real, not cosmetic, bottleneck on this graph.
  parent_edge <- list()
  frontier <- out_edges_by_vertex[[as.character(E[from[1], 2])]]
  for (g in frontier) parent_edge[[as.character(g)]] <- from[1]
  visited_edges <- frontier
  while (!(to[1] %in% visited_edges) && length(frontier) > 0) {
    next_frontier <- integer(0)
    for (g in frontier) {
      children <- setdiff(out_edges_by_vertex[[as.character(E[g, 2])]], visited_edges)
      for (child in children) parent_edge[[as.character(child)]] <- g
      next_frontier <- c(next_frontier, children)
    }
    visited_edges <- c(visited_edges, next_frontier)
    frontier <- next_frontier
  }
  if (!(to[1] %in% visited_edges)) {
    stop("directional_transfer_walk(): to's edge is not reachable forward from from's edge; from must be upstream of to.")
  }

  # Reconstruct the edge sequence from from's edge to to's edge by following
  # parent pointers backward from to's edge, then reversing.
  edge_seq <- to[1]
  while (edge_seq[1] != from[1]) {
    edge_seq <- c(parent_edge[[as.character(edge_seq[1])]], edge_seq)
  }
  edge_seq <- edge_seq[-1]  # drop from[1] itself; keep only descendants

  log_abs <- -kappa * (edge_lengths[from[1]] - from[2])  # remainder of from's edge
  sign_acc <- 1
  prev <- from[1]
  for (g in edge_seq) {
    beta <- -setup$in_by_edge[prev] / setup$out_by_edge[g]
    log_abs <- log_abs + log(abs(beta))
    sign_acc <- sign_acc * sign(beta)
    traversal <- if (g != to[1]) edge_lengths[g] else to[2]
    log_abs <- log_abs - kappa * traversal
    prev <- g
  }
  list(sign = sign_acc, log_abs = log_abs)
}

#' Log-transfer factor A(from, to) -- fast oriented-tree path or explicit walk
#'
#' Dispatches to the O(1) root-normalised accumulators (`logG_tail`/
#' `signG_tail` from `directional_ou_setup()`) for in/out trees. In-trees
#' use outlet-normalised ratios; out-trees use source-normalised ratios in
#' the opposite subtraction order. Irregular trees retain the explicit walk.
#'
#' @param setup A setup list from `directional_ou_setup()`.
#' @param from,to Length-2 `c(edge_number, absolute_distance)` vectors, with
#'  `from` upstream of (or equal to) `to` on the unique directed tree path
#'  (caller's responsibility, not re-verified here).
#' @return `list(sign, log_abs)` with `sign * exp(log_abs) == A(from, to)`.
#' @noRd
directional_log_transfer <- function(setup, from, to) {
  if (identical(setup$tree_orientation, "out")) {
    logG_at <- function(point) {
      setup$logG_tail[point[1]] - setup$kappa * point[2]
    }
    log_abs <- logG_at(to) - logG_at(from)
    sign_val <- setup$signG_tail[from[1]] * setup$signG_tail[to[1]]
    return(list(sign = sign_val, log_abs = log_abs))
  }

  if (!setup$is_dendritic) {
    return(directional_transfer_walk(setup, from, to))
  }

  # logG_tail[e] + kappa*t equals logG_tail[e] at t=0 and logG_head[e] at
  # t=edge_lengths[e], since logG_tail[e] = logG_head[e] - kappa*length(e).
  logG_at <- function(point) setup$logG_tail[point[1]] + setup$kappa * point[2]

  log_abs <- logG_at(from) - logG_at(to)
  sign_val <- setup$signG_tail[from[1]] * setup$signG_tail[to[1]]
  list(sign = sign_val, log_abs = log_abs)
}
