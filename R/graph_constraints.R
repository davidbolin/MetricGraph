# Internal, read-only construction for integer bases and fractional templates.
# Columns are endpoint derivative blocks in edge order: start(e), end(e).
# Non-directional setups use unweighted constraints, including degree-one
# Neumann rows for alpha >= 2 (edge_constraint = TRUE in buildC()).
# No existing C / CoB cache or graph methods are used here.

validate_constraint_case <- function(alpha, directional = FALSE) {
  if (!is.numeric(alpha) || is.complex(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0.5 || alpha > 3) {
    stop("alpha must be finite, greater than 0.5, and at most 3")
  }
  if (!is.logical(directional) || length(directional) != 1L || is.na(directional)) {
    stop("directional must be TRUE or FALSE")
  }
  if (directional && alpha != 1) {
    stop("directional constraints are only supported for alpha = 1")
  }
  invisible(NULL)
}

validate_constraint_topology <- function(E, nV) {
  if (!is.numeric(nV) || is.complex(nV) || length(nV) != 1L || !is.finite(nV) ||
      nV < 0 || nV > .Machine$integer.max || nV != floor(nV)) {
    stop("nV must be a nonnegative integer")
  }
  if (!is.matrix(E) || !is.numeric(E) || is.complex(E) || ncol(E) != 2L ||
      any(!is.finite(E)) || any(E != floor(E)) ||
      any(E < 1 | E > nV)) {
    stop("E must be a two-column matrix of vertex indices in 1:nV")
  }
  invisible(NULL)
}

build_kirchhoff_constraint_matrix <- function(E, nV, alpha) {
  validate_constraint_case(alpha)
  if (alpha != floor(alpha)) {
    stop("integer constraint matrices require alpha = 1, 2, or 3")
  }
  validate_constraint_topology(E, nV)
  alpha <- as.integer(alpha)
  nV <- as.integer(nV)

  # Stable radix ordering groups incidences while retaining their endpoint
  # order within each vertex. A loop contributes two distinct incidences.
  endpoints <- as.integer(c(t(E)))
  ordered <- order(endpoints, method = "radix")
  left <- head(ordered, -1L)
  right <- tail(ordered, -1L)
  same_vertex <- endpoints[left] == endpoints[right]
  left <- left[same_vertex]
  right <- right[same_vertex]
  n_continuity <- length(left)
  continuity_rows <- rep(seq_len(n_continuity), 2L)
  continuity_values <- rep(c(1, -1), each = n_continuity)

  # Equality between consecutive incidences gives d_v - 1 independent
  # continuity rows at a vertex of degree d_v, with no rows at degree one.
  i <- continuity_rows
  j <- alpha * (c(left, right) - 1L) + 1L
  x <- continuity_values
  n_constraints <- n_continuity

  if (alpha >= 2L) {
    # Kirchhoff rows include degree-one vertices (Neumann conditions).
    # Isolated vertices contribute no rows because they have no edge ends.
    active_vertices <- which(tabulate(endpoints, nbins = nV) > 0L)
    vertex_rows <- integer(nV)
    vertex_rows[active_vertices] <- seq_along(active_vertices)
    i <- c(i, n_constraints + vertex_rows[endpoints])
    j <- c(j, alpha * (seq_along(endpoints) - 1L) + 2L)
    x <- c(x, rep(c(1, -1), times = nrow(E)))
    n_constraints <- n_constraints + length(active_vertices)
  }

  if (alpha == 3L) {
    # Even derivatives have the same orientation at both edge ends.
    i <- c(i, n_constraints + continuity_rows)
    j <- c(j, 3L * (c(left, right) - 1L) + 3L)
    x <- c(x, continuity_values)
    n_constraints <- n_constraints + n_continuity
  }

  return(Matrix::sparseMatrix(
    i = i, j = j, x = x,
    dims = c(n_constraints, alpha * length(endpoints))
  ))
}


# Common row-basis representation. No vertex selector is needed here.
constraint_basis <- function(K) {
  if (nrow(K) == 0L) {
    T <- Matrix::Diagonal(ncol(K))
    S <- numeric(0)
  } else {
    CoB <- c_basis2(K)
    T <- t(CoB$T)
    S <- CoB$S
  }
  keep <- setdiff(seq_len(nrow(T)), seq_along(S))
  Tc <- T[keep, , drop = FALSE]
  list(K = K, T = T, S = S, keep = keep, Tc = Tc)
}

# Build exactly one non-directional order, including its PtV observation map.
build_constraint_matrices <- function(E, nV, PtV, alpha) {
  K <- build_kirchhoff_constraint_matrix(E, nV, alpha)
  if (is.null(PtV)) {
    stop("PtV is NULL; run observation_to_vertex() before build_constraints()")
  }
  if (!is.numeric(PtV) || is.complex(PtV) ||
      !is.null(dim(PtV)) || any(!is.finite(PtV)) ||
      any(PtV != floor(PtV)) || any(PtV < 1 | PtV > nV)) {
    stop("PtV must be a vector of vertex indices in 1:nV")
  }

  # Preserve gives.indices(): prefer the first start, then the first end.
  start_edges <- match(PtV, E[, 1])
  endpoints <- 2L * start_edges - 1L
  use_end <- is.na(start_edges)
  endpoints[use_end] <- 2L * match(PtV[use_end], E[, 2])
  if (anyNA(endpoints)) {
    stop("Every vertex in PtV must be incident to an edge")
  }

  setup <- constraint_basis(K)
  setup$indices <- as.integer(alpha) * (endpoints - 1L) + 1L
  setup$A <- t(setup$Tc)[setup$indices, , drop = FALSE]
  setup
}

# Store the first component and the repeated-component templates separately.
# No rational order, rational coefficients or large joint matrix is needed.
build_fractional_constraint_templates <- function(E, nV, PtV, alpha,
                                                  integer_setups = NULL) {
  validate_constraint_case(alpha)
  if (alpha == floor(alpha)) {
    stop("fractional constraint templates require noninteger alpha")
  }
  get_setup <- function(order) {
    setup <- integer_setups[[paste0("alpha", order)]]
    if (is.null(setup)) setup <- build_constraint_matrices(E, nV, PtV, order)
    setup
  }
  q <- floor(alpha)
  # This also validates topology, PtV and the selected endpoint incidences
  # when the upper-order setup is not already present in the valid cache.
  upper <- get_setup(q + 1L)
  if (q == 0) {
    first_Tc <- Matrix::Diagonal(nV)
    first_A <- Matrix::sparseMatrix(
      i = seq_along(PtV), j = PtV, x = 1,
      dims = c(length(PtV), nV)
    )
  } else {
    lower <- get_setup(q)
    first_Tc <- lower$Tc
    first_A <- lower$A
  }
  list(
    Tc_joint = list(Tc1 = first_Tc, Tc2 = upper$Tc),
    A = list(A1 = first_A, A2 = upper$A)
  )
}

# Same construction as buildDirectionalConstraints(alpha = 1), with local
# variables only. Keep the legacy method and its C / CoB side effects separate.
build_directional_constraint_matrices <- function(
    E, nV, weight, DirectionalWeightFunction_out, DirectionalWeightFunction_in) {
  validate_constraint_topology(E, nV)
  nE <- nrow(E)
  if (!is.numeric(weight) || is.complex(weight) || length(weight) != nE ||
      any(!is.finite(weight))) {
    stop("directional weights must be a finite numeric vector with one value per edge")
  }
  if (!is.function(DirectionalWeightFunction_out) ||
      !is.function(DirectionalWeightFunction_in)) {
    stop("directional weighting functions must be functions")
  }
  weight_vectors <- directional_weight_vectors(
    E = E, nE = nE, weight = weight,
    DirectionalWeightFunction_out = DirectionalWeightFunction_out,
    DirectionalWeightFunction_in = DirectionalWeightFunction_in
  )
  if (any(!is.finite(weight_vectors$out_by_edge)) ||
      any(!is.finite(weight_vectors$in_by_edge))) {
    stop("directional weighting functions must produce finite weights")
  }
  temp_E <- matrix(as.integer(E), ncol = 2)
  K <- construct_directional_constraint_matrix_fast(
    temp_E, as.integer(nV), as.integer(nE), 1L,
    as.integer(tabulate(E[, 2], nbins = nV)),
    as.integer(tabulate(E[, 1], nbins = nV)),
    as.numeric(weight_vectors$out_by_edge), as.numeric(weight_vectors$in_by_edge)
  )
  constraint_basis(K)
}

# Prevent R's partial matching from treating alpha1 as alpha1_directional.
#' @export
`$.metric_graph_constraints` <- function(x, name) {
  x[[name, exact = TRUE]]
}
