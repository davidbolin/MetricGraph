# Internal, read-only construction for integer-order graph constraints.
# Columns are endpoint derivative blocks in edge order: start(e), end(e).
# These are the unweighted constraints, including degree-one Neumann rows
# for alpha >= 2 (the edge_constraint = TRUE convention of buildC()).
# No existing C / CoB cache or graph methods are used here.

build_kirchhoff_constraint_matrix <- function(E, nV, alpha) {
  if (!is.numeric(alpha) || is.complex(alpha) || length(alpha) != 1L ||
      is.na(alpha) || !(alpha %in% 1:3)) {
    stop("alpha must be one of 1, 2, or 3")
  }
  alpha <- as.integer(alpha)
  if (!is.numeric(nV) || is.complex(nV) || length(nV) != 1L || !is.finite(nV) ||
      nV < 0 || nV > .Machine$integer.max || nV != floor(nV)) {
    stop("nV must be a nonnegative integer")
  }
  nV <- as.integer(nV)
  if (!is.matrix(E) || !is.numeric(E) || is.complex(E) || ncol(E) != 2L ||
      any(!is.finite(E)) || any(E != floor(E)) ||
      any(E < 1 | E > nV)) {
    stop("E must be a two-column matrix of vertex indices in 1:nV")
  }

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


build_constraint_matrices <- function(E, nV, PtV) {
  # Also validates the topology before using it to select endpoint values.
  K1 <- build_kirchhoff_constraint_matrix(E, nV, alpha = 1)
  K2 <- build_kirchhoff_constraint_matrix(E, nV, alpha = 2)
  K3 <- build_kirchhoff_constraint_matrix(E, nV, alpha = 3)

  vertices <- PtV
  if (is.null(vertices)) {
    stop("PtV is NULL; run observation_to_vertex() before build_constraints()")
  }
  if (!is.numeric(vertices) || is.complex(vertices) ||
      !is.null(dim(vertices)) || any(!is.finite(vertices)) ||
      any(vertices != floor(vertices)) ||
      any(vertices < 1 | vertices > nV)) {
    stop("PtV must be a vector of vertex indices in 1:nV")
  }

  # Match the existing vertex ordering: prefer a start endpoint and fall
  # back to an end endpoint. Compute the incidences once for all alphas.
  start_edges <- match(vertices, E[, 1])
  endpoints <- 2L * start_edges - 1L
  use_end <- is.na(start_edges)
  endpoints[use_end] <- 2L * match(vertices[use_end], E[, 2])
  if (anyNA(endpoints)) {
    stop("Every vertex in PtV must be incident to an edge")
  }

  build_case <- function(K, alpha) {
    if (nrow(K) == 0L) {
      T <- Matrix::Diagonal(ncol(K))
      S <- numeric(0)
    } else {
      CoB <- c_basis2(K)
      # The native routine returns basis columns; T uses basis rows.
      T <- t(CoB$T)
      S <- CoB$S
    }
    keep <- setdiff(seq_len(nrow(T)), seq_along(S))
    Tc <- T[keep, , drop = FALSE]
    indices <- alpha * (endpoints - 1L) + 1L
    A <- t(Tc)[indices, , drop = FALSE]

    return(list(K = K, T = T, S = S, keep = keep,
                Tc = Tc, indices = indices, A = A))
  }

  return(list(
    alpha1 = build_case(K1, alpha = 1),
    alpha2 = build_case(K2, alpha = 2),
    alpha3 = build_case(K3, alpha = 3)
  ))
}
