# Fractional backend, following the same dispatch and build conventions as
# Qalpha1 / Qalpha2. The public spde_precision() validates common arguments.
Qfractional <- function(theta, alpha, graph, build = TRUE, m = 2,
                        type_rational_approx = "chebfun", type_interp = "spline") {
  if (inherits(graph, "graph_components")) {
    if (!build) {
      stop("build = FALSE is not supported for 'graph_components'.")
    }
    Qs <- lapply(graph$graphs, function(g) {
      Qfractional(theta = theta, alpha = alpha, graph = g, build = TRUE,
                  m = m, type_rational_approx = type_rational_approx,
                  type_interp = type_interp)
    })
    # Compatibility path for the deprecated graph_components class only.
    # Keep Qi$Gj matrices separate, with Gj following graph$graphs order.
    # A disconnected metric_graph uses the full-graph path below instead.
    graph_names <- paste0("G", seq_along(Qs))
    components <- lapply(names(Qs[[1L]]), function(component) {
      stats::setNames(lapply(Qs, function(Q) Q[[component]]), graph_names)
    })
    names(components) <- names(Qs[[1L]])
    return(components)
  }

  # Connected and disconnected metric_graph objects share this path. Each
  # Qi retains all vertices or edge endpoints in the original graph order,
  # matching constraint templates built on the same full graph. Do not split
  # or reorder connected pieces here; get_components() has its own ordering.
  Qraw <- fractional_spde_precision_components(
    kappa = theta[2], tau = theta[1], alpha = alpha, graph = graph, m = m,
    type_rational_approx = type_rational_approx, type_interp = type_interp
  )
  if (!build) return(Qraw)
  return(lapply(Qraw, function(entries) {
    Matrix::sparseMatrix(i = entries$i, j = entries$j,
                         x = entries$x, dims = entries$dims)
  }))
}

# Scaled boundaryless component precisions, before any graph conditioning.
# Each component is a list(i, j, x, dims), with both matrix triangles.
# Keeping components separate fixes their order and allows different
# constraint bases, without assembling sparse matrices during construction.
fractional_spde_precision_components <- function(
    kappa, tau, alpha, graph, m,
    type_rational_approx = "chebfun", type_interp = "spline") {
  if (!is.numeric(m) || is.complex(m) || length(m) != 1L ||
      !is.finite(m) || m < 1 || m != floor(m)) {
    stop("m must be a positive integer")
  }
  if (!is.character(type_rational_approx) || length(type_rational_approx) != 1L ||
      is.na(type_rational_approx) ||
      !type_rational_approx %in% c("chebfun", "brasil", "chebfunLB")) {
    stop("type_rational_approx must be 'chebfun', 'brasil', or 'chebfunLB'")
  }
  if (!is.character(type_interp) || length(type_interp) != 1L ||
      is.na(type_interp) ||
      !type_interp %in% c("linear", "log", "spline", "logspline")) {
    stop("type_interp must be 'linear', 'log', 'spline', or 'logspline'")
  }
  lengths <- validate_spde_precision_graph(graph)
  coeff <- tryCatch(
    rSPDE:::interp_rational_coefficients(
      order = m, type_rational_approx = type_rational_approx,
      type_interp = type_interp, alpha = alpha
    ),
    error = function(e) {
      stop("Could not obtain rational coefficients: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  # The m+1-component representation requires a positive constant and
  # residues, with negative poles. Reject unsupported or degenerate tables
  # before dividing by coefficients or constructing endpoint covariances.
  if (length(coeff$k) != 1L || !is.finite(coeff$k) || coeff$k <= 0 ||
      length(coeff$r) != m || length(coeff$p) != m ||
      any(!is.finite(coeff$r)) || any(coeff$r <= 0) ||
      any(!is.finite(coeff$p)) || any(coeff$p >= 0)) {
    stop("Rational coefficients must have positive k and residues, and negative poles")
  }
  nu <- alpha - 1/2
  sigma2 <- gamma(nu) /
    (tau^2 * kappa^(2*nu) * sqrt(4*pi) * gamma(nu + 1/2))
  if (!is.finite(sigma2) || sigma2 <= 0) {
    stop("The supplied kappa and tau produce a non-finite or zero variance scale")
  }
  c_alpha <- gamma(alpha) / gamma(alpha - 1/2)
  q <- floor(alpha)

  if (q == 0) {
    indices <- seq_len(graph$nV)
    value <- 1 / (coeff$k * c_alpha * sqrt(4*pi) * sigma2 / kappa)
    Q0 <- list(i = indices, j = indices, x = rep(value, graph$nV),
               dims = rep(as.integer(graph$nV), 2L))
  } else {
    c_q <- gamma(q) / gamma(q - 1/2)
    Q0 <- spde_endpoint_precision_triplets(lengths, kappa, p = 0, alpha = q)
    Q0$x <- Q0$x / (coeff$k * c_alpha * sigma2 / c_q)
  }
  Qi <- lapply(seq_along(coeff$p), function(i) {
    entries <- spde_endpoint_precision_triplets(lengths, kappa, coeff$p[i], alpha)
    entries$x <- entries$x / (coeff$r[i] * sigma2)
    entries
  })
  Qraw <- c(list(Q0), Qi)
  names(Qraw) <- paste0("Q", seq_along(Qraw) - 1L)
  if (any(vapply(Qraw, function(Q) any(!is.finite(Q$x)), logical(1)))) {
    stop("Non-finite fractional precision; check the parameters and edge lengths")
  }
  Qraw
}

# Invert in standardized derivative coordinates, then restore the original
# precision units: C = D R D implies C^{-1} = D^{-1} R^{-1} D^{-1}.
# This removes coordinate-scale imbalance, not dependence between endpoints.
# Keep solve()'s default singularity check for intrinsically singular inputs.
invert_spde_covariance <- function(covariance) {
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      is.complex(covariance) || nrow(covariance) == 0L ||
      nrow(covariance) != ncol(covariance) || any(!is.finite(covariance))) {
    stop("Covariance must be a finite, nonempty, real square matrix", call. = FALSE)
  }
  variances <- diag(covariance)
  if (any(variances <= 0)) {
    stop("Covariance must have strictly positive diagonal entries", call. = FALSE)
  }
  scales <- outer(sqrt(variances), sqrt(variances))
  if (any(!is.finite(scales)) || any(scales <= 0)) {
    stop("Covariance scales must be finite and strictly positive", call. = FALSE)
  }
  precision <- solve(covariance / scales) / scales
  if (any(!is.finite(precision))) {
    stop("Non-finite covariance precision; check covariance scaling", call. = FALSE)
  }
  return(precision)
}

# Work only with small dense endpoint covariances (at most 6 x 6). The
# stationary covariance comes from rSPDE, and half the same-point precision
# is removed at each endpoint to give the boundaryless edge contribution.
spde_endpoint_precision_triplets <- function(lengths, kappa, p, alpha) {
  R00 <- rSPDE:::matern.p.joint(s = 0, t = 0, kappa = kappa, p = p, alpha = alpha)
  n_state <- nrow(R00)
  indices <- seq_len(n_state)
  R00_inverse <- invert_spde_covariance(R00)
  correction <- matrix(0, 2L * n_state, 2L * n_state)
  correction[indices, indices] <- R00_inverse / 2
  correction[n_state + indices, n_state + indices] <- R00_inverse / 2

  # Reuse identical edge lengths within this parameter evaluation.
  distinct_lengths <- unique(lengths)
  blocks <- lapply(distinct_lengths, function(l_e) {
    if (n_state == 1L) {
      # The scalar state has an exponential covariance. expm1 avoids loss
      # of accuracy when the endpoints are close relative to the range.
      rate <- kappa * sqrt(1 - p)
      rho <- exp(-rate * l_e)
      denominator <- -expm1(-2 * rate * l_e)
      block <- matrix(c(1, -rho, -rho, 1), 2L, 2L) /
        (R00[1, 1] * denominator)
    } else {
      R0l <- rSPDE:::matern.p.joint(
        s = 0, t = l_e, kappa = kappa, p = p, alpha = alpha
      )
      covariance <- rbind(cbind(R00, R0l), cbind(t(R0l), R00))
      block <- invert_spde_covariance(covariance)
    }
    block <- block - correction
    # Retain the same upper-triangle convention as the stationary precision
    # assembler and emit its mirror explicitly for general triplet storage.
    block[lower.tri(block)] <- t(block)[lower.tri(block)]
    nonzero <- which(block != 0, arr.ind = TRUE)
    list(i = nonzero[, 1], j = nonzero[, 2], x = block[nonzero],
         dims = rep(2L * n_state, 2L))
  })
  combine_spde_precision_triplets(blocks[match(lengths, distinct_lengths)])
}

# Concatenate block-diagonal triplets by offsetting their local indices.
# Used for edge blocks within each component; rational components stay separate.
# No sparse matrix is created here, including for the build = FALSE path.
combine_spde_precision_triplets <- function(blocks) {
  nrow_total <- sum(vapply(blocks, function(block) block$dims[1], numeric(1)))
  ncol_total <- sum(vapply(blocks, function(block) block$dims[2], numeric(1)))
  n_entries <- sum(vapply(blocks, function(block) length(block$x), numeric(1)))
  if (any(c(nrow_total, ncol_total, n_entries) > .Machine$integer.max)) {
    stop("SPDE precision exceeds the supported sparse matrix size")
  }
  i <- integer(n_entries)
  j <- integer(n_entries)
  x <- numeric(n_entries)
  row_offset <- col_offset <- entry_offset <- 0L
  for (block in blocks) {
    n <- length(block$x)
    if (n > 0L) {
      positions <- entry_offset + seq_len(n)
      i[positions] <- block$i + row_offset
      j[positions] <- block$j + col_offset
      x[positions] <- block$x
    }
    entry_offset <- entry_offset + n
    row_offset <- row_offset + block$dims[1]
    col_offset <- col_offset + block$dims[2]
  }
  list(i = i, j = j, x = x, dims = as.integer(c(nrow_total, ncol_total)))
}

# Validate public arguments independently of the chosen backend. The existing
# integer assemblers still receive the same parameter values and defaults.
validate_spde_precision_arguments <- function(kappa, tau, alpha, graph, BC, build) {
  if (!inherits(graph, c("metric_graph", "graph_components"))) {
    stop("graph must be a 'metric_graph' or 'graph_components' object.")
  }
  if (!is.numeric(alpha) || is.complex(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0.5 || alpha > 3) {
    stop("alpha must be strictly greater than 0.5 and at most 3")
  }
  if (!is.numeric(kappa) || is.complex(kappa) || length(kappa) != 1L ||
      !is.finite(kappa) || kappa < 0 || (kappa == 0 && alpha != 1)) {
    stop("kappa must be positive; zero is supported only for alpha = 1")
  }
  if (!is.numeric(tau) || is.complex(tau) || length(tau) != 1L ||
      !is.finite(tau) || tau <= 0) {
    stop("tau must be a finite positive number")
  }
  # Only alpha one and two inspect BC; accept the package's boundary codes.
  if (alpha %in% c(1, 2) &&
      ((!is.numeric(BC) && !is.logical(BC)) || is.complex(BC) ||
       length(BC) != 1L || !is.finite(BC) || !BC %in% 0:3)) {
    stop("BC must be one of 0, 1, 2, or 3 for alpha = 1 or 2")
  }
  if (!is.logical(build) || length(build) != 1L || is.na(build)) {
    stop("build must be TRUE or FALSE")
  }
  invisible(NULL)
}

# Shared graph validation for raw endpoint construction. No graph mutation
# or observation data is required for either fractional alpha or alpha 3.
validate_spde_precision_graph <- function(graph) {
  for (name in c("nV", "nE")) {
    value <- graph[[name]]
    if (!is.numeric(value) || is.complex(value) || length(value) != 1L ||
        !is.finite(value) || value < 0 || value != floor(value)) {
      stop("graph$", name, " must be a nonnegative integer")
    }
  }
  E <- graph$E
  if (!is.matrix(E) || !is.numeric(E) || is.complex(E) ||
      !identical(as.integer(dim(E)), as.integer(c(graph$nE, 2L))) ||
      any(!is.finite(E)) || any(E != floor(E)) ||
      any(E < 1 | E > graph$nV)) {
    stop("graph$E must contain one pair of valid vertex indices per edge")
  }
  lengths <- graph$edge_lengths
  if (!is.numeric(lengths) || is.complex(lengths) ||
      any(!is.finite(lengths)) || any(lengths <= 0) ||
      length(lengths) != graph$nE) {
    stop("graph$edge_lengths must contain a finite positive length for each edge")
  }
  lengths
}
