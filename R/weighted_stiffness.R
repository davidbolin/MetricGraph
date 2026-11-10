#' Assemble a weighted stiffness matrix on a graph mesh
#'
#' Uses a piecewise linear positive coefficient on each mesh edge. The
#' coefficient integral on an edge is its length times the mean of its two
#' endpoint values.
#'
#' @param edges Two-column integer matrix of mesh edge endpoints.
#' @param lengths Positive mesh edge lengths, in the same order as rows of
#'   `edges`.
#' @param weights Positive coefficient values at mesh nodes.
#' @return A symmetric sparse positive-semidefinite stiffness matrix.
#' @export
weighted_stiffness <- function(edges, lengths, weights) {
  if (!is.matrix(edges) || ncol(edges) != 2L ||
      anyNA(edges) || any(edges != as.integer(edges))) {
    stop("edges must be a two-column integer matrix.")
  }
  n <- length(weights)
  if (any(edges < 1L) || any(edges > n) || any(edges[, 1] == edges[, 2])) {
    stop("edges must connect distinct valid mesh nodes.")
  }
  if (length(lengths) != nrow(edges) || anyNA(lengths) ||
      any(!is.finite(lengths)) || any(lengths <= 0)) {
    stop("lengths must contain one positive finite value per edge.")
  }
  if (anyNA(weights) || any(!is.finite(weights)) || any(weights <= 0)) {
    stop("weights must be positive finite nodal values.")
  }
  g <- (weights[edges[, 1]] + weights[edges[, 2]]) / (2 * lengths)
  Matrix::sparseMatrix(
    i = c(edges[, 1], edges[, 2], edges[, 1], edges[, 2]),
    j = c(edges[, 1], edges[, 2], edges[, 2], edges[, 1]),
    x = c(g, g, -g, -g), dims = c(n, n)
  )
}
