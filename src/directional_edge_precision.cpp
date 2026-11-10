#include <Rcpp.h>

using namespace Rcpp;

// Assemble the independent-edge alpha=1 precision in triplet form. Boundary
// validation and source-to-edge resolution stay in R, where the metric_graph
// object is available; this kernel only performs the numeric O(nE) work.
//
// [[Rcpp::export]]
List directional_edge_precision_triplets_cpp(
    const NumericVector& edge_lengths,
    const IntegerVector& stationary_edges,
    double tau,
    double kappa,
    double w) {
  const int n_edges = edge_lengths.size();
  const int n_stationary = stationary_edges.size();
  const int n_triplets = 4 * n_edges + n_stationary;

  IntegerVector i(n_triplets);
  IntegerVector j(n_triplets);
  NumericVector x(n_triplets);
  const double scale = 2.0 * kappa * tau * tau;

  int offset = 0;
  for (int edge = 0; edge < n_edges; ++edge) {
    const double decay = std::exp(-kappa * edge_lengths[edge]);
    const double decay_squared = decay * decay;
    const double one_minus_decay_squared = 1.0 - decay_squared;
    const double upper = w + decay_squared / one_minus_decay_squared;
    const double lower = (1.0 - w) + decay_squared / one_minus_decay_squared;
    const double off_diagonal = -decay / one_minus_decay_squared;
    const int tail = 2 * edge + 1;
    const int head = tail + 1;

    i[offset] = tail;
    j[offset] = tail;
    x[offset++] = scale * upper;

    i[offset] = head;
    j[offset] = head;
    x[offset++] = scale * lower;

    i[offset] = tail;
    j[offset] = head;
    x[offset++] = scale * off_diagonal;

    i[offset] = head;
    j[offset] = tail;
    x[offset++] = scale * off_diagonal;
  }

  for (int index = 0; index < n_stationary; ++index) {
    const int edge = stationary_edges[index];
    if (edge < 1 || edge > n_edges) {
      stop("stationary_edges contains an invalid edge index.");
    }
    const int tail = 2 * (edge - 1) + 1;
    i[offset] = tail;
    j[offset] = tail;
    x[offset++] = scale * (1.0 - w);
  }

  return List::create(
      _["i"] = i,
      _["j"] = j,
      _["x"] = x,
      _["dims"] = IntegerVector::create(2 * n_edges, 2 * n_edges));
}
