// directional_ou_covariance.cpp
//
// C++ port of directional_ou_setup_numeric() (R/covariance_directional_ou.R):
// the kappa/tau/sigma_source-dependent per-edge recursion of the directional
// OU covariance setup (var_tail/var_head, logG_head/logG_tail/signG_head/
// signG_tail). Direct 1:1 translation -- see the R function's comments for
// the math; this file only translates, it does not redesign.
//
// Takes a "skeleton" List (directional_ou_skeleton()'s output, extended with
// topo_order) plus a flat, R-resolved source_var_by_vertex vector (see
// directional_ou_resolve_source_var() in R -- the sigma_source name-matching
// logic is resolved in R, once, up front, exactly like
// directional_weight_vectors() resolves the DirectionalWeightFunction_*
// closures before crossing into C++).
//
// Note: E's vertex entries and topo_order's edge entries are 1-indexed
// (R convention); they are converted to 0-indexed C++ array/vector indices
// on use below, but topo_order's *iteration order* is used as-is.

#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
List directional_ou_setup_numeric_cpp(List skeleton, double kappa, double tau,
                                       NumericVector source_var_by_vertex) {
  int nE = as<int>(skeleton["nE"]);
  NumericVector edge_lengths = as<NumericVector>(skeleton["edge_lengths"]);
  IntegerMatrix E = as<IntegerMatrix>(skeleton["E"]);
  IntegerVector V_indegree = as<IntegerVector>(skeleton["V_indegree"]);
  IntegerVector V_outdegree = as<IntegerVector>(skeleton["V_outdegree"]);
  NumericVector out_by_edge = as<NumericVector>(skeleton["out_by_edge"]);
  NumericVector in_by_edge = as<NumericVector>(skeleton["in_by_edge"]);
  IntegerVector topo_order = as<IntegerVector>(skeleton["topo_order"]);
  std::string tree_orientation = as<std::string>(skeleton["tree_orientation"]);

  int nV = V_indegree.size();

  // Bucket in-edges and out-edges by vertex, mirroring the R
  // split(seq_len(nE), E[, 2]) / split(seq_len(nE), E[, 1]) buckets --
  // built here, self-contained, rather than assuming the R-side
  // in_edges_by_vertex/out_edges_by_vertex lists are passed in.
  std::vector<std::vector<int>> in_edges_by_vertex(nV);
  std::vector<std::vector<int>> out_edges_by_vertex(nV);
  for (int e = 0; e < nE; ++e) {
    int tail = E(e, 0) - 1;  // 0-indexed tail vertex
    int head = E(e, 1) - 1;  // 0-indexed head vertex
    out_edges_by_vertex[tail].push_back(e);
    in_edges_by_vertex[head].push_back(e);
  }

  double sigma_stationary = 1.0 / (2.0 * kappa * tau * tau);

  // First pass, forward in topo_order: var_tail/var_head.
  NumericVector var_tail(nE), var_head(nE);
  for (int idx = 0; idx < nE; ++idx) {
    int e = topo_order[idx] - 1;  // 0-indexed edge
    int v = E(e, 0) - 1;          // tail vertex of e, 0-indexed
    if (V_indegree[v] == 0) {
      var_tail[e] = source_var_by_vertex[v];
    } else {
      const std::vector<int>& E_in = in_edges_by_vertex[v];
      double sum = 0.0;
      for (std::size_t k = 0; k < E_in.size(); ++k) {
        int f = E_in[k];
        double beta = in_by_edge[f] / out_by_edge[e];
        sum += beta * beta * var_head[f];
      }
      var_tail[e] = sum;
    }
    double c2 = std::exp(-2.0 * kappa * edge_lengths[e]);
    var_head[e] = c2 * var_tail[e] + (1.0 - c2) * sigma_stationary;
  }

  // Second pass: root-normalised log-transfer accumulators. In-trees are
  // normalised toward their outlet and processed in reverse topological
  // order. Out-trees are normalised from their source and processed forward.
  NumericVector logG_head(nE), logG_tail(nE), signG_head(nE), signG_tail(nE);
  if (tree_orientation == "out") {
    for (int idx = 0; idx < nE; ++idx) {
      int e = topo_order[idx] - 1;  // 0-indexed edge
      int v = E(e, 0) - 1;          // tail vertex of e, 0-indexed
      if (V_indegree[v] == 0) {
        logG_tail[e] = 0.0;
        signG_tail[e] = 1.0;
      } else {
        int f = in_edges_by_vertex[v][0];  // unique parent under indegree <= 1
        double beta_ef = -in_by_edge[f] / out_by_edge[e];
        logG_tail[e] = logG_head[f] + std::log(std::fabs(beta_ef));
        double sign_beta_ef = (beta_ef > 0.0) - (beta_ef < 0.0);
        signG_tail[e] = signG_head[f] * sign_beta_ef;
      }
      logG_head[e] = logG_tail[e] - kappa * edge_lengths[e];
      signG_head[e] = signG_tail[e];
    }
  } else {
    for (int idx = nE - 1; idx >= 0; --idx) {
      int e = topo_order[idx] - 1;  // 0-indexed edge
      int v = E(e, 1) - 1;          // head vertex of e, 0-indexed
      if (V_outdegree[v] == 0) {
        logG_head[e] = 0.0;
        signG_head[e] = 1.0;
      } else {
        int g = out_edges_by_vertex[v][0];  // first out-edge at v, 0-indexed
        double beta_ge = -in_by_edge[e] / out_by_edge[g];
        logG_head[e] = logG_tail[g] + std::log(std::fabs(beta_ge));
        double sign_beta_ge = (beta_ge > 0.0) - (beta_ge < 0.0);
        signG_head[e] = signG_tail[g] * sign_beta_ge;
      }
      logG_tail[e] = logG_head[e] - kappa * edge_lengths[e];
      signG_tail[e] = signG_head[e];
    }
  }

  return List::create(
      Named("kappa") = kappa,
      Named("tau") = tau,
      Named("sigma_stationary") = sigma_stationary,
      Named("var_tail") = var_tail,
      Named("var_head") = var_head,
      Named("logG_head") = logG_head,
      Named("logG_tail") = logG_tail,
      Named("signG_head") = signG_head,
      Named("signG_tail") = signG_tail);
}

// Build a graph-only O(1)-query LCA index for an out-tree (or forest). The
// Euler tour and sparse-table RMQ are cached by
// directional_ou_setup_structure(), alongside enter/exit/depth/parent_edge,
// so repeated covariance evaluations at new parameter values do not rebuild
// them. Edge identifiers stored in the returned R object are 1-indexed.
// [[Rcpp::export]]
List directional_ou_out_tree_lca_index_cpp(IntegerVector parent_edge,
                                            IntegerVector depth) {
  int nE = parent_edge.size();
  if (depth.size() != nE) {
    stop("parent_edge and depth must have the same length");
  }

  std::vector<std::vector<int>> children(nE);
  std::vector<int> roots;
  for (int e = 0; e < nE; ++e) {
    if (IntegerVector::is_na(parent_edge[e])) {
      roots.push_back(e);
    } else {
      int parent = parent_edge[e] - 1;
      if (parent < 0 || parent >= nE) {
        stop("parent_edge contains an invalid edge number");
      }
      children[parent].push_back(e);
    }
  }

  IntegerVector first(nE, NA_INTEGER);
  IntegerVector root_edge(nE, NA_INTEGER);
  std::vector<int> euler;
  std::vector<int> euler_depth;
  if (nE > 0) {
    euler.reserve(2 * nE - roots.size());
    euler_depth.reserve(2 * nE - roots.size());
  }

  for (std::size_t root_idx = 0; root_idx < roots.size(); ++root_idx) {
    int root = roots[root_idx];
    first[root] = euler.size() + 1;
    root_edge[root] = root + 1;
    euler.push_back(root);
    euler_depth.push_back(depth[root]);

    std::vector<int> edge_stack(1, root);
    std::vector<std::size_t> next_child(1, 0);
    while (!edge_stack.empty()) {
      int edge = edge_stack.back();
      std::size_t &child_idx = next_child.back();
      if (child_idx < children[edge].size()) {
        int child = children[edge][child_idx++];
        first[child] = euler.size() + 1;
        root_edge[child] = root + 1;
        euler.push_back(child);
        euler_depth.push_back(depth[child]);
        edge_stack.push_back(child);
        next_child.push_back(0);
      } else {
        edge_stack.pop_back();
        next_child.pop_back();
        if (!edge_stack.empty()) {
          int parent = edge_stack.back();
          euler.push_back(parent);
          euler_depth.push_back(depth[parent]);
        }
      }
    }
  }

  for (int e = 0; e < nE; ++e) {
    if (IntegerVector::is_na(first[e])) {
      stop("parent_edge does not describe an acyclic rooted forest");
    }
  }

  int tour_size = euler.size();
  IntegerVector log2_floor(tour_size + 1);
  for (int i = 2; i <= tour_size; ++i) {
    log2_floor[i] = log2_floor[i / 2] + 1;
  }
  int n_levels = tour_size == 0 ? 0 : log2_floor[tour_size] + 1;
  IntegerMatrix rmq(n_levels, tour_size);
  for (int i = 0; i < tour_size; ++i) {
    rmq(0, i) = euler[i] + 1;
  }
  for (int level = 1; level < n_levels; ++level) {
    int half_span = 1 << (level - 1);
    int span = half_span << 1;
    for (int i = 0; i + span <= tour_size; ++i) {
      int left = rmq(level - 1, i) - 1;
      int right = rmq(level - 1, i + half_span) - 1;
      rmq(level, i) = depth[left] <= depth[right] ? left + 1 : right + 1;
    }
  }

  return List::create(
      Named("first") = first,
      Named("root_edge") = root_edge,
      Named("log2_floor") = log2_floor,
      Named("rmq") = rmq);
}

// C++ O(n^2) pairwise covariance fill shared by both oriented-tree cases:
//   * converging in-trees (the K1/K2 Columbia paths), where two incomparable
//     edges have no common ancestor and hence zero covariance;
//   * diverging out-trees (the reversed-continuity path), where incomparable
//     edges use the cached Euler/RMQ index above to find their LCA in O(1).
//
// This is the compiled counterpart of directional_ou_pair_covariance() /
// directional_lca() / directional_point_variance() /
// directional_log_transfer(). Irregularly oriented trees remain on the R
// fallback in directional_ou_covariance_from_setup_cpp().
//
// E's edge-number columns (E, PtE_abs/PtE2_abs col 1) and enter/exit are
// 1-indexed (R convention); converted to 0-indexed C++ array access
// consistently below.
// [[Rcpp::export]]
NumericMatrix directional_ou_covariance_oriented_tree_cpp(
    NumericMatrix E, NumericVector edge_lengths,
    std::string tree_orientation, double kappa, double sigma_stationary,
    NumericVector var_tail, NumericVector logG_tail, NumericVector signG_tail,
    IntegerVector enter, IntegerVector exit, IntegerVector depth,
    List out_tree_lca_index,
    NumericMatrix PtE_abs, NumericMatrix PtE2_abs, bool same_point_set) {

  int n1 = PtE_abs.nrow();
  int n2 = PtE2_abs.nrow();
  NumericMatrix Sigma(n1, n2);
  bool is_out_tree = tree_orientation == "out";
  if (!is_out_tree && tree_orientation != "in") {
    stop("tree_orientation must be 'in' or 'out'");
  }

  IntegerVector first;
  IntegerVector root_edge;
  IntegerVector log2_floor;
  IntegerMatrix rmq;
  if (is_out_tree) {
    first = as<IntegerVector>(out_tree_lca_index["first"]);
    root_edge = as<IntegerVector>(out_tree_lca_index["root_edge"]);
    log2_floor = as<IntegerVector>(out_tree_lca_index["log2_floor"]);
    rmq = as<IntegerMatrix>(out_tree_lca_index["rmq"]);
  }

  // point_variance(e, t): e is a 0-indexed edge number here (caller passes
  // the already-decremented edge index).
  auto point_variance = [&](int e0, double t) -> double {
    double decay = std::exp(-2.0 * kappa * t);
    return decay * var_tail[e0] + (1.0 - decay) * sigma_stationary;
  };

  // is_downstream(from_e0, to_e0): is to_e0 reachable forward from from_e0
  // (i.e. from_e0 is upstream of or equal to to_e0)? 0-indexed edges.
  auto is_downstream = [&](int from_e0, int to_e0) -> bool {
    if (is_out_tree) {
      return enter[from_e0] <= enter[to_e0] &&
             exit[to_e0] <= exit[from_e0];
    }
    return enter[to_e0] <= enter[from_e0] && exit[from_e0] <= exit[to_e0];
  };

  auto logG_at = [&](int e0, double t) -> double {
    return is_out_tree ? logG_tail[e0] - kappa * t
                       : logG_tail[e0] + kappa * t;
  };

  auto out_tree_lca = [&](int edge_a0, int edge_b0) -> int {
    int left = first[edge_a0] - 1;
    int right = first[edge_b0] - 1;
    if (left > right) {
      std::swap(left, right);
    }
    int interval_length = right - left + 1;
    int level = log2_floor[interval_length];
    int block_length = 1 << level;
    int candidate_a0 = rmq(level, left) - 1;
    int candidate_b0 = rmq(level, right - block_length + 1) - 1;
    return depth[candidate_a0] <= depth[candidate_b0]
        ? candidate_a0 : candidate_b0;
  };

  for (int i = 0; i < n1; ++i) {
    int j_start = same_point_set ? i : 0;
    int e_s0 = static_cast<int>(PtE_abs(i, 0)) - 1;
    double t_s = PtE_abs(i, 1);

    for (int j = j_start; j < n2; ++j) {
      int e_t0 = static_cast<int>(PtE2_abs(j, 0)) - 1;
      double t_t = PtE2_abs(j, 1);

      // LCA search specialised to the two oriented-tree cases.
      int lca_e0;
      double lca_t;
      if (e_s0 == e_t0) {
        // Case 1: same edge -- LCA is the upstream (smaller-distance) point.
        lca_e0 = e_s0;
        lca_t = std::min(t_s, t_t);
      } else if (is_downstream(e_s0, e_t0)) {
        // Case 2: point_s's edge is upstream of point_t's edge.
        lca_e0 = e_s0;
        lca_t = t_s;
      } else if (is_downstream(e_t0, e_s0)) {
        lca_e0 = e_t0;
        lca_t = t_t;
      } else if (is_out_tree) {
        int root_s0 = root_edge[e_s0] - 1;
        int root_t0 = root_edge[e_t0] - 1;
        if (root_s0 != root_t0) {
          if (E(root_s0, 0) == E(root_t0, 0)) {
            stop("directional_lca(): LCA at a branching source vertex (indegree 0, outdegree > 1) is not yet supported -- vertex %d",
                 static_cast<int>(E(root_s0, 0)));
          }
          Sigma(i, j) = 0.0;
          continue;
        }
        lca_e0 = out_tree_lca(e_s0, e_t0);
        lca_t = edge_lengths[lca_e0];
      } else {
        // Case 4: no common ancestor -- covariance 0.
        Sigma(i, j) = 0.0;
        continue;
      }

      double r_aa = point_variance(lca_e0, lca_t);

      double logG_lca = logG_at(lca_e0, lca_t);
      // transfer(lca, point_s): out-trees are source-normalised and therefore
      // use the opposite subtraction order from outlet-normalised in-trees.
      double log_abs_s = is_out_tree
          ? logG_at(e_s0, t_s) - logG_lca
          : logG_lca - logG_at(e_s0, t_s);
      double sign_s = signG_tail[lca_e0] * signG_tail[e_s0];

      // transfer(lca, point_t)
      double log_abs_t = is_out_tree
          ? logG_at(e_t0, t_t) - logG_lca
          : logG_lca - logG_at(e_t0, t_t);
      double sign_t = signG_tail[lca_e0] * signG_tail[e_t0];

      Sigma(i, j) = r_aa * sign_s * sign_t * std::exp(log_abs_s + log_abs_t);
    }
  }

  if (same_point_set) {
    for (int i = 0; i < n1; ++i) {
      for (int j = 0; j < i; ++j) {
        Sigma(i, j) = Sigma(j, i);
      }
    }
  }

  return Sigma;
}
