#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static const double METRIC_GRAPH_R_EARTH = 6371008.8; // WGS84 mean radius (m)
static const double METRIC_GRAPH_D2R = M_PI / 180.0;

//' @name compute_PtE_edges_cpp
//' @title Per-edge cumulative relative positions
//' @description Given a list of edges (each a 2-column numeric matrix of
//' vertex coordinates), compute, for every edge, the cumulative arc-length
//' normalized to lie in \[0, 1\]. Returns a list of numeric vectors, one per
//' edge, each starting at 0 and ending at 1.
//'
//' Degenerate (zero-length) edges return a vector of NaN values.
//'
//' @param edges List of two-column numeric matrices.
//' @param longlat Logical. If TRUE, use haversine on (lon, lat) in degrees;
//' otherwise use planar Euclidean distances.
//' @return A list of numeric vectors, the same length as \code{edges}.
//' @noRd
// [[Rcpp::export]]
List compute_PtE_edges_cpp(List edges, bool longlat) {
  int nE = edges.size();
  List out(nE);

  for (int i = 0; i < nE; i++) {
    NumericMatrix e = edges[i];
    int n = e.nrow();
    NumericVector pte(n);

    if (n < 2) {
      if (n == 1) pte[0] = 0.0;
      out[i] = pte;
      continue;
    }

    pte[0] = 0.0;
    double total = 0.0;

    if (longlat) {
      double lon1 = e(0, 0);
      double lat1 = e(0, 1);
      double lat1r = lat1 * METRIC_GRAPH_D2R;
      for (int k = 1; k < n; k++) {
        double lon2 = e(k, 0);
        double lat2 = e(k, 1);
        double lat2r = lat2 * METRIC_GRAPH_D2R;
        double dlat = (lat2 - lat1) * METRIC_GRAPH_D2R;
        double dlon = (lon2 - lon1) * METRIC_GRAPH_D2R;
        double sd_lat = std::sin(dlat * 0.5);
        double sd_lon = std::sin(dlon * 0.5);
        double a = sd_lat * sd_lat
                 + std::cos(lat1r) * std::cos(lat2r) * sd_lon * sd_lon;
        if (a > 1.0) a = 1.0;
        total += 2.0 * METRIC_GRAPH_R_EARTH * std::asin(std::sqrt(a));
        pte[k] = total;
        lon1 = lon2;
        lat1 = lat2;
        lat1r = lat2r;
      }
    } else {
      for (int k = 1; k < n; k++) {
        double dx = e(k, 0) - e(k - 1, 0);
        double dy = e(k, 1) - e(k - 1, 1);
        total += std::sqrt(dx * dx + dy * dy);
        pte[k] = total;
      }
    }

    if (total > 0.0) {
      for (int k = 0; k < n; k++) pte[k] /= total;
    } else {
      for (int k = 0; k < n; k++) pte[k] = R_NaN;
    }

    out[i] = pte;
  }

  return out;
}

//' @name compute_edge_lengths_cpp
//' @title Per-edge total polyline lengths
//' @description Given a list of edges (each a 2-column numeric matrix of
//' vertex coordinates), compute the total polyline length of each edge.
//' Returns a numeric vector of length \code{length(edges)}.
//'
//' @param edges List of two-column numeric matrices.
//' @param longlat Logical. If TRUE, treat columns as (lon, lat) in degrees
//' and return lengths in metres (haversine on a sphere of radius
//' 6371008.8 m). If FALSE, return planar Euclidean lengths in whatever
//' coordinate units the edges happen to be in.
//' @return A numeric vector of edge lengths.
//' @noRd
// [[Rcpp::export]]
NumericVector compute_edge_lengths_cpp(List edges, bool longlat) {
  int nE = edges.size();
  NumericVector out(nE);

  for (int i = 0; i < nE; i++) {
    NumericMatrix e = edges[i];
    int n = e.nrow();
    if (n < 2) { out[i] = 0.0; continue; }

    double total = 0.0;

    if (longlat) {
      double lon1 = e(0, 0);
      double lat1 = e(0, 1);
      double lat1r = lat1 * METRIC_GRAPH_D2R;
      for (int k = 1; k < n; k++) {
        double lon2 = e(k, 0);
        double lat2 = e(k, 1);
        double lat2r = lat2 * METRIC_GRAPH_D2R;
        double dlat = (lat2 - lat1) * METRIC_GRAPH_D2R;
        double dlon = (lon2 - lon1) * METRIC_GRAPH_D2R;
        double sd_lat = std::sin(dlat * 0.5);
        double sd_lon = std::sin(dlon * 0.5);
        double a = sd_lat * sd_lat
                 + std::cos(lat1r) * std::cos(lat2r) * sd_lon * sd_lon;
        if (a > 1.0) a = 1.0;
        total += 2.0 * METRIC_GRAPH_R_EARTH * std::asin(std::sqrt(a));
        lon1 = lon2;
        lat1 = lat2;
        lat1r = lat2r;
      }
    } else {
      for (int k = 1; k < n; k++) {
        double dx = e(k, 0) - e(k - 1, 0);
        double dy = e(k, 1) - e(k - 1, 1);
        total += std::sqrt(dx * dx + dy * dy);
      }
    }

    out[i] = total;
  }

  return out;
}

//' @name postprocess_edges_cpp
//' @title Dedupe interior polyline points while preserving endpoints
//' @description Bit-faithful translation of the post-merge "clean edges"
//' \code{lapply} that the metric_graph constructor runs after the
//' vertex-merge step. For each edge, deduplicates rows by exact
//' \code{(x, y)} equality, while forcing the original first row to remain
//' at position 1 and the original last row to be appended at the end.
//'
//' @param edges List of two-column numeric matrices.
//' @return A list of two-column numeric matrices with deduped interior
//' rows. May have fewer rows than the input.
//' @noRd
// [[Rcpp::export]]
List postprocess_edges_cpp(List edges) {
  int nE = edges.size();
  List out(nE);

  // scratch buffers reused per edge
  std::vector<int> kept;  kept.reserve(64);
  std::vector<int> kept2; kept2.reserve(64);

  for (int i = 0; i < nE; i++) {
    NumericMatrix e = edges[i];
    int n = e.nrow();

    if (n <= 2) {
      // For n <= 2 the original code skips the if(nrow > 2) branch
      // entirely; just copy through (rebuild to drop any rownames attr).
      NumericMatrix copy(n, 2);
      for (int k = 0; k < n; k++) {
        copy(k, 0) = e(k, 0);
        copy(k, 1) = e(k, 1);
      }
      out[i] = copy;
      continue;
    }

    // Step 1+2: unique() on rows 0..n-2 (R indices 1..n-1)
    kept.clear();
    for (int k = 0; k < n - 1; k++) {
      double xk = e(k, 0), yk = e(k, 1);
      bool seen = false;
      for (size_t m = 0; m < kept.size(); m++) {
        if (e(kept[m], 0) == xk && e(kept[m], 1) == yk) {
          seen = true;
          break;
        }
      }
      if (!seen) kept.push_back(k);
    }

    // Step 3: append the last row (R index n, C++ index n-1)
    kept.push_back(n - 1);
    int nk = (int) kept.size();

    if (nk > 2) {
      // Step 4: drop kept[0], unique() the remainder, prepend the
      // ORIGINAL row 0 (not whatever kept[0] was -- if the original
      // row 0 appeared anywhere in the interior, kept[0] is still it).
      kept2.clear();
      for (int m = 1; m < nk; m++) {
        int row = kept[m];
        double xr = e(row, 0), yr = e(row, 1);
        bool seen = false;
        for (size_t p = 0; p < kept2.size(); p++) {
          if (e(kept2[p], 0) == xr && e(kept2[p], 1) == yr) {
            seen = true;
            break;
          }
        }
        if (!seen) kept2.push_back(row);
      }
      int n_out = 1 + (int) kept2.size();
      NumericMatrix out_e(n_out, 2);
      out_e(0, 0) = e(0, 0);
      out_e(0, 1) = e(0, 1);
      for (int k = 0; k < (int) kept2.size(); k++) {
        out_e(k + 1, 0) = e(kept2[k], 0);
        out_e(k + 1, 1) = e(kept2[k], 1);
      }
      out[i] = out_e;
    } else {
      NumericMatrix out_e(nk, 2);
      for (int k = 0; k < nk; k++) {
        out_e(k, 0) = e(kept[k], 0);
        out_e(k, 1) = e(kept[k], 1);
      }
      out[i] = out_e;
    }
  }

  return out;
}

//' @name aeqd_project_cpp
//' @title Closed-form spherical Azimuthal Equidistant projection
//' @description Forward projects (lon, lat) coordinates to local AEQD-style
//' planar coordinates centered on a given (lon0, lat0). Output is in metres
//' on a sphere of radius 6371008.8 m (WGS84 mean Earth radius).
//'
//' Forward formula:
//' \preformatted{
//'   c = acos(sin(lat0)*sin(lat) + cos(lat0)*cos(lat)*cos(lon - lon0))
//'   k = c / sin(c)             (-> 1 as c -> 0)
//'   x = R * k * cos(lat) * sin(lon - lon0)
//'   y = R * k * (cos(lat0)*sin(lat) - sin(lat0)*cos(lat)*cos(lon - lon0))
//' }
//'
//' @param pts Two-column numeric matrix of (lon, lat) coordinates in degrees.
//' @param lon0 Longitude of the projection center, in degrees.
//' @param lat0 Latitude of the projection center, in degrees.
//' @return A two-column numeric matrix of planar (x, y) coordinates in metres.
//' @noRd
// [[Rcpp::export]]
NumericMatrix aeqd_project_cpp(NumericMatrix pts, double lon0, double lat0) {
  int n = pts.nrow();
  NumericMatrix out(n, 2);
  const double R = METRIC_GRAPH_R_EARTH;
  const double D2R = METRIC_GRAPH_D2R;

  double lat0r = lat0 * D2R;
  double lon0r = lon0 * D2R;
  double cos_lat0 = std::cos(lat0r);
  double sin_lat0 = std::sin(lat0r);

  for (int i = 0; i < n; i++) {
    double lonr = pts(i, 0) * D2R;
    double latr = pts(i, 1) * D2R;
    double cos_lat = std::cos(latr);
    double sin_lat = std::sin(latr);
    double dlon = lonr - lon0r;
    double cos_dlon = std::cos(dlon);
    double sin_dlon = std::sin(dlon);

    double cos_c = sin_lat0 * sin_lat + cos_lat0 * cos_lat * cos_dlon;
    if (cos_c > 1.0) cos_c = 1.0;
    if (cos_c < -1.0) cos_c = -1.0;
    double c_ang = std::acos(cos_c);
    double k = (c_ang < 1e-12) ? 1.0 : (c_ang / std::sin(c_ang));

    out(i, 0) = R * k * cos_lat * sin_dlon;
    out(i, 1) = R * k * (cos_lat0 * sin_lat - sin_lat0 * cos_lat * cos_dlon);
  }

  return out;
}

//' @name split_one_edge_cpp
//' @title Inner edge-splitting work for split_edge_batch (one edge)
//' @description Processes a single edge: builds coords_list1 (the first
//' segment of the original edge ending at the first split point), coords_list2
//' (a List of subsequent segments between split points and from the last
//' split to the end of the edge), segment_lengths (a NumericVector of length
//' n_split + 1), and the aux_matrix (an IntegerMatrix of (n_split + 1) x 2
//' giving the new E rows). Each segment carries a PtE attribute renormalized
//' to \[0, 1\].
//' @param edge nx2 matrix: the polyline
//' @param PtE_edge n vector: PtE values at each polyline point
//' @param t_values n_split vector: sorted positions in \[0, 1\] where to split
//' @param edge_len double: original edge length
//' @param E_row IntegerVector(2): current E\[Ei,\] = (start_v, end_v)
//' @param first_new_v int: ID of the first new vertex for this edge
//' @noRd
// [[Rcpp::export]]
List split_one_edge_cpp(
    NumericMatrix edge,
    NumericVector PtE_edge,
    NumericVector t_values,
    double edge_len,
    IntegerVector E_row,
    int first_new_v
) {
  int n_edge  = edge.nrow();
  int n_split = t_values.size();

  // Cumulative arc length normalized to [0, 1]
  std::vector<double> dist_vec(n_edge, 0.0);
  for (int i = 0; i < n_edge - 1; i++) {
    double dx = edge(i + 1, 0) - edge(i, 0);
    double dy = edge(i + 1, 1) - edge(i, 1);
    dist_vec[i + 1] = dist_vec[i] + std::sqrt(dx * dx + dy * dy);
  }
  double total_d = dist_vec[n_edge - 1];
  std::vector<double> dist_norm(n_edge);
  if (total_d > 0.0) {
    for (int i = 0; i < n_edge; i++) dist_norm[i] = dist_vec[i] / total_d;
  } else {
    for (int i = 0; i < n_edge; i++) dist_norm[i] = 0.0;
    dist_norm[n_edge - 1] = 1.0;
  }

  // Interpolation: val_lines (n_split x 2) and idx_positions (n_split, 1-indexed)
  NumericMatrix val_lines(n_split, 2);
  IntegerVector idx_positions(n_split);
  for (int i = 0; i < n_split; i++) {
    double pr = t_values[i];
    if (pr < 0) pr = 0;
    if (pr > 1) pr = 1;
    int tmp_ind = -1;
    for (int j = 0; j < n_edge - 1; j++) {
      if (pr >= dist_norm[j] && pr <= dist_norm[j + 1]) tmp_ind = j;
    }
    if (tmp_ind < 0) tmp_ind = n_edge - 2;
    double denom = dist_norm[tmp_ind + 1] - dist_norm[tmp_ind];
    double dist_pos = (denom == 0.0) ? 0.0 : (pr - dist_norm[tmp_ind]) / denom;
    val_lines(i, 0) = edge(tmp_ind, 0) + (edge(tmp_ind + 1, 0) - edge(tmp_ind, 0)) * dist_pos;
    val_lines(i, 1) = edge(tmp_ind, 1) + (edge(tmp_ind + 1, 1) - edge(tmp_ind, 1)) * dist_pos;
    idx_positions[i] = tmp_ind + 1;
  }

  // coords_list1: edge[1..idx_positions[0]] + val_lines[0], with dedup
  int idx0 = idx_positions[0];
  int cap1 = idx0 + 1;
  NumericMatrix tmp1(cap1, 2);
  NumericVector pte1(cap1);
  int n1 = 0;
  for (int k = 0; k < idx0; k++) {
    double x = edge(k, 0);
    double y = edge(k, 1);
    if (n1 == 0 || x != tmp1(n1 - 1, 0) || y != tmp1(n1 - 1, 1)) {
      tmp1(n1, 0) = x; tmp1(n1, 1) = y; pte1[n1] = PtE_edge[k]; n1++;
    }
  }
  {
    double x = val_lines(0, 0);
    double y = val_lines(0, 1);
    if (n1 == 0 || x != tmp1(n1 - 1, 0) || y != tmp1(n1 - 1, 1)) {
      tmp1(n1, 0) = x; tmp1(n1, 1) = y; pte1[n1] = t_values[0]; n1++;
    }
  }
  NumericMatrix coords_list1(n1, 2);
  NumericVector pte1_out(n1);
  for (int k = 0; k < n1; k++) {
    coords_list1(k, 0) = tmp1(k, 0);
    coords_list1(k, 1) = tmp1(k, 1);
    pte1_out[k] = pte1[k];
  }
  double base1 = pte1_out[0];
  double norm1 = (n1 > 0) ? (pte1_out[n1 - 1] - base1) : 0.0;
  if (norm1 == 0.0) {
    for (int k = 0; k < n1; k++) pte1_out[k] = 0.0;
  } else {
    for (int k = 0; k < n1; k++) pte1_out[k] = (pte1_out[k] - base1) / norm1;
  }
  coords_list1.attr("PtE") = pte1_out;

  // coords_list2[[j]] for j = 0..n_split-1
  List coords_list2(n_split);
  int max_seg_size = n_edge + 2;
  NumericMatrix tmp2(max_seg_size, 2);
  NumericVector pte2(max_seg_size);
  for (int j = 0; j < n_split; j++) {
    int n2 = 0;
    tmp2(n2, 0) = val_lines(j, 0);
    tmp2(n2, 1) = val_lines(j, 1);
    pte2[n2]    = t_values[j];
    n2++;
    if (j < n_split - 1) {
      int ij  = idx_positions[j];
      int ij1 = idx_positions[j + 1];
      if (ij != ij1) {
        for (int k = ij; k < ij1; k++) {
          double x = edge(k, 0), y = edge(k, 1);
          if (x != tmp2(n2 - 1, 0) || y != tmp2(n2 - 1, 1)) {
            tmp2(n2, 0) = x; tmp2(n2, 1) = y; pte2[n2] = PtE_edge[k]; n2++;
          }
        }
      }
      double xe = val_lines(j + 1, 0), ye = val_lines(j + 1, 1);
      if (xe != tmp2(n2 - 1, 0) || ye != tmp2(n2 - 1, 1)) {
        tmp2(n2, 0) = xe; tmp2(n2, 1) = ye; pte2[n2] = t_values[j + 1]; n2++;
      }
    } else {
      int ij = idx_positions[j];
      for (int k = ij; k < n_edge; k++) {
        double x = edge(k, 0), y = edge(k, 1);
        if (x != tmp2(n2 - 1, 0) || y != tmp2(n2 - 1, 1)) {
          tmp2(n2, 0) = x; tmp2(n2, 1) = y; pte2[n2] = PtE_edge[k]; n2++;
        }
      }
    }
    NumericMatrix cl2(n2, 2);
    NumericVector p2_out(n2);
    for (int k = 0; k < n2; k++) {
      cl2(k, 0) = tmp2(k, 0);
      cl2(k, 1) = tmp2(k, 1);
      p2_out[k] = pte2[k];
    }
    double base2 = p2_out[0];
    double norm2 = (n2 > 0) ? (p2_out[n2 - 1] - base2) : 0.0;
    if (norm2 == 0.0) {
      for (int k = 0; k < n2; k++) p2_out[k] = 0.0;
    } else {
      for (int k = 0; k < n2; k++) p2_out[k] = (p2_out[k] - base2) / norm2;
    }
    cl2.attr("PtE") = p2_out;
    coords_list2[j] = cl2;
  }

  // aux_matrix: ((n_split + 1) x 2) rows for E
  IntegerMatrix aux_matrix(n_split + 1, 2);
  aux_matrix(0, 0) = E_row[0];
  aux_matrix(0, 1) = first_new_v;
  for (int j = 0; j < n_split - 1; j++) {
    aux_matrix(j + 1, 0) = first_new_v + j;
    aux_matrix(j + 1, 1) = first_new_v + j + 1;
  }
  aux_matrix(n_split, 0) = first_new_v + n_split - 1;
  aux_matrix(n_split, 1) = E_row[1];

  // segment_lengths
  NumericVector segment_lengths(n_split + 1);
  segment_lengths[0] = t_values[0] * edge_len;
  for (int j = 1; j < n_split; j++) {
    segment_lengths[j] = (t_values[j] - t_values[j - 1]) * edge_len;
  }
  segment_lengths[n_split] = (1.0 - t_values[n_split - 1]) * edge_len;

  return List::create(
    Named("coords_list1")    = coords_list1,
    Named("coords_list2")    = coords_list2,
    Named("segment_lengths") = segment_lengths,
    Named("aux_matrix")      = aux_matrix,
    Named("val_lines")       = val_lines
  );
}

//' @name split_edges_batch_cpp
//' @title Batched inner edge-splitting work for split_edge_batch
//' @description Processes all edge groups in one call. Returns a List of
//' per-group results, each matching split_one_edge_cpp's output structure.
//' @param edges_full List: full self$edges, looked up by Ei (1-indexed)
//' @param PtE_full List: precomputed PtE attributes for each edge in edges_full
//' (caller extracts these once to avoid attr() lookups inside C++)
//' @param E_full IntegerMatrix: full self$E
//' @param edge_lens NumericVector: full self$edge_lengths
//' @param edge_ids IntegerVector: which edge each group refers to (1-indexed)
//' @param t_values_list List: per-group t_values
//' @param first_new_vs IntegerVector: ID of the first new vertex for each group
//' @noRd
// [[Rcpp::export]]
List split_edges_batch_cpp(
    List edges_full,
    List PtE_full,
    IntegerMatrix E_full,
    NumericVector edge_lens,
    IntegerVector edge_ids,
    List t_values_list,
    IntegerVector first_new_vs
) {
  int n_groups = edge_ids.size();
  List out(n_groups);
  IntegerVector E_row(2);
  for (int g = 0; g < n_groups; g++) {
    int Ei = edge_ids[g] - 1;
    NumericMatrix edge      = as<NumericMatrix>(edges_full[Ei]);
    NumericVector PtE_edge  = as<NumericVector>(PtE_full[Ei]);
    NumericVector t_values  = as<NumericVector>(t_values_list[g]);
    double edge_len         = edge_lens[Ei];
    E_row[0] = E_full(Ei, 0);
    E_row[1] = E_full(Ei, 1);
    int first_new_v = first_new_vs[g];
    out[g] = split_one_edge_cpp(edge, PtE_edge, t_values, edge_len,
                                E_row, first_new_v);
  }
  return out;
}
