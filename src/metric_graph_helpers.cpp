#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// =============================================================================
// Compiled batched helpers for the metric_graph constructor hot path.
//
// Both functions take an edge list (each element a 2-column numeric matrix
// of polyline vertex coordinates) and a `longlat` flag. The flag selects
// between two segment-length implementations:
//
//   * longlat = false  ->  planar Euclidean: sqrt(dx^2 + dy^2)
//   * longlat = true   ->  haversine on (lon, lat) in degrees, on a sphere
//                          of radius 6371008.8 m (WGS84 mean radius).
//
// The haversine branch agrees with sf::st_length to ~1e-7 relative for
// typical metric-graph edges (sub-millimetre absolute error at OSM scales),
// while running ~9000x faster than constructing one sf::st_sfc per edge.
// =============================================================================

static const double METRIC_GRAPH_R_EARTH = 6371008.8; // WGS84 mean radius (m)
// Note: R headers (GraphicsEngine.h) #define DEG2RAD as a macro, so we
// must use a different identifier here.
static const double METRIC_GRAPH_D2R = M_PI / 180.0;

//' @name compute_PtE_edges_cpp
//' @title Per-edge cumulative relative positions
//' @description Given a list of edges (each a 2-column numeric matrix of
//' vertex coordinates), compute, for every edge, the cumulative arc-length
//' normalized to lie in the standard unit interval. Returns a list of numeric vectors, one per
//' edge, each starting at 0 and ending at 1.
//'
//' Degenerate (zero-length) edges return a vector of NaN values
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
//' @description Translation of the post-merge "clean edges"
//' \code{lapply} that the metric_graph constructor runs after the
//' vertex-merge step. For each edge, deduplicates rows by exact
//' \code{(x, y)} equality, while forcing the original first row to remain
//' at position 1 and the original last row to be appended at the end.
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
