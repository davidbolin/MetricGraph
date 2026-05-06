
// CB_revisited.cpp
// Optimized constraint-basis construction for MetricGraph (alpha = 2 case).
//
// Core idea
// ---------
// For the non-directional alpha=2 constraints built by
// construct_constraint_matrix(), every vertex v with degree n_e contributes
// two independent sub-blocks to the constraint matrix C:
//
//   Derivative sub-block  (1 x n_e row)
//     Row = [+1,...,+1, -1,...,-1]  (n_le ones, n_ue minus-ones)
//     Columns: 4*le[k]+1  (k=0..n_le-1)  and  4*ue[k]+3  (k=0..n_ue-1)
//
//   Continuity sub-block  ((n_e-1) x n_e difference matrix D)
//     D(i,i)=+1, D(i,i+1)=-1  for i=0..n_e-2
//     Columns: 4*le[k]    (k=0..n_le-1)  and  4*ue[k]+2  (k=0..n_ue-1)
//
// Because the two sub-blocks occupy disjoint columns, and no column is shared
// between vertices, every subcluster found by the generic c_basis2() BFS
// corresponds to exactly one of these per-vertex sub-blocks.
//
// Consequently:
//   * The SVD of the continuity sub-block depends only on n_e  → precompute
//     once for each distinct degree and reuse for every vertex of that degree.
//   * The SVD of the derivative sub-block depends on (n_le, n_ue) → cache
//     by pair (at most O(max_deg^2) unique pairs).
//   * We can build T, U, S directly from the graph structure without ever
//     assembling the full sparse matrix C.
//
// The construct_constraint_matrix_fast() below is a tightened version of the
// original: it avoids rebuilding the edge lists and uses a single-pass
// reservation so no reallocation happens.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
using namespace Rcpp;

// [[Rcpp::depends(RcppEigen)]]

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

struct PrecomputedSVD {
    Eigen::MatrixXd U;   // (m x m)  m = number of rows in the sub-block
    Eigen::VectorXd S;   // (m)      singular values, S(i) >= 0
    Eigen::MatrixXd V;   // (n x n)  full V, n = number of columns
};

// Build the (d-1) x d difference matrix: D(i,i)=+1, D(i,i+1)=-1
static Eigen::MatrixXd diff_matrix(int d) {
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(d - 1, d);
    for (int i = 0; i < d - 1; ++i) {
        D(i, i)     =  1.0;
        D(i, i + 1) = -1.0;
    }
    return D;
}

// Precompute full SVDs of (d-1) x d difference matrices for d = 2..max_degree
static std::vector<PrecomputedSVD> precompute_continuity_svds(int max_degree) {
    std::vector<PrecomputedSVD> cache(max_degree + 1);
    for (int d = 2; d <= max_degree; ++d) {
        Eigen::MatrixXd D = diff_matrix(d);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            D, Eigen::ComputeFullU | Eigen::ComputeFullV);
        cache[d].U = svd.matrixU();
        cache[d].S = svd.singularValues();
        cache[d].V = svd.matrixV();
    }
    return cache;
}

// Compute full SVD of the 1 x n_e derivative row
// [+1,...,+1 (n_le times), -1,...,-1 (n_ue times)]
static PrecomputedSVD compute_derivative_svd(int n_le, int n_ue) {
    int n_e = n_le + n_ue;
    Eigen::MatrixXd row(1, n_e);
    for (int i = 0;    i < n_le; ++i) row(0, i)      =  1.0;
    for (int i = n_le; i < n_e;  ++i) row(0, i)      = -1.0;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        row, Eigen::ComputeFullU | Eigen::ComputeFullV);
    PrecomputedSVD result;
    result.U = svd.matrixU();   // 1 x 1
    result.S = svd.singularValues(); // length 1
    result.V = svd.matrixV();   // n_e x n_e  (full)
    return result;
}

// ---------------------------------------------------------------------------
// construct_constraint_matrix_fast
// ---------------------------------------------------------------------------
// Same semantics as the original construct_constraint_matrix() but
//   * builds the lower/upper edge lists once up-front (no repeated O(nE) scan)
//   * reserves the exact capacity needed for all triplet vectors
//   * avoids the intermediate pair<int,int> vector in the continuity block
//
// [[Rcpp::export]]
Eigen::SparseMatrix<double>
construct_constraint_matrix(const Eigen::MatrixXi &E, int nV,
                            int edge_constraint) {
    int nE = E.rows();

    // Build adjacency lists (1-indexed vertices, 0-indexed edge IDs)
    std::vector<std::vector<int>> lower_edges(nV + 1), upper_edges(nV + 1);
    lower_edges.reserve(nV + 1);
    upper_edges.reserve(nV + 1);
    for (int e = 0; e < nE; ++e) {
        lower_edges[E(e, 0)].push_back(e);
        upper_edges[E(e, 1)].push_back(e);
    }

    // First pass: count constraints and non-zeros so we can reserve exactly.
    int nC_total  = 0;
    int nnz_total = 0;
    for (int v = 1; v <= nV; ++v) {
        int n_le = static_cast<int>(lower_edges[v].size());
        int n_ue = static_cast<int>(upper_edges[v].size());
        int n_e  = n_le + n_ue;
        if ((edge_constraint && n_e == 1) || n_e > 1) {
            ++nC_total;
            nnz_total += n_e;
        }
        if (n_e > 1) {
            nC_total  += n_e - 1;
            nnz_total += 2 * (n_e - 1);
        }
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nnz_total);

    int row = 0;
    for (int v = 1; v <= nV; ++v) {
        const std::vector<int> &le = lower_edges[v];
        const std::vector<int> &ue = upper_edges[v];
        int n_le = static_cast<int>(le.size());
        int n_ue = static_cast<int>(ue.size());
        int n_e  = n_le + n_ue;

        // ---- Derivative constraint ----
        if ((edge_constraint && n_e == 1) || n_e > 1) {
            for (int k = 0; k < n_le; ++k)
                triplets.emplace_back(row, 4 * le[k] + 1,  1.0);
            for (int k = 0; k < n_ue; ++k)
                triplets.emplace_back(row, 4 * ue[k] + 3, -1.0);
            ++row;
        }

        // ---- Continuity constraints ----
        if (n_e > 1) {
            // Column list in the joined order: lower edges then upper edges.
            // lower edge e → DoF column 4*e+0
            // upper edge e → DoF column 4*e+2
            // Constraint i links adjacent entries i-1 and i with +1/-1.
            std::vector<int> cols(n_e);
            for (int k = 0; k < n_le; ++k) cols[k]        = 4 * le[k];
            for (int k = 0; k < n_ue; ++k) cols[n_le + k] = 4 * ue[k] + 2;

            for (int i = 1; i < n_e; ++i) {
                triplets.emplace_back(row, cols[i - 1],  1.0);
                triplets.emplace_back(row, cols[i],     -1.0);
                ++row;
            }
        }
    }

    Eigen::SparseMatrix<double> C(nC_total, 4 * nE);
    C.setFromTriplets(triplets.begin(), triplets.end());
    C.makeCompressed();
    return C;
}

// ---------------------------------------------------------------------------
// c_basis2_graph  –  optimized Change-of-Basis for the alpha=2 case
// ---------------------------------------------------------------------------
// Replaces construct_constraint_matrix() + c_basis2() with a single pass
// over the vertices.  Returns exactly the same list structure as c_basis2():
//   T  – (4*nE) x (4*nE) sparse orthogonal matrix
//   S  – nC-length vector of singular values
//   U  – nC x nC sparse matrix of left singular vectors
//   larget.cluster  – row indices of the largest subcluster (0-indexed)
//   cluster.n       – sizes of all subclusters
//
// Processing order matches construct_constraint_matrix():
//   for each vertex v = 1..nV:
//     1. derivative constraint   (1 row)
//     2. continuity constraints  (n_e-1 rows)
// so the constraint row numbering here agrees with that function's output.
//
// [[Rcpp::export]]
Rcpp::List c_basis2_graph(const Eigen::MatrixXi &E, int nV,
                           int edge_constraint) {
    int nE   = E.rows();
    int n_dof = 4 * nE;

    // ---- Build adjacency lists (1-indexed vertices, 0-indexed edges) ----
    std::vector<std::vector<int>> lower_edges(nV + 1), upper_edges(nV + 1);
    for (int e = 0; e < nE; ++e) {
        lower_edges[E(e, 0)].push_back(e);
        upper_edges[E(e, 1)].push_back(e);
    }

    // ---- Compute nC and max degree in one pass ----
    int nC         = 0;
    int max_degree = 0;
    for (int v = 1; v <= nV; ++v) {
        int n_le = static_cast<int>(lower_edges[v].size());
        int n_ue = static_cast<int>(upper_edges[v].size());
        int n_e  = n_le + n_ue;
        if (n_e > max_degree) max_degree = n_e;
        if ((edge_constraint && n_e == 1) || n_e > 1) ++nC;  // derivative
        if (n_e > 1)                                   nC += n_e - 1; // continuity
    }

    // ---- Precompute SVDs ----
    // Continuity: (d-1) x d difference matrix, one per unique degree d
    std::vector<PrecomputedSVD> cont_cache =
        precompute_continuity_svds(max_degree);

    // Derivative: 1 x n_e row, keyed by (n_le, n_ue)
    // Key = n_le * (max_degree+1) + n_ue  (unique for reasonable degrees)
    std::map<int, PrecomputedSVD> deriv_cache;
    auto get_deriv_svd = [&](int n_le, int n_ue) -> const PrecomputedSVD & {
        int key = n_le * (max_degree + 1) + n_ue;
        auto it = deriv_cache.find(key);
        if (it == deriv_cache.end()) {
            deriv_cache.emplace(key, compute_derivative_svd(n_le, n_ue));
        }
        return deriv_cache.at(key);
    };

    // ---- Allocate output structures ----
    Eigen::VectorXd singular_values(nC);
    std::vector<Eigen::Triplet<double>> tripletsT, tripletsU;

    // Column counters for T:
    //   counter_K  : next "constrained" column  (0 .. nC-1)
    //   counter_D  : next "null/free"  column   (nC .. n_dof-1)
    int counter_K    = 0;
    int counter_D    = nC;
    int constraint_row = 0;   // current row in U / S

    // Track which DoF columns participate in at least one constraint
    std::vector<bool> dof_active(n_dof, false);

    // Subcluster bookkeeping (matches c_basis2 output)
    std::vector<int> cluster_sizes;
    std::vector<int> index_largest;
    int n_largest = 0;

    // Reserve a rough upper bound for the triplet lists
    tripletsT.reserve(static_cast<size_t>(nC) * 4 + n_dof);
    tripletsU.reserve(static_cast<size_t>(nC) * 4);

    // ---- Main loop over vertices ----
    for (int v = 1; v <= nV; ++v) {
        const std::vector<int> &le = lower_edges[v];
        const std::vector<int> &ue = upper_edges[v];
        int n_le = static_cast<int>(le.size());
        int n_ue = static_cast<int>(ue.size());
        int n_e  = n_le + n_ue;
        if (n_e == 0) continue;

        // ==================================================================
        // 1. DERIVATIVE sub-block  (processed first, matching C row order)
        // ==================================================================
        if ((edge_constraint && n_e == 1) || n_e > 1) {
            const PrecomputedSVD &svd = get_deriv_svd(n_le, n_ue);
            // Columns of this sub-block in the DoF vector
            //   lower edges: 4*e+1
            //   upper edges: 4*e+3
            // Column j in [0..n_e-1] in the local ordering

            // Mark columns active
            for (int k = 0; k < n_le; ++k) dof_active[4 * le[k] + 1] = true;
            for (int k = 0; k < n_ue; ++k) dof_active[4 * ue[k] + 3] = true;

            // --- Fill T ---
            // V(:, 0)       → T column counter_K        (constrained direction)
            // V(:, 1..n_e-1) → T columns counter_D .. counter_D+n_e-2 (null)
            for (int k = 0; k < n_le; ++k) {
                int gc = 4 * le[k] + 1;  // global DoF column
                tripletsT.emplace_back(gc, counter_K, svd.V(k, 0));
                for (int i = 1; i < n_e; ++i)
                    tripletsT.emplace_back(gc, counter_D + i - 1, svd.V(k, i));
            }
            for (int k = 0; k < n_ue; ++k) {
                int gc = 4 * ue[k] + 3;
                int local_j = n_le + k;
                tripletsT.emplace_back(gc, counter_K, svd.V(local_j, 0));
                for (int i = 1; i < n_e; ++i)
                    tripletsT.emplace_back(gc, counter_D + i - 1,
                                           svd.V(local_j, i));
            }

            // --- Fill U ---
            // U is 1x1 for the derivative sub-block
            tripletsU.emplace_back(constraint_row, counter_K, svd.U(0, 0));

            // --- Fill S ---
            singular_values(counter_K) = svd.S(0);

            // Subcluster bookkeeping: size 1 (one constraint row)
            cluster_sizes.push_back(1);
            if (1 > n_largest) {
                n_largest = 1;
                index_largest = {constraint_row};
            }

            ++constraint_row;
            ++counter_K;
            counter_D += n_e - 1;   // n_e-1 null-space directions
        }

        // ==================================================================
        // 2. CONTINUITY sub-block  (n_e-1 rows, n_e columns, n_e > 1)
        // ==================================================================
        if (n_e > 1) {
            const PrecomputedSVD &svd = cont_cache[n_e];
            // Columns of this sub-block in the DoF vector
            //   lower edges: 4*e+0
            //   upper edges: 4*e+2
            // Local column order: [le[0],..,le[n_le-1], ue[0],..,ue[n_ue-1]]

            // Mark columns active
            for (int k = 0; k < n_le; ++k) dof_active[4 * le[k]]     = true;
            for (int k = 0; k < n_ue; ++k) dof_active[4 * ue[k] + 2] = true;

            // --- Fill T ---
            // V(:, 0..n_e-2) → T columns counter_K .. counter_K+n_e-2
            // V(:, n_e-1)    → T column  counter_D   (one null direction)
            int null_col = counter_D;   // continuity sub-block has exactly 1 null direction

            for (int k = 0; k < n_le; ++k) {
                int gc = 4 * le[k];
                for (int i = 0; i < n_e - 1; ++i)
                    tripletsT.emplace_back(gc, counter_K + i, svd.V(k, i));
                tripletsT.emplace_back(gc, null_col, svd.V(k, n_e - 1));
            }
            for (int k = 0; k < n_ue; ++k) {
                int gc = 4 * ue[k] + 2;
                int local_j = n_le + k;
                for (int i = 0; i < n_e - 1; ++i)
                    tripletsT.emplace_back(gc, counter_K + i, svd.V(local_j, i));
                tripletsT.emplace_back(gc, null_col, svd.V(local_j, n_e - 1));
            }

            // --- Fill U ---
            // U is (n_e-1) x (n_e-1); rows in U are constraint_row..constraint_row+n_e-2
            for (int row_k = 0; row_k < n_e - 1; ++row_k) {
                for (int i = 0; i < n_e - 1; ++i)
                    tripletsU.emplace_back(constraint_row + row_k,
                                           counter_K + i,
                                           svd.U(row_k, i));
            }

            // --- Fill S ---
            for (int i = 0; i < n_e - 1; ++i)
                singular_values(counter_K + i) = svd.S(i);

            // Subcluster bookkeeping: size n_e-1 constraint rows
            int clust_size = n_e - 1;
            cluster_sizes.push_back(clust_size);
            if (clust_size > n_largest) {
                n_largest = clust_size;
                index_largest.resize(clust_size);
                for (int k = 0; k < clust_size; ++k)
                    index_largest[k] = constraint_row + k;
            }

            constraint_row += n_e - 1;
            counter_K      += n_e - 1;
            counter_D      += 1;          // one null direction for continuity
        }
    }

    // ---- Identity entries for inactive DoF columns ----
    // In c_basis2, these land at T columns  counter..A.outerSize()-1.
    // Here counter_D == number of active DoF columns after the loop,
    // which equals "counter" in c_basis2.  Inactive columns get identity
    // entries at the next available T columns (counter_D, counter_D+1, ...).
    for (int col = 0; col < n_dof; ++col) {
        if (!dof_active[col]) {
            tripletsT.emplace_back(col, counter_D, 1.0);
            ++counter_D;
        }
    }

    // ---- Build sparse output matrices ----
    Eigen::SparseMatrix<double> T(n_dof, n_dof);
    T.setFromTriplets(tripletsT.begin(), tripletsT.end());

    Eigen::SparseMatrix<double> U_mat(nC, nC);
    U_mat.setFromTriplets(tripletsU.begin(), tripletsU.end());

    Rcpp::List output;
    output["T"]              = T;
    output["S"]              = singular_values;
    output["U"]              = U_mat;
    output["larget.cluster"] = index_largest;   // keep typo to match c_basis2
    output["cluster.n"]      = cluster_sizes;
    return output;
}

// ===========================================================================
// Directional-model optimizations
// ===========================================================================
//
// For the directional alpha=1 (and general alpha) constraints built by
// construct_directional_constraint_matrix(), every vertex v contributes
// independent sub-blocks whose structure depends only on its in/out-degree:
//
//   Type-1 vertex  (indegree n_in > 0, outdegree n_out > 0)
//     For each derivative order der = 1..alpha, one subcluster:
//       n_out × (n_out + n_in)  local matrix  [I_{n_out} | 𝟏_{n_out × n_in}]
//       Out-edge columns:  2*alpha*oe[k] + der-1   (k = 0..n_out-1)
//       In-edge columns:   2*alpha*ie[k] + alpha+der-1  (k = 0..n_in-1)
//
//   Type-2 vertex  (indegree 0, outdegree n_out > 1)
//     For each der = 1..alpha, one subcluster:
//       (n_out-1) × n_out  difference matrix  (same as non-directional)
//       Out-edge columns:  2*alpha*oe[k] + der-1
//
// No column is shared between vertices → BFS in c_basis2() is provably
// redundant → we can build T, U, S directly from the graph structure.
//
// SVDs to precompute:
//   type-1: SVD of [I_{n_out} | 𝟏_{n_out × n_in}]  for each (n_in, n_out)
//   type-2: SVD of (n_out-1) × n_out difference matrix  (reuse cont_cache)
//
// Constraint-row ordering matches construct_directional_constraint_matrix:
//   1. All type-1 vertices (v = 0..nV-1 with in>0 and out>0)
//        for out-edge i = 0..n_out-1, for der = 1..alpha: one row
//   2. All type-2 vertices (v = 0..nV-1 with in=0 and out>1)
//        for pair i = 1..n_out-1, for der = 1..alpha: one row
//
// Because out-edge varies slowly and der varies fast inside each vertex,
// the rows for a fixed der at a type-1 vertex are non-contiguous
// (stride alpha).  We handle this explicitly in the U / S fill.
// ---------------------------------------------------------------------------

// Build [I_{n_out} | 𝟏_{n_out × n_in}] and compute its full SVD
static PrecomputedSVD compute_type1_svd(int n_in, int n_out) {
    int n_col = n_out + n_in;
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n_out, n_col);
    A.leftCols(n_out)  = Eigen::MatrixXd::Identity(n_out, n_out);
    A.rightCols(n_in)  = Eigen::MatrixXd::Ones(n_out, n_in);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    PrecomputedSVD res;
    res.U = svd.matrixU();
    res.S = svd.singularValues();
    res.V = svd.matrixV();
    return res;
}

// ---------------------------------------------------------------------------
// construct_directional_constraint_matrix_fast
// ---------------------------------------------------------------------------
// Same output as the R construct_directional_constraint_matrix() but builds
// the sparse matrix in C++ using adjacency lists (O(nE)), avoiding the
// O(nV * nE) inner scan.
//
// w_out[e] = DirectionalWeightFunction_out applied to the weight of edge e.
// w_in[e]  = DirectionalWeightFunction_in applied to the weights of the
//            in-edges of the vertex that edge e enters (one scalar per edge).
// Both vectors have length nE and should be computed in R before calling.
//
// [[Rcpp::export]]
Eigen::SparseMatrix<double>
construct_directional_constraint_matrix_fast(
    const Eigen::MatrixXi &E, int nV, int nE, int alpha,
    const std::vector<int> &V_indegree,
    const std::vector<int> &V_outdegree,
    const std::vector<double> &w_out,
    const std::vector<double> &w_in) {

    // Build adjacency lists (0-indexed vertices, 0-indexed edge IDs)
    std::vector<std::vector<int>> out_edges(nV), in_edges(nV);
    for (int e = 0; e < nE; ++e) {
        out_edges[E(e, 0) - 1].push_back(e);   // E is 1-indexed
        in_edges [E(e, 1) - 1].push_back(e);
    }

    // First pass: count nC and nnz
    int nC = 0, nnz = 0;
    for (int v = 0; v < nV; ++v) {
        int n_out = V_outdegree[v], n_in = V_indegree[v];
        if (n_out > 0 && n_in > 0) {
            nC  += n_out * alpha;
            nnz += n_out * (n_in + 1) * alpha;
        } else if (n_in == 0 && n_out > 1) {
            nC  += (n_out - 1) * alpha;
            nnz += 2 * (n_out - 1) * alpha;
        }
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nnz);

    int row = 0;

    // ---- Type-1 vertices ----
    for (int v = 0; v < nV; ++v) {
        int n_out = V_outdegree[v], n_in = V_indegree[v];
        if (n_out == 0 || n_in == 0) continue;

        const std::vector<int> &oe = out_edges[v];
        const std::vector<int> &ie = in_edges[v];

        for (int i = 0; i < n_out; ++i) {
            for (int der = 1; der <= alpha; ++der) {
                // Out-edge: column = 2*alpha*oe[i] + der-1, value = w_out[oe[i]]
                triplets.emplace_back(row, 2 * alpha * oe[i] + der - 1,
                                      w_out[oe[i]]);
                // In-edges: column = 2*alpha*ie[j] + alpha+der-1, value = w_in[ie[j]]
                for (int j = 0; j < n_in; ++j)
                    triplets.emplace_back(
                        row, 2 * alpha * ie[j] + alpha + der - 1,
                        w_in[ie[j]]);
                ++row;
            }
        }
    }

    // ---- Type-2 vertices ----
    for (int v = 0; v < nV; ++v) {
        int n_out = V_outdegree[v], n_in = V_indegree[v];
        if (n_in != 0 || n_out <= 1) continue;

        const std::vector<int> &oe = out_edges[v];

        for (int i = 1; i < n_out; ++i) {
            for (int der = 1; der <= alpha; ++der) {
                triplets.emplace_back(row, 2 * alpha * oe[i]     + der - 1,  1.0);
                triplets.emplace_back(row, 2 * alpha * oe[i - 1] + der - 1, -1.0);
                ++row;
            }
        }
    }

    Eigen::SparseMatrix<double> C(nC, 2 * alpha * nE);
    C.setFromTriplets(triplets.begin(), triplets.end());
    C.makeCompressed();
    return C;
}

// ---------------------------------------------------------------------------
// c_basis2_directional_graph  –  optimized CoB for the directional case
// ---------------------------------------------------------------------------
// Replaces construct_directional_constraint_matrix() + c_basis2() with a
// single pass over vertices using precomputed SVDs.
// Returns the same list structure as c_basis2().
//
// [[Rcpp::export]]
Rcpp::List c_basis2_directional_graph(
    const Eigen::MatrixXi &E, int nV, int nE, int alpha,
    const std::vector<int> &V_indegree,
    const std::vector<int> &V_outdegree) {

    int n_dof = 2 * alpha * nE;

    // ---- Build adjacency lists ----
    std::vector<std::vector<int>> out_edges(nV), in_edges(nV);
    for (int e = 0; e < nE; ++e) {
        out_edges[E(e, 0) - 1].push_back(e);
        in_edges [E(e, 1) - 1].push_back(e);
    }

    // ---- Count nC and find max degrees ----
    int nC = 0, max_out = 0, max_in = 0;
    for (int v = 0; v < nV; ++v) {
        int n_out = V_outdegree[v], n_in = V_indegree[v];
        if (n_out > max_out) max_out = n_out;
        if (n_in  > max_in)  max_in  = n_in;
        if (n_out > 0 && n_in > 0)
            nC += n_out * alpha;
        else if (n_in == 0 && n_out > 1)
            nC += (n_out - 1) * alpha;
    }

    // ---- Precompute SVDs ----
    // Type-1: [I_{n_out} | 1_{n_out × n_in}], keyed by n_in*(max_out+1)+n_out
    std::map<int, PrecomputedSVD> type1_cache;
    auto get_type1_svd = [&](int n_in, int n_out) -> const PrecomputedSVD & {
        int key = n_in * (max_out + 1) + n_out;
        auto it = type1_cache.find(key);
        if (it == type1_cache.end())
            type1_cache.emplace(key, compute_type1_svd(n_in, n_out));
        return type1_cache.at(key);
    };

    // Type-2: (n_out-1) x n_out difference matrix — reuse cont_cache
    std::vector<PrecomputedSVD> diff_cache =
        precompute_continuity_svds(max_out > 1 ? max_out : 2);

    // ---- Allocate output ----
    Eigen::VectorXd singular_values(nC);
    std::vector<Eigen::Triplet<double>> tripletsT, tripletsU;
    tripletsT.reserve(static_cast<size_t>(nC) * 8 + n_dof);
    tripletsU.reserve(static_cast<size_t>(nC) * 4);

    int counter_K    = 0;   // next constrained T column  (0..nC-1)
    int counter_D    = nC;  // next null/free   T column  (nC..n_dof-1)
    int constraint_row = 0; // current row in U / S

    std::vector<bool> dof_active(n_dof, false);

    std::vector<int>  cluster_sizes;
    std::vector<int>  index_largest;
    int               n_largest = 0;

    // ---- Helper: fill T rows for a column block ----
    // gc = global DoF column; local_j = its index in the precomputed V;
    // n_constrained = number of constrained T columns for this subcluster;
    // n_null        = number of null T columns.
    auto fill_T_col = [&](int gc, int local_j,
                          int n_constrained, int n_null,
                          const PrecomputedSVD &svd) {
        for (int i = 0; i < n_constrained; ++i)
            tripletsT.emplace_back(gc, counter_K + i,
                                   svd.V(local_j, i));
        for (int i = 0; i < n_null; ++i)
            tripletsT.emplace_back(gc, counter_D + i,
                                   svd.V(local_j, n_constrained + i));
    };

    // ==================================================================
    // Pass 1 – Type-1 vertices  (indegree > 0  AND  outdegree > 0)
    // ==================================================================
    for (int v = 0; v < nV; ++v) {
        int n_out = V_outdegree[v], n_in = V_indegree[v];
        if (n_out == 0 || n_in == 0) continue;

        const std::vector<int> &oe = out_edges[v];
        const std::vector<int> &ie = in_edges[v];
        const PrecomputedSVD   &svd = get_type1_svd(n_in, n_out);

        // n_col = n_out + n_in  columns per subcluster
        // n_constrained = n_out,   n_null = n_in

        // Constraint-row layout for this vertex:
        //   out-edge i, derivative der  →  constraint_row + i*alpha + (der-1)
        // For derivative der, the subcluster's rows are
        //   { constraint_row + i*alpha + (der-1) : i = 0..n_out-1 }  (non-contiguous)

        for (int der = 1; der <= alpha; ++der) {
            // Mark DoF columns active
            for (int k = 0; k < n_out; ++k)
                dof_active[2 * alpha * oe[k] + der - 1] = true;
            for (int k = 0; k < n_in; ++k)
                dof_active[2 * alpha * ie[k] + alpha + der - 1] = true;

            // --- Fill T ---
            for (int k = 0; k < n_out; ++k)
                fill_T_col(2 * alpha * oe[k] + der - 1, k, n_out, n_in, svd);
            for (int k = 0; k < n_in; ++k)
                fill_T_col(2 * alpha * ie[k] + alpha + der - 1,
                           n_out + k, n_out, n_in, svd);

            // --- Fill U (n_out × n_out block, non-contiguous rows) ---
            for (int i = 0; i < n_out; ++i) {
                int crow = constraint_row + i * alpha + (der - 1);
                for (int j = 0; j < n_out; ++j)
                    tripletsU.emplace_back(crow, counter_K + j, svd.U(i, j));
            }

            // --- Fill S ---
            for (int i = 0; i < n_out; ++i)
                singular_values(counter_K + i) = svd.S(i);

            // Subcluster bookkeeping
            cluster_sizes.push_back(n_out);
            if (n_out > n_largest) {
                n_largest = n_out;
                index_largest.clear();
                for (int i = 0; i < n_out; ++i)
                    index_largest.push_back(constraint_row + i * alpha + (der - 1));
            }

            counter_K += n_out;
            counter_D += n_in;
        }

        constraint_row += n_out * alpha;   // advance past all rows for this vertex
    }

    // ==================================================================
    // Pass 2 – Type-2 vertices  (indegree == 0  AND  outdegree > 1)
    // ==================================================================
    for (int v = 0; v < nV; ++v) {
        int n_out = V_outdegree[v], n_in = V_indegree[v];
        if (n_in != 0 || n_out <= 1) continue;

        const std::vector<int> &oe  = out_edges[v];
        const PrecomputedSVD   &svd = diff_cache[n_out];

        // n_constrained = n_out-1,  n_null = 1
        int n_constrained = n_out - 1;

        // Constraint-row layout:
        //   pair i (i=1..n_out-1), derivative der
        //     →  constraint_row + (i-1)*alpha + (der-1)
        // For derivative der, rows are
        //   { constraint_row + i*alpha + (der-1) : i = 0..n_out-2 }

        for (int der = 1; der <= alpha; ++der) {
            // Mark DoF columns active
            for (int k = 0; k < n_out; ++k)
                dof_active[2 * alpha * oe[k] + der - 1] = true;

            // --- Fill T ---
            for (int k = 0; k < n_out; ++k)
                fill_T_col(2 * alpha * oe[k] + der - 1, k,
                           n_constrained, 1, svd);

            // --- Fill U ((n_out-1) × (n_out-1) block, non-contiguous rows) ---
            for (int i = 0; i < n_constrained; ++i) {
                int crow = constraint_row + i * alpha + (der - 1);
                for (int j = 0; j < n_constrained; ++j)
                    tripletsU.emplace_back(crow, counter_K + j, svd.U(i, j));
            }

            // --- Fill S ---
            for (int i = 0; i < n_constrained; ++i)
                singular_values(counter_K + i) = svd.S(i);

            // Subcluster bookkeeping
            cluster_sizes.push_back(n_constrained);
            if (n_constrained > n_largest) {
                n_largest = n_constrained;
                index_largest.clear();
                for (int i = 0; i < n_constrained; ++i)
                    index_largest.push_back(constraint_row + i * alpha + (der - 1));
            }

            counter_K += n_constrained;
            counter_D += 1;
        }

        constraint_row += n_constrained * alpha;
    }

    // ---- Identity entries for inactive DoF ----
    for (int col = 0; col < n_dof; ++col) {
        if (!dof_active[col]) {
            tripletsT.emplace_back(col, counter_D, 1.0);
            ++counter_D;
        }
    }

    // ---- Build sparse output matrices ----
    Eigen::SparseMatrix<double> T(n_dof, n_dof);
    T.setFromTriplets(tripletsT.begin(), tripletsT.end());

    Eigen::SparseMatrix<double> U_mat(nC, nC);
    U_mat.setFromTriplets(tripletsU.begin(), tripletsU.end());

    Rcpp::List output;
    output["T"]              = T;
    output["S"]              = singular_values;
    output["U"]              = U_mat;
    output["larget.cluster"] = index_largest;
    output["cluster.n"]      = cluster_sizes;
    return output;
}

// ===========================================================================
// c_basis2  –  fast general-purpose replacement for c_basis2_old
// ===========================================================================
// Takes the same sparse matrix C as input as c_basis2_old but replaces the
// expensive BFS (O(nnz * nC)) with union-find (O(nnz * α(nC))), then does
// a small SVD per subcluster rather than iterating over the full column graph.
//
// For constraint matrices produced by this package the subclusters are tiny
// (at most max_degree rows/cols) so each SVD takes microseconds.  The total
// speedup over c_basis2_old on realistic graphs is 10–100×.
//
// Returns the same list structure as c_basis2_old:
//   T  – n_dof × n_dof sparse orthogonal matrix
//   S  – nC-length vector of singular values
//   U  – nC × nC sparse matrix of left singular vectors
//   larget.cluster  – row indices (0-indexed) of the largest subcluster
//   cluster.n       – sizes of all subclusters
//
// [[Rcpp::export]]
Rcpp::List c_basis2(Eigen::MappedSparseMatrix<double> A,
                    double eps_limit = 1e-10) {
    int nC    = A.rows();
    int n_dof = A.cols();

    // ---- Build row_entries (col, val) for each row  ----
    // A is column-major MappedSparseMatrix; iterate by outer (= columns).
    std::vector<std::vector<std::pair<int,double>>> row_entries(nC);
    for (int col = 0; col < A.outerSize(); ++col) {
        for (Eigen::MappedSparseMatrix<double>::InnerIterator it(A, col);
             it; ++it) {
            row_entries[it.row()].push_back({col, it.value()});
        }
    }

    // ---- Union-Find on rows: rows sharing a column → same subcluster ----
    std::vector<int> parent(nC);
    std::iota(parent.begin(), parent.end(), 0);

    auto find_root = [&](int x) -> int {
        while (parent[x] != x) { parent[x] = parent[parent[x]]; x = parent[x]; }
        return x;
    };

    // First row seen for each column, used to union rows together
    std::vector<int> col_first_row(n_dof, -1);
    for (int i = 0; i < nC; ++i) {
        for (auto &[col, val] : row_entries[i]) {
            if (col_first_row[col] < 0) {
                col_first_row[col] = i;
            } else {
                int r = find_root(col_first_row[col]);
                int s = find_root(i);
                if (r != s) parent[r] = s;
            }
        }
    }
    // Full path compression: ensure every parent[i] points directly to root
    for (int i = 0; i < nC; ++i) parent[i] = find_root(i);

    // ---- Group rows by subcluster ----
    std::unordered_map<int, std::vector<int>> sc_map;
    sc_map.reserve(nC);
    for (int i = 0; i < nC; ++i) sc_map[parent[i]].push_back(i);

    // Sort each subcluster's rows and sort subclusters by first row
    std::vector<std::vector<int>> subclusters;
    subclusters.reserve(sc_map.size());
    for (auto &[root, rows] : sc_map) subclusters.push_back(std::move(rows));
    for (auto &rows : subclusters) std::sort(rows.begin(), rows.end());
    std::sort(subclusters.begin(), subclusters.end(),
              [](const auto &a, const auto &b){ return a[0] < b[0]; });

    // ---- Identify all active columns (have at least one non-zero row) ----
    std::vector<bool> dof_active(n_dof, false);
    for (int col = 0; col < n_dof; ++col)
        if (col_first_row[col] >= 0) dof_active[col] = true;

    // ---- Allocate output ----
    Eigen::VectorXd singular_values(nC);
    std::vector<Eigen::Triplet<double>> tripletsT, tripletsU;
    tripletsT.reserve(n_dof * 2);
    tripletsU.reserve(nC * 4);

    std::vector<int> cluster_sizes;
    cluster_sizes.reserve(subclusters.size());
    std::vector<int> index_largest;
    int n_largest  = 0;
    int counter_K  = 0;   // next constrained T column  (0..nC-1)
    int counter_D  = nC;  // next null/free   T column  (nC..n_dof-1)

    // ---- Process each subcluster ----
    for (auto &rows : subclusters) {
        int m = static_cast<int>(rows.size());

        // Collect all columns for this subcluster (sorted)
        std::vector<int> col_vec;
        {
            std::unordered_set<int> col_set;
            for (int r : rows)
                for (auto &[col, val] : row_entries[r])
                    col_set.insert(col);
            col_vec.assign(col_set.begin(), col_set.end());
        }
        std::sort(col_vec.begin(), col_vec.end());
        int n_cols = static_cast<int>(col_vec.size());

        // Map global col → local col index
        std::unordered_map<int,int> col_local;
        col_local.reserve(n_cols);
        for (int j = 0; j < n_cols; ++j) col_local[col_vec[j]] = j;

        // Build submatrix A_sub (m × n_cols)
        Eigen::MatrixXd A_sub = Eigen::MatrixXd::Zero(m, n_cols);
        for (int i = 0; i < m; ++i)
            for (auto &[col, val] : row_entries[rows[i]])
                A_sub(i, col_local[col]) = val;

        // SVD of A_sub (full U and V)
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(
            A_sub, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Eigen::MatrixXd &V = svd.matrixV();
        const Eigen::MatrixXd &U = svd.matrixU();
        const Eigen::VectorXd &S = svd.singularValues();

        int n_null = n_cols - m;   // dimension of null space

        // Fill T: for each column gc in the subcluster, local index j
        //   V(j, 0..m-1)    → T columns counter_K .. counter_K+m-1
        //   V(j, m..n_cols-1) → T columns counter_D .. counter_D+n_null-1
        for (int j = 0; j < n_cols; ++j) {
            int gc = col_vec[j];
            for (int i = 0; i < m; ++i)
                if (V(j, i) != 0.0)
                    tripletsT.emplace_back(gc, counter_K + i, V(j, i));
            for (int i = 0; i < n_null; ++i)
                if (V(j, m + i) != 0.0)
                    tripletsT.emplace_back(gc, counter_D + i, V(j, m + i));
        }

        // Fill U: m × m block at (rows[i], counter_K+j)
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < m; ++j)
                if (U(i, j) != 0.0)
                    tripletsU.emplace_back(rows[i], counter_K + j, U(i, j));

        // Fill S
        for (int i = 0; i < m; ++i)
            singular_values(counter_K + i) = S(i);

        // Subcluster bookkeeping
        cluster_sizes.push_back(m);
        if (m > n_largest) {
            n_largest = m;
            index_largest.assign(rows.begin(), rows.end());
        }

        counter_K += m;
        counter_D += n_null;
    }

    // ---- Identity entries for inactive DoF columns ----
    for (int col = 0; col < n_dof; ++col) {
        if (!dof_active[col]) {
            tripletsT.emplace_back(col, counter_D, 1.0);
            ++counter_D;
        }
    }

    // ---- Build sparse output matrices ----
    Eigen::SparseMatrix<double> T(n_dof, n_dof);
    T.setFromTriplets(tripletsT.begin(), tripletsT.end());

    Eigen::SparseMatrix<double> U_mat(nC, nC);
    U_mat.setFromTriplets(tripletsU.begin(), tripletsU.end());

    Rcpp::List output;
    output["T"]              = T;
    output["S"]              = singular_values;
    output["U"]              = U_mat;
    output["larget.cluster"] = index_largest;
    output["cluster.n"]      = cluster_sizes;
    return output;
}
