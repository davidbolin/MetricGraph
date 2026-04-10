
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

//' @name unique_vector
//' @noRd
//'
// sorts integervector
void unique_vector(std::vector<int> &veci) {
  std::sort(veci.begin(), veci.end());
  std::vector<int>::iterator it_unique = std::unique(veci.begin(), veci.end());
  veci.resize(std::distance(veci.begin(), it_unique));
}
//' @name set_diff
//' @noRd
//'
// setdiff(A,B) stored in C
void set_diff(std::vector<int> &A, std::vector<int> &B, std::vector<int> &C) {
  C.resize(0);
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  std::set_difference(A.begin(), A.end(), B.begin(), B.end(),
                      std::inserter(C, C.begin()));
}

//' @name c_basis2
//' @title CB construction
//' @description The SVD-based constraint basis construction for non-overlapping
//' subsets of constraints. Algorithm 1 from the reference.
//' Creating a basis from the matrix A
//' @param A `nxk matrix` must have rank k
//' @param eps_limit `double` used as a limit of small value
//' @return T (n x n) the basis matrix
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List c_basis2(Eigen::MappedSparseMatrix<double> A,
                    double eps_limit = 1e-10) {

  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P(A.cols());
  P.setIdentity();
  std::vector<int> index(0);
  int counter = 0;
  for (int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
      index.push_back(k);
      std::swap(P.indices()[counter], P.indices()[k]);
      counter++;
      break;
    }
  }
  std::vector<Eigen::Triplet<double>> tripletListT;
  // tripletListT.reserve(static_cast<size_t>(counter) *
  //                          static_cast<size_t>(counter) +
  //                      (A.cols() - counter));
  Eigen::SparseMatrix<double> A_ID = A * P;
  std::vector<int> index_A(A.rows());
  for (int i = 0; i < A.rows(); i++)
    index_A[i] = i;

  // creating a indexing putting relevant columns first so the first
  // k columns spanns A
  int counter_K = 0;
  int counter_D = A.rows();
  //
  std::vector<Eigen::Triplet<double>> tripletListU;
  Eigen::VectorXd singular_values(A.rows());
  //
  Eigen::SparseMatrix<double, Eigen::RowMajor> A_ID_rowm(A_ID);
  int count_subcluster = 0;
  std::vector<int> n_subcluster(0);
  std::vector<int> index_largest;
  int n_largest_cluster = 0;
  while (index_A.size() != 0) {
    count_subcluster++;
    std::vector<int> index_temp(1);
    std::vector<int> index_new(1);
    std::vector<int> col_index_A(0);
    index_temp[0] = index_A[0];
    index_new[0] = index_A[0];

    while (index_new.size() > 0) {
      std::vector<int> col_index(0);
      for (std::vector<int>::iterator it = index_new.begin();
           it != index_new.end(); ++it) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_A(
                 A_ID_rowm, *it);
             it_A; ++it_A) {
          col_index.push_back(it_A.col());
        }
      }
      unique_vector(col_index);

      std::vector<int> row_index(0);
      for (std::vector<int>::iterator it = col_index.begin();
           it != col_index.end(); ++it) {
        for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it_A(
                 A_ID, *it);
             it_A; ++it_A)
          row_index.push_back(it_A.row());
      }

      unique_vector(row_index);
      set_diff(row_index, index_temp, index_new);
      for (std::vector<int>::iterator it = index_new.begin();
           it != index_new.end(); ++it) {
        index_temp.push_back(*it);
      }
      for (int i = 0; i < col_index.size(); i++)
        col_index_A.push_back(col_index[i]);
    }
    unique_vector(col_index_A);
    std::vector<int> col_index_full(counter + 1);
    for (int i = 0; i < col_index_A.size(); i++)
      col_index_full[col_index_A[i]] = i;

    Eigen::MatrixXd A_id_temp =
        Eigen::MatrixXd::Zero(index_temp.size(), col_index_A.size());
    int i_temp = 0;
    if (index_temp.size() > n_largest_cluster) {
      n_largest_cluster = index_temp.size();
      index_largest = index_temp;
    }

    n_subcluster.push_back(index_temp.size());
    for (std::vector<int>::iterator it = index_temp.begin();
         it != index_temp.end(); ++it) {
      for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_A(
               A_ID_rowm, *it);
           it_A; ++it_A) {
        A_id_temp(i_temp, col_index_full[it_A.col()]) = it_A.value();
      }
      i_temp++;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A_id_temp, Eigen::ComputeFullV |
                                                         Eigen::ComputeFullU);
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::VectorXd s = svd.singularValues();

    for (int i = 0; i < s.size(); i++)
      singular_values(counter_K + i) = s(i);
    for (int i = 0; i < s.size(); i++) {
      int i_j = counter_K + i;
      for (int j = 0; j < s.size(); j++)
        tripletListU.push_back(
            Eigen::Triplet<double>(index_temp[j], i_j, U(j, i)));
    }

    for (int i = 0; i < col_index_A.size(); i++) {

      int i_j;
      if (i < index_temp.size()) {
        i_j = counter_K;
        counter_K++;
      } else {
        i_j = counter_D;
        counter_D++;
      }
      for (int j = 0; j < col_index_A.size(); j++) {
        tripletListT.push_back(
            Eigen::Triplet<double>(index[col_index_A[j]], i_j, V(j, i)));
      }
    }

    std::vector<int> index_A_copy(index_A.size());
    std::copy(index_A.data(), index_A.data() + index_A.size(),
              index_A_copy.begin());
    set_diff(index_A_copy, index_temp, index_A);
  }

  //
  // building the basis
  //
  for (int i = counter; i < A.outerSize(); i++) {
    int i_i = P.indices()[i];
    tripletListT.push_back(Eigen::Triplet<double>(i_i, i, 1.));
  }
  Eigen::SparseMatrix<double> T;
  T.resize(A.cols(), A.cols());
  T.setFromTriplets(tripletListT.begin(), tripletListT.end());
  Eigen::SparseMatrix<double> U;
  U.resize(A.rows(), A.rows());
  U.setFromTriplets(tripletListU.begin(), tripletListU.end());
  Rcpp::List output;
  output["T"] = T;
  output["S"] = singular_values;
  output["U"] = U;
  output["larget.cluster"] = index_largest;
  // output["number cluster"] = count_subcluster;
  output["cluster.n"] = n_subcluster;
  return (output);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double>
construct_constraint_matrix(const Eigen::MatrixXi &E, int nV,
                            int edge_constraint) {
  int nE = E.rows(); // Number of edges inferred from the rows of E

  // Reserve memory based on expected constraints
  std::vector<int> i_;
  std::vector<int> j_;
  std::vector<double> x_;
  i_.reserve(2 * nE);
  j_.reserve(2 * nE);
  x_.reserve(2 * nE);

  int count_constraint = 0;
  int count = 0;

  // Precompute lower and upper edge indices with memory consideration
  std::vector<std::vector<int>> lower_edges(nV + 1), upper_edges(nV + 1);
  for (int e = 0; e < nE; ++e) {
    lower_edges[E(e, 0)].push_back(e);
    upper_edges[E(e, 1)].push_back(e);
  }

  // Loop over each vertex
  for (int v = 1; v <= nV; ++v) {
    const std::vector<int> &le = lower_edges[v];
    const std::vector<int> &ue = upper_edges[v];
    int n_e = le.size() + ue.size();

    // Derivative constraint
    if ((edge_constraint && n_e == 1) || n_e > 1) {
      for (int k = 0; k < n_e; ++k) {
        i_.push_back(count_constraint); // Use zero-based indexing
        if (k < le.size()) {
          j_.push_back(4 * (le[k]) + 1); // Corrected to zero-based
          x_.push_back(1.0);
        } else {
          j_.push_back(4 * (ue[k - le.size()]) + 3); // Corrected to zero-based
          x_.push_back(-1.0);
        }
        count++;
      }
      count_constraint++;
    }

    // Internal constraints for nodes with more than one edge
    if (n_e > 1) {
      std::vector<std::pair<int, int>> edges;
      for (int e : le)
        edges.emplace_back(e, 1); // Adjusted to zero-based
      for (int e : ue)
        edges.emplace_back(e, 3); // Adjusted to zero-based

      for (int i = 1; i < n_e; ++i) {
        i_.push_back(count_constraint); // Zero-based indexing
        j_.push_back(4 * (edges[i - 1].first) + edges[i - 1].second -
                     1); // Adjusted to zero-based
        x_.push_back(1.0);

        i_.push_back(count_constraint); // Zero-based indexing
        j_.push_back(4 * (edges[i].first) + edges[i].second -
                     1); // Adjusted to zero-based
        x_.push_back(-1.0);

        count_constraint++;
        count += 2;
      }
    }
  }

  // Populate triplet list and construct the sparse matrix
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(count);
  for (int k = 0; k < count; ++k) {
    tripletList.emplace_back(i_[k], j_[k], x_[k]);
  }

  // Create the sparse matrix `C` with dimensions based on the constraints and
  // edges
  Eigen::SparseMatrix<double> C(count_constraint, 4 * nE);
  C.setFromTriplets(tripletList.begin(), tripletList.end());
  C.makeCompressed(); // Ensure it is in compressed column storage

  return C;
}

// Gets the correct sparsity however cannot fill correct x because we cannot
// pass the directional weight functions to C.

// [[Rcpp::export]]
Eigen::SparseMatrix<double> construct_directional_constraint_matrix(
    const Eigen::MatrixXi &E, int nV, int nE, int alpha,
    const std::vector<int> &V_indegree, const std::vector<int> &V_outdegree) {

  // Determine index conditions
  std::vector<bool> index_outdegree(nV, false);
  std::vector<bool> index_in0(nV, false);
  for (int v = 0; v < nV; ++v) {
    index_outdegree[v] = V_outdegree[v] > 0 && V_indegree[v] > 0;
    index_in0[v] = V_indegree[v] == 0;
  }

  // Calculate nC based on conditions
  int nC = 0;
  for (int v = 0; v < nV; ++v) {
    if (index_outdegree[v]) {
      nC += V_outdegree[v] * (1 + V_indegree[v]);
    } else if (index_in0[v]) {
      nC += V_outdegree[v] - 1;
    }
  }
  nC *= alpha;

  // Reserve memory based on expected constraints
  std::vector<int> i_;
  std::vector<int> j_;
  std::vector<double> x_;
  i_.reserve(nC);
  j_.reserve(nC);
  x_.reserve(nC);

  int count_constraint = 0;
  int count = 0;

  // Process vertices with outdegree and indegree
  for (int v = 0; v < nV; ++v) {
    if (index_outdegree[v]) {
      std::vector<int> out_edges;
      std::vector<int> in_edges;

      // Find out_edges and in_edges for vertex v
      for (int e = 0; e < nE; ++e) {
        if (E(e, 0) == v + 1)
          out_edges.push_back(e);
        if (E(e, 1) == v + 1)
          in_edges.push_back(e);
      }

      int n_in = in_edges.size();
      for (int i = 0; i < out_edges.size(); ++i) {
        for (int der = 1; der <= alpha; ++der) {
          i_.insert(i_.end(), n_in + 1, count_constraint);
          j_.push_back(2 * alpha * (out_edges[i]) + der - 1);

          // Add indices for each in_edge, and set all x_ values to 1
          for (int j = 0; j < n_in; ++j) {
            j_.push_back(2 * alpha * (in_edges[j]) + alpha + der - 1);
          }

          // Set x_ values to 1
          x_.insert(x_.end(), n_in + 1, 1.0);

          count += (n_in + 1);
          count_constraint++;
        }
      }
    }
  }

  // Process vertices with indegree == 0
  for (int v = 0; v < nV; ++v) {
    if (index_in0[v]) {
      std::vector<int> out_edges;
      for (int e = 0; e < nE; ++e) {
        if (E(e, 0) == v + 1)
          out_edges.push_back(e);
      }

      if (out_edges.size() > 1) {
        for (int i = 1; i < out_edges.size(); ++i) {
          for (int der = 1; der <= alpha; ++der) {
            i_.push_back(count_constraint);
            j_.push_back(2 * alpha * (out_edges[i]) + der - 1);
            x_.push_back(1.0);

            i_.push_back(count_constraint);
            j_.push_back(2 * alpha * (out_edges[i - 1]) + der - 1);
            x_.push_back(-1.0);

            count += 2;
            count_constraint++;
          }
        }
      }
    }
  }

  // Populate triplet list and construct the sparse matrix
  std::vector<Eigen::Triplet<double>> tripletList;
  tripletList.reserve(count);
  for (int k = 0; k < count; ++k) {
    tripletList.emplace_back(i_[k], j_[k], x_[k]);
  }

  Eigen::SparseMatrix<double> C(count_constraint, 2 * alpha * nE);
  C.setFromTriplets(tripletList.begin(), tripletList.end());
  C.makeCompressed();

  return C;
}
