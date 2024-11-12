#include <Rcpp.h>
#include <RcppEigen.h>
#include <unordered_set>
#include <algorithm>
#include <iterator>
#include <Eigen/Sparse>
#include <vector>

using namespace Rcpp;


typedef Eigen::SparseMatrix<double> SparseMatrix;
typedef Eigen::Triplet<double> Triplet;

// [[Rcpp::depends(RcppEigen)]]

void unique_vector_v2(std::vector<int> &veci) {
  std::unordered_set<int> set(veci.begin(), veci.end());
  veci.assign(set.begin(), set.end());
  std::sort(veci.begin(), veci.end());
}

void set_diff_v2(const std::vector<int> &A, const std::vector<int> &B, std::vector<int> &C) {
  std::unordered_set<int> setA(A.begin(), A.end());
  for (const int &b : B) {
    setA.erase(b);
  }
  C.assign(setA.begin(), setA.end());
  std::sort(C.begin(), C.end());
}

// [[Rcpp::export]]
Rcpp::List c_basis2_v2(Eigen::MappedSparseMatrix<double> A, double eps_limit = 1e-10) {
  Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> P(A.cols());
  P.setIdentity();
  std::vector<int> index;
  index.reserve(A.cols()); // Reserve memory to avoid reallocations
  
  int counter = 0;
  for (int k = 0; k < A.outerSize(); ++k) {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
      index.push_back(k);
      std::swap(P.indices()[counter], P.indices()[k]);
      counter++;
      break;
    }
  }
  
  Eigen::SparseMatrix<double> A_ID = A * P;
  std::vector<int> index_A(A.rows());
  std::iota(index_A.begin(), index_A.end(), 0); // Fill with 0,1,...,A.rows()-1

  std::vector<Eigen::Triplet<double>> tripletListT, tripletListU;
  tripletListT.reserve(counter * counter + (A.cols() - counter));
  tripletListU.reserve(A.rows() * A.rows());

  Eigen::SparseMatrix<double, Eigen::RowMajor> A_ID_rowm(A_ID);
  Eigen::VectorXd singular_values(A.rows());
  
  int count_subcluster = 0, n_largest_cluster = 0;
  std::vector<int> index_largest, n_subcluster;
  int counter_K = 0, counter_D = A.rows();

  while (!index_A.empty()) {
    count_subcluster++;
    
    std::vector<int> index_temp(1, index_A[0]), index_new(1, index_A[0]), col_index_A;
    
    while (!index_new.empty()) {
      std::vector<int> col_index;
      for (const int &idx : index_new) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_A(A_ID_rowm, idx); it_A; ++it_A) {
          col_index.push_back(it_A.col());
        }
      }
      
      unique_vector_v2(col_index);
      
      std::vector<int> row_index;
      for (const int &col : col_index) {
        for (Eigen::SparseMatrix<double, Eigen::ColMajor>::InnerIterator it_A(A_ID, col); it_A; ++it_A) {
          row_index.push_back(it_A.row());
        }
      }
      
      unique_vector_v2(row_index);
      set_diff_v2(row_index, index_temp, index_new);
      index_temp.insert(index_temp.end(), index_new.begin(), index_new.end());
      col_index_A.insert(col_index_A.end(), col_index.begin(), col_index.end());
    }
    
    unique_vector_v2(col_index_A);
    std::vector<int> col_index_full(counter + 1, -1); // Initialize with -1 for unassigned elements
    for (size_t i = 0; i < col_index_A.size(); ++i) {
      col_index_full[col_index_A[i]] = i;
    }
    
    Eigen::MatrixXd A_id_temp = Eigen::MatrixXd::Zero(index_temp.size(), col_index_A.size());
    if (index_temp.size() > n_largest_cluster) {
      n_largest_cluster = index_temp.size();
      index_largest = index_temp;
    }
    
    n_subcluster.push_back(index_temp.size());
    for (size_t i = 0; i < index_temp.size(); ++i) {
      for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it_A(A_ID_rowm, index_temp[i]); it_A; ++it_A) {
        A_id_temp(i, col_index_full[it_A.col()]) = it_A.value();
      }
    }
    
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(A_id_temp, Eigen::ComputeFullV | Eigen::ComputeFullU);
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::VectorXd s = svd.singularValues();
    
    for (int i = 0; i < s.size(); ++i) {
      singular_values(counter_K + i) = s(i);
      for (int j = 0; j < index_temp.size(); ++j) {
        tripletListU.emplace_back(index_temp[j], counter_K + i, U(j, i));
      }
    }
    
    for (int i = 0; i < col_index_A.size(); ++i) {
      int i_j = (i < index_temp.size()) ? counter_K++ : counter_D++;
      for (int j = 0; j < col_index_A.size(); ++j) {
        tripletListT.emplace_back(index[col_index_A[j]], i_j, V(j, i));
      }
    }
    
    set_diff_v2(index_A, index_temp, index_A);
  }
  
  // Build the basis matrix T
  for (int i = counter; i < A.outerSize(); ++i) {
    int i_i = P.indices()[i];
    tripletListT.emplace_back(i_i, i, 1.0);
  }
  
  Eigen::SparseMatrix<double> T(A.cols(), A.cols());
  T.setFromTriplets(tripletListT.begin(), tripletListT.end());
  
  Eigen::SparseMatrix<double> U(A.rows(), A.rows());
  U.setFromTriplets(tripletListU.begin(), tripletListU.end());
  
  return Rcpp::List::create(
    Rcpp::Named("T") = T,
    Rcpp::Named("S") = singular_values,
    Rcpp::Named("U") = U,
    Rcpp::Named("largest.cluster") = index_largest,
    Rcpp::Named("cluster.n") = n_subcluster
  );
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> construct_constraint_matrix(const Eigen::MatrixXi& E, int nV, int edge_constraint) {
    int nE = E.rows();  // Number of edges inferred from the rows of E

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
        const std::vector<int>& le = lower_edges[v];
        const std::vector<int>& ue = upper_edges[v];
        int n_e = le.size() + ue.size();

        // Derivative constraint
        if ((edge_constraint && n_e == 1) || n_e > 1) {
            for (int k = 0; k < n_e; ++k) {
                i_.push_back(count_constraint);  // Use zero-based indexing
                if (k < le.size()) {
                    j_.push_back(4 * (le[k]) + 1);  // Corrected to zero-based
                    x_.push_back(1.0);
                } else {
                    j_.push_back(4 * (ue[k - le.size()]) + 3);  // Corrected to zero-based
                    x_.push_back(-1.0);
                }
                count++;
            }
            count_constraint++;
        }

        // Internal constraints for nodes with more than one edge
        if (n_e > 1) {
            std::vector<std::pair<int, int>> edges;
            for (int e : le) edges.emplace_back(e, 1);  // Adjusted to zero-based
            for (int e : ue) edges.emplace_back(e, 3);  // Adjusted to zero-based

            for (int i = 1; i < n_e; ++i) {
                i_.push_back(count_constraint);  // Zero-based indexing
                j_.push_back(4 * (edges[i - 1].first) + edges[i - 1].second - 1);  // Adjusted to zero-based
                x_.push_back(1.0);

                i_.push_back(count_constraint);  // Zero-based indexing
                j_.push_back(4 * (edges[i].first) + edges[i].second - 1);  // Adjusted to zero-based
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

    // Create the sparse matrix `C` with dimensions based on the constraints and edges
    Eigen::SparseMatrix<double> C(count_constraint, 4 * nE);
    C.setFromTriplets(tripletList.begin(), tripletList.end());
    C.makeCompressed();  // Ensure it is in compressed column storage

    return C;
}