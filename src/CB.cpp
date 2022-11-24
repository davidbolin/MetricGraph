
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]


//' @name unique_vector
//' @noRd
//'
// sorts integervector
void unique_vector(std::vector<int> & veci){
  std::sort(veci.begin(), veci.end());
  std::vector<int>::iterator  it_unique = std::unique(veci.begin(), veci.end());
  veci.resize( std::distance(veci.begin(),it_unique) );

}
//' @name set_diff
//' @noRd
//'
// setdiff(A,B) stored in C
void set_diff(std::vector<int> & A,
              std::vector<int> & B,
              std::vector<int> & C){
  C.resize(0);
  std::sort(A.begin(), A.end());
  std::sort(B.begin(), B.end());
  std::set_difference(A.begin(), A.end(),
                      B.begin(), B.end(),
                      std::inserter(C, C.begin()));
}

//' @name c_basis2
//' @title CB construction
//' @description The SVD-based constraint basis construction for non-overlapping
//' subsets of constraints. Algorithm 1 from the reference.
//' Creating a basis from the matrix A
//' @param A [nxk matrix] must have rank k
//' @param eps_limit [double] used as a limit of small value
//' @return T (n x n) the basis matrix
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List  c_basis2(Eigen::MappedSparseMatrix<double> A,
                             double eps_limit = 1e-10) {

  Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> P(A.cols());
  P.setIdentity();
  std::vector<int> index(0);
  int counter = 0;
  for(int k = 0; k < A.outerSize(); ++k) {
    for(Eigen::MappedSparseMatrix<double>::InnerIterator it(A,k); it; ++it) {
      index.push_back(k);
      std::swap(P.indices()[counter],P.indices()[k]);
      counter++;
      break;
    }
  }
  std::vector< Eigen::Triplet<double> > tripletListT;
  tripletListT.reserve(counter*counter + (A.cols()-counter));
  Eigen::SparseMatrix<double>  A_ID = A *P;
  std::vector<int> index_A(A.rows());
  for(int i =0; i < A.rows(); i++)
    index_A[i]  = i;


  // creating a indexing putting relevant columns first so the first
  // k columns spanns A
  int counter_K = 0;
  int counter_D = A.rows();
  //
  std::vector< Eigen::Triplet<double> > tripletListU;
  Eigen::VectorXd singular_values(A.rows());
  //
  Eigen::SparseMatrix<double,Eigen::RowMajor> A_ID_rowm(A_ID);
  int count_subcluster = 0;
  std::vector<int> n_subcluster(0);
  std::vector<int> index_largest;
  int n_largest_cluster =0 ;
  while(index_A.size()!=0){
    count_subcluster++;
    std::vector<int> index_temp(1);
    std::vector<int> index_new(1);
    std::vector<int> col_index_A(0);
    index_temp[0] = index_A[0];
    index_new[0]  = index_A[0];

    while(index_new.size()>0){
      std::vector<int> col_index(0);
      for (std::vector<int>::iterator it = index_new.begin() ; it != index_new.end(); ++it){
        for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it_A(A_ID_rowm, *it); it_A; ++it_A){
          col_index.push_back(it_A.col());
        }
      }
      unique_vector(col_index);

      std::vector<int> row_index(0);
      for (std::vector<int>::iterator it = col_index.begin() ; it != col_index.end(); ++it){
        for (Eigen::SparseMatrix<double,Eigen::ColMajor>::InnerIterator it_A(A_ID, *it); it_A; ++it_A)
          row_index.push_back(it_A.row());
      }

      unique_vector(row_index);
      set_diff(row_index, index_temp, index_new);
      for(std::vector<int>::iterator  it  = index_new.begin();
                                      it != index_new.end();
                                      ++it){
        index_temp.push_back(*it);
      }
      for(int i=0; i < col_index.size(); i++)
        col_index_A.push_back(col_index[i]);
    }
    unique_vector(col_index_A);
    std::vector<int> col_index_full(counter+1);
    for(int i = 0; i < col_index_A.size(); i++)
      col_index_full[col_index_A[i]] = i;

    Eigen::MatrixXd A_id_temp = Eigen::MatrixXd::Zero(index_temp.size(),
                                                      col_index_A.size());
    int i_temp = 0;
    if( index_temp.size() > n_largest_cluster){
      n_largest_cluster = index_temp.size();
      index_largest = index_temp;

    }

    n_subcluster.push_back(index_temp.size());
    for (std::vector<int>::iterator it = index_temp.begin() ; it != index_temp.end(); ++it){
      for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it_A(A_ID_rowm, *it); it_A; ++it_A){
        A_id_temp(i_temp, col_index_full[it_A.col()]) = it_A.value();
      }
      i_temp++;
    }
    Eigen::JacobiSVD<Eigen::MatrixXd> svd( A_id_temp, Eigen::ComputeFullV | Eigen::ComputeFullU );
    Eigen::MatrixXd V = svd.matrixV();
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::VectorXd s = svd.singularValues();


    for(int i=0; i < s.size();i++)
      singular_values(counter_K + i) = s(i);
    for(int i=0; i < s.size();i++){
      int i_j = counter_K + i;
      for(int j=0; j< s.size(); j++)
        tripletListU.push_back(Eigen::Triplet<double>(index_temp[j], i_j, U(j,i)));
    }


    for(int i=0; i< col_index_A.size(); i++){

      int i_j ;
      if(i < index_temp.size()){
        i_j = counter_K;
        counter_K++;
      }else{
        i_j = counter_D;
        counter_D++;
      }
      for(int j=0; j< col_index_A.size(); j++){
        tripletListT.push_back(Eigen::Triplet<double>(index[col_index_A[j]] ,i_j, V(j,i)));
      }

    }

    std::vector<int> index_A_copy(index_A.size());
    std::copy ( index_A.data(),
                index_A.data()+index_A.size(),
                index_A_copy.begin() );
    set_diff(index_A_copy,
             index_temp,
             index_A);
  }

  //
  // building the basis
  //
  for(int i = counter; i < A.outerSize(); i++){
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
  //output["number cluster"] = count_subcluster;
  output["cluster.n"] =  n_subcluster;
  return(output);
}
