#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// 

//' @name assemble_fem
//' @title Construction of FEM matrices
//' @description Function used to construct FEM matrices on metric graphs.
//' @param E [nx2 matrix] Matrix of edges
//' @param h_e [n vector] Vector of h's
//' @param nV [int] Number of vertices
//' @noRd
//'
// [[Rcpp::export]]

Rcpp::List assemble_fem(Eigen::MatrixXd E, Eigen::VectorXd h_e, int nV){

    typedef Eigen::Triplet<double> Trip;
    std::vector<Trip> trp_C, trp_G, trp_B;
    int i;
    int v1,v2;

    Eigen::SparseMatrix<double> C(nV,nV), G(nV,nV), B(nV,nV);

    for(i=0; i<E.rows(); i++){
        v1 = E(i, 0)-1;
        v2 = E(i, 1)-1;

        // Assembling C
        trp_C.push_back(Trip(v1,v1,h_e(i)/3));
        trp_C.push_back(Trip(v2,v2,h_e(i)/3));
        trp_C.push_back(Trip(v1,v2,h_e(i)/6));
        trp_C.push_back(Trip(v2,v1,h_e(i)/6));

        // Assembling G
        
        trp_G.push_back(Trip(v1,v1,1/h_e(i)));
        trp_G.push_back(Trip(v2,v2,1/h_e(i)));
        trp_G.push_back(Trip(v1,v2,-1/h_e(i)));
        trp_G.push_back(Trip(v2,v1,-1/h_e(i)));

        // Assembling B
        
        trp_B.push_back(Trip(v1,v1,-0.5));
        trp_B.push_back(Trip(v1,v2,-0.5));
        trp_B.push_back(Trip(v2,v2,0.5));
        trp_B.push_back(Trip(v2,v1,0.5));   
    }

    C.setFromTriplets(trp_C.begin(), trp_C.end());
    G.setFromTriplets(trp_G.begin(), trp_G.end());          
    B.setFromTriplets(trp_B.begin(), trp_B.end());      

    Rcpp::List out;
    out["C"] = C;
    out["G"] = G;
    out["B"] = B;
    return(out);
}