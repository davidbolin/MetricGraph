#include <Rcpp.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

namespace {
    Eigen::SparseMatrix<double> selected_inv_aux(const Eigen::SparseMatrix<double>& Lq) {
        Eigen::SparseMatrix<double> S = Lq.selfadjointView<Eigen::Lower>();

        int n = Lq.rows();

        for (int i = n - 1; i >= 0; --i) {
            Eigen::SparseMatrix<double>::ReverseInnerIterator Si(S, i);
            for (Eigen::SparseMatrix<double>::ReverseInnerIterator ij(Lq, i); ij; --ij) {

                Eigen::SparseMatrix<double>::ReverseInnerIterator iL(Lq, i);
                Eigen::SparseMatrix<double>::ReverseInnerIterator iS(S, ij.row());

                Si.valueRef() = 0.0;
                while (iL && iL.row() > i) {
                    while (iS && (iL.row() < iS.row())) {
                        --iS;
                    }
                    if (iS && (iL.row() == iS.row())) {
                        Si.valueRef() -= iL.value() * iS.value();
                        --iS;
                    }
                    --iL;
                }

                if (i == ij.row()) {
                    Si.valueRef() += 1 / iL.value();
                    Si.valueRef() /= iL.value();
                } else {
                    Si.valueRef() /= iL.value();
                    while (iS && iS.row() > i) {
                        --iS;
                    }
                    iS.valueRef() = Si.value();
                }
                --Si;
            }
        }

        return S;
    }
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> selected_inv_cpp(const Eigen::SparseMatrix<double>& Q) {
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(Q);
    solver.factorize(Q);

    if (solver.info() != Eigen::Success) {
        Rcpp::stop("Cholesky decomposition failed");
    }

    Eigen::SparseMatrix<double> L = solver.matrixL();
    Eigen::SparseMatrix<double> result_perm = selected_inv_aux(L);

    // Correct permutation using solver.permutationPinv()
    return solver.permutationPinv() * result_perm * solver.permutationPinv().transpose();
}