// Fast full-graph directional alpha=1 plug-in LOO.
//
// cv_loo_selinv_cpp uses a Takahashi selected inverse (O(fill) memory).
//   Q~ is built sparse; SimplicialLLT factorises once; selinv_aux_local
//   (ported from src/selected_inv.cpp) gives S = Q~^{-1} on fill pattern;
//   d_i and Pr_i are gathered without ever forming Z.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;


// =========================================================================
// Takahashi selected-inverse — ported from src/selected_inv.cpp.
// Named selinv_aux_local to avoid ODR conflict with that file's
// anonymous-namespace version. Do NOT modify src/selected_inv.cpp.
// =========================================================================
static Eigen::SparseMatrix<double>
selinv_aux_local(const Eigen::SparseMatrix<double>& Lq) {
  typedef Eigen::SparseMatrix<double> SpMat;
  SpMat S = Lq.selfadjointView<Eigen::Lower>();
  int n = Lq.rows();
  for (int i = n - 1; i >= 0; --i) {
    SpMat::ReverseInnerIterator Si(S, i);
    for (SpMat::ReverseInnerIterator ij(Lq, i); ij; --ij) {
      SpMat::ReverseInnerIterator iL(Lq, i);
      SpMat::ReverseInnerIterator iS(S, ij.row());
      Si.valueRef() = 0.0;
      while (iL && iL.row() > i) {
        while (iS && (iL.row() < iS.row())) --iS;
        if (iS && (iL.row() == iS.row())) {
          Si.valueRef() -= iL.value() * iS.value();
          --iS;
        }
        --iL;
      }
      if (i == ij.row()) {
        Si.valueRef() += 1.0 / iL.value();
        Si.valueRef() /= iL.value();
      } else {
        Si.valueRef() /= iL.value();
        while (iS && iS.row() > i) --iS;
        iS.valueRef() = Si.value();
      }
      --Si;
    }
  }
  return S;
}

// Compute q_k^T S q_k  where  q_k = Tc[:,dof0]*sib0 + Tc[:,dof1]*sib1.
// For self-loops: dof1 is ignored, sib0 = SiB[k,0]+SiB[k,1], sib1 = 0.
// S.coeff(r,s) returns 0 for structural zeros (safe; design doc guarantees
// all needed (r,s) are in the fill pattern for standard river networks).
static double selinv_gather(
    const Eigen::SparseMatrix<double>& S,
    const Eigen::SparseMatrix<double>& Tc_sp,
    int dof0, int dof1,
    double sib0, double sib1, bool is_loop
) {
  typedef Eigen::SparseMatrix<double> SpMat;
  double T0ST0 = 0.0;
  for (SpMat::InnerIterator ir(Tc_sp, dof0); ir; ++ir) {
    double tc0r = ir.value();
    int r = ir.row();
    for (SpMat::InnerIterator ic(Tc_sp, dof0); ic; ++ic)
      T0ST0 += tc0r * S.coeff(r, ic.row()) * ic.value();
  }
  if (is_loop) return sib0 * sib0 * T0ST0;

  double T1ST1 = 0.0;
  for (SpMat::InnerIterator ir(Tc_sp, dof1); ir; ++ir) {
    double tc1r = ir.value();
    int r = ir.row();
    for (SpMat::InnerIterator ic(Tc_sp, dof1); ic; ++ic)
      T1ST1 += tc1r * S.coeff(r, ic.row()) * ic.value();
  }
  double T0ST1 = 0.0;
  for (SpMat::InnerIterator ir(Tc_sp, dof0); ir; ++ir) {
    double tc0r = ir.value();
    int r = ir.row();
    for (SpMat::InnerIterator ic(Tc_sp, dof1); ic; ++ic)
      T0ST1 += tc0r * S.coeff(r, ic.row()) * ic.value();
  }
  return sib0*sib0*T0ST0 + 2.0*sib0*sib1*T0ST1 + sib1*sib1*T1ST1;
}


// r_1 covariance matrix: S(i,j) = 1/(2*kappa*tau^2) * exp(-kappa * |D(i,j)|)
// Matches R's r_1(D, kappa, tau) exactly.
// tau = 1 / reciprocal_tau  (passed in natural scale)
static MatrixXd build_r1_cov(const NumericMatrix& D,
                              double kappa, double tau) {
  int m = D.nrow();
  double scale = 1.0 / (2.0 * kappa * tau * tau);
  MatrixXd S(m, m);
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < m; ++j)
      S(i, j) = scale * std::exp(-kappa * std::abs(D(i, j)));
  return S;
}
// =========================================================================
// cv_loo_selinv_cpp — LOO via Takahashi selected inverse.
//
// Replaces the dense Z matrix (nFree × n_obs, O(GB) for large networks)
// with the sparse selected inverse S = Q~^{-1} on the fill pattern.
// Memory: O(nnz(L)) ≈ O(fill) instead of O(nFree × n_obs).
//
// Algorithm:
//   1. Build Q~ = Tc (Q + B^T Σ^{-1} B) Tc^T  as SPARSE.
//   2. SimplicialLLT factorisation of Q~.
//   3. Takahashi recursion → S = selected_inv(Q~) on fill pattern.
//   4. For LOO obs i on edge e:
//        d_i  = [Σ^{-1}]_ii − q_i^T S q_i          (sparse gather)
//        Pr_i = sinv_r_adj_i − q_i^T v_g             (dot with w = Tc^T Q~^{-1} g)
//      where  q_i = Tc[:,dof0]*sib0 + Tc[:,dof1]*sib1,
//             g   = Tc Qpmu.
//   5. mu_i = y_i − Pr_i/d_i,  var_i = 1/d_i.
// =========================================================================
//' @noRd
// [[Rcpp::export]]
List cv_loo_selinv_cpp(
  List        precomputed_data,
  IntegerMatrix edge_endpoints,
  List        Q_list,
  double      sigma_e,
  double      reciprocal_tau,
  double      kappa
) {
  typedef Eigen::SparseMatrix<double> SpMat;

  // ---- 1. Extract Tc as sparse -------------------------------------------
  SEXP Tc_sexp = precomputed_data["Tc"];
  int nFree = 0, n_dof = 0;
  SpMat Tc_sp;
  if (Rf_isS4(Tc_sexp)) {
    S4 Tc_s4 = as<S4>(Tc_sexp);
    IntegerVector p_v = Tc_s4.slot("p");
    IntegerVector i_v = Tc_s4.slot("i");
    NumericVector x_v = Tc_s4.slot("x");
    IntegerVector d_v = Tc_s4.slot("Dim");
    nFree = d_v[0]; n_dof = d_v[1];
    std::vector<Eigen::Triplet<double>> t;
    t.reserve(x_v.size());
    for (int col = 0; col < n_dof; ++col)
      for (int idx = p_v[col]; idx < p_v[col+1]; ++idx)
        t.emplace_back(i_v[idx], col, x_v[idx]);
    Tc_sp.resize(nFree, n_dof);
    Tc_sp.setFromTriplets(t.begin(), t.end());
  } else {
    Function as_mat("as.matrix");
    NumericMatrix TcR = as<NumericMatrix>(as_mat(Tc_sexp));
    nFree = TcR.nrow(); n_dof = TcR.ncol();
    std::vector<Eigen::Triplet<double>> t;
    for (int col = 0; col < n_dof; ++col)
      for (int row = 0; row < nFree; ++row)
        if (std::abs(TcR(row, col)) > 1e-14)
          t.emplace_back(row, col, TcR(row, col));
    Tc_sp.resize(nFree, n_dof);
    Tc_sp.setFromTriplets(t.begin(), t.end());
  }
  Tc_sp.makeCompressed();

  // ---- 2. Scalars and precomputed arrays ----------------------------------
  int n_cov    = as<int>(precomputed_data["n_cov"]);
  CharacterVector u_repl  = as<CharacterVector>(precomputed_data["u_repl"]);
  IntegerVector obs_edges = precomputed_data["obs.edges"];
  IntegerMatrix E_mat     = edge_endpoints;
  int nE     = obs_edges.size();
  int n_repl = u_repl.size();
  double tau = 1.0 / reciprocal_tau;

  // ---- 3. Build sparse Q from Q_list triplets ----------------------------
  SpMat Q_sp(n_dof, n_dof);
  {
    IntegerVector Qi = Q_list["i"];
    IntegerVector Qj = Q_list["j"];
    NumericVector Qx = Q_list["x"];
    std::vector<Eigen::Triplet<double>> t;
    t.reserve(Qx.size());
    for (int k = 0; k < (int)Qx.size(); ++k)
      t.emplace_back(Qi[k]-1, Qj[k]-1, Qx[k]);
    Q_sp.setFromTriplets(t.begin(), t.end());
  }

  // ---- 4. Per-obs cache struct -------------------------------------------
  struct ObsEntry {
    double y, diag_sinv, sinv_r, sib0, sib1;
    int edge;      // 1-based
    int local;     // 0-based within edge
    int repl_idx;  // 0-based replicate index
    bool is_loop;
    std::vector<double> sinv_x;  // n_cov entries (always allocated)
  };
  std::vector<ObsEntry> obs_cache;
  obs_cache.reserve(1024);

  // ---- 5. Global H/h accumulators ----------------------------------------
  MatrixXd H_mat = MatrixXd::Zero(n_cov, n_cov);
  VectorXd h_vec = VectorXd::Zero(n_cov);

  // ---- 6. Per-replicate Qpmu, QpmuX, XtSX, XtSy -------------------------
  std::vector<VectorXd> Qpmu_store(n_repl, VectorXd::Zero(n_dof));
  std::vector<MatrixXd> QpmuX_store(n_repl,
                          MatrixXd::Zero(n_dof, std::max(n_cov, 1)));
  std::vector<MatrixXd> XtSX_store(n_repl,
                          MatrixXd::Zero(std::max(n_cov,1), std::max(n_cov,1)));
  std::vector<VectorXd> XtSy_store(n_repl, VectorXd::Zero(std::max(n_cov,1)));

  // BtSinvB sparse triplets — built from first replicate only
  std::vector<Eigen::Triplet<double>> btsb_trips;

  List y_data   = precomputed_data["y_data"];
  List D_data   = precomputed_data["D_data"];
  List x_data_list;
  if (n_cov > 0) x_data_list = as<List>(precomputed_data["x_data"]);

  // =========================================================================
  // PHASE 1 — replicate × edge loop
  //   First replicate: also build BtSinvB triplets and collect obs structure.
  //   All replicates: accumulate Qpmu_r, XtSX_r, XtSy_r, QpmuX_r and store
  //   per-obs (y, sinv_r, sinv_x).
  // =========================================================================
  bool btsb_done = false;  // BtSinvB built after first replicate's edges

  for (int repl_y = 0; repl_y < n_repl; ++repl_y) {
    std::string curr_repl = as<std::string>(u_repl[repl_y]);
    std::string repl_name = "repl_" + curr_repl;

    if (!y_data.containsElementNamed(repl_name.c_str())) continue;
    List y_repl = as<List>(y_data[repl_name]);
    List D_repl = as<List>(D_data[repl_name]);
    List x_repl;
    if (n_cov > 0) x_repl = as<List>(x_data_list[repl_name]);

    VectorXd& Qpmu_r = Qpmu_store[repl_y];

    for (int j = 0; j < nE; ++j) {
      int e = obs_edges[j];
      std::string ename = "edge_" + std::to_string(e);
      if (!y_repl.containsElementNamed(ename.c_str())) continue;
      SEXP y_sexp = y_repl[ename];
      if (Rf_isNull(y_sexp)) continue;
      NumericVector y_i_r = as<NumericVector>(y_sexp);
      int n_i = y_i_r.size();
      if (n_i == 0) continue;

      NumericMatrix D_r = as<NumericMatrix>(D_repl[ename]);
      MatrixXd Scov = build_r1_cov(D_r, kappa, tau);

      MatrixXd S_EE = Scov.topLeftCorner(2, 2);
      MatrixXd S_EO = Scov.topRightCorner(2, n_i);
      MatrixXd S_OE = Scov.bottomLeftCorner(n_i, 2);
      MatrixXd S_OO = Scov.bottomRightCorner(n_i, n_i);

      Eigen::LLT<MatrixXd> llt_see(S_EE);
      MatrixXd Bt = llt_see.solve(S_EO);

      MatrixXd Sigma_i = S_OO - S_OE * Bt;
      for (int d = 0; d < n_i; ++d) Sigma_i(d, d) += sigma_e * sigma_e;

      Eigen::LLT<MatrixXd> llt_sig(Sigma_i);
      if (llt_sig.info() != Eigen::Success)
        stop("cv_loo_selinv_cpp: Cholesky of Sigma_i failed (edge " +
             std::to_string(e) + ")");

      MatrixXd Sinv_e   = llt_sig.solve(MatrixXd::Identity(n_i, n_i));
      MatrixXd Sigma_iB = llt_sig.solve(Bt.transpose());  // n_i x 2

      VectorXd y_eigen(n_i);
      for (int ii = 0; ii < n_i; ++ii) y_eigen[ii] = y_i_r[ii];
      VectorXd v_i = llt_sig.solve(y_eigen);

      int e0   = e - 1;
      int dof0 = 2 * e0;
      int dof1 = 2 * e0 + 1;
      bool is_loop = (E_mat(e0, 0) == E_mat(e0, 1));

      // BtSinvB and Qpmu accumulation
      MatrixXd BtSinvB_e = Bt * Sigma_iB;
      if (!btsb_done) {
        // BtSinvB is the same for every replicate (same obs locations);
        // build once from first replicate.
        if (is_loop) {
          btsb_trips.emplace_back(dof0, dof0, BtSinvB_e.sum());
        } else {
          btsb_trips.emplace_back(dof0, dof0, BtSinvB_e(0, 0));
          btsb_trips.emplace_back(dof0, dof1, BtSinvB_e(0, 1));
          btsb_trips.emplace_back(dof1, dof0, BtSinvB_e(1, 0));
          btsb_trips.emplace_back(dof1, dof1, BtSinvB_e(1, 1));
        }
      }
      if (is_loop) {
        Qpmu_r[dof0] += (Sigma_iB.transpose() * y_eigen).sum();
      } else {
        VectorXd tmp = Sigma_iB.transpose() * y_eigen;
        Qpmu_r[dof0] += tmp[0];
        Qpmu_r[dof1] += tmp[1];
      }

      // Covariate accumulators
      MatrixXd X_i;
      MatrixXd SinvX_i;
      if (n_cov > 0) {
        NumericMatrix X_r = as<NumericMatrix>(x_repl[ename]);
        X_i.resize(n_i, n_cov);
        for (int rr = 0; rr < n_i; ++rr)
          for (int cc = 0; cc < n_cov; ++cc)
            X_i(rr, cc) = X_r(rr, cc);
        SinvX_i = llt_sig.solve(X_i);
        XtSX_store[repl_y] += X_i.transpose() * SinvX_i;
        XtSy_store[repl_y] += X_i.transpose() * v_i;
        MatrixXd dQX = Sigma_iB.transpose() * X_i;
        if (is_loop) {
          QpmuX_store[repl_y].row(dof0) += dQX.row(0) + dQX.row(1);
        } else {
          QpmuX_store[repl_y].row(dof0) += dQX.row(0);
          QpmuX_store[repl_y].row(dof1) += dQX.row(1);
        }
      }

      // Per-obs cache entries
      VectorXd sinv_diag = Sinv_e.diagonal();
      for (int k = 0; k < n_i; ++k) {
        ObsEntry oe;
        oe.y        = y_i_r[k];
        oe.diag_sinv = sinv_diag[k];
        oe.sinv_r   = v_i[k];
        oe.edge     = e;
        oe.local    = k;
        oe.repl_idx = repl_y;
        oe.is_loop  = is_loop;
        if (is_loop) {
          oe.sib0 = Sigma_iB(k, 0) + Sigma_iB(k, 1);
          oe.sib1 = 0.0;
        } else {
          oe.sib0 = Sigma_iB(k, 0);
          oe.sib1 = Sigma_iB(k, 1);
        }
        if (n_cov > 0) {
          oe.sinv_x.resize(n_cov);
          for (int cc = 0; cc < n_cov; ++cc) oe.sinv_x[cc] = SinvX_i(k, cc);
        }
        obs_cache.push_back(std::move(oe));
      }
    }  // end edge loop

    if (!btsb_done) btsb_done = true;
  }  // end replicate loop

  int n_cache = (int)obs_cache.size();

  // =========================================================================
  // PHASE 2 — build sparse Q~, factorize, Takahashi selinv
  // =========================================================================
  SpMat BtSinvB_sp(n_dof, n_dof);
  BtSinvB_sp.setFromTriplets(btsb_trips.begin(), btsb_trips.end());

  SpMat Qp_sp = Q_sp + BtSinvB_sp;

  // Q~ = Tc (Q + BtSinvB) Tc^T  (sparse triple product)
  SpMat Qtilde_sp = Tc_sp * Qp_sp * Tc_sp.transpose();
  // Symmetrize for numerical safety
  SpMat Qt_sym = SpMat(Qtilde_sp.transpose());
  Qtilde_sp = 0.5 * Qtilde_sp + 0.5 * Qt_sym;
  Qtilde_sp.makeCompressed();

  Eigen::SimplicialLLT<SpMat> solver;
  solver.analyzePattern(Qtilde_sp);
  solver.factorize(Qtilde_sp);
  if (solver.info() != Eigen::Success)
    stop("cv_loo_selinv_cpp: SimplicialLLT of Q_tilde failed");

  SpMat L_factor = solver.matrixL();
  SpMat S_perm   = selinv_aux_local(L_factor);
  // Permute back: S_orig = P^{-1} S_perm (P^{-1})^T
  SpMat S = solver.permutationPinv() * S_perm *
            solver.permutationPinv().transpose();
  S.makeCompressed();

  // =========================================================================
  // PHASE 3 — per-replicate solves: v_g_r, w_r, VX_r, wX_r, H/h
  // =========================================================================
  std::vector<VectorXd> w_store(n_repl, VectorXd::Zero(n_dof));
  std::vector<MatrixXd> wX_store(n_repl,
                          MatrixXd::Zero(n_dof, std::max(n_cov, 1)));

  for (int repl_y = 0; repl_y < n_repl; ++repl_y) {
    VectorXd TcQpmu_r = Tc_sp * Qpmu_store[repl_y];
    VectorXd v_g_r    = solver.solve(TcQpmu_r);
    w_store[repl_y]   = Tc_sp.transpose() * v_g_r;

    if (n_cov > 0) {
      MatrixXd TcQpmuX_r = Tc_sp * QpmuX_store[repl_y];  // nFree x n_cov
      MatrixXd VX_r      = solver.solve(TcQpmuX_r);       // nFree x n_cov
      // wX_r = Tc^T VX_r  (n_dof x n_cov), used for Pr and REML adjustments
      MatrixXd wX_r(n_dof, n_cov);
      for (int cc = 0; cc < n_cov; ++cc)
        wX_r.col(cc) = Tc_sp.transpose() * VX_r.col(cc);
      wX_store[repl_y] = wX_r;
      // H = XtSX - QpmuX^T Tc^T Q~^{-1} Tc QpmuX = XtSX - TcQpmuX^T VX_full
      // h = XtSy - QpmuX^T Tc^T Q~^{-1} Tc Qpmu   = XtSy - TcQpmuX^T v_g
      // NOT VX_r^T VX_r / VX_r^T v_g (those would be Q~^{-2} based).
      H_mat += XtSX_store[repl_y] - TcQpmuX_r.transpose() * VX_r;
      h_vec += XtSy_store[repl_y] - TcQpmuX_r.transpose() * v_g_r;
    }
  }

  // =========================================================================
  // PHASE 4 — beta_hat = H^{-1} h
  // =========================================================================
  VectorXd beta_hat_vec = VectorXd::Zero(n_cov);
  if (n_cov > 0) {
    Eigen::LLT<MatrixXd> llt_H(H_mat);
    if (llt_H.info() != Eigen::Success)
      stop("cv_loo_selinv_cpp: LLT of H failed");
    beta_hat_vec = llt_H.solve(h_vec);
  }

  // =========================================================================
  // PHASE 5 — LOO predictions from obs_cache
  // =========================================================================
  NumericVector mu_out(n_cache, NA_REAL);
  NumericVector var_out(n_cache, NA_REAL);

  for (int gi = 0; gi < n_cache; ++gi) {
    const ObsEntry& oe = obs_cache[gi];
    int e    = oe.edge;
    int e0   = e - 1;
    int dof0 = 2 * e0;
    int dof1 = 2 * e0 + 1;
    double sib0 = oe.sib0, sib1 = oe.sib1;

    // d_i = diag_sinv_i - q_i^T S q_i  (Takahashi gather)
    double d_i = oe.diag_sinv -
                 selinv_gather(S, Tc_sp, dof0, dof1, sib0, sib1, oe.is_loop);

    // sinv_r_adj = (Sigma^{-1} y)_i - (Sigma^{-1} X)_i beta_hat
    double sinv_r_adj = oe.sinv_r;
    if (n_cov > 0)
      for (int cc = 0; cc < n_cov; ++cc)
        sinv_r_adj -= oe.sinv_x[cc] * beta_hat_vec[cc];

    // Pr_i = sinv_r_adj - q_i^T Tc^T Q~^{-1} Tc (Qpmu - QpmuX beta_hat)
    // w_r = Tc^T Q~^{-1} Tc Qpmu  (y-based);  w_r_adj subtracts beta contribution.
    const VectorXd& w_r = w_store[oe.repl_idx];
    double w0 = w_r[dof0], w1 = w_r[dof1];
    if (n_cov > 0) {
      const MatrixXd& wX_loc = wX_store[oe.repl_idx];
      for (int cc = 0; cc < n_cov; ++cc) {
        double bcc = beta_hat_vec[cc];
        w0 -= wX_loc(dof0, cc) * bcc;
        if (!oe.is_loop) w1 -= wX_loc(dof1, cc) * bcc;
      }
    }
    double Pr_i = sinv_r_adj;
    Pr_i -= sib0 * w0;
    if (!oe.is_loop) Pr_i -= sib1 * w1;

    mu_out[gi]  = oe.y - Pr_i / d_i;
    var_out[gi] = 1.0 / d_i;
  }

  // ---- Return value -------------------------------------------------------
  NumericVector beta_hat_ret(n_cov);
  for (int cc = 0; cc < n_cov; ++cc) beta_hat_ret[cc] = beta_hat_vec[cc];

  NumericMatrix H_ret(std::max(n_cov,1), std::max(n_cov,1));
  H_ret.fill(0.0);
  for (int r = 0; r < n_cov; ++r)
    for (int c = 0; c < n_cov; ++c)
      H_ret(r, c) = H_mat(r, c);

  return List::create(
    Named("mu")       = mu_out,
    Named("var")      = var_out,
    Named("beta_hat") = beta_hat_ret,
    Named("H")        = H_ret
  );
}
