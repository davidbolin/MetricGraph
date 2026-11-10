#include <Rcpp.h>
#include <RcppEigen.h>
#include <map>
#include <cmath>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

// ============================================================
// Kernel scalars
// ============================================================

// r_1(d) = (1/(2*kappa*tau^2)) * exp(-kappa*|d|)
static inline double r1_s(double d, double kappa, double inv2ktau2) {
    return inv2ktau2 * std::exp(-kappa * std::abs(d));
}

// r_2 and its derivatives
//   c  = 1 / (4*kappa^3*tau^2)
//   deriv=0: c*(1+kappa*|d|)*exp(-kappa*|d|)
//   deriv=1: -kappa^2*c*d*exp(-kappa*|d|)
//   deriv=2: kappa^2*c*(kappa*|d|-1)*exp(-kappa*|d|)
static inline double r2_s(double d, double kappa, double c_r2, int deriv) {
    double ad = std::abs(d);
    double R0 = std::exp(-kappa * ad);
    if (deriv == 0) return c_r2 * (1.0 + kappa * ad) * R0;
    if (deriv == 1) return -kappa * kappa * c_r2 * d * R0;
    return kappa * kappa * c_r2 * (kappa * ad - 1.0) * R0;   // deriv == 2
}


// ============================================================
// Bridge predictors: S_e (m x b) and Sigma*_e (m x m)
// ============================================================

// alpha = 1: S only — O(m), used by kriging (does not need Sigma*_e)
static Eigen::MatrixXd build_S1(const Eigen::VectorXd& t_abs, double l_e,
                                  double kappa, double tau)
{
    int m = t_abs.size();
    double c = 1.0 / (2.0 * kappa * tau * tau);

    Eigen::Matrix2d Sbb;
    Sbb(0,0) = r1_s(0.0,  kappa, c);
    Sbb(0,1) = r1_s(-l_e, kappa, c);
    Sbb(1,0) = Sbb(0,1);
    Sbb(1,1) = Sbb(0,0);

    Eigen::MatrixXd Sxb(m, 2);
    for (int i = 0; i < m; i++) {
        Sxb(i,0) = r1_s(t_abs[i],       kappa, c);
        Sxb(i,1) = r1_s(t_abs[i] - l_e, kappa, c);
    }
    return Sxb * Sbb.inverse();
}

// alpha = 1: S and Sigma*_e — O(m^2), used by direct
static void build_bridge1(const Eigen::VectorXd& t_abs, double l_e,
                           double kappa, double tau,
                           Eigen::MatrixXd& S, Eigen::MatrixXd& Ss)
{
    int m = t_abs.size();
    double c = 1.0 / (2.0 * kappa * tau * tau);

    Eigen::Matrix2d Sbb;
    Sbb(0,0) = r1_s(0.0,   kappa, c);
    Sbb(0,1) = r1_s(-l_e,  kappa, c);
    Sbb(1,0) = Sbb(0,1);
    Sbb(1,1) = Sbb(0,0);

    Eigen::MatrixXd Sxb(m, 2);
    for (int i = 0; i < m; i++) {
        Sxb(i,0) = r1_s(t_abs[i],        kappa, c);
        Sxb(i,1) = r1_s(t_abs[i] - l_e,  kappa, c);
    }

    Eigen::MatrixXd Sxx(m, m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            Sxx(i,j) = r1_s(t_abs[i] - t_abs[j], kappa, c);

    S  = Sxb * Sbb.inverse();
    Eigen::MatrixXd tmp1 = Sxx - S * Sxb.transpose();
    Ss = 0.5 * (tmp1 + tmp1.transpose());
}

// alpha = 2: S only — O(m), used by kriging
static Eigen::MatrixXd build_S2(const Eigen::VectorXd& t_abs, double l_e,
                                  double kappa, double tau)
{
    int m = t_abs.size();
    double c = 1.0 / (4.0 * kappa * kappa * kappa * tau * tau);

    Eigen::Matrix4d Sbb;
    Sbb.setZero();
    Sbb(0,0) = r2_s(0.0,  kappa, c, 0);
    Sbb(0,2) = r2_s(-l_e, kappa, c, 0);
    Sbb(2,0) = Sbb(0,2);
    Sbb(2,2) = Sbb(0,0);
    Sbb(1,1) = -r2_s(0.0, kappa, c, 2);
    Sbb(1,3) = -r2_s(l_e, kappa, c, 2);
    Sbb(3,1) =  Sbb(1,3);
    Sbb(3,3) =  Sbb(1,1);
    Sbb(1,0) = r2_s(0.0 - 0.0, kappa, c, 1);
    Sbb(1,2) = r2_s(0.0 - l_e, kappa, c, 1);
    Sbb(3,0) = r2_s(l_e - 0.0, kappa, c, 1);
    Sbb(3,2) = r2_s(l_e - l_e, kappa, c, 1);
    Sbb(0,1) = Sbb(1,0); Sbb(2,1) = Sbb(1,2);
    Sbb(0,3) = Sbb(3,0); Sbb(2,3) = Sbb(3,2);

    Eigen::MatrixXd Sxb(m, 4);
    for (int i = 0; i < m; i++) {
        double d0  = t_abs[i];
        double dle = t_abs[i] - l_e;
        Sxb(i,0) =  r2_s( d0,  kappa, c, 0);
        Sxb(i,2) =  r2_s( dle, kappa, c, 0);
        Sxb(i,1) =  r2_s(-d0,  kappa, c, 1);
        Sxb(i,3) =  r2_s(-dle, kappa, c, 1);
    }
    return Sxb * Sbb.inverse();
}

// alpha = 2: S and Sigma*_e — O(m^2), used by direct
static void build_bridge2(const Eigen::VectorXd& t_abs, double l_e,
                           double kappa, double tau,
                           Eigen::MatrixXd& S, Eigen::MatrixXd& Ss)
{
    int m = t_abs.size();
    double c = 1.0 / (4.0 * kappa * kappa * kappa * tau * tau);

    // Sigma_bb (4x4)
    // Layout: rows/cols = [u(0), u'(0), u(l_e), u'(l_e)]  (0-indexed: 0,1,2,3)
    // R code indices c(1,3) → 0-based {0,2}; c(2,4) → 0-based {1,3}
    Eigen::Matrix4d Sbb;
    Sbb.setZero();

    // {0,2} x {0,2}: r_2(D_vv, ., 0), D_vv=outer({0,l_e},{0,l_e},-)
    Sbb(0,0) = r2_s(0.0,   kappa, c, 0);
    Sbb(0,2) = r2_s(-l_e,  kappa, c, 0);   // r_2 uses |d|, same as +l_e
    Sbb(2,0) = Sbb(0,2);
    Sbb(2,2) = Sbb(0,0);

    // {1,3} x {1,3}: -r_2(D_dd, ., 2), D_dd = absolute distances
    Sbb(1,1) = -r2_s(0.0,  kappa, c, 2);
    Sbb(1,3) = -r2_s(l_e,  kappa, c, 2);
    Sbb(3,1) =  Sbb(1,3);
    Sbb(3,3) =  Sbb(1,1);

    // {1,3} x {0,2}: r_2(D_dv, ., 1), D_dv=outer({0,l_e},{0,l_e},-)
    Sbb(1,0) = r2_s(0.0 - 0.0, kappa, c, 1);   // = 0
    Sbb(1,2) = r2_s(0.0 - l_e, kappa, c, 1);
    Sbb(3,0) = r2_s(l_e - 0.0, kappa, c, 1);
    Sbb(3,2) = r2_s(l_e - l_e, kappa, c, 1);   // = 0

    // {0,2} x {1,3}: transpose of {1,3} x {0,2}
    Sbb(0,1) = Sbb(1,0);
    Sbb(2,1) = Sbb(1,2);
    Sbb(0,3) = Sbb(3,0);
    Sbb(2,3) = Sbb(3,2);

    // Sigma_xb (m x 4)
    // D_xv = outer(t_abs, {0,l_e}, -)
    // cols {0,2}: r_2(D_xv, ., 0)
    // cols {1,3}: r_2(-D_xv, ., 1)
    Eigen::MatrixXd Sxb(m, 4);
    for (int i = 0; i < m; i++) {
        double d0  = t_abs[i];
        double dle = t_abs[i] - l_e;
        Sxb(i,0) =  r2_s( d0,  kappa, c, 0);
        Sxb(i,2) =  r2_s( dle, kappa, c, 0);
        Sxb(i,1) =  r2_s(-d0,  kappa, c, 1);
        Sxb(i,3) =  r2_s(-dle, kappa, c, 1);
    }

    // Sigma_xx (m x m)
    Eigen::MatrixXd Sxx(m, m);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            Sxx(i,j) = r2_s(t_abs[i] - t_abs[j], kappa, c, 0);

    S  = Sxb * Sbb.inverse();
    Eigen::MatrixXd tmp2 = Sxx - S * Sxb.transpose();
    Ss = 0.5 * (tmp2 + tmp2.transpose());
}


// ============================================================
// draw_edge_direct_cpp  (Method A — O(m^3) per edge)
// ============================================================

//' @name draw_edge_direct_cpp
//' @title C++ direct bridge draw (Method A)
//' @description
//' Same algorithm as the R reference `draw_edge_direct`.
//' Builds S_e, Sigma*_e via the closed-form kernel, then draws
//' S_e b_e + chol(Sigma*_e)^T z.  Uses Eigen::LLT and R::norm_rand().
//' @param kappa,tau SPDE parameters.
//' @param b_e Boundary state: length-2 (alpha=1) or length-4 (alpha=2).
//' @param l_e Edge length.
//' @param t_abs Interior locations in absolute coordinates (sorted).
//' @param alpha Smoothness: 1 or 2.
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd draw_edge_direct_cpp(
        double kappa, double tau,
        const Eigen::VectorXd& b_e,
        double l_e,
        const Eigen::VectorXd& t_abs,
        int alpha)
{
    int m = t_abs.size();
    if (m == 0) return Eigen::VectorXd(0);

    Eigen::MatrixXd S, Ss;
    if (alpha == 1) build_bridge1(t_abs, l_e, kappa, tau, S, Ss);
    else            build_bridge2(t_abs, l_e, kappa, tau, S, Ss);

    Eigen::VectorXd mu = S * b_e;

    // Sigma*_e = L * L^T  (LLT gives lower Cholesky L)
    // Matches R: mu + t(chol(Sigma*_e)) %*% z  since t(upper_R) == lower_L
    Eigen::LLT<Eigen::MatrixXd> llt(Ss);

    Eigen::VectorXd z(m);
    for (int i = 0; i < m; i++) z[i] = R::norm_rand();

    if (llt.info() == Eigen::Success) {
        return mu + llt.matrixL() * z;
    }

    // Cholesky failed: Sigma* is near-singular due to floating-point error.
    // Fall back to eigendecomposition, clamping small negative eigenvalues to zero.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Ss);
    if (es.info() != Eigen::Success) {
        Rcpp::stop("Eigen decomposition failed for bridge covariance");
    }
    Eigen::VectorXd vals = es.eigenvalues();
    double max_eval = vals.cwiseAbs().maxCoeff();
    double tol = std::max(1e-14, 1e-6 * max_eval);
    if (vals.minCoeff() < -tol) {
        Rcpp::warning("Bridge covariance has materially negative eigenvalues");
    }
    for (int i = 0; i < vals.size(); ++i) {
        vals[i] = vals[i] > 0.0 ? std::sqrt(vals[i]) : 0.0;
    }
    return mu + es.eigenvectors() * vals.asDiagonal() * z;
}


// ============================================================
// Markov chain helpers (Method B internals)
// ============================================================

// Alpha = 1: OU process
// Draws n normals: 1 for initial state + (n-1) for transitions.
// t_aug must be sorted (0 first, l_e last).
static Eigen::VectorXd markov1(const std::vector<double>& t_aug,
                                double kappa, double tau)
{
    int n     = (int)t_aug.size();
    double r0 = 1.0 / (2.0 * kappa * tau * tau);

    Eigen::VectorXd x(n);
    x[0] = R::norm_rand() * std::sqrt(r0);
    for (int i = 0; i < n - 1; i++) {
        double h   = t_aug[i+1] - t_aug[i];
        double phi = std::exp(-kappa * h);
        double sd  = std::sqrt(r0 * (1.0 - phi * phi));
        x[i+1]     = phi * x[i] + sd * R::norm_rand();
    }
    return x;
}

// Per-lag cache entry for the 2-D Markov recursion (alpha=2)
struct MCEntry {
    Eigen::Matrix2d A;
    Eigen::Matrix2d Lom;   // lower Cholesky of innovation covariance Omega
};

// Alpha = 2: 2-D Matérn-3/2 process
// State: [u, u']
// Draws 2*(n) normals: 2 for initial state + 2*(n-1) for transitions.
// Unique-lag cache: builds (A, chol(Omega)) once per distinct h.
static Eigen::MatrixXd markov2(const std::vector<double>& t_aug,
                                double kappa, double tau)
{
    int n       = (int)t_aug.size();
    double kap2 = kappa * kappa;
    double c    = 1.0 / (4.0 * kap2 * kappa * tau * tau);

    double c_val = r2_s(0.0, kappa, c, 0);
    double r0_du = -r2_s(0.0, kappa, c, 2);

    // Build unique-lag cache
    std::map<double, MCEntry> cache;
    for (int i = 0; i < n - 1; i++) {
        double h = t_aug[i+1] - t_aug[i];
        if (cache.count(h)) continue;

        double phi = std::exp(-kappa * h);
        double kh  = kappa * h;

        Eigen::Matrix2d A;
        A(0,0) =  phi * (1.0 + kh);
        A(0,1) =  phi * h;
        A(1,0) = -phi * kap2 * h;
        A(1,1) =  phi * (1.0 - kh);

        // Cross-covariance Ch:  [[r2(h,0), -r2(h,1)], [r2(h,1), -r2(h,2)]]
        Eigen::Matrix2d Ch;
        Ch(0,0) =  r2_s(h, kappa, c, 0);
        Ch(0,1) = -r2_s(h, kappa, c, 1);
        Ch(1,0) =  r2_s(h, kappa, c, 1);
        Ch(1,1) = -r2_s(h, kappa, c, 2);

        Eigen::Matrix2d R0d;
        R0d.setZero();
        R0d(0,0) = c_val;
        R0d(1,1) = r0_du;

        // Innovation covariance: Omega = R0_diag - A * Ch^T
        Eigen::Matrix2d Omega = R0d - A * Ch.transpose();
        Omega = 0.5 * (Omega + Omega.transpose());

        Eigen::LLT<Eigen::Matrix2d> llt(Omega);
        MCEntry mc;
        mc.A   = A;
        mc.Lom = llt.matrixL();
        cache[h] = mc;
    }

    // Simulate
    Eigen::MatrixXd X(2, n);
    // Initial draw: [sqrt(c_val)*z1, sqrt(r0_du)*z2]  (same order as R's rnorm(2))
    X(0,0) = std::sqrt(c_val)  * R::norm_rand();
    X(1,0) = std::sqrt(r0_du)  * R::norm_rand();

    for (int i = 0; i < n - 1; i++) {
        double h        = t_aug[i+1] - t_aug[i];
        const MCEntry& mc = cache.at(h);
        Eigen::Vector2d z;
        z[0] = R::norm_rand();
        z[1] = R::norm_rand();
        X.col(i+1) = mc.A * X.col(i) + mc.Lom * z;
    }
    return X;
}


// ============================================================
// draw_edge_kriging_cpp  (Method B — O(m) per edge)
// ============================================================

//' @name draw_edge_kriging_cpp
//' @title C++ kriging-corrected Markov draw (Method B)
//' @description
//' Same algorithm as the R reference `draw_edge_kriging`.
//' Simulates a full Markov path on the edge, then applies the kriging
//' correction to enforce the vertex boundary conditions.
//' Uses a unique-lag cache for (A, chol(Omega)) and R::norm_rand().
//' @param kappa,tau SPDE parameters.
//' @param b_e Boundary state: length-2 (alpha=1) or length-4 (alpha=2).
//' @param l_e Edge length.
//' @param t_abs Interior locations in absolute coordinates (must be sorted).
//' @param alpha Smoothness: 1 or 2.
//' @noRd
// [[Rcpp::export]]
Eigen::VectorXd draw_edge_kriging_cpp(
        double kappa, double tau,
        const Eigen::VectorXd& b_e,
        double l_e,
        const Eigen::VectorXd& t_abs,
        int alpha)
{
    int m = t_abs.size();
    if (m == 0) return Eigen::VectorXd(0);

    // Build augmented grid: 0, t_abs[0..m-1], l_e  (t_abs assumed sorted)
    std::vector<double> t_aug(m + 2);
    t_aug[0] = 0.0;
    for (int i = 0; i < m; i++) t_aug[i+1] = t_abs[i];
    t_aug[m+1] = l_e;

    if (alpha == 1) {
        Eigen::VectorXd x_aug = markov1(t_aug, kappa, tau);

        Eigen::VectorXd x_int(m);
        for (int i = 0; i < m; i++) x_int[i] = x_aug[i+1];

        Eigen::Vector2d x_ends(x_aug[0], x_aug[m+1]);
        Eigen::MatrixXd S = build_S1(t_abs, l_e, kappa, tau);

        return x_int + S * (b_e - x_ends);

    } else {
        Eigen::MatrixXd X_aug = markov2(t_aug, kappa, tau);

        Eigen::VectorXd x_int(m);
        for (int i = 0; i < m; i++) x_int[i] = X_aug(0, i+1);

        Eigen::Vector4d x_ends_state;
        x_ends_state[0] = X_aug(0, 0);
        x_ends_state[1] = X_aug(1, 0);
        x_ends_state[2] = X_aug(0, m+1);
        x_ends_state[3] = X_aug(1, m+1);

        Eigen::MatrixXd S = build_S2(t_abs, l_e, kappa, tau);

        return x_int + S * (b_e - x_ends_state);
    }
}
