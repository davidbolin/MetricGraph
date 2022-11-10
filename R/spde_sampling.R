#' Sample Whittle-Matérn field on graph
#' @details Samples a Gaussian Whittle-Matérn field on a graph, either from the prior
#' or conditionally on observations
#' \deqn{y_i = u(t_i) + sigma_e e_i}{y_i = u(t_i) + sigma_e e_i}
#' in the graph,  where \eqn{e_i} are independent standard Gaussian variables.
#' @param kappa parameter kappa
#' @param theta parameter theta
#' @param sigma_e parameter sigma_e
#' @param posterior sample conditionally on the observations?
#' @param u sample evaluated at the mesh nodes in the graph
#' @export
sample_spde_mesh <- function(kappa, sigma, sigma_e = 0, alpha = 1, graph,
                             posterior = FALSE) {

  if (!posterior) {
    if (alpha == 1) {
      Q <- spde_precision(kappa = kappa, sigma = sigma, alpha = 1, graph = graph)
      R <- Cholesky(Q,LDL = FALSE, perm = TRUE)
      V0 <- as.vector(solve(R, solve(R,rnorm(graph$nV), system = 'Lt'), system = 'Pt'))
      u <- NULL

      inds_PtE <- sort(unique(graph$mesh$PtE[,1])) #inds
      for (i in inds_PtE) {
        t <- graph$mesh$PtE[graph$mesh$PtE[,1] == i,2]
        samp <- sample_alpha1_line(kappa = kappa, sigma = sigma,
                                   u_e = V0[graph$E[i, ]], t = t,
                                   l_e = graph$edge_lengths[i])
        u <- c(u, samp[,2])
      }
    } else if (alpha == 2) {
      Q <- spde_precision(kappa = kappa, sigma = sigma,
                          alpha = 2, graph = graph, BC = 1)
      if(is.null(graph$C))
        graph$buildC(2)

      Qmod <- (graph$CBobj$T) %*% Q %*% t(graph$CBobj$T)
      Qtilde <- Qmod[-c(1:dim(graph$CBobj$U)[1]),-c(1:dim(graph$CBobj$U)[1])]
      R <- Cholesky(forceSymmetric(Qtilde),LDL = FALSE, perm = TRUE)
      V0 <- as.vector(solve(R, solve(R,rnorm(4*graph$nV - dim(graph$CBobj$U)[1]),
                                     system = 'Lt')
                            , system = 'Pt'))
      u_e <- t(graph$CBobj$T) %*% c(rep(0, dim(graph$CBobj$U)[1]), V0)
      u <- NULL
      inds_PtE <- sort(unique(graph$mesh$PtE[,1])) #inds
      for (i in inds_PtE) {
        t <- graph$mesh$PtE[graph$mesh$PtE[,1] == i,2]
        samp <- sample_alpha2_line(kappa = kappa, sigma = sigma,
                                   sigma_e = sigma_e,
                                   u_e = u_e[4*(i-1) +1:4],
                                   t = t,
                                   l_e = graph$edge_lengths[i])
        u <- c(u, samp[,2])
      }
    } else {
      stop("only alpha = 1 and alpha = 2 implemented.")
    }
  } else {
    stop("TODO: implement posterior sampling")
  }
  return(u)
}
#' Sample Gaussian process with exponential covariance on a line given end points
#' @details Samples a Gaussian process \eqn{u(t)} with an exponential covariance function
#' \deqn{r(h) = \sigma^2\exp(-\kappa h)/(2\kappa)}{r(h) = sigma^2*(exp(-kappa*h)/(2*kappa)}
#' on an interval \eqn{(0,l_e)} conditionally on \eqn{u(0), u(l_e)}.
#' If `y` and `py` are supplied, the sampling is done conditionally on observations
#' \deqn{y_i = u(t_i) + sigma_e e_i}{y_i = u(t_i) + sigma_e e_i}
#' where \eqn{e_i} are independent standard Gaussian variables.
#' @param kappa parameter kappa
#' @param theta parameter theta
#' @param sigma_e parameter sigma_e
#' @param u_e  (2 x 1) the two end points
#' @param l_e (1 x 1) line length
#' @param  t (n x 1) distance on the line to be sampled from (not end points)
#' @param  nt (1 x 1) number of equidistance points to sample from if t is  null
#' @param  py  (n x 1) observation locations
#' @param  y (n x 1) observations
#' @param  sample (bool) if true sample else return posterior mean
#' @export
sample_alpha1_line <- function(kappa, sigma, sigma_e,
                               u_e, l_e, t = NULL,
                               nt = 100,  py = NULL,
                               y = NULL, sample = TRUE) {

  if (is.null(t)) {
    t  <- seq(0, 1, length.out = nt)
    t <- t * l_e
  }
  t_end <- c(0, l_e)
  t <- unique(t)
  t0 <- t
  if (is.null(y) == FALSE) {
    ind_remove <- which(py %in% t_end)
    if (length(ind_remove) > 0) {
      y  <- y[-ind_remove]
      py <- py[-ind_remove]
    }

    ind_remove_t <- which(t %in% c(t_end, py))
    if (length(ind_remove_t) > 0)
      t <- t[-ind_remove_t]
  } else {
    ind_remove_t <- which(t %in% t_end)
    if(length(ind_remove_t) > 0)
      t <- t[-ind_remove_t]
  }

  t <- c(t_end, t)
  if (is.null(py) == FALSE) {
    t <- c(py, t)
  }

  Q <- precision_exp_line(kappa = kappa, sigma = sigma, t = t)

  index_E <- length(py) + 1:2
  Q_X <- Q[-index_E,-index_E, drop=F]
  mu_X <- as.vector(Matrix::solve(Q_X, -Q[-index_E, index_E] %*% u_e))
  if (is.null(py) == FALSE) {
    Matrix::diag(Q_X)[1:length(py)] <- Matrix::diag(Q_X)[1:length(py)] + 1/sigma_e^2
    AtY <- rep(0,dim(Q_X)[1])
    AtY[1:length(py)] <- (y - mu_X[1:length(py)]) / sigma_e^2
    mu_X <- mu_X + as.vector(Matrix::solve(Q_X, AtY))
  }

  x <- rep(0, length(t))

  if(sample){
    R_X <- Matrix::Cholesky(Q_X, LDL = FALSE, perm = TRUE)
    z <- rnorm(dim(R_X)[1])
    x[-index_E] <- mu_X + as.vector(Matrix::solve(R_X, Matrix::solve(R_X,z,system = 'Lt'),
                                                  system='Pt'))
    x[index_E] <- u_e
  }else{
    x[-index_E] <- mu_X
    x[index_E] <- u_e

  }

  x_out <- matrix(0, nrow=length(t0), 2)
  x_out[, 1] <- t0
  for (i in 1:length(t0))
  {
    ind <- which(t == t0[i])
    x_out[i,2] <- x[ind]
  }
  return(x_out)
}


#' Sample Gaussian process with alpha = 2 on a line given end points
#' @details Samples a Gaussian process \eqn{u(t)} with alpha = 2 on an
#' interval \eqn{(0,l_e)} conditionally on \eqn{u(0), u(l_e)}.
#' If `y` and `py` are supplied, the sampling is done conditionally on observations
#' \deqn{y_i = u(t_i) + sigma_e e_i}{y_i = u(t_i) + sigma_e e_i}
#' where \eqn{e_i} are independent standard Gaussian variables.
#' @param kappa parameter kappa
#' @param theta parameter theta
#' @param sigma_e parameter sigma_e
#' @param u_e  (4 x 1) process and derivative at the two end points
#' @param l_e (1 x 1) line length
#' @param  t (n x 1) distance on the line to be sampled from (not end points)
#' @param  nt (1 x 1) number of equidistance points to sample from if t is  null
#' @param  py  (n x 1) observation locations
#' @param  y (n x 1) observations
#' @param  sample (bool) if true sample else return posterior mean
#' @export
sample_alpha2_line <-function(kappa, sigma, sigma_e,
                              u_e, l_e, t=NULL, Line=NULL,
                              nt=100,  py=NULL, y=NULL, sample=TRUE) {

  if(is.null(t)){
    t  = seq(0, 1, length.out = nt)
    t <- t * l_e
  }
  t_end <- c(0, l_e)
  t <- unique(t)
  t0 <- t
  if (is.null(y) == FALSE) {
    ind_remove = which(py %in% t_end)
    if (length(ind_remove) > 0) {
      y  <- y[-ind_remove]
      py <- py[-ind_remove]
    }

    ind_remove_t <- which(t %in% c(t_end,py))
    if(length(ind_remove_t) > 0)
      t <- t[-ind_remove_t]

    if (length(py) == 0)
      py <- NULL
  } else {
    ind_remove_t <- which(t %in% t_end)
    if (length(ind_remove_t)>0)
      t <- t[-ind_remove_t]
  }

  t <- c(t_end, t)
  if (is.null(py) == FALSE) {
    t <- c(py, t)
  }
  Sigma <- matrix(0, length(t) + 2, length(t) + 2)
  d.index <- c(1, 2)
  index_E <- 2 + length(py) + 1:2
  D <- outer (t, t, `-`)
  Sigma[-d.index, -d.index] <- r_2(D, kappa = kappa,
                                   sigma = sigma, deriv = 0)
  Sigma[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l_e))),
                                  kappa = kappa, sigma = sigma, deriv = 2)
  Sigma[d.index, -d.index] <- -r_2(D[index_E-2,],kappa = kappa,
                                   sigma = sigma, deriv = 1)
  Sigma[-d.index,  d.index] <- t(Sigma[d.index,  -d.index])

  index_boundary <- c(d.index,index_E)
  u_e <- u_e[c(2, 4, 1, 3)]
  SinvS <- solve(Sigma[index_boundary, index_boundary],
                 Sigma[index_boundary, -index_boundary])
  Sigma_X <- Sigma[-index_boundary, -index_boundary] -
    Sigma[-index_boundary, index_boundary] %*% SinvS
  mu_X <- - t(SinvS) %*% (0-u_e)

  if(is.null(py) == FALSE){
    index_y <- 1:length(py)
    Sigma_Y <- Sigma_X[index_y, index_y, drop = FALSE]
    Matrix::diag(Sigma_Y) <-  Matrix::diag(Sigma_Y) + sigma_e^2

    SinvS <- solve(Sigma_Y, Sigma_X[index_y,,drop = FALSE])
    Sigma_X <- Sigma_X - Sigma_X[,index_y, drop = FALSE] %*% SinvS
    mu_X <- mu_X + t(SinvS) %*% (y- mu_X[index_y])
  }

  x <- rep(0, length(t))

  if(sample){
    R_X <- chol(Sigma_X)
    z <- rnorm(dim(Sigma_X)[1])
    x[-c(1:2)] <- mu_X + t(R_X) %*% z
    x[c(1:2)] <- u_e[c(3, 4)]
  }else{
    x[-c(1:2)] <- mu_X
    x[1:2] <- u_e[c(3, 4)]

  }
  x_out <- matrix(0, nrow = length(t0), 2)
  x_out[, 1] <- t0
  for(i in 1:length(t0))
  {
    ind <- which(t == t0[i])
    x_out[i, 2] <- x[ind]
  }
  return(x_out)
}
