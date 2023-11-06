#' Samples a Whittle-Matérn field on a metric graph
#' 
#' Obtains samples of a Whittle-Matérn field on a metric graph.
#' 
#' @details Samples a Gaussian Whittle-Matérn field on a metric graph, either
#' from the prior or conditionally on observations
#' \deqn{y_i = u(t_i) + \sigma_e e_i}{y_i = u(t_i) + \sigma_e e_i}
#' on the graph,  where \eqn{e_i} are independent standard Gaussian variables.
#' The parameters for the field can either be specified in terms of tau and kappa
#' or practical correlation range and marginal standard deviation.
#' @param kappa Range parameter.
#' @param tau Precision parameter.
#' @param sigma Marginal standard deviation parameter.
#' @param range Practical correlation range parameter.
#' @param sigma_e Standard deviation of the measurement noise.
#' @param alpha Smoothness parameter.
#' @param graph A `metric_graph` object.
#' @param PtE Matrix with locations (edge, normalized distance on edge) where
#' the samples should be generated.
#' @param type If "manual" is set, then sampling is done at the locations
#' specified in `PtE`. Set to "mesh" for simulation at mesh nodes, and to "obs"
#' for simulation at observation locations.
#' @param posterior Sample conditionally on the observations?
#' @param nsim Number of samples to be generated.
#' @param method Which method to use for the sampling? The options are
#' "conditional" and "Q". Here, "Q" is more stable but takes longer.
#' @param BC Boundary conditions for degree 1 vertices. BC = 0 gives Neumann
#' boundary conditions and BC = 1 gives stationary boundary conditions.
#' @return Matrix or vector with the samples.
#' @export
sample_spde <- function(kappa, tau, range, sigma, sigma_e = 0, alpha = 1, graph,
                        PtE = NULL,
                        type = "manual", posterior = FALSE,
                        nsim = 1,
                        method = c("conditional", "Q"),
                        BC = 1) {

  check <- check_graph(graph)
  method <- method[[1]]

  if((missing(kappa) || missing(tau)) && (missing(sigma) || missing(range))){
    stop("You should either provide either kappa and tau, or sigma and range.")
  } else if(!missing(sigma) && !missing(range)){
    nu <- alpha - 0.5
    kappa <- sqrt(8 * nu) / range
    tau <- sqrt(gamma(nu) / (sigma^2 * kappa^(2 * nu) *
    (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))
  }

  if (!(type %in% c("manual","mesh", "obs"))) {
    stop("Type must be 'manual', 'mesh' or 'obs'.")
  }
  if( type == "mesh" && !check$has.mesh) {
    stop("mesh must be provided")
  }

  if(posterior && !check$has.data){
    stop("The graph contains no data.")
  }

  if (!(type %in% c("manual","obs", "mesh"))) {
    stop("Type must be 'manual', 'obs' or 'mesh'.")
  }
  if(type == "mesh" && !check$has.mesh) {
    stop("No mesh provided in the graph object.")
  }
  if (type == "obs" && !check$has.data) {
    stop("no observation locations in mesh object.")
  }
  if(is.null(PtE) && type == "manual") {
    stop("must provide PtE for manual mode.")
  }
  if(!is.null(PtE) && !(type == "manual")) {
    warning("PtE provided but mode is not manual.")
  }
  if((nsim == 1) || (method == "Q")){
  if (!posterior) {
    if (alpha == 1) {
        if(method == "conditional"){
              Q <- spde_precision(kappa = kappa, tau = tau,
                                  alpha = 1, graph = graph, BC=BC)
              R <- Cholesky(Q,LDL = FALSE, perm = TRUE)
              V0 <- as.vector(solve(R, solve(R,rnorm(graph$nV),
                                             system = 'Lt'), system = 'Pt'))


              if(type == "mesh") {
                u <- V0
                inds_PtE <- unique(graph$mesh$PtE[,1])
              } else if (type == "obs") {
                u <- NULL
                inds_PtE <- unique(graph$PtE[,1])
              } else {
                order_PtE <- order(PtE[,1], PtE[,2])
                ordered_PtE <- PtE[order_PtE,]
                u <- NULL
                inds_PtE <- unique(ordered_PtE[,1])
              }

              for (i in inds_PtE) {
                if(type == "mesh") {
                  t <- graph$mesh$PtE[graph$mesh$PtE[,1] == i, 2]
                } else if (type == "obs") {
                  t <- graph$PtE[graph$PtE[,1] == i, 2]
                } else {
                  t <- ordered_PtE[ordered_PtE[,1] == i, 2]
                }

                samp <- sample_alpha1_line(kappa = kappa, tau = tau,
                                           u_e = V0[graph$E[i, ]], t = t,
                                           l_e = graph$edge_lengths[i])
                u <- c(u, samp[,2])
              }
                if(type == "manual"){
                  u[order_PtE] <- u
                }
    } else if(method == "Q"){
        if(type == "manual"){
          graph_tmp <- graph$get_initial_graph()
          order_PtE <- order(PtE[,1], PtE[,2])
          n_obs_add <- nrow(PtE)
          n_obs_tmp <- n_obs_add
          y_tmp <- rep(NA, n_obs_add)
          if(max(PtE[,2])>1){
            stop("You should provide normalized locations!")
          }
          df_graph <- data.frame(y = y_tmp, edge_number = PtE[,1],
                      distance_on_edge = PtE[,2])
          graph_tmp$add_observations(data = df_graph, normalized=TRUE)
          graph_tmp$observation_to_vertex()
          Q_tmp <- Qalpha1(theta = c(tau, kappa), graph_tmp, BC=BC)
        } else if(type == "obs"){
          Q_tmp <- Qalpha1(theta = c(tau, kappa), graph_tmp, BC=BC)
          n_obs_tmp <- length(graph$data[[".group"]])
          order_PtE <- 1:n_obs_tmp
        } else if(type == "mesh"){
          graph_tmp <- graph$get_initial_graph()
          n_obs_mesh <- nrow(graph$mesh$PtE)
          y_tmp <- rep(NA, n_obs_mesh)
          df_graph <- data.frame(y = y_tmp, edge_number = graph$mesh$PtE[,1],
                      distance_on_edge = graph$mesh$PtE[,2])
          graph_tmp$add_observations(data = df_graph, normalized=TRUE)
          graph_tmp$observation_to_vertex()
          Q_tmp <- Qalpha1(theta = c(tau, kappa), graph_tmp, BC=BC)
          n_obs_tmp <- dim(Q_tmp)[1]
          order_PtE <- 1:n_obs_tmp
        }
          sizeQ <- nrow(Q_tmp)
          Z <- rnorm(sizeQ * nsim)
          dim(Z) <- c(sizeQ, nsim)
          LQ <- chol(forceSymmetric(Q_tmp))
          u <- solve(LQ, Z)
          gap <- sizeQ - n_obs_tmp
          u <- as.vector(u[(gap+1):sizeQ])
          u[order_PtE] <- u

    } else{
      stop("Method should be either 'conditional' or 'Q'!")
    }

    } else if (alpha == 2) {

      Q <- spde_precision(kappa = kappa, tau = tau,
                          alpha = 2, graph = graph, BC = BC)
      if(is.null(graph$C))
        graph$buildC(2)

      Qmod <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
      Qtilde <- Qmod[-c(1:dim(graph$CoB$U)[1]),-c(1:dim(graph$CoB$U)[1])]
      R <- Cholesky(forceSymmetric(Qtilde),LDL = FALSE, perm = TRUE)
      V0 <- as.vector(solve(R, solve(R,rnorm(4*graph$nE - dim(graph$CoB$U)[1]),
                                     system = 'Lt'), system = 'Pt'))
      u_e <- t(graph$CoB$T) %*% c(rep(0, dim(graph$CoB$U)[1]), V0)
      VtE <- graph$VtEfirst()


      if(type == "mesh") {
        initial_graph <- graph$get_initial_graph()
        u_s <- u_e[seq(from=1, by = 2, to = length(u_e))]
        u <- u_s[which(!duplicated(c(t(initial_graph$E))))]
        inds_PtE <- unique(graph$mesh$PtE[,1])
      } else if (type == "obs") {
        u <- NULL
        inds_PtE <- unique(graph$PtE[,1])
      } else {
        order_PtE <- order(PtE[,1], PtE[,2])
        ordered_PtE <- PtE[order_PtE,]
        u <- NULL
        inds_PtE <- unique(ordered_PtE[,1])
      }

      for (i in inds_PtE) {
        if(type == "mesh") {
          t <- graph$mesh$PtE[graph$mesh$PtE[,1] == i, 2]
        } else if (type == "obs") {
          t <- graph$PtE[graph$PtE[,1] == i, 2]
        } else {
          t <- ordered_PtE[ordered_PtE[,1] == i, 2]
        }
        samp <- sample_alpha2_line(kappa = kappa, tau = tau,
                                   sigma_e = sigma_e,
                                   u_e = u_e[4*(i-1) +1:4],
                                   t = t,
                                   l_e = graph$edge_lengths[i])
        u <- c(u, samp[,2])
      }
      if(type == "manual"){
        u[order_PtE] <- u
      }
    } else {
      stop("only alpha = 1 and alpha = 2 implemented.")
    }
  } else {
    stop("TODO: implement posterior sampling")
  }
  return(u)
  } else if ((nsim%%1 == 0) && nsim>1 && method != "Q"){
    u_rep <- unlist(lapply(1:nsim, function(i){
      sample_spde(kappa=kappa, tau=tau, range=range, sigma=sigma, sigma_e = sigma_e,
      alpha = alpha, graph = graph,
                        PtE = PtE,
                        type = type,
                        posterior = posterior,
                        nsim = 1)
    }))
    return(matrix(u_rep, ncol = nsim))
  } else{
    stop("The number of simulations must be an integer greater than zero!")
  }
}
#' Samples a Gaussian process with exponential covariance on an interval given
#' the values at the end points.
#' @details Samples a Gaussian process \eqn{u(t)} with an exponential covariance
#' function
#' \deqn{r(h) = \sigma^2\exp(-\kappa h)/(2\kappa)}{r(h) = sigma^2*(exp(-kappa*h)/(2*kappa)}
#' on an interval \eqn{(0,l_e)} conditionally on \eqn{u(0), u(l_e)}.
#' If `y` and `py` are supplied, the sampling is done conditionally on
#' observations
#' \deqn{y_i = u(t_i) + sigma_e e_i}{y_i = u(t_i) + sigma_e e_i}
#' where \eqn{e_i} are independent standard Gaussian variables.
#' @param kappa parameter kappa
#' @param tau parameter tau
#' @param sigma_e parameter sigma_e
#' @param u_e  (2 x 1) the two end points
#' @param l_e (1 x 1) line length
#' @param  t (n x 1) distance on the line to be sampled from (not end points)
#' @param  nt (1 x 1) number of equidistance points to sample from if t is  null
#' @param  py  (n x 1) observation locations
#' @param  y (n x 1) observations
#' @param  sample (bool) if true sample else return posterior mean
#' @noRd
sample_alpha1_line <- function(kappa, tau, sigma_e,
                               u_e, l_e, t = NULL,
                               nt = 100,  py = NULL,
                               y = NULL, sample = TRUE) {

  if (is.null(t)) {
    t  <- seq(0, 1, length.out = nt)
  }
  t <- t * l_e

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

  Q <- precision_exp_line(kappa = kappa, tau = tau, t = t)

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
#' @param tau parameter tau
#' @param sigma_e parameter sigma_e
#' @param u_e  (4 x 1) process and derivative at the two end points
#' @param l_e (1 x 1) line length
#' @param  t (n x 1) distance on the line to be sampled from (not end points)
#' @param  nt (1 x 1) number of equidistance points to sample from if t is  null
#' @param  py  (n x 1) observation locations
#' @param  y (n x 1) observations
#' @param  sample (bool) if true sample else return posterior mean
#' @noRd
sample_alpha2_line <-function(kappa, tau, sigma_e,
                              u_e, l_e, t=NULL, Line=NULL,
                              nt=100,  py=NULL, y=NULL, sample=TRUE) {

  if(is.null(t)){
    t  = seq(0, 1, length.out = nt)
  }
    t <- t * l_e
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
                                   tau = tau, deriv = 0)
  Sigma[d.index, d.index] <- -r_2(as.matrix(dist(c(0,l_e))),
                                  kappa = kappa, tau = tau, deriv = 2)
  Sigma[d.index, -d.index] <- -r_2(D[index_E-2,],kappa = kappa,
                                   tau = tau, deriv = 1)
  Sigma[-d.index,  d.index] <- t(Sigma[d.index,  -d.index])

  index_boundary <- c(d.index,index_E)
  u_e <- u_e[c(2, 4, 1, 3)]
  if(length(Sigma[index_boundary, -index_boundary])>0){
    SinvS <- solve(Sigma[index_boundary, index_boundary],
                   Sigma[index_boundary, -index_boundary])
    Sigma_X <- Sigma[-index_boundary, -index_boundary] -
      Sigma[-index_boundary, index_boundary] %*% SinvS
    mu_X <- - t(SinvS) %*% (0-u_e)
  } else{
    Sigma_X <- Sigma
    mu_X <- 0
  }


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
