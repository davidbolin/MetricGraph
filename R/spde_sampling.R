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
#' @param directional should we use directional model currently only for alpha=1
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
sample_spde <- function(kappa, tau, range, sigma, sigma_e = 0, alpha = 1,
                        directional = FALSE,
                        graph,
                        PtE = NULL,
                        type = "manual", posterior = FALSE,
                        nsim = 1,
                        method = c("conditional", "Q"),
                        BC = 1) {

  if (inherits(graph, "graph_components")) {
    method <- match.arg(method, c("conditional", "Q"))
    has_kappa <- !missing(kappa)
    has_tau   <- !missing(tau)
    has_range <- !missing(range)
    has_sigma <- !missing(sigma)
    build_args <- function(extra) {
      a <- c(extra, list(posterior = posterior, nsim = nsim,
                         method = method, BC = BC, sigma_e = sigma_e,
                         alpha = alpha, directional = directional))
      if (has_kappa) a$kappa <- kappa
      if (has_tau)   a$tau   <- tau
      if (has_range) a$range <- range
      if (has_sigma) a$sigma <- sigma
      a
    }

    if (type == "manual") {
      if (is.null(PtE)) stop("must provide PtE for manual mode.")
      if (NCOL(PtE) != 3) {
        stop("For 'graph_components', PtE must have 3 columns: (component, edge, distance).")
      }
      PtE <- as.matrix(PtE)
      comp_id <- as.integer(PtE[, 1])
      orig_order <- seq_len(nrow(PtE))
      out_pieces <- vector("list", graph$n)
      out_idx    <- vector("list", graph$n)
      for (k in seq_len(graph$n)) {
        sel <- which(comp_id == k)
        if (length(sel) == 0L) next
        sub_PtE <- PtE[sel, 2:3, drop = FALSE]
        out_pieces[[k]] <- do.call(sample_spde, build_args(list(
          graph = graph$graphs[[k]], PtE = sub_PtE, type = "manual"
        )))
        out_idx[[k]] <- orig_order[sel]
      }
      flat_idx <- unlist(out_idx, use.names = FALSE)
      if (nsim == 1L) {
        u <- numeric(nrow(PtE))
        u[flat_idx] <- unlist(out_pieces, use.names = FALSE)
        return(u)
      }
      u <- matrix(0, nrow = nrow(PtE), ncol = nsim)
      for (k in seq_len(graph$n)) {
        if (is.null(out_pieces[[k]])) next
        u[out_idx[[k]], ] <- out_pieces[[k]]
      }
      return(u)
    }

    # type = "mesh" or "obs": per-component sampling, concat
    if (type == "mesh") {
      if (any(vapply(graph$graphs, function(g) is.null(g$mesh),
                     logical(1)))) {
        stop("Every component must have a mesh; call graph$build_mesh() first.")
      }
    }
    out_pieces <- vector("list", graph$n)
    for (k in seq_len(graph$n)) {
      g <- graph$graphs[[k]]
      if (type == "obs" && is.null(g$.__enclos_env__$private$data)) next
      out_pieces[[k]] <- do.call(sample_spde, build_args(list(
        graph = g, type = type
      )))
    }
    if (nsim == 1L) {
      return(unlist(out_pieces, use.names = FALSE))
    }
    return(do.call(rbind, out_pieces))
  }

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
    if (alpha == 1 && directional == F) {
        if(method == "conditional"){
              Q <- spde_precision(kappa = kappa, tau = tau,
                                  alpha = 1, graph = graph, BC=BC)
              R <- Cholesky(Q,LDL = FALSE, perm = TRUE)
              V0 <- as.vector(solve(R, solve(R,rnorm(graph$nV),
                                             system = 'Lt'), system = 'Pt'))


              if(type == "mesh") {
                u <- V0
                inds_PtE <- unique(graph$mesh$PtE[,1])
                t_by_edge <- split(graph$mesh$PtE[,2], graph$mesh$PtE[,1])
              } else if (type == "obs") {
                inds_PtE <- unique(graph$PtE[,1])
                t_by_edge <- split(graph$PtE[,2], graph$PtE[,1])
              } else {
                order_PtE <- order(PtE[,1], PtE[,2])
                ordered_PtE <- PtE[order_PtE, , drop = FALSE]
                inds_PtE <- unique(ordered_PtE[,1])
                t_by_edge <- split(ordered_PtE[,2], ordered_PtE[,1])
              }

              E_loc <- graph$E
              el_loc <- graph$edge_lengths
              u_list <- vector("list", length(inds_PtE))
              for (k in seq_along(inds_PtE)) {
                i <- inds_PtE[k]
                t <- t_by_edge[[as.character(i)]]
                samp <- sample_alpha1_line(kappa = kappa, tau = tau,
                                           u_e = V0[E_loc[i, ]], t = t,
                                           l_e = el_loc[i])
                u_list[[k]] <- samp[,2]
              }
              if(type == "mesh") {
                u <- c(u, unlist(u_list, use.names = FALSE))
              } else {
                u <- unlist(u_list, use.names = FALSE)
              }
              if(type == "manual"){
                u[order_PtE] <- u
              }
    }else if(method == "Q"){
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
          graph_tmp$add_observations(data = df_graph, normalized=TRUE,
                  suppress_warnings = TRUE, verbose=0)
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
          graph_tmp$add_observations(data = df_graph, normalized=TRUE,
                  suppress_warnings = TRUE, verbose=0)
          graph_tmp$observation_to_vertex()
          Q_tmp <- Qalpha1(theta = c(tau, kappa), graph_tmp, BC=BC)
          n_obs_tmp <- dim(Q_tmp)[1]
          order_PtE <- 1:n_obs_tmp
        }
          sizeQ <- nrow(Q_tmp)
          Z <- rnorm(sizeQ * nsim)
          dim(Z) <- c(sizeQ, nsim)
          LQ <- Cholesky(forceSymmetric(Q_tmp), LDL = FALSE, perm = TRUE)
          u <- solve(LQ, solve(LQ, Z, system = "Lt"), system = "Pt")
          gap <- sizeQ - n_obs_tmp
          u <- as.matrix(u[(gap+1):sizeQ, , drop = FALSE])
          u[order_PtE, ] <- u
          if (nsim == 1) u <- as.vector(u)

    } else{
      stop("Method should be either 'conditional' or 'Q'!")
    }

    }else if(alpha == 1 && directional == T) {

      Q <- Qalpha1_edges(c( tau,kappa),
                         graph,
                         w = 0,
                         BC=1, build=T)
      if(is.null(graph$C) || graph$CoB$alpha != 1){
        graph$buildDirectionalConstraints(alpha = 1)
      }
      n_const <- length(graph$CoB$S)
      ind.const <- c(1:n_const)
      Tc <- graph$CoB$T[-ind.const,]
      Q <- Tc %*% Q %*% t(Tc)
      R <- Cholesky(Q, LDL = FALSE, perm = TRUE)
      V0 <- as.vector(solve(R,
                      solve(R,rnorm(dim(R)[1]),system = 'Lt'),
                      system = 'Pt'))
      u_e <- t(graph$CoB$T) %*% c(rep(0, dim(graph$CoB$U)[1]), V0)
      VtE <- graph$VtEfirst()
      if(type == "mesh") {
        initial_graph <- graph$get_initial_graph()
        u_s <- u_e
        u <- u_s[which(!duplicated(c(t(initial_graph$E))))]
        inds_PtE <- unique(graph$mesh$PtE[,1])
        t_by_edge <- split(graph$mesh$PtE[,2], graph$mesh$PtE[,1])
      } else if (type == "obs") {
        inds_PtE <- unique(graph$PtE[,1])
        t_by_edge <- split(graph$PtE[,2], graph$PtE[,1])
      } else {
        order_PtE <- order(PtE[,1], PtE[,2])
        ordered_PtE <- PtE[order_PtE,]
        inds_PtE <- unique(ordered_PtE[,1])
        t_by_edge <- split(ordered_PtE[,2], ordered_PtE[,1])
      }

      el_loc <- graph$edge_lengths
      u_list <- vector("list", length(inds_PtE))
      for (k in seq_along(inds_PtE)) {
        i <- inds_PtE[k]
        t <- t_by_edge[[as.character(i)]]
        samp <- sample_alpha1_line(kappa = kappa, tau = tau,
                                   sigma_e = sigma_e,
                                   u_e = u_e[2*(i-1) +1:2],
                                   t = t,
                                   l_e = el_loc[i])
        u_list[[k]] <- samp[,2]
      }
      if(type == "mesh") {
        u <- c(u, unlist(u_list, use.names = FALSE))
      } else {
        u <- unlist(u_list, use.names = FALSE)
      }
      if(type == "manual"){
        u[order_PtE] <- u
      }

    }else if (alpha == 2) {

      Q <- spde_precision(kappa = kappa, tau = tau,
                          alpha = 2, graph = graph, BC = BC)
      if(is.null(graph$CoB)){
        graph$buildC(2)
      } else if(graph$CoB$alpha == 1){
        graph$buildC(2)
      }

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
        t_by_edge <- split(graph$mesh$PtE[,2], graph$mesh$PtE[,1])
      } else if (type == "obs") {
        inds_PtE <- unique(graph$PtE[,1])
        t_by_edge <- split(graph$PtE[,2], graph$PtE[,1])
      } else {
        order_PtE <- order(PtE[,1], PtE[,2])
        ordered_PtE <- PtE[order_PtE,]
        inds_PtE <- unique(ordered_PtE[,1])
        t_by_edge <- split(ordered_PtE[,2], ordered_PtE[,1])
      }

      el_loc <- graph$edge_lengths
      u_list <- vector("list", length(inds_PtE))
      for (k in seq_along(inds_PtE)) {
        i <- inds_PtE[k]
        t <- t_by_edge[[as.character(i)]]
        samp <- sample_alpha2_line(kappa = kappa, tau = tau,
                                   sigma_e = sigma_e,
                                   u_e = u_e[4*(i-1) +1:4],
                                   t = t,
                                   l_e = el_loc[i])
        u_list[[k]] <- samp[,2]
      }
      if(type == "mesh") {
        u <- c(u, unlist(u_list, use.names = FALSE))
      } else {
        u <- unlist(u_list, use.names = FALSE)
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
#' @param  py  (m x 1) observation locations
#' @param  y (m x 1) observations
#' @param  sample (bool) if true sample else return posterior mean
#' @return x (n x 2)  1- position on the edge, 2- value of the simulations
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

  # Build tridiagonal precision matrix as dense (small per-edge matrices;
  # dense base-R chol/solve avoids S4 dispatch overhead from Matrix package)
  l_t <- length(t)
  order_t <- order(t)
  t_sorted <- t[order_t]
  # inv_order_t maps sorted indices back to original positions
  inv_order_t <- integer(l_t)
  inv_order_t[order_t] <- seq_len(l_t)

  Q_sorted <- matrix(0, l_t, l_t)
  scale <- 2 * kappa * tau^2
  for (i in 2:l_t) {
    c1 <- exp(-kappa * (t_sorted[i] - t_sorted[i - 1]))
    c2 <- c1 * c1
    one_m_c2 <- 1 - c2
    c_1 <- scale * (0.5 + c2 / one_m_c2)
    c_2 <- scale * (-c1 / one_m_c2)
    Q_sorted[i, i]         <- Q_sorted[i, i]         + c_1
    Q_sorted[i - 1, i - 1] <- Q_sorted[i - 1, i - 1] + c_1
    Q_sorted[i, i - 1]     <- Q_sorted[i, i - 1]     + c_2
    Q_sorted[i - 1, i]     <- Q_sorted[i - 1, i]     + c_2
  }
  Q_sorted[1, 1]     <- Q_sorted[1, 1]     + scale * 0.5
  Q_sorted[l_t, l_t] <- Q_sorted[l_t, l_t] + scale * 0.5

  # Unsort to match original t ordering
  Q <- Q_sorted[inv_order_t, inv_order_t]

  index_E <- length(py) + 1:2
  Q_X <- Q[-index_E, -index_E, drop = FALSE]
  rhs <- -Q[-index_E, index_E, drop = FALSE] %*% u_e
  mu_X <- as.vector(solve(Q_X, rhs))

  if (!is.null(py)) {
    diag(Q_X)[1:length(py)] <- diag(Q_X)[1:length(py)] + 1 / sigma_e^2
    AtY <- rep(0, nrow(Q_X))
    AtY[1:length(py)] <- (y - mu_X[1:length(py)]) / sigma_e^2
    mu_X <- mu_X + as.vector(solve(Q_X, AtY))
  }

  x <- rep(0, l_t)

  if (sample) {
    R_X <- chol(Q_X)
    z <- rnorm(nrow(Q_X))
    x[-index_E] <- mu_X + backsolve(R_X, z)
    x[index_E] <- u_e
  } else {
    x[-index_E] <- mu_X
    x[index_E] <- u_e
  }

  x_out <- matrix(0, nrow = length(t0), 2)
  x_out[, 1] <- t0
  x_out[, 2] <- x[match(t0, t)]
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
    diag(Sigma_Y) <- diag(Sigma_Y) + sigma_e^2

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
  x_out[, 2] <- x[match(t0, t)]
  return(x_out)
}
