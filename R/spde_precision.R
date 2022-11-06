#' The precision matrix for all vertices for a Whittle-Matérn field 
#' @param kappa range parameter kappa
#' @param sigma variance parameter
#' @param alpha smoothness parameter (1 or 2)
#' @param graph metric_graph object
#' @param BC boundary conditions for degree=1 vertices. BC =0 gives Neumann
#' boundary conditions and BC=1 gives stationary boundary conditions
#' @param build (bool) if TRUE return the precision matrix otherwise return
#' a list(i,j,x, nv)
#' @return Precision matrix or list
#' @export
spde_precision <- function(kappa, sigma, alpha, graph, BC = 1, build = TRUE) {
    if (alpha == 1) {
        return(Qalpha1(theta = c(kappa, sigma),
                       graph = graph,
                       BC = BC,
                       build = build))
    } else if (alpha == 2) {
        return(Qalpha2(theta = c(kappa, sigma),
                       graph = graph,
                       BC = BC,
                       build = build))
    }
}

#' The precision matrix for all vertices in the alpha=1 case
#' @param theta - kappa, sigma
#' @param graph metric_graph object
#' @param BC boundary conditions for degree=1 vertices. BC =0 gives Neumann
#' boundary conditions and BC=1 gives stationary boundary conditions
#' @param build (bool) if TRUE return the precision matrix otherwise return
#' a list(i,j,x, nv)
#' @return Precision matrix or list
Qalpha1 <- function(theta, graph, BC = 1, build = TRUE) {

  kappa <- theta[1]
  sigma <- theta[2]
  i_ <- j_ <- x_ <- rep(0, dim(graph$V)[1]*4)
  count <- 0
  for(i in 1:graph$nE){
    l_e <- graph$edge_lengths[i]
    c1 <- exp(-kappa*l_e)
    c2 <- c1^2
    one_m_c2 = 1-c2
    c_1 = 0.5 + c2/one_m_c2
    c_2 = -c1/one_m_c2

    if (graph$E[i, 1] != graph$E[i, 2]) {

      i_[count + 1] <- graph$E[i, 1]
      j_[count + 1] <- graph$E[i, 1]
      x_[count + 1] <- c_1

      i_[count + 2] <- graph$E[i, 2]
      j_[count + 2] <- graph$E[i, 2]
      x_[count + 2] <- c_1


      i_[count + 3] <- graph$E[i, 1]
      j_[count + 3] <- graph$E[i, 2]
      x_[count + 3] <- c_2

      i_[count + 4] <- graph$E[i, 2]
      j_[count + 4] <- graph$E[i, 1]
      x_[count + 4] <- c_2
      count <- count + 4
    }else{
      i_[count + 1] <- graph$E[i, 1]
      j_[count + 1] <- graph$E[i, 1]
      x_[count + 1] <- tanh(0.5 * theta * l_e)
      count <- count + 1
    }
  }
  if(BC == 1){
    #does this work for circle?
    i.table <- table(i_[1:count])
    index = as.integer(names(which(i.table < 3)))
    i_ <- c(i_[1:count], index)
    j_ <- c(j_[1:count], index)
    x_ <- c(x_[1:count], rep(0.5, length(index)))
    count <- count + length(index)
  }
  if(build){
    Q <- Matrix::sparseMatrix(i = i_[1:count],
                              j = j_[1:count],
                              x = (2 * kappa / sigma^2) * x_[1:count],
                              dims = c(graph$nV, graph$nV))


    return(Q)
  } else {
    return(list(i = i_[1:count],
                j = j_[1:count],
                x = (2 * kappa / sigma^2) * x_[1:count],
                dims = c(graph$nV, graph$nV)))
  }
}


#' The precision matrix for all vertices in the alpha=2 case
#' @param theta - kappa, sigma
#' @param graph metric_graph object
#' @param BC boundary conditions for degree=1 vertices. BC =0 gives Neumann
#' boundary conditions and BC=1 gives stationary boundary conditions
#' @param build (bool) if TRUE return the precision matrix otherwise return
#' a list(i,j,x, nv)
#' @details This is the unconstrained precision matrix of the process and its 
#' derivatives. The ordering of the variables is acording to graph$E, where for
#' each edge there are four random variables: processes and derivate for
#' lower and upper edge end points
#' @return Precision matrix or list
Qalpha2 <- function(theta, graph, BC = 1, build = TRUE) {

  kappa <- theta[1]
  sigma <- theta[2]
  i_ <- j_ <- x_ <- rep(0, graph$nE * 16)
  count <- 0

  R_00 <- matrix(c( r_2(0, c(kappa, sigma), 0),
                   -r_2(0, c(kappa, sigma), 1), 
                   -r_2(0, c(kappa, sigma), 1), 
                   -r_2(0, c(kappa, sigma), 2)), 2, 2)
  R_node <- rbind(cbind(R_00, matrix(0, 2, 2)),
                  cbind(matrix(0, 2, 2), R_00))
  Ajd <- -0.5 * solve(rbind(cbind(R_00, matrix(0, 2, 2)),
                            cbind(matrix(0, 2, 2), R_00)))
  for (i in 1:graph$nE) {

    l_e <- graph$edge_lengths[i]
    #lots of redundant caculations
    d_ <- c(0, l_e)
    D <- outer(d_, d_, "-")
    r_0l <-   r_2(l_e, c(kappa, sigma))
    r_11 <- - r_2(l_e, c(kappa, sigma),2)
    # order by node not derivative
    R_01 <- matrix(c(r_0l, r_2(-l_e, c(kappa, sigma),1),
                     r_2(l_e, c(kappa, sigma), 1), r_11), 2, 2)
    
    R_node[1:2, 3:4] <- R_01
    R_node[3:4, 1:2] <- t(R_01)
    Q_adj <- solve(R_node) + Ajd

    if (graph$E[i, 1] != graph$E[i, 2]) {

      #lower edge precision u
      i_[count + 1] <- 4 * (i - 1) + 1
      j_[count + 1] <- 4 * (i - 1) + 1
      x_[count + 1] <- Q_adj[1, 1]

      #lower edge  u'
      i_[count + 2] <- 4 * (i - 1) + 2
      j_[count + 2] <- 4 * (i - 1) + 2
      x_[count + 2] <- Q_adj[2, 2]

      #upper edge  u
      i_[count + 3] <- 4 * (i - 1) + 3
      j_[count + 3] <- 4 * (i - 1) + 3
      x_[count + 3] <- Q_adj[3, 3]

      #upper edge  u'
      i_[count + 4] <- 4 * (i - 1) + 4
      j_[count + 4] <- 4 * (i - 1) + 4
      x_[count + 4] <- Q_adj[4, 4]

      #lower edge  u, u'
      i_[count + 5] <- 4 * (i - 1) + 1
      j_[count + 5] <- 4 * (i - 1) + 2
      x_[count + 5] <- Q_adj[1, 2]
      i_[count + 6] <- 4 * (i - 1) + 2
      j_[count + 6] <- 4 * (i - 1) + 1
      x_[count + 6] <- Q_adj[1, 2]

      #upper edge  u, u'
      i_[count + 7] <- 4 * (i - 1) + 3
      j_[count + 7] <- 4 * (i - 1) + 4
      x_[count + 7] <- Q_adj[3, 4]
      i_[count + 8] <- 4 * (i - 1) + 4
      j_[count + 8] <- 4 * (i - 1) + 3
      x_[count + 8] <- Q_adj[3, 4]

      #lower edge  u, upper edge  u,
      i_[count + 9]  <- 4 * (i - 1) + 1
      j_[count + 9]  <- 4 * (i - 1) + 3
      x_[count + 9]  <- Q_adj[1, 3]
      i_[count + 10] <- 4 * (i - 1) + 3
      j_[count + 10] <- 4 * (i - 1) + 1
      x_[count + 10] <- Q_adj[1, 3]

      #lower edge  u, upper edge  u',
      i_[count + 11] <- 4 * (i - 1) + 1
      j_[count + 11] <- 4 * (i - 1) + 4
      x_[count + 11] <- Q_adj[1, 4]
      i_[count + 12] <- 4 * (i - 1) + 4
      j_[count + 12] <- 4 * (i - 1) + 1
      x_[count + 12] <- Q_adj[1, 4]

      #lower edge  u', upper edge  u,
      i_[count + 13] <- 4 * (i - 1) + 2
      j_[count + 13] <- 4 * (i - 1) + 3
      x_[count + 13] <- Q_adj[2, 3]
      i_[count + 14] <- 4 * (i - 1) + 3
      j_[count + 14] <- 4 * (i - 1) + 2
      x_[count + 14] <- Q_adj[2, 3]

      #lower edge  u', upper edge  u,
      i_[count + 13] <- 4 * (i - 1) + 2
      j_[count + 13] <- 4 * (i - 1) + 3
      x_[count + 13] <- Q_adj[2, 3]
      i_[count + 14] <- 4 * (i - 1) + 3
      j_[count + 14] <- 4 * (i - 1) + 2
      x_[count + 14] <- Q_adj[2, 3]

      #lower edge  u', upper edge  u',
      i_[count + 15] <- 4 * (i - 1) + 2
      j_[count + 15] <- 4 * (i - 1) + 4
      x_[count + 15] <- Q_adj[2, 4]
      i_[count + 16] <- 4 * (i - 1) + 4
      j_[count + 16] <- 4 * (i - 1) + 2
      x_[count + 16] <- Q_adj[2, 4]

      count <- count + 16
    } else {
        error("circlular edges are not implemented")
    }
  }

  if(BC == 1){
    #Vertices with of degree 1
    i.table <- table(c(graph$E))
    index <- as.integer(names(which(i.table == 1)))
    #for this vertices locate position
    lower.edges <- which(graph$E[, 1] %in% index)
    upper.edges <- which(graph$E[, 2] %in% index)
    for (le in lower.edges) {
      ind <- c(4 * (le - 1) + 1, 4 * (le-1) + 2)

      i_ <- c(i_, ind)
      j_ <- c(j_, ind)
      x_ <- c(x_, 0.5*c(1 / R_00[1, 1], 1 / R_00[2, 2]))
      count <- count + 2
    }
    for (ue in upper.edges) {
      ind <- c(4 * (ue - 1) + 3, 4 * (ue - 1) + 4)
      i_ <- c(i_, ind)
      j_ <- c(j_, ind)
      x_ <- c(x_, 0.5 * c(1 / R_00[1, 1], 1 / R_00[2, 2]))
      count <- count + 2
    }
  }
  if (build) {
    Q <- Matrix::sparseMatrix(i    = i_[1:count],
                              j    = j_[1:count],
                              x    = x_[1:count],
                              dims = c(4 * graph$nE, 4 * graph$nE))

    return(Q)
  }else{
    return(list(i = i_[1:count],
                j = j_[1:count],
                x  = x_[1:count],
                dims=c(4 * graph$nE, 4 * graph$nE)))
  }
}


#' The exponential covariance
#' @param D vector or matrix with distances
#' @param theta parameters kappa and sigma
#' @export
r_1 <- function(D, theta) {
  return((theta[2]^2 / (2 * theta[1])) * exp(-theta[1] * abs(D)))
}

#' plot the exponential covariance for parameter set
#' @param theta parameters sigma_e (nugget), kappa and sigma
#' @export
plot_r_1 <- function(theta, t = NULL) {
  if (is.null(t)) {
    r_0 <- r_1(0, theta[2:3])
    r_p <- -log(0.05)/theta[2]
    t <- seq(0, r_p, length.out = 100)
  }
  r_ <- r_1(t, theta[2:3])
  if (t[1] == 0)
    r_[1] = r_[1] + theta[1]^2
  plot(t, r_, type = "l")

}

#' the Whittle--Matern covariance with alpha=1 on a circle
#' @param t locations
#' @param l_e circle perimeter
#' @param theta parameters kappa and sigma
#' @export
r_1_circle <- function(t, l_e, theta) {

  kappa <- theta[1]
  c <- theta[2]^2 / (2 * kappa * sinh(kappa *l_e/2))
  r <- matrix(0, length(t), length(t))
  r_0 <- cosh(-kappa * l_e / 2)
  if (length(t) == 1) {
    r[1, 1] <- c * r_0
    return(r)
  }

  for (i in 1:(length(t) - 1)) {

    for (ii in (i + 1):length(t)){
      r[i, ii] <- cosh(kappa * (abs(t[i] - t[ii]) - l_e / 2))
    }
  }
  r <- r + t(r)
  diag(r) <- r_0
  return(c * r)
}