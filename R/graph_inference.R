
#' Likelihood evaluation not using sparsity
#'
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph.obj      - graph object
#' @param model Model for Gaussian process, supported options are alpha1 (SPDE with alpha=1),
#' alpha2 (SPDE with alpha=2), GL (graph Laplacian model with alpha=1), GL2 (graph Laplacian
#' model with alpha=2) and isoExp (model with isotropic exponential covariance)
#' @return The log-likelihood
#' @export
likelihood.graph.covariance <- function(theta, graph, model = "alpha1"){

  n.o <- length(graph$y)
  n.v <- dim(graph$V)[1]
  #build covariance matrix
  if(model == "alpha1"){

    Q <- Q.exp(theta[2:3], graph)
    Sigma <- as.matrix(solve(Q))[graph$PtV,graph$PtV]

  } else if (model == "alpha2"){

    n.c <- 1:length(graph$CBobj$S)
    Q <- Q.matern2(c(theta[2],theta[3]), graph.obj$V, graph.obj$EtV, graph.obj$edge_lengths, BC = 1)
    Qtilde <- (graph.obj$CBobj$T)%*%Q%*%t(graph.obj$CBobj$T)
    Qtilde <- Qtilde[-n.c,-n.c]
    Sigma.overdetermined  = t(graph.obj$CBobj$T[-n.c,])%*%solve(Qtilde)%*%(graph.obj$CBobj$T[-n.c,])
    index.obs <-  4*(graph.obj$PtE[,1]-1) + (1 * (abs(graph.obj$PtE[,2])<1e-14)) + (3 * (abs(graph.obj$PtE[,2])>1e-14))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])

  } else if (model == "GL"){
    Q <- theta[2]*(theta[3]*Diagonal(graph.obj$nV,1) + graph.obj$Laplacian)
    Sigma <- as.matrix(solve(Q))[graph.obj$PtV,graph.obj$PtV]
  } else if (model == "GL2"){
    Q <- theta[2]*(theta[3]*Diagonal(graph.obj$nV,1) + graph.obj$Laplacian)
    Q <- Q%*%Q
    Sigma <- as.matrix(solve(Q))[graph.obj$PtV,graph.obj$PtV]
  } else if (model == "isoExp"){
    Sigma <- as.matrix(theta[3]^2*exp(-theta[2]*graph.obj$res.dist[graph.obj$PtV,graph.obj$PtV]))
  } else {
    stop("wrong model choice.")
  }

  diag(Sigma) <- diag(Sigma)  +  theta[1]^2

  R <- chol(Sigma)

  return(as.double(-sum(log(diag(R))) - 0.5*t(graph.obj$y)%*%solve(Sigma,graph.obj$y) - n.o*log(2*pi)/2))
}

#' Prediction for models assuming observations at vertices
#'
#' @param theta Estimated model parameters (sigma_e, sigma, kappa)
#' @param graph.obj Graph object
#' @param model Type of model: "alpha1" gives SPDE with alpha=1, "GL" gives
#' graph Laplacian, and "isoExp" gives isotropic exponential
#' @param ind Indices for cross validation. It should be a vector of the same length as the
#' data, with integer values representing each group in the cross-validation.
#' If NULL, leave-one-out cross validation is performed.
#' @details This function does not use sparsity for any model.
#'
#' @return Vector with all predictions
#' @export
posterior.crossvalidation.covariance <- function(theta,
                                                 graph.obj,
                                                 model = "alpha1",
                                                 ind = NULL)
{
  n.o <- length(graph.obj$y)
  n.v <- dim(graph.obj$V)[1]
  #build covariance matrix
  if(model == "alpha1"){
    Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$edge_lengths)
    Sigma <- as.matrix(solve(Q))[graph.obj$PtV,graph.obj$PtV]
  } else if (model == "alpha2"){
    n.c <- 1:length(graph.obj$CBobj$S)
    Q <- Q.matern2(c(theta[3],theta[2]), graph.obj$V, graph.obj$EtV, graph.obj$edge_lengths, BC = 1)
    Qtilde <- (graph.obj$CBobj$T)%*%Q%*%t(graph.obj$CBobj$T)
    Qtilde <- Qtilde[-n.c,-n.c]
    Sigma.overdetermined  = t(graph.obj$CBobj$T[-n.c,])%*%solve(Qtilde)%*%(graph.obj$CBobj$T[-n.c,])
    index.obs <-  4*(graph.obj$PtE[,1]-1) + (1 * (graph.obj$PtE[,2]==0)) + (3 * (graph.obj$PtE[,2]!= 0))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
  } else if (model == "GL"){
    Q <- theta[2]*(theta[3]*Diagonal(graph.obj$nV,1) + graph.obj$Laplacian)
    Sigma <- as.matrix(solve(Q))[graph.obj$PtV,graph.obj$PtV]
  } else if (model == "GL2"){
    Q <- theta[2]*(theta[3]*Diagonal(graph.obj$nV,1) + graph.obj$Laplacian)
    Q <- Q%*%Q
    Sigma <- as.matrix(solve(Q))[graph.obj$PtV,graph.obj$PtV]
  } else if (model == "isoExp"){
    Sigma <- as.matrix(theta[3]^2*exp(-theta[2]*graph.obj$res.dist[graph.obj$PtV,graph.obj$PtV]))
  } else {
    stop("wrong model choice.")
  }
  Sigma.o <- Sigma
  diag(Sigma.o) <- diag(Sigma.o)  +  theta[1]^2
  if(is.null(ind)){
    ind <- 1:length(graph.obj$y)
  }
  mu.p <- var.p <- logscore <- crps <- scrps <- rep(0,length(graph.obj$y))
  mae <- rmse <- rep(0,length(graph.obj$y))
  for(j in 1:length(unique(ind))){
    i <- which(ind == j)
    mu.p[i] <-Sigma[i,-i]%*%solve(Sigma.o[-i,-i],graph.obj$y[-i])
    Sigma.p <- Sigma.o[i,i] - Sigma.o[i,-i]%*%solve(Sigma.o[-i,-i],Sigma.o[-i,i])
    var.p[i] <- diag(Sigma.p)
    logscore[i] <- LS(graph.obj$y[i], mu.p[i], sqrt(var.p[i]))
    crps[i] <- CRPS(graph.obj$y[i], mu.p[i], sqrt(var.p[i]))
    scrps[i] <- SCRPS(graph.obj$y[i], mu.p[i], sqrt(var.p[i]))
    mae[i] <- abs(graph.obj$y[i] - mu.p[i])
    rmse[i] <- (graph.obj$y[i] - mu.p[i])^2
  }
  return(list(mu = mu.p,
              var = var.p,
              logscore = mean(logscore),
              crps = mean(crps),
              scrps = mean(scrps),
              mae = mean(mae),
              rmse = mean(rmse)))
}


#' posterior mean calculation not using sparsity
#'
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph      - metric_graph object
#'
#' @return The posterior mean
#' @export
posterior.mean.stupid <- function(theta, graph){
  sigma_e <- theta[1]
  #build Q
  Q <- Q.exp(theta[2:3], graph)
  Sigma <- as.matrix(solve(Q))
  SigmaO <- Sigma[graph$PtV,graph$PtV]
  diag(SigmaO) <- diag(SigmaO)  +  sigma_e^2

  return( Sigma[,graph$PtV]%*%solve(SigmaO,graph$y))
}

#' prediction not using sparsity
#'
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph      - metric_graph object
#'
#' @return Leave-one-out predictions
#' @export
posterior.leave.stupid <- function(theta, graph){
  sigma_e <- theta[1]
  #build Q
  Q <- Q.exp(theta[2:3], graph)
  Sigma <- as.matrix(solve(Q))
  SigmaO <- Sigma[graph$PtV,graph$PtV]
  diag(SigmaO) <- diag(SigmaO)  +  sigma_e^2
  y_p <- rep(0,length(graph$y))
  for(i in 1:length(graph$y)){
    mu_p <-  Sigma[,graph$PtV[-i]]%*%solve(SigmaO[-i,-i],graph$y[-i])
    y_p[i] <- mu_p[graph$PtV[i]]
  }

  return( y_p)
}


#' Prediction for models assuming observations at vertices
#'
#' @param theta Estimated model parameters (sigma_e, sigma, kappa)
#' @param graph.obj Graph object
#' @param model Type of model: "alpha1" gives SPDE with alpha=1, "GL" gives
#' graph Laplacian, and "isoExp" gives isotropic exponential
#' @param ind Indices for cross validation. It should be a vector of the same length as the
#' data, with integer values representing each group in the cross-validation.
#' If NULL, leave-one-out cross validation is performed.
#'
#' @return Vector with all predictions
#' @export
posterior.crossvalidation <- function(theta,
                                      graph.obj,
                                      model = "alpha1",
                                      ind = NULL)
{
  #setup matrices for prediction
  sigma_e <- theta[1]
  if(model == "isoExp"){
    Sigma <- theta[3]^2*exp(-theta[2]*graph.obj$res.dist[graph.obj$PtV,graph.obj$PtV])
    Sigma.o <- Sigma
    diag(Sigma.o) <- diag(Sigma.o) + theta[1]^2
  } else if(model == "alpha2"){
    n.c <- 1:length(graph.obj$CBobj$S)
    Q <- Q.matern2(c(theta[3],theta[2]), graph.obj$V, graph.obj$EtV, graph.obj$edge_lengths, BC = 1)
    Qtilde <- (graph.obj$CBobj$T)%*%Q%*%t(graph.obj$CBobj$T)
    Qtilde <- Qtilde[-n.c,-n.c]
    Sigma.overdetermined  = t(graph.obj$CBobj$T[-n.c,])%*%solve(Qtilde)%*%(graph.obj$CBobj$T[-n.c,])
    index.obs <-  4*(graph.obj$PtE[,1]-1) + (1 * (graph.obj$PtE[,2]==0)) + (3 * (graph.obj$PtE[,2]!= 0))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
    Sigma.o <- Sigma
    diag(Sigma.o) <- diag(Sigma.o) + theta[1]^2
  } else if(model == "alpha1" || model == "GL"){
    if(model == "alpha1"){
      Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$edge_lengths)
    } else if(model == "GL"){
      Q <- theta[2]*(theta[3]*Diagonal(graph.obj$nV,1) + graph.obj$Laplacian)
    }
    A <- Diagonal(graph.obj$nV,rep(1,graph.obj$nV))[graph.obj$PtV,]
    Q.p <- Q  + t(A)%*%A/sigma_e^2
  } else {
    stop("Wrong model choice.")
  }

  if(is.null(ind)){
    ind <- 1:length(graph.obj$y)
  }
  y_p <- rep(0,length(graph.obj$y))
  for(j in 1:length(unique(ind))){
    i <- which(ind == j)
    if(model == "isoExp" || model == "alpha2"){
      y_p[i] <-Sigma[i,-i]%*%solve(Sigma.o[-i,-i],graph.obj$y[-i])
    } else {
      A <- Diagonal(graph.obj$nV,rep(1,graph.obj$nV))[graph.obj$PtV[-i],]
      mu.p <- solve(Q  + t(A)%*%A/sigma_e^2, as.vector(t(A)%*%graph.obj$y[-i]/sigma_e^2))
      y_p[i] <- mu.p[graph.obj$PtV[i]]
    }
  }
  return( y_p)
}

#' Log-likelihood calculation for graph Laplacian model
#'
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph.obj      - graph object
#'
#' @return The log-likelihood
#' @export
likelihood.graph_laplacian <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q

  Q <- theta[2]*(theta[3]*Diagonal(graph.obj$nV,1) + graph.obj$Laplacian)

  A <- Diagonal(graph.obj$nV,rep(1,graph.obj$nV))[graph.obj$PtV,]
  Q.p <- Q  + t(A)%*%A/sigma_e^2
  mu.p <- solve(Q.p,as.vector(t(A)%*%graph.obj$y/sigma_e^2))

  R <- chol(Q)
  R.p <- chol(Q.p)

  n.o <- length(graph.obj$y)
  l <- sum(log(diag(R))) - sum(log(diag(R.p))) - n.o*log(sigma_e)
  v <- graph$y - A%*%mu.p
  l <- l - 0.5*(t(mu.p)%*%Q%*%mu.p + t(v)%*%v/sigma_e^2) - 0.5 * n.o*log(2*pi)
  return(as.double(l))
}

#'
#' Computes the posterior expectation for each node in the graph
#' @param theta     - (sigma_e, sigma, kappa)
#' @param graph - metric_graph object
#' @param rem.edge  - remove edge
#' @export
posterior.mean.exp <- function(theta, graph, rem.edge=NULL){
  sigma_e <- theta[1]
  #build Q
  Qp.list <- Q.exp(theta[2:3], graph, build=FALSE)
  #build BSIGMAB
  Qpmu <- rep(0, graph$nV)

  obs.edges <- graph$nE
  if(is.null(rem.edge)==FALSE)
    obs.edges <- setdiff(obs.edges, rem.edge)
  i_ <- j_ <- x_ <- rep(0,4*length(obs.edges))
  count <- 0
  for(e in obs.edges){
    y_i <- graph.obj$y[e]
    l <- graph.obj$edge_lengths[e]
    D_matrix <- as.matrix(dist(c(0,l,graph.obj$PtE[e,2])))
    S <- r_1(D_matrix,theta[2:3])

    #covariance update see Art p.17
    E.ind         <- c(1:2)
    Obs.ind       <- -E.ind
    Bt            <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind])
    Sigma_i       <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
    R <- chol(Sigma_i)
    Sigma_iB      <- solve(Sigma_i, t(Bt))
    BtSinvB       <- Bt %*% Sigma_iB

    E <- graph$E[e,]
    if(E[1]==E[2]){
      Qpmu[E[1]]    <- Qpmu[E[1]]     +  sum(t(Sigma_iB)%*%y_i)
      Qp[E[1],E[1]] <-  Qp[E[1],E[1]] +  sum(Bt %*% Sigma_iB)
      i_[count+1] <- E[1]
      j_[count+1] <- E[1]
      x_[count+1] <- sum(Bt %*% Sigma_iB)
      count <- count + 1
    }else{
      i_[count+(1:4)] <- c(E[1],E[1],E[2],E[2])
      j_[count+(1:4)] <- c(E[1],E[2],E[1],E[2])
      x_[count+(1:4)] <- c(BtSinvB[1,1], BtSinvB[1,2], BtSinvB[1,2], BtSinvB[2,2] )
      count <- count + 4
      Qpmu[E]    <- Qpmu[E] + t(Sigma_iB)%*%y_i
    }
  }
  i_ <- c(Qp.list$i, i_[1:count])
  j_ <- c(Qp.list$j, j_[1:count])
  x_ <- c(Qp.list$x, x_[1:count])
  Qp <- Matrix::sparseMatrix(i= i_,
                             j= j_,
                             x= x_,
                             dims=Qp.list$dims)


  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)

  v <- c(as.matrix(Matrix::solve(R,Matrix::solve(R, Qpmu,system = 'P'), system='L')))
  Qpmu <- as.vector(Matrix::solve(R,Matrix::solve(R, v,system = 'Lt'), system='Pt'))

  return(Qpmu)

}


#'
#' Computes the posterior expectation for each observation in the graph
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph      - metric_graph object
#' @param leave.edge.out - compute the expectation of the graph if the observatrions are not on the edge
#' @export
posterior.mean.obs.exp <- function(theta, graph, leave.edge.out = FALSE){

  sigma_e = theta[1]

  Qp <- Q.exp(theta[2:3], graph)
  if(leave.edge.out==FALSE)
    V.post <- posterior.mean.exp(theta, graph)

  y_hat <- rep(0, length(graph$y))
  Qpmu <- rep(0, nrow(graph.obj$V))
  obs.edges <- unique(graph.obj$PtE[,1])

  for(e in obs.edges){

    if(leave.edge.out==TRUE)
      V.post <- posterior.mean.exp(theta, graph, rem.edge = e)

    obs.id <- graph.obj$PtE[,1] == e
    y_i <- graph.obj$y[obs.id]
    l <- graph.obj$edge_lengths[e]
    D_matrix <- as.matrix(dist(c(0,l,graph.obj$PtE[obs.id,2])))
    S <- r_1(D_matrix,theta[2:3])

    #covariance update see Art p.17
    E.ind         <- c(1:2)
    Obs.ind       <- -E.ind
    Bt            <- solve(S[E.ind, E.ind],S[E.ind, Obs.ind])

    E <- graph$E[e,]
    y_hat[obs.id] <- t(Bt)%*%V.post[E]
    if(leave.edge.out==F){
      Sigma_i       <- S[Obs.ind,Obs.ind] - S[Obs.ind, E.ind] %*% Bt
      Sigma_noise  <- Sigma_i
      diag(Sigma_noise) <- diag(Sigma_noise) + sigma_e^2

      y_hat[obs.id] <- y_hat[obs.id] +  Sigma_i%*%solve(Sigma_noise,y_i-y_hat[obs.id] )
    }
  }
  return(y_hat)
  #()
}

#' Log-likelihood calculation for alpha=1 without integration
#'
#' @param theta          - (sigma_e, sigma, kappa)
#' @param graph.obj      - graph object
#'
#' @return The log-likelihood
#' @export
likelihood.exp.graph.v2 <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q
  Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$edge_lengths)
  A <- Diagonal(graph.obj$nV,rep(1,graph.obj$nV))[graph.obj$PtV,]
  Q.p <- Q  + t(A)%*%A/sigma_e^2
  mu.p <- solve(Q.p,as.vector(t(A)%*%graph.obj$y/sigma_e^2))

  R <- chol(Q)
  R.p <- chol(Q.p)

  n.o <- length(graph.obj$y)
  l <- sum(log(diag(R))) - sum(log(diag(R.p))) - n.o*log(sigma_e)
  v <- graph.obj$y  - A%*%mu.p
  l <- l - 0.5*(t(mu.p)%*%Q%*%mu.p + t(v)%*%v/sigma_e^2) - 0.5 * n.o*log(2*pi)
  return(as.double(l))
}

#'
#' Computes the log likelihood function fo theta for the graph object
#' @param theta     - (sigma_e, sigma, kappa)
#' @param graph.obj - graph object
#' @export
likelihood.exp.graph <- function(theta, graph.obj){
  sigma_e <- theta[1]
  #build Q
  #Q <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$edge_lengths)

  Q.list <- Q.exp(theta[2:3], graph.obj$V, graph.obj$EtV, graph.obj$edge_lengths, build=F)

  Qp <- Matrix::sparseMatrix(i = Q.list$i,
                             j = Q.list$j,
                             x = Q.list$x,
                             dims=Q.list$dims)
  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)
  loglik <- Matrix::determinant(R)$modulus[1]
  #build BSIGMAB
  Qpmu <- rep(0, nrow(graph.obj$V))
  obs.edges <- unique(graph.obj$PtE[,1])

  i_ <- j_ <- x_ <- rep(0,4*length(obs.edges))

  count <- 0
  for(e in obs.edges){
    obs.id <- graph.obj$PtE[,1] == e
    y_i <- graph.obj$y[obs.id]
    l <- graph.obj$edge_lengths[e]
    D_matrix <- as.matrix(dist(c(0,l,graph.obj$PtE[obs.id,2])))
    S <- r_1(D_matrix,theta[2:3])

    #covariance update see Art p.17
    E.ind         <- c(1:2)
    Obs.ind       <- -E.ind
    Bt            <- solve(S[E.ind, E.ind,drop=F],S[E.ind, Obs.ind,drop=F])
    Sigma_i       <- S[Obs.ind,Obs.ind,drop=F] - S[Obs.ind, E.ind,drop=F] %*% Bt
    diag(Sigma_i) <- diag(Sigma_i) + sigma_e^2
    R <- chol(Sigma_i)
    Sigma_iB      <- solve(Sigma_i, t(Bt))
    BtSinvB       <- Bt %*% Sigma_iB

    E <- graph.obj$EtV[e,2:3]
    if(E[1]==E[2]){
      Qpmu[E[1]]    <- Qpmu[E[1]]     +  sum(t(Sigma_iB)%*%y_i)
      i_[count+1] <- E[1]
      j_[count+1] <- E[1]
      x_[count+1] <- sum(Bt %*% Sigma_iB)
    }else{
      Qpmu[E]    <- Qpmu[E] + t(Sigma_iB)%*%y_i
      i_[count+(1:4)] <- c(E[1],E[1],E[2],E[2])
      j_[count+(1:4)] <- c(E[1],E[2],E[1],E[2])
      x_[count+(1:4)] <- c(BtSinvB[1,1], BtSinvB[1,2], BtSinvB[1,2], BtSinvB[2,2] )
      count <- count + 4
    }

    loglik <- loglik - 0.5  * t(y_i)%*%solve(Sigma_i,y_i)
    loglik <- loglik - sum(log(diag(R)))
  }

  i_ <- c(Q.list$i, i_[1:count])
  j_ <- c(Q.list$j, j_[1:count])
  x_ <- c(Q.list$x, x_[1:count])
  Qp <- Matrix::sparseMatrix(i= i_,
                             j= j_,
                             x= x_,
                             dims=Q.list$dims)

  R <- Matrix::Cholesky(Qp, LDL = FALSE, perm = TRUE)

  loglik <- loglik - Matrix::determinant(R)$modulus[1]
  n.o <- length(graph.obj$y)
  v <- c(as.matrix(Matrix::solve(R,Matrix::solve(R, Qpmu,system = 'P'), system='L')))
  #Qpmu <- as.vector(solve(R,solve(R, v,system = 'Lt'), system='Pt'))

  loglik <- loglik + 0.5  * t(v)%*%v - 0.5 * n.o*log(2*pi)
  return(loglik[1])
}
