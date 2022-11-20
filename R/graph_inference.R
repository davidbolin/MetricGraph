
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
  check <- gpgraph_check_graph(graph)

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
  check <- gpgraph_check_graph(graph)

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


