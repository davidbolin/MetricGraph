
#' Posterior mean for models assuming observations at vertices
#'
#' @param theta estimated model parameters (sigma_e, sigma, kappa)
#' @param graph metric graph  object
#' @param model Type of model: "alpha1" gives SPDE with alpha=1, "GL1" gives
#' the model based on the graph Laplacian with smoothness 1, "GL2" gives the
#' model based on the graph Laplacian with smoothness 2, and "isoExp" gives a
#' model with isotropic exponential covariance
#' @details This function does not use sparsity for any model.
#'
#' @return Vector with the posterior mean evaluated at the observation locations
#' @export
posterior_mean_covariance <- function(theta, graph, model = "alpha1")
{
  check <- gpgraph_check_graph(graph)

  if(is.null(graph$PtV) && (model != "isoExp")){
    stop("You must run graph$observation_to_vertex() first.")
  }

  n.o <- length(graph$y)
  n.v <- dim(graph$V)[1]

  sigma_e <- theta[1]
  sigma <- theta[2]
  kappa <- theta[3]

  #build covariance matrix
  if (model == "alpha1") {
    Q <- Qalpha1(c(sigma, kappa), graph)
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "alpha2") {
    n.c <- 1:length(graph$CBobj$S)
    Q <- Qalpha2(c(sigma, kappa), graph, BC = 1)
    Qtilde <- (graph$CBobj$T) %*% Q %*% t(graph$CBobj$T)
    Qtilde <- Qtilde[-n.c, -n.c]
    Sigma.overdetermined = t(graph$CBobj$T[-n.c, ]) %*%
      solve(Qtilde) %*% (graph$CBobj$T[-n.c, ])
    index.obs <- 4*(graph$PtE[, 1] - 1) + (1 * (graph$PtE[, 2] == 0)) +
      (3 * (graph$PtE[, 2] != 0))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
  } else if (model == "GL1"){
    Q <- (kappa^2 * Matrix::Diagonal(graph$nV, 1) + graph$Laplacian) / sigma^2
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "GL2"){
    K <- kappa^2*Matrix::Diagonal(graph$nV, 1) + graph$Laplacian
    Q <- K %*% K / sigma^2
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "isoExp"){
    Sigma <- as.matrix(sigma^2 * exp(-kappa*graph$res.dist))
  } else {
    stop("wrong model choice")
  }

  Sigma.o <- Sigma

  diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2

  return(as.vector(Sigma %*% solve(Sigma.o, graph$y)))
}

#' Prediction for models assuming observations at vertices
#'
#' @param theta Estimated model parameters (sigma_e, sigma, kappa)
#' @param graph metric graph object
#' @param model Type of model: "alpha1" gives SPDE with alpha=1, "GL1" gives
#' the model based on the graph Laplacian with smoothness 1, "GL2" gives the
#' model based on the graph Laplacian with smoothness 2, and "isoExp" gives a
#' model with isotropic exponential covariance
#' @param ind Indices for cross validation. It should be a vector of the same
#' length as the data, with integer values representing each group in the
#' cross-validation. If NULL, leave-one-out cross validation is performed.
#' @details This function does not use sparsity for any model.
#'
#' @return Vector with all predictions
#' @export
posterior.crossvalidation.covariance <- function(theta,
                                                 graph,
                                                 model = "alpha1",
                                                 ind = NULL)
{
  check <- gpgraph_check_graph(graph)

  if(is.null(graph$PtV)){
    stop("You must run graph$observation_to_vertex() first.")
  }

  sigma_e <- theta[1]
  sigma <- theta[2]
  kappa <- theta[3]

  #build covariance matrix
  if (model == "alpha1") {
    Q <- Qalpha1(c(sigma, kappa), graph)
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "alpha2") {
    n.c <- 1:length(graph$CBobj$S)
    Q <- Qalpha2(c(sigma, kappa), graph, BC = 1)
    Qtilde <- (graph$CBobj$T) %*% Q %*% t(graph$CBobj$T)
    Qtilde <- Qtilde[-n.c, -n.c]
    Sigma.overdetermined = t(graph$CBobj$T[-n.c, ]) %*%
      solve(Qtilde) %*% (graph$CBobj$T[-n.c, ])
    index.obs <- 4*(graph$PtE[, 1] - 1) + (1 * (graph$PtE[, 2] == 0)) +
      (3 * (graph$PtE[, 2] != 0))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
  } else if (model == "GL1"){
    Q <- (kappa^2 * Matrix::Diagonal(graph$nV, 1) + graph$Laplacian) / sigma^2
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "GL2"){
    K <- kappa^2*Matrix::Diagonal(graph$nV, 1) + graph$Laplacian
    Q <- K %*% K / sigma^2
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "isoExp"){
    Sigma <- as.matrix(sigma^2 * exp(-kappa*graph$res.dist))
  } else {
    stop("wrong model choice")
  }
  Sigma.o <- Sigma

  diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2

  if(is.null(ind)){
    ind <- 1:length(graph$y)
  }
  mu.p <- var.p <- logscore <- crps <- scrps <- rep(0, length(graph$y))
  mae <- rmse <- rep(0, length(graph$y))
  for (j in 1:length(unique(ind))) {
    i <- which(ind == j)
    mu.p[i] <-Sigma[i, -i] %*% solve(Sigma.o[-i, -i], graph$y[-i])
    Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                        Sigma.o[-i, i])
    var.p[i] <- diag(Sigma.p)
    logscore[i] <- LS(graph$y[i], mu.p[i], sqrt(var.p[i]))
    crps[i] <- CRPS(graph$y[i], mu.p[i], sqrt(var.p[i]))
    scrps[i] <- SCRPS(graph$y[i], mu.p[i], sqrt(var.p[i]))
    mae[i] <- abs(graph$y[i] - mu.p[i])
    rmse[i] <- (graph$y[i] - mu.p[i])^2
  }
  return(list(mu = mu.p,
              var = var.p,
              logscore = -mean(logscore),
              crps = -mean(crps),
              scrps = -mean(scrps),
              mae = mean(mae),
              rmse = mean(rmse)))
}


#' Prediction for models assuming observations at vertices
#'
#' @param theta Estimated model parameters (sigma_e, sigma, kappa)
#' @param graph metric graph object
#' @param model Type of model: "alpha1" gives SPDE with alpha=1, "GL1" gives
#' the model based on the graph Laplacian with smoothness 1, "GL2" gives the
#' model based on the graph Laplacian with smoothness 2, and "isoExp" gives a
#' model with isotropic exponential covariance
#' @param ind Indices for cross validation. It should be a vector of the same
#' length as the data, with integer values representing each group in the
#' cross-validation. If NULL, leave-one-out cross validation is performed.
#' @return Vector with the posterior expectations and variances as well as
#' mean absolute error (MAE), root mean squared errors (RMSE), and three
#' negatively oriented proper scoring rules: log-score, CRPS, and scaled
#' CRPS.
#' @export
posterior.crossvalidation <- function(theta,
                                      graph,
                                      model = "alpha1",
                                      ind = NULL)
{
  check <- gpgraph_check_graph(graph)
  if(is.null(graph$PtV)){
    stop("You must run graph$observation_to_vertex() first.")
  }
  sigma_e <- theta[1]
  sigma <- theta[2]
  kappa <- theta[3]

  #setup matrices for prediction
  if(model == "isoExp"){
    Sigma <- as.matrix(sigma^2 * exp(-kappa*graph$res.dist))
    Sigma.o <- Sigma
    diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2
  } else if(model == "alpha2"){
    n.c <- 1:length(graph$CBobj$S)
    Q <- Qalpha2(c(sigma, kappa), graph, BC = 1)
    Qtilde <- (graph$CBobj$T) %*% Q %*% t(graph$CBobj$T)
    Qtilde <- Qtilde[-n.c, -n.c]
    Sigma.overdetermined = t(graph$CBobj$T[-n.c, ]) %*%
      solve(Qtilde) %*% (graph$CBobj$T[-n.c, ])
    index.obs <- 4*(graph$PtE[, 1] - 1) + (1 * (graph$PtE[, 2] == 0)) +
      (3 * (graph$PtE[, 2] != 0))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
    Sigma.o <- Sigma
    diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2
  } else if (model %in% c("alpha1", "GL1", "GL2")) {
    if(model == "alpha1"){
      Q <- Qalpha1(c(sigma, kappa), graph)
    } else if(model == "GL1"){
      Q <- (kappa^2*Matrix::Diagonal(graph$nV,1) + graph$Laplacian) / sigma^2
    } else if (model == "GL2") {
      K <- (kappa^2*Matrix::Diagonal(graph$nV,1) + graph$Laplacian)
      Q <- K %*% K / sigma^2
    }
    A <- graph$A
    Q.p <- Q  + t(A)%*%A/sigma_e^2
  } else {
    stop("Wrong model choice.")
  }

  if(is.null(ind)){
    ind <- 1:length(graph$y)
  }

  mu.p <- var.p <- logscore <- crps <- scrps <- rep(0, length(graph$y))
  mae <- rmse <- rep(0, length(graph$y))

  for (j in 1:length(unique(ind))) {
    i <- which(ind == j)
    if(model == "isoExp" || model == "alpha2"){
      mu.p[i] <-Sigma[i,-i] %*% solve(Sigma.o[-i,-i], graph$y[-i])
      Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                          Sigma.o[-i, i])
      var.p[i] <- diag(Sigma.p)
    } else {
      A <- Matrix::Diagonal(graph$nV, rep(1, graph$nV))[graph$PtV[-i], ]
      Q.p <- Q + t(A) %*% A / sigma_e^2
      mu.p[i] <- solve(Q.p,
                       as.vector(t(A) %*% graph$y[-i] / sigma_e^2))[graph$PtV[i]]
      v <- rep(0,dim(Q.p)[1])
      v[graph$PtV[i]] <- 1
      var.p[i] <- solve(Q.p, v)[graph$PtV[i]] + sigma_e^2
    }
    logscore[i] <- LS(graph$y[i], mu.p[i], sqrt(var.p[i]))
    crps[i] <- CRPS(graph$y[i], mu.p[i], sqrt(var.p[i]))
    scrps[i] <- SCRPS(graph$y[i], mu.p[i], sqrt(var.p[i]))
    mae[i] <- abs(graph$y[i] - mu.p[i])
    rmse[i] <- (graph$y[i] - mu.p[i])^2
  }
  return(list(mu = mu.p,
              var = var.p,
              logscore = -mean(logscore),
              crps = -mean(crps),
              scrps = -mean(scrps),
              mae = mean(mae),
              rmse = mean(rmse)))
}



CRPS <- function(y, mu, sigma)
{
  return(-Exy(mu, sigma, y) + 0.5 * Exx(mu, sigma))
}

SCRPS <- function(y, mu, sigma)
{
  return(-Exy(mu, sigma, y) / Exx(mu, sigma) - 0.5 * log(Exx(mu, sigma)))
}

LS <- function(y, mu, sigma)
{
  return(dnorm(y, mean = mu, sd = sigma, log = TRUE))
}

Exx <- function(mu, sigma) {
  #X-X' = N(0,2*sigma^2)
  return(Efnorm(0, sqrt(2) * sigma))
}
#compute E[|X-y|] when X is N(mu,sigma^2)
Exy <- function(mu, sigma, y) {
  #X-y = N(mu-y,sigma^2)
  return(Efnorm(mu - y, sigma))
}

Efnorm <- function(mu, sigma) {
  return(sigma * sqrt(2 / pi) * exp(-(mu ^ 2) / (2 * sigma ^ 2)) + mu * (1 -
                                                                           2 * pnorm(-mu / sigma)))
}
