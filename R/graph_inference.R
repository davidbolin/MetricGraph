
#' Posterior mean for Gaussian random field models on metric graphs assuming
#' observations at vertices
#'
#' @param theta Estimated model parameters (sigma_e, tau, kappa).
#' @param graph A `metric_graph`  object.
#' @param model Type of model: "alpha1" gives SPDE with alpha=1, "GL1" gives
#' the model based on the graph Laplacian with smoothness 1, "GL2" gives the
#' model based on the graph Laplacian with smoothness 2, and "isoExp" gives a
#' model with isotropic exponential covariance
#' @details This function does not use sparsity for any model.
#'
#' @return Vector with the posterior mean evaluated at the observation locations
#' @noRd
posterior_mean_covariance <- function(theta, graph, model = "alpha1")
{
  check <- check_graph(graph)

  if(is.null(graph$PtV) && (model != "isoExp")){
    stop("You must run graph$observation_to_vertex() first.")
  }

  n.o <- length(graph$y)
  n.v <- dim(graph$V)[1]

  sigma_e <- theta[1]
  tau <- theta[2]
  kappa <- theta[3]

  #build covariance matrix
  if (model == "alpha1") {
    Q <- Qalpha1(c(tau, kappa), graph)
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "alpha2") {
    n.c <- 1:length(graph$CoB$S)
    Q <- Qalpha2(c(tau, kappa), graph, BC = 1)
    Qtilde <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
    Qtilde <- Qtilde[-n.c, -n.c]
    Sigma.overdetermined = t(graph$CoB$T[-n.c, ]) %*%
      solve(Qtilde) %*% (graph$CoB$T[-n.c, ])
    index.obs <- 4*(graph$PtE[, 1] - 1) + (1 * (graph$PtE[, 2] == 0)) +
      (3 * (graph$PtE[, 2] != 0))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
  } else if (model == "GL1"){
    Q <- (kappa^2 * Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]]) * tau^2
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "GL2"){
    K <- kappa^2*Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]]
    Q <- K %*% K * tau^2
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "isoExp"){
    Sigma <- as.matrix(tau^(-2) * exp(-kappa*graph$res_dist))
  } else {
    stop("wrong model choice")
  }

  Sigma.o <- Sigma

  diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2

  return(as.vector(Sigma %*% solve(Sigma.o, graph$y)))
}

#' Prediction for models assuming observations at vertices
#'
#' @param theta Estimated model parameters (sigma_e, tau, kappa).
#' @param graph A `metric_graph` object.
#' @param model Type of model: "alpha1" gives SPDE with alpha=1, "GL1" gives
#' the model based on the graph Laplacian with smoothness 1, "GL2" gives the
#' model based on the graph Laplacian with smoothness 2, and "isoExp" gives a
#' model with isotropic exponential covariance
#' @param data_name Name of the column of the response variable.
#' @param ind Indices for cross validation. It should be a vector of the same
#' length as the data, with integer values representing each group in the
#' cross-validation. If NULL, leave-one-out cross validation is performed.
#' @param BC Which boundary condition to use for the Whittle-Matern models (0,1).
#' Here 0 denotes Neumann boundary conditions and 1 stationary conditions.
#' @details This function does not use sparsity for any model.
#'
#' @return Vector with all predictions
#' @noRd
posterior_crossvalidation_covariance_manual <- function(theta,
                                                 graph,
                                                 data_name,
                                                 model = "alpha1",
                                                 ind = NULL,
                                                 BC=1)
{
  check <- check_graph(graph)

  if(is.null(graph$PtV)){
    stop("You must run graph$observation_to_vertex() first.")
  }

  sigma_e <- theta[1]
  tau <- theta[2]
  kappa <- theta[3]

  #build covariance matrix
  if (model == "alpha1") {
    Q <- Qalpha1(c(tau, kappa), graph)
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "alpha2") {

    graph$buildC(2,BC==0)
    n.c <- 1:length(graph$CoB$S)
    Q <- Qalpha2(c(tau, kappa), graph, BC = BC)
    Qtilde <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
    Qtilde <- Qtilde[-n.c, -n.c]
    Sigma.overdetermined = t(graph$CoB$T[-n.c, ]) %*%
      solve(Qtilde) %*% (graph$CoB$T[-n.c, ])
    PtE = graph$get_PtE()
    index.obs <- 4*(PtE[, 1] - 1) + (1 * (PtE[, 2] == 0)) +
      (3 * (PtE[, 2] != 0))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
  } else if (model == "GL1"){
    Q <- (kappa^2 * Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]]) * tau^2 # DOES NOT WORK FOR REPLICATES
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "GL2"){
    K <- kappa^2*Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]] # DOES NOT WORK FOR REPLICATES
    Q <- K %*% K * tau^2
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "isoExp"){
    Sigma <- as.matrix(tau^(-2) * exp(-kappa*graph$res_dist))
  } else {
    stop("wrong model choice")
  }
  Sigma.o <- Sigma

  diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2

  if(is.null(ind)){
    ind <- 1:length(graph$data[[data_name]])
  }
  mu.p <- var.p <- logscore <- crps <- scrps <- rep(0, length(graph$data[[data_name]]))
  mae <- rmse <- rep(0, length(graph$data[[data_name]]))
  for (j in 1:length(unique(ind))) {
    i <- which(ind == j)
    mu.p[i] <-Sigma[i, -i] %*% solve(Sigma.o[-i, -i], graph$data[[data_name]][-i])
    Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                        Sigma.o[-i, i])
    var.p[i] <- diag(Sigma.p)
    logscore[i] <- LS(graph$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    crps[i] <- CRPS(graph$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    scrps[i] <- SCRPS(graph$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    mae[i] <- abs(graph$data[[data_name]][i] - mu.p[i])
    rmse[i] <- (graph$data[[data_name]][i] - mu.p[i])^2
  }
  return(list(mu = mu.p,
              var = var.p,
              logscore = -mean(logscore),
              crps = -mean(crps),
              scrps = -mean(scrps),
              mae = mean(mae),
              rmse = sqrt(mean(rmse))))
}


#' Prediction for models assuming observations at vertices
#'
#' @param theta Estimated model parameters (sigma_e, tau, kappa).
#' @param graph A `metric_graph` object.
#' @param data_name Name of the data.
#' @param model Type of model: "alpha1" gives SPDE with alpha=1, "GL1" gives
#' the model based on the graph Laplacian with smoothness 1, "GL2" gives the
#' model based on the graph Laplacian with smoothness 2, and "isoExp" gives a
#' model with isotropic exponential covariance.
#' @param ind Indices for cross validation. It should be a vector of the same
#' length as the data, with integer values representing each group in the
#' cross-validation. If NULL, leave-one-out cross validation is performed.
#' @param BC Which boundary condition to use for the Whittle-Matern models (0,1).
#' Here 0 denotes Neumann boundary conditions and 1 stationary conditions.
#' @return Vector with the posterior expectations and variances as well as
#' mean absolute error (MAE), root mean squared errors (RMSE), and three
#' negatively oriented proper scoring rules: log-score, CRPS, and scaled
#' CRPS.
#' @noRd
posterior_crossvalidation_manual <- function(theta,
                                      graph,
                                      data_name,
                                      model = "alpha1",
                                      ind = NULL,
                                      BC = 1)
{
  check <- check_graph(graph)
  if(is.null(graph$PtV)){
    stop("You must run graph$observation_to_vertex() first.")
  }
  sigma_e <- theta[1]
  tau <- theta[2]
  kappa <- theta[3]

  #setup matrices for prediction
  if(model == "isoExp"){
    graph$compute_resdist()
    Sigma <- as.matrix(tau^(-2) * exp(-kappa*graph$res_dist[[1]]))
    Sigma.o <- Sigma
    diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2
  } else if(model == "alpha2"){
    n.c <- 1:length(graph$CoB$S)
    Q <- Qalpha2(c(tau, kappa), graph, BC = 1)
    Qtilde <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
    Qtilde <- Qtilde[-n.c, -n.c]
    Sigma.overdetermined = t(graph$CoB$T[-n.c, ]) %*%
      solve(Qtilde) %*% (graph$CoB$T[-n.c, ])
    PtE = graph$get_PtE()
    index.obs <- 4*(PtE[, 1] - 1) + (1 * (PtE[, 2] == 0)) +
      (3 * (PtE[, 2] != 0))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
    Sigma.o <- Sigma
    diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2
  } else if (model %in% c("alpha1", "GL1", "GL2")) {
    if(model == "alpha1"){
      Q <- Qalpha1(c(tau, kappa), graph, BC = BC)
    } else if(model == "GL1"){
      Q <- (kappa^2*Matrix::Diagonal(graph$nV,1) + graph$Laplacian[[1]]) * tau^2
    } else if (model == "GL2") {
      K <- (kappa^2*Matrix::Diagonal(graph$nV,1) + graph$Laplacian[[1]])
      Q <- K %*% K * tau^2
    }
    A <- graph$A()
    Q.p <- Q  + t(A)%*%A/sigma_e^2
  } else {
    stop("Wrong model choice.")
  }

  if(is.null(ind)){
    ind <- 1:length(graph$data[[data_name]])
  }

  mu.p <- var.p <- logscore <- crps <- scrps <- rep(0, length(graph$data[[data_name]]))
  mae <- rmse <- rep(0, length(graph$data[[data_name]]))

  for (j in 1:length(unique(ind))) {
    i <- which(ind == j)
    if(model == "isoExp" || model == "alpha2"){
      mu.p[i] <-Sigma[i,-i] %*% solve(Sigma.o[-i,-i], graph$data[[data_name]][-i])
      Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                          Sigma.o[-i, i])
      var.p[i] <- diag(Sigma.p)
    } else {
      A <- Matrix::Diagonal(graph$nV, rep(1, graph$nV))[graph$PtV[-i], ]
      Q.p <- Q + t(A) %*% A / sigma_e^2
      mu.p[i] <- solve(Q.p,
                       as.vector(t(A) %*% graph$data[[data_name]][-i] / sigma_e^2))[graph$PtV[i]]
      v <- rep(0,dim(Q.p)[1])
      v[graph$PtV[i]] <- 1
      var.p[i] <- solve(Q.p, v)[graph$PtV[i]] + sigma_e^2
    }
    logscore[i] <- LS(graph$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    crps[i] <- CRPS(graph$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    scrps[i] <- SCRPS(graph$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    mae[i] <- abs(graph$data[[data_name]][i] - mu.p[i])
    rmse[i] <- (graph$data[[data_name]][i] - mu.p[i])^2
  }
  return(list(mu = mu.p,
              var = var.p,
              logscore = -mean(logscore),
              crps = -mean(crps),
              scrps = -mean(scrps),
              mae = mean(mae),
              rmse = sqrt(mean(rmse))))
}


#' Leave-one-out crossvalidation for `graph_lme` models assuming observations at
#' the vertices of metric graphs
#'
#' @param object A fitted model using the `graph_lme()` function.
#' @return Vector with the posterior expectations and variances as well as
#' mean absolute error (MAE), root mean squared errors (RMSE), and three
#' negatively oriented proper scoring rules: log-score, CRPS, and scaled
#' CRPS.
#' @export
posterior_crossvalidation <- function(object)
{
  if(!inherits(object,"graph_lme")){
    stop("object should be of class graph_lme.")
  }

  graph <- object$graph$clone()

  graph$observation_to_vertex()

  beta_cov <- object$coeff$fixed_effects
  sigma_e <- object$coeff$measurement_error
  tau <- object$coeff$random_effects[1]
  kappa <- object$coeff$random_effects[2]

  # if(tolower(object$latent_model$type) == "whittlematern"){
  #   if(object$parameterization_latent == "matern"){
  #     kappa <- ifelse(object$latent_model$alpha == 1,
  #                 sqrt(8 * 0.5) / kappa, sqrt(8 * (1.5)) / kappa)
  #   }
  # }

  if(tolower(object$latent_model$type) == "isocov"){
    if(object$latent_model$cov_function_name == "other"){
      stop("Currently the cross-validation is only implemented for the exponential covariance function.")
    }
    model <- "isoExp"
    sigma <- object$coeff$random_effects[1]
  } else if(tolower(object$latent_model$type) == "whittlematern"){
    if(object$latent_model$alpha == 1){
      model <- "alpha1"
    } else {
      model <- "alpha2"
    }
  } else{
    graph$compute_laplacian(full=TRUE)
    if(object$latent_model$alpha == 1){
      model <- "GL1"
    } else {
      model <- "GL2"
    }
  }

  BC <- object$BC

  #setup matrices for prediction
  if(model == "isoExp"){
    graph$compute_resdist(full = TRUE)
    Sigma <- as.matrix(sigma^2 * exp(-kappa*graph$res_dist[[1]]))
    Sigma.o <- Sigma
    diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2
  } else if(model == "alpha2"){
    n.c <- 1:length(graph$CoB$S)
    Q <- Qalpha2(c(tau, kappa), graph, BC = BC)
    Qtilde <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
    Qtilde <- Qtilde[-n.c, -n.c]
    Sigma.overdetermined = t(graph$CoB$T[-n.c, ]) %*%
      solve(Qtilde) %*% (graph$CoB$T[-n.c, ])
    PtE = graph$get_PtE()
    index.obs <- 4*(PtE[, 1] - 1) + (1 * (PtE[, 2] == 0)) +
      (3 * (PtE[, 2] != 0))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
    Sigma.o <- Sigma
    diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2
  } else if (model %in% c("alpha1", "GL1", "GL2")) {
    if(model == "alpha1"){
      Q <- Qalpha1(c(tau, kappa), graph, BC = BC)
    } else if(model == "GL1"){
      Q <- (kappa^2*Matrix::Diagonal(graph$nV,1) + graph$Laplacian[[1]]) * tau^2
    } else if (model == "GL2") {
      K <- (kappa^2*Matrix::Diagonal(graph$nV,1) + graph$Laplacian[[1]])
      Q <- K %*% K * tau^2
    }
    A <- graph$A()
    Q.p <- Q  + t(A)%*%A/sigma_e^2
  } else {
    stop("Wrong model choice.")
  }

  if(!is.matrix(object$model_matrix)){
    object$model_matrix <- matrix(object$model_matrix, ncol=1)
  }

  y_graph <- object$model_matrix[,1]

  ind <- 1:length(y_graph)

  if(ncol(object$model_matrix) > 1){
    X_cov <- object$model_matrix[,-1]
  } else{
    X_cov <- NULL
  }


  mu.p <- var.p <- logscore <- crps <- scrps <- rep(0, length(y_graph))
  mae <- rmse <- rep(0, length(y_graph))

  for (j in 1:length(unique(ind))) {
    i <- which(ind == j)

    y_cv <- y_graph[-i]
    v_cv <- y_cv
    if(!is.null(X_cov)){
      v_cv <- v_cv - as.vector(X_cov[-i, ] %*% beta_cov)
      mu_fe <- as.vector(X_cov[i, ] %*% beta_cov)
    } else {
      mu_fe <- 0
    }

    if(model == "isoExp" || model == "alpha2"){
      mu.p[i] <-Sigma[i,-i] %*% solve(Sigma.o[-i,-i], v_cv) + mu_fe
      Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                          Sigma.o[-i, i])
      var.p[i] <- diag(Sigma.p)
    } else {
      A <- Matrix::Diagonal(graph$nV, rep(1, graph$nV))[graph$PtV[-i], ]
      Q.p <- Q + t(A) %*% A / sigma_e^2
      mu.p[i] <- solve(Q.p,
                       as.vector(t(A) %*% v_cv / sigma_e^2))[graph$PtV[i]] + mu_fe
      v <- rep(0,dim(Q.p)[1])
      v[graph$PtV[i]] <- 1
      var.p[i] <- solve(Q.p, v)[graph$PtV[i]] + sigma_e^2
    }
    logscore[i] <- LS(y_graph[i], mu.p[i], sqrt(var.p[i]))
    crps[i] <- CRPS(y_graph[i], mu.p[i], sqrt(var.p[i]))
    scrps[i] <- SCRPS(y_graph[i], mu.p[i], sqrt(var.p[i]))
    mae[i] <- abs(y_graph[i] - mu.p[i])
    rmse[i] <- (y_graph[i] - mu.p[i])^2
  }
  return(list(mu = mu.p,
              var = var.p,
              logscore = -mean(logscore),
              crps = -mean(crps),
              scrps = -mean(scrps),
              mae = mean(mae),
              rmse = sqrt(mean(rmse))))
}

#' Leave-one-out crossvalidation for `graph_lme` models assuming observations at
#' the vertices of metric graphs
#'
#' @param object A fitted model using the `graph_lme()` function.
#' @details This function does not use sparsity for any model.
#'
#' @return Vector with all predictions
#' @noRd
posterior_crossvalidation_covariance <- function(object)
{

  if(!inherits(object,"graph_lme")){
    stop("object should be of class graph_lme.")
  }

  graph <- object$graph$clone()

  graph$observation_to_vertex()

  check <- check_graph(graph)

  beta_cov <- object$coeff$fixed_effects
  sigma_e <- object$coeff$measurement_error
  tau <- object$coeff$random_effects[1]
  kappa <- object$coeff$random_effects[2]

  # if(tolower(object$latent_model$type) == "whittlematern"){
  #   if(object$parameterization_latent == "matern"){
  #     kappa <- ifelse(object$latent_model$alpha == 1,
  #                 sqrt(8 * 0.5) / kappa, sqrt(8 * (1.5)) / kappa)
  #   }
  # }

  if(tolower(object$latent_model$type) == "isocov"){
    if(object$latent_model$cov_function_name == "other"){
      stop("Currently the cross-validation is only implemented for the exponential covariance function.")
    }
    model <- "isoExp"
    sigma <- object$coeff$random_effects[1]
  } else if(tolower(object$latent_model$type) == "whittlematern"){
    if(object$latent_model$alpha == 1){
      model <- "alpha1"
    } else {
      model <- "alpha2"
    }
  } else{
    graph$compute_laplacian(full=TRUE)
    if(object$latent_model$alpha == 1){
      model <- "GL1"
    } else {
      model <- "GL2"
    }
  }

  BC <- object$BC

  if(!is.matrix(object$model_matrix)){
    object$model_matrix <- matrix(object$model_matrix, ncol=1)
  }

  y_graph <- object$model_matrix[,1]

  ind <- 1:length(y_graph)

  if(ncol(object$model_matrix) > 1){
    X_cov <- object$model_matrix[,-1]
  } else{
    X_cov <- NULL
  }

  #build covariance matrix
  if (model == "alpha1") {
    Q <- Qalpha1(c(tau, kappa), graph)
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "alpha2") {

    graph$buildC(2,BC==0)
    n.c <- 1:length(graph$CoB$S)
    Q <- Qalpha2(c(tau, kappa), graph, BC = BC)
    Qtilde <- (graph$CoB$T) %*% Q %*% t(graph$CoB$T)
    Qtilde <- Qtilde[-n.c, -n.c]
    Sigma.overdetermined = t(graph$CoB$T[-n.c, ]) %*%
      solve(Qtilde) %*% (graph$CoB$T[-n.c, ])
    PtE = graph$get_PtE()
    index.obs <- 4*(PtE[, 1] - 1) + (1 * (PtE[, 2] == 0)) +
      (3 * (PtE[, 2] != 0))
    Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
  } else if (model == "GL1"){
    Q <- (kappa^2 * Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]]) * tau^2 # DOES NOT WORK FOR REPLICATES
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "GL2"){
    K <- kappa^2*Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]] # DOES NOT WORK FOR REPLICATES
    Q <- K %*% K * tau^2
    Sigma <- as.matrix(solve(Q))[graph$PtV, graph$PtV]
  } else if (model == "isoExp"){
    graph$compute_resdist(full = TRUE)
    Sigma <- as.matrix(sigma^2 * exp(-kappa*graph$res_dist[[1]]))
  } else {
    stop("wrong model choice")
  }
  Sigma.o <- Sigma

  diag(Sigma.o) <- diag(Sigma.o) + sigma_e^2

  if(is.null(ind)){
    ind <- 1:length(y_graph)
  }
  mu.p <- var.p <- logscore <- crps <- scrps <- rep(0, length(y_graph))
  mae <- rmse <- rep(0, length(y_graph))
  for (j in 1:length(unique(ind))) {
    i <- which(ind == j)

    y_cv <- y_graph[-i]
    v_cv <- y_cv
    if(!is.null(X_cov)){
      v_cv <- v_cv - as.vector(X_cov[-i, ] %*% beta_cov)
      mu_fe <- as.vector(X_cov[i, ] %*% beta_cov)
    } else {
      mu_fe <- 0
    }

    mu.p[i] <-Sigma[i, -i] %*% solve(Sigma.o[-i, -i], v_cv) + mu_fe
    Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                        Sigma.o[-i, i])
    var.p[i] <- diag(Sigma.p)
    logscore[i] <- LS(y_graph[i], mu.p[i], sqrt(var.p[i]))
    crps[i] <- CRPS(y_graph[i], mu.p[i], sqrt(var.p[i]))
    scrps[i] <- SCRPS(y_graph[i], mu.p[i], sqrt(var.p[i]))
    mae[i] <- abs(y_graph[i] - mu.p[i])
    rmse[i] <- (y_graph[i] - mu.p[i])^2
  }
  return(list(mu = mu.p,
              var = var.p,
              logscore = -mean(logscore),
              crps = -mean(crps),
              scrps = -mean(scrps),
              mae = mean(mae),
              rmse = sqrt(mean(rmse))))
}


#' @noRd
CRPS <- function(y, mu, sigma)
{
  return(-Exy(mu, sigma, y) + 0.5 * Exx(mu, sigma))
}

#' @noRd
SCRPS <- function(y, mu, sigma)
{
  return(-Exy(mu, sigma, y) / Exx(mu, sigma) - 0.5 * log(Exx(mu, sigma)))
}


#' @noRd
LS <- function(y, mu, sigma)
{
  return(dnorm(y, mean = mu, sd = sigma, log = TRUE))
}


#' @noRd
Exx <- function(mu, sigma) {
  #X-X' = N(0,2*sigma^2)
  return(Efnorm(0, sqrt(2) * sigma))
}



#compute E[|X-y|] when X is N(mu,sigma^2)

#' @noRd
Exy <- function(mu, sigma, y) {
  #X-y = N(mu-y,sigma^2)
  return(Efnorm(mu - y, sigma))
}

#' @noRd
Efnorm <- function(mu, sigma) {
  return(sigma * sqrt(2 / pi) * exp(-(mu ^ 2) / (2 * sigma ^ 2)) + mu * (1 -
                                                                           2 * pnorm(-mu / sigma)))
}
