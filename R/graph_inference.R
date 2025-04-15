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
    ind <- 1:length(graph$.__enclos_env__$private$data[[data_name]])
  }
  mu.p <- var.p <- logscore <- crps <- scrps <- rep(0, length(graph$.__enclos_env__$private$data[[data_name]]))
  mae <- rmse <- rep(0, length(graph$.__enclos_env__$private$data[[data_name]]))
  for (j in 1:length(unique(ind))) {
    i <- which(ind == j)
    mu.p[i] <-Sigma[i, -i] %*% solve(Sigma.o[-i, -i], graph$.__enclos_env__$private$data[[data_name]][-i])
    Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                        Sigma.o[-i, i])
    var.p[i] <- diag(Sigma.p)
    logscore[i] <- LS(graph$.__enclos_env__$private$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    crps[i] <- CRPS(graph$.__enclos_env__$private$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    scrps[i] <- SCRPS(graph$.__enclos_env__$private$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    mae[i] <- abs(graph$.__enclos_env__$private$data[[data_name]][i] - mu.p[i])
    rmse[i] <- (graph$.__enclos_env__$private$data[[data_name]][i] - mu.p[i])^2
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
    Sigma <- as.matrix(tau^(-2) * exp(-kappa*graph$res_dist[[".complete"]]))
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
  } else {
    stop("Wrong model choice.")
  }

  if(is.null(ind)){
    ind <- 1:length(graph$.__enclos_env__$private$data[[data_name]])
  }

  mu.p <- var.p <- logscore <- crps <- scrps <- rep(0, length(graph$.__enclos_env__$private$data[[data_name]]))
  mae <- rmse <- rep(0, length(graph$.__enclos_env__$private$data[[data_name]]))

  for (j in 1:length(unique(ind))) {
    i <- which(ind == j)
    if(model == "isoExp" || model == "alpha2"){
      mu.p[i] <-Sigma[i,-i] %*% solve(Sigma.o[-i,-i], graph$.__enclos_env__$private$data[[data_name]][-i])
      Sigma.p <- Sigma.o[i, i] - Sigma.o[i, -i] %*% solve(Sigma.o[-i, -i],
                                                          Sigma.o[-i, i])
      var.p[i] <- diag(Sigma.p)
    } else {
      A <- Matrix::Diagonal(graph$nV, rep(1, graph$nV))[graph$PtV[-i], ]
      Q.p <- Q + t(A) %*% A / sigma_e^2
      mu.p[i] <- solve(Q.p,
                       as.vector(t(A) %*% graph$.__enclos_env__$private$data[[data_name]][-i] / sigma_e^2))[graph$PtV[i]]
      v <- rep(0,dim(Q.p)[1])
      v[graph$PtV[i]] <- 1
      var.p[i] <- solve(Q.p, v)[graph$PtV[i]] + sigma_e^2
    }
    logscore[i] <- LS(graph$.__enclos_env__$private$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    crps[i] <- CRPS(graph$.__enclos_env__$private$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    scrps[i] <- SCRPS(graph$.__enclos_env__$private$data[[data_name]][i], mu.p[i], sqrt(var.p[i]))
    mae[i] <- abs(graph$.__enclos_env__$private$data[[data_name]][i] - mu.p[i])
    rmse[i] <- (graph$.__enclos_env__$private$data[[data_name]][i] - mu.p[i])^2
  }
  return(list(mu = mu.p,
              var = var.p,
              logscore = -mean(logscore),
              crps = -mean(crps),
              scrps = -mean(scrps),
              mae = mean(mae),
              rmse = sqrt(mean(rmse))))
}

#' Cross-validation for `graph_lme` models assuming observations at
#' the vertices of metric graphs
#'
#' This function performs cross-validation by computing predictions for test data
#' using either the posterior distribution from a fitted model (pseudo-CV) or by
#' refitting the model for each fold (true CV).
#'
#' @param object A fitted model using the `graph_lme()` function or a named list of fitted objects using the `graph_lme()` function.
#' @param scores A vector of scores to compute. The options are "logscore", "crps", "scrps", "mae", and "rmse". By default, all scores are computed.
#' @param mode Cross-validation mode. Options are "k-fold", "loo" (leave-one-out), or "lpo" (leave-percentage-out). Default is "k-fold".
#' @param k Number of folds for k-fold cross-validation. Default is 10.
#' @param percentage The percentage (from 1 to 99) of the data to be used to train the model. Will only be used if `mode` is "lpo". Default is 80.
#' @param number_folds Number of folds to be done if `mode` is "lpo". Default is 10.
#' @param train_test_indices Optional list containing train and test indices for each fold. If provided, k, mode, and percentage are ignored.
#' @param true_CV Logical indicating whether to refit the model for each fold (TRUE) or use the posterior distribution from the fitted model (FALSE). Default is FALSE.
#' @param factor Which factor to multiply the scores. The default is 1.
#' @param tibble Return the scores as a `tidyr::tibble()`
#' @param parallel_folds Logical indicating whether to run computations in parallel across folds. Default is FALSE.
#' @param parallel_fitting Logical indicating whether to run model fitting in parallel. Default is FALSE.
#' @param n_cores Number of cores to use for parallel computation. Default is parallel::detectCores() - 1.
#' @param print Logical indicating whether to print progress of which fold is being processed. Default is FALSE.
#' @param seed Random seed for reproducibility in fold creation. Default is NULL.
#' @param return_indices Logical indicating whether to return the train/test indices used. Default is FALSE.
#' @param use_precomputed Logical indicating whether to use precomputation for faster CV. Default is TRUE.
#' @return Vector with the posterior expectations and variances as well as
#' mean absolute error (MAE), root mean squared errors (RMSE), and three
#' negatively oriented proper scoring rules: log-score, CRPS, and scaled
#' CRPS.
#' @export
posterior_crossvalidation <- function(object, scores = c("logscore", "crps", "scrps", "mae", "rmse"), mode = "k-fold", k = 10, percentage = 20, number_folds = 10,
                                     train_test_indices = NULL, true_CV = FALSE, factor = 1, tibble = TRUE, 
                                     parallel_folds = FALSE, parallel_fitting = FALSE, n_cores = parallel::detectCores() - 1, 
                                     print = FALSE, seed = NULL, return_indices = FALSE, use_precomputed = TRUE)
{
  
  if(!inherits(object,"graph_lme") && !is.list(object)){
    stop("object should be of class graph_lme or a list of objects of class graph_lme.")
  }
  
  # Convert scores to lowercase for robustness
  scores <- tolower(scores)
  
  # Check if scores are valid
  valid_scores <- c("logscore", "crps", "scrps", "mae", "rmse")
  invalid_scores <- setdiff(scores, valid_scores)

  # Pre-clean the data to only include non-NA values for the response variable
  if(inherits(object, "graph_lme")) {
    # Get the response variable name
    response_var <- as.character(object$response_var)
    
    # Access the private data from the graph object
    data <- object$graph$.__enclos_env__$private$data
    
    # Check if the response variable exists in the data
    if(response_var %in% names(data)) {
      # Filter out NA values in the response variable
      non_na_indices <- !is.na(data[[response_var]])
      
      # Update all columns in the data to only include non-NA response values
      data <- lapply(data, function(col) col[non_na_indices])
      
      # Update the private data in the graph object
      object$graph$.__enclos_env__$private$data <- data
      
      # Update the number of observations in the object
      object$nobs <- sum(non_na_indices)
    } else {
      warning(paste("Response variable", response_var, "not found in the data."))
    }
  }
  
  if(length(invalid_scores) > 0) {
    warning(paste("Invalid scores:", paste(invalid_scores, collapse=", "), 
                  "- will be ignored. Valid options are:", paste(valid_scores, collapse=", ")))
    scores <- intersect(scores, valid_scores)
  }

  if(parallel_folds && parallel_fitting){
    stop("parallel_folds and parallel_fitting cannot both be TRUE.")
  }

  if(!inherits(object,"graph_lme")){
    if(is.null(names(object))){
      warning("The list with fitted models does not contain names for the models, thus the results will not be properly named.")
    }
    results_list <- lapply(object, function(obj){
      posterior_crossvalidation(obj, k = k, mode = mode, percentage = percentage, 
                               number_folds = number_folds,
                               train_test_indices = train_test_indices, 
                               true_CV = true_CV, factor = factor, 
                               tibble = FALSE, parallel_folds = parallel_folds, 
                               parallel_fitting = parallel_fitting,
                               n_cores = n_cores, print = print, 
                               seed = seed, return_indices = return_indices,
                               scores = scores)
    })
    res <- list()
    res[["mu"]] <- lapply(results_list, function(dat){dat[["mu"]]})
    res[["var"]] <- lapply(results_list, function(dat){dat[["var"]]})
    res[["scores"]] <- lapply(results_list, function(dat){dat[["scores"]]})
    scores <- do.call(rbind, res[["scores"]])
    row.names(scores) <- names(res[["scores"]])
    if(tibble){
      scores[["Model"]] <- row.names(scores)
      scores <- tidyr::as_tibble(scores)
      scores <- scores[, c("Model", "logscore", "crps", "scrps", "mae", "rmse")]
    }
    res[["scores"]] <- scores
    if(return_indices && !is.null(results_list[[1]][["indices"]])) {
      res[["indices"]] <- results_list[[1]][["indices"]]
    }
    return(res)
  }

  # Create or use provided train/test indices
  if(is.null(train_test_indices)) {
    if(!is.null(seed)) {
      set.seed(seed)
    }
    
    # Get number of observations
    n_obs <- object$nobs
    
    # Create indices for cross-validation
    if(mode == "loo") {
      # Leave-one-out: each observation is a fold
      train_test_indices <- lapply(1:n_obs, function(i) {
        list(train = setdiff(1:n_obs, i), test = i)
      })
      k <- n_obs
    } else if(mode == "k-fold") {
      # k-fold cross-validation
      fold_indices <- sample(rep(1:k, length.out = n_obs))
      train_test_indices <- lapply(1:k, function(i) {
        test <- which(fold_indices == i)
        train <- setdiff(1:n_obs, test)
        list(train = train, test = test)
      })
    } else if(mode == "lpo") {
      # Leave-percentage-out cross-validation
      if(percentage <= 0 || percentage >= 100) {
        stop("percentage must be a number between 1 and 99!")
      }
      
      # Calculate number of observations for training
      n_train <- round(n_obs * percentage / 100)
      
      # Create multiple folds
      train_test_indices <- lapply(1:number_folds, function(i) {
        # Randomly select indices for training
        train <- sample(1:n_obs, n_train)
        test <- setdiff(1:n_obs, train)
        list(train = train, test = test)
      })
      k <- number_folds
    } else {
      stop("Invalid mode. Choose 'k-fold', 'loo', or 'lpo'")
    }
  } else {
    # Use provided indices
    k <- length(train_test_indices)
  }

  # Initialize result vectors
  n_obs <- object$nobs
  mu.p <- var.p <- rep(NA, n_obs)
  logscore <- crps <- scrps <- mae <- rmse <- rep(NA, n_obs)
  
  # Different precomputation strategy based on CV type
  need_variances <- any(scores %in% c("logscore", "crps", "scrps"))
  precomputed_data <- NULL
  precomputed_graph <- NULL
  
  if(print) {
    cat("Pre-computing data structures for faster cross-validation...\n")
  }
  
  if(use_precomputed) {
    if(!true_CV) {
      # For pseudo-CV: precompute parameter-dependent structures (Q, Sigma, etc.)
      # Get a dummy test index to create precomputed structures
      dummy_test_idx <- train_test_indices[[1]]$test[1]
      
      # Create dummy prediction data for a single point
      dummy_data <- object$graph$.__enclos_env__$private$data
      dummy_data <- lapply(dummy_data, function(x){x[dummy_test_idx]})
      
      # Make initial prediction to get precomputed data
      precomp_pred <- predict(object, 
                            newdata = dummy_data, 
                            edge_number = ".edge_number", 
                            distance_on_edge = ".distance_on_edge", 
                            normalized = TRUE, 
                            compute_variances = need_variances,
                            advanced_options = list(precompute_data = TRUE, precompute_type = "full"))
                            
      # Extract precomputed data
      precomputed_data <- precomp_pred$precomputed_data

      precomputed_graph <- precomputed_data$graph_bkp
      
      if(is.null(precomputed_data)) {
        warning("Precomputation for pseudo-CV failed, falling back to standard prediction")
      } else if(print) {
        cat("Precomputation for pseudo-CV successful.\n")
      }
    } else {
      # For true-CV: just precompute the graph structure
      # Create a dummy test index to create precomputed structures
      dummy_test_idx <- train_test_indices[[1]]$test[1]
      
      # Create dummy prediction data for a single point
      dummy_data <- object$graph$.__enclos_env__$private$data
      dummy_data <- lapply(dummy_data, function(x){x[dummy_test_idx]})
      
      # Make initial prediction to get structure-only precomputed data
      precomp_pred <- predict(object, 
                            newdata = dummy_data, 
                            edge_number = ".edge_number", 
                            distance_on_edge = ".distance_on_edge", 
                            normalized = TRUE, 
                            compute_variances = need_variances,
                            advanced_options = list(precompute_data = TRUE, precompute_type = "structure"))
                            
      # Extract precomputed data
      precomputed_data <- precomp_pred$precomputed_data
      
      if(is.null(precomputed_data)) {
        warning("Precomputation for true-CV failed, falling back to standard prediction")
      } else if(print) {
        cat("Precomputation for true-CV successful.\n")
      }
      
      # Extract just the graph for use
      precomputed_graph <- precomputed_data$graph_bkp
    }
  # Reorder train and test indices based on the ordering from precomputed_graph
  if (!is.null(precomputed_graph) && !is.null(precomputed_graph$.__enclos_env__$private$data[[".dummy_order_var"]])) {
    # Get the ordering vector
    order_var <- precomputed_graph$.__enclos_env__$private$data[[".dummy_order_var"]]
    
    # Apply the reordering to train_test_indices
    for (i in 1:length(train_test_indices)) {
      # Create temporary copies
      train_orig <- train_test_indices[[i]]$train
      test_orig <- train_test_indices[[i]]$test
            
      # Reorder using the ordering vector
      # Find the positions in order_var that match the original indices
      train_test_indices[[i]]$reordered_train <- match(train_orig, order_var)
      train_test_indices[[i]]$reordered_test <- match(test_orig, order_var)

    }
    
    if (print) {
        cat("Train and test indices reordered according to graph structure.\n")
      }
    }
  }

  # Function to process a single fold
  process_fold <- function(fold_idx) {
    if(print && !parallel_folds) {
      cat("Processing fold", fold_idx, "of", k, "\n")
    }
    
    train_indices <- train_test_indices[[fold_idx]]$train
    test_indices <- train_test_indices[[fold_idx]]$test

    if(use_precomputed){
      reordered_train_indices <- train_test_indices[[fold_idx]]$reordered_train
      reordered_test_indices <- train_test_indices[[fold_idx]]$reordered_test
    } 
    
    local_results <- list(
      test_indices = test_indices,
      mu.p = rep(NA, length(test_indices)),
      var.p = rep(NA, length(test_indices)),
      logscore = rep(NA, length(test_indices)),
      crps = rep(NA, length(test_indices)),
      scrps = rep(NA, length(test_indices)),
      mae = rep(NA, length(test_indices)),
      rmse = rep(NA, length(test_indices))
    )

    response_name <- as.character(object$response_var)
    
    if(true_CV) {
      # True cross-validation: refit the model for each fold
      # Create a copy of the graph with NA values for test indices
      cv_graph <- NULL
      
      cv_graph <- object$graph$clone()
      cv_graph$.__enclos_env__$private$data[[response_name]][test_indices] <- NA
      
      model_options <- object$options_model
      # Set starting values based on the original model parameters
      if(tolower(object$latent_model$type) %in% c("whittlematern", "graphlaplacian")) {
        # For WhittleMatern or graphLaplacian models
        if(!is.null(object$coeff$random_effects)) {
          model_options$start_kappa <- object$coeff$random_effects[2]
          model_options$start_tau <- object$coeff$random_effects[1]
        }
        if(!is.null(object$coeff$measurement_error)) {
          model_options$start_sigma_e <- object$coeff$measurement_error
        }
      } else if(tolower(object$latent_model$type) == "isocov") {
        # For isoCov models
        if(!is.null(object$coeff$random_effects)) {
          model_options$start_par_vec <- object$coeff$random_effects
        }
        if(!is.null(object$coeff$measurement_error)) {
          model_options$start_sigma_e <- object$coeff$measurement_error
        }
      }

      # Refit the model using the same optimization method
      cv_model <- graph_lme(
        formula = object$formula,
        graph = cv_graph,
        model = object$latent_model,
        BC = object$BC,
        optim_method = object$optim_method,
        optim_controls = object$optim_controls,
        model_options = model_options,
        parallel = parallel_fitting,
        n_cores = n_cores
      )

    } else {
      # Pseudo cross-validation: use the original model parameters
      # Create a copy of the object with NA values for test indices
      cv_model <- update_graph_lme_with_na(object, test_indices)
    }

    if(!use_precomputed) {
      new_data <- object$graph$.__enclos_env__$private$data
      new_data <- lapply(new_data, function(x){x[test_indices]})
    }


    # Make predictions for test indices
    if(!is.null(precomputed_data)) {
      # Use precomputed data structures
      if(need_variances){
        pred <- predict(cv_model, 
                      edge_number = ".edge_number", 
                      distance_on_edge = ".distance_on_edge", 
                      normalized = TRUE, 
                      compute_variances = TRUE, 
                      advanced_options = list(precompute_data = precomputed_data, precompute_type = "full", test_idx = reordered_test_indices, train_idx = reordered_train_indices))
      } else {
        pred <- predict(cv_model, 
                      edge_number = ".edge_number", 
                      distance_on_edge = ".distance_on_edge", 
                      normalized = TRUE, 
                      compute_variances = FALSE,
                      advanced_options = list(precompute_data = precomputed_data, precompute_type = "full", test_idx = reordered_test_indices, train_idx = reordered_train_indices))
      }
    } else {
      # Fallback if precomputation failed
      if(need_variances){
        pred <- predict(cv_model, 
                      newdata = new_data, 
                      edge_number = ".edge_number", 
                      distance_on_edge = ".distance_on_edge", 
                      normalized = TRUE, 
                      compute_variances = TRUE)
      } else {
        pred <- predict(cv_model, 
                      newdata = new_data, 
                      edge_number = ".edge_number", 
                      distance_on_edge = ".distance_on_edge", 
                      normalized = TRUE, 
                      compute_variances = FALSE)
      }
    }
      
    # Extract predictions for test indices
    for(i in seq_along(test_indices)) {
      idx <- test_indices[i]
      local_results$mu.p[i] <- pred$mean[i]
      if(need_variances) {
        local_results$var.p[i] <- pred$variance[i] + object$coeff$measurement_error^2
      }
      y_test <- object$graph$.__enclos_env__$private$data[[response_name]][idx]
      
      # Calculate scores
      if("logscore" %in% scores) {
        local_results$logscore[i] <- LS(y_test, local_results$mu.p[i], sqrt(local_results$var.p[i]))
      }
      if("crps" %in% scores) {
        local_results$crps[i] <- CRPS(y_test, local_results$mu.p[i], sqrt(local_results$var.p[i]))
      }
      if("scrps" %in% scores) {
        local_results$scrps[i] <- SCRPS(y_test, local_results$mu.p[i], sqrt(local_results$var.p[i]))
      }
      if("mae" %in% scores) {
        local_results$mae[i] <- abs(y_test - local_results$mu.p[i])
      }
      if("rmse" %in% scores) {
        local_results$rmse[i] <- (y_test - local_results$mu.p[i])^2
      }
    }
    
    return(local_results)
  }

  # Process folds in parallel or sequentially
  if(parallel_folds) {
    if(!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' is not available. Using sequential processing instead.")
      parallel_folds <- FALSE
    }
  }
  
  if(parallel_folds) {
    cl <- parallel::makeCluster(n_cores)
    
    # Load necessary packages and functions in each worker
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(Matrix))
      
      # Define helper functions inside each worker
      LS <- function(y, mu, sigma) {
        return(dnorm(y, mean = mu, sd = sigma, log = TRUE))
      }
      
      CRPS <- function(y, mu, sigma) {
        return(-Exy(mu, sigma, y) + 0.5 * Exx(mu, sigma))
      }
      
      SCRPS <- function(y, mu, sigma) {
        return(-Exy(mu, sigma, y) / Exx(mu, sigma) - 0.5 * log(Exx(mu, sigma)))
      }
      
      Exx <- function(mu, sigma) {
        return(Efnorm(0, sqrt(2) * sigma))
      }
      
      Exy <- function(mu, sigma, y) {
        return(Efnorm(mu - y, sigma))
      }
      
      Efnorm <- function(mu, sigma) {
        return(sigma * sqrt(2 / pi) * exp(-(mu ^ 2) / (2 * sigma ^ 2)) + mu * (1 - 2 * pnorm(-mu / sigma)))
      }
      
      NULL
    })
    
    # Export necessary objects to the workers
    parallel::clusterExport(cl, 
                          varlist = c("object", "train_test_indices", 
                                     "update_graph_lme_with_na", "graph_lme", "predict", 
                                     "print", "true_CV", "precomputed_data", "precomputed_graph",
                                     "need_variances", "scores", "use_precomputed"),
                          envir = environment())
    
    if(print) {
      cat("Starting parallel processing with", n_cores, "cores\n")
    }
    
    results <- parallel::parLapply(cl, 1:k, process_fold)
    parallel::stopCluster(cl)
  } else {
    results <- lapply(1:k, process_fold)
  }
  
  # Combine results from all folds
  for(i in 1:k) {
    test_indices <- results[[i]]$test_indices
    for(j in seq_along(test_indices)) {
      idx <- test_indices[j]
      mu.p[idx] <- results[[i]]$mu.p[j]
      var.p[idx] <- results[[i]]$var.p[j]
      logscore[idx] <- results[[i]]$logscore[j]
      crps[idx] <- results[[i]]$crps[j]
      scrps[idx] <- results[[i]]$scrps[j]
      mae[idx] <- results[[i]]$mae[j]
      rmse[idx] <- results[[i]]$rmse[j]
    }
  }
  
  # Compile final results
  res <- list(
    mu = mu.p,
    var = var.p
  )
  
  # Create a data frame with only the requested scores
  # Initialize scores data frame with at least one row to avoid the error
  # "replacement has 1 row, data has 0"
  res[["scores"]] <- data.frame(dummy = NA)
  
  if("logscore" %in% scores) {
    res[["scores"]]$logscore <- -factor * mean(logscore, na.rm = TRUE)
  }
  if("crps" %in% scores) {
    res[["scores"]]$crps <- -factor * mean(crps, na.rm = TRUE)
  }
  if("scrps" %in% scores) {
    res[["scores"]]$scrps <- -factor * mean(scrps, na.rm = TRUE)
  }
  if("mae" %in% scores) {
    res[["scores"]]$mae <- factor * mean(mae, na.rm = TRUE)
  }
  if("rmse" %in% scores) {
    res[["scores"]]$rmse <- factor * sqrt(mean(rmse, na.rm = TRUE))
  }
  
  # Remove the dummy column after adding the real scores
  res[["scores"]]$dummy <- NULL
  
  attr(res[["scores"]], "factor") <- factor
  
  if(return_indices) {
    res[["indices"]] <- train_test_indices
  }
  
  if(tibble){
    res[["scores"]] <- tidyr::as_tibble(res[["scores"]])
  }
  
  return(res)
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
      model <- "WM alpha2"
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
  } else if (model == "WM alpha2") {

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
    Sigma <- as.matrix(sigma^2 * exp(-kappa*graph$res_dist[[".complete"]]))
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
      v_cv <- v_cv - as.vector(X_cov[-i, , drop = FALSE] %*% beta_cov)
      mu_fe <- as.vector(X_cov[i, , drop = FALSE] %*% beta_cov)
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
