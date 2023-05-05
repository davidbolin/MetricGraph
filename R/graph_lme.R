#' Metric graph linear mixed effects models
#'
#' Fitting linear mixed effects model in metric graphs. The random effects can be
#' Whittle-Matern fields on metric graphs, Whittle-Matern fields based on the
#' graph Laplacian as well as fields with isotropic covariance structure.
#'
#' @param formula Formula object describing the relation between the response variables and the fixed effects.
#' @param graph A `metric_graph` object.
#' @param model The random effects model that will be used. A list containing the elements `type`, which can be
#' `WhittleMatern`, `graphLaplacian` or `isoCov`. For `Whittle-Matern` and `graph-Laplacian` models, the list must also contain a parameter `alpha` (which is 1 by default). For `isoCov` models, the list must 
#' contain a parameter `cov_function`, containing the covariance function. The function accepts a string input for the following covariance functions: 'exp_covariance', 'alpha1', 'alpha2', 'GL1', 'GL2'. For another covariance function, the function itself must be provided as the `cov_function` argument. The default is 'exp_covariance', the
#' exponential covariance. We also have covariance-based versions of the Whittle-Matern and graph Laplacian models, however they are much slower, they are the following (string) values for 'cov_function': 'alpha1' and 'alpha2' for Whittle-Matern fields, and 'GL1' and 'GL2' for graph Laplacian models. Finally, for `Whittle-Matern` models, there is an additional parameter
#' `version`, which can be either 1 or 2, to tell which version of the likelihood should be used. Version is 1 by default. 
#' @param repl vector or list containing which replicates to consider in the model.
#' If `NULL` all replicates will be considered.
#' @param optim_method The method to be used with `optim` function.
#' @param starting_values_latent A vector containing the starting values for the latent model. If the latent model is `WhittleMatern` or `graphLaplacian`, then the starting values should be provided as a vector of the form c(sigma,kappa) or c(sigma,range) depending on the parameterization. If the model is `isoCov`, then the starting values should be provided as a vector containing the parameters of the covariance function.
#' @param start_sigma_e Starting value for the standard deviation of the measurament error.
#' @param parameterization_latent The parameterization for `WhittleMatern` and `graphLaplacian` models. The options are 'matern' and 'spde'. The 'matern' parameterizes as 'sigma' and 'range', whereas the 'spde' parameterization is given in terms of 'sigma' and 'kappa'.
#' @param BC For `WhittleMatern` models. Which boundary condition to use (0,1) 0 is no adjustment on boundary point
#'        1 is making the boundary condition stationary.
#' @param model_matrix logical indicating whether the model matrix should be returned as component of the returned value.
#' @param optim_controls Additional controls to be passed to `optim` or `optimParallel`.
#' @return A list containing the fitted model.
#' @rdname graph_lme
#' @export
#' 

graph_lme <- function(formula, graph, 
                model = list(type = "WhittleMatern", alpha = 1, version = 1), 
                repl = NULL,
                optim_method = "L-BFGS-B", 
                starting_values_latent = NULL,
                start_sigma_e = NULL,
                parameterization_latent = c("matern", "spde"),
                BC = 1, 
                model_matrix = TRUE,
                optim_controls = list()) {
  model_type <- model[["type"]]
  model_type <- tolower(model_type)

  parameterization_latent <- parameterization_latent[[1]]

  if(!(parameterization_latent%in%c("matern", "spde"))){
    stop("The possible values for 'parameterization_latent' are 'matern' and 'spde'!")
  }

  if(!(BC%in%c(0,1))){
    stop("The possible values for 'BC' are 0 and 1!")
  }

  if(!(model_type%in% c("whittlematern", "graphlaplacian", "isocov"))){
    stop("The possible models are 'WhittleMatern', 'graphLaplacian', 'isoCov')!")
  }

  if(model_type%in% c("whittlematern", "graphlaplacian")){
    if(is.null(model[["alpha"]])){
      model[["alpha"]] <- 1
    }
    if(!(model[["alpha"]]%in%c(1,2))){
      stop("alpha should be either 1 or 2.")
    }
    if(parameterization_latent == "spde"){
      par_names <- c("sigma", "kappa")
    } else{
      par_names <- c("sigma", "range")
    }
  }

  if(model_type =="whittlematern"){
      if(is.null(model[["version"]])){
        model[["version"]] <- 1
      }
      if(!(model[["version"]]%in%c(1,2))){
        stop("version should be either 1 or 2.")
      }
  }

  if(model_type =="isocov"){
    if(is.null(model[["cov_function"]])){
      model[["cov_function"]] <- "exp_covariance"
    }
  }

  if (is.null(formula)) {
    stop("No formula provided!")
  }

  if (is.null(graph)) {
    stop("No graph provided!")
  }

  if (!inherits(graph, "metric_graph")) {
    stop("The graph must be of class 'metric_graph'!")
  }

  call_graph_lme <- match.call()

  if(is.null(graph$data)){
    stop("No data found in the graph. Please, add observations before fitting the model.")
  }

  nV_orig <- NULL

  if(model_type %in% c("graphlaplacian", "isocov")){
    graph_bkp <- graph$clone()
    graph_bkp$observation_to_vertex()
    nV_orig <- graph_bkp$nV
    data <- graph_bkp$data
  } else if((model_type == "whittlematern") && (model[["alpha"]] == 1) && (model[["version"]]==2)){
    graph_bkp <- graph$clone()
    graph_bkp$observation_to_vertex()
    data <- graph_bkp$data
  } else{
    data <- graph$data
    graph_bkp <- graph
  }

  y_term <- stats::terms(formula)[[2]]

  y_graph <- eval(y_term, envir = data, enclos = parent.frame())
  y_graph <- as.numeric(y_graph)

  cov_term <- stats::delete.response(terms(formula))

  X_cov <- stats::model.matrix(cov_term, data)

  if(all(dim(X_cov) == c(0,1))){
    names_temp <- colnames(X_cov)
    X_cov <- matrix(1, nrow = length(y_graph))
    colnames(X_cov) <- names_temp
  }

  if(is.null(starting_values_latent)){
    if(model_type == "whittlematern"){
      if(model[["alpha"]] == 1){
        model_start <- "alpha1"
      } else{
        model_start <- "alpha2"
      }
    } else if(model_type == "graphlaplacian"){
      if(model[["alpha"]] == 1){
        model_start <- "GL1"
      } else{
        model_start <- "GL2"
      }
    } else{
      graph_bkp$res_dist <- NULL
      if(is.character(model[["cov_function"]])){
        if(model[["cov_function"]] == "exp_covariance"){
          model_start <- "isoExp"
          par_names <- c("sigma", "kappa")
        } else if(model[["cov_function"]] %in% c("alpha1","alpha2", "GL1", "GL2")){
          model_start <- model[["cov_function"]]
        } 
      } else{
        stop("For 'isoCov' models with a non-exponential covariance, that are not 'alpha1', 'alpha2', 'GL1' or 'GL2', you should provide the starting values!")
      }
    }

      range_par <- FALSE
      if(model_type == "whittlematern"){
        range_par <- ifelse(parameterization_latent == "matern",TRUE,FALSE)
      }
      start_values <- graph_starting_values(graph = graph_bkp,
                    model = model_start, 
                    manual_data = unlist(y_graph),
                    log_scale = TRUE,
                    range_par = range_par)
  } else {
    start_values <- c(log(0.1*sd(y_graph),log(starting_values_latent)))
    par_names <- names(starting_values_latent)
  }

  if(!is.null(start_sigma_e)){
        start_values[1] <- log(start_sigma_e)
  }

  if(ncol(X_cov)>0){
    names_tmp <- colnames(X_cov)
    data_tmp <- cbind(y_graph, X_cov)
    data_tmp <- na.omit(data_tmp)
    temp_coeff <- lm(data_tmp[,1] ~ data_tmp[,-1] - 1)$coeff
    names(temp_coeff) <- names_tmp
    start_values <- c(start_values, temp_coeff)
    rm(data_tmp)
  }
  
  if(model_type == "whittlematern"){
    if(model[["alpha"]] == 1){
      if(model[["version"]] == 2){
        likelihood <- function(theta){
          return(-likelihood_alpha1_v2(theta = theta, graph = graph_bkp, 
              X_cov = X_cov, y = y_graph, repl = repl, BC = BC, parameterization = parameterization_latent))
        }
      } else {
        likelihood <- function(theta){
          return(-likelihood_alpha1(theta = theta, graph = graph_bkp, data_name = NULL, manual_y = y_graph,
                             X_cov = X_cov, repl = repl, BC = BC, parameterization = parameterization_latent))
        }
      }
    } else{
      likelihood <- function(theta){
          return(-likelihood_alpha2(theta = theta, graph = graph_bkp, data_name = NULL, manual_y = y_graph,
                             X_cov = X_cov, repl = repl, BC = BC, parameterization = parameterization_latent))
        }
    }
  } else if (model_type == "graphlaplacian"){
      likelihood <- likelihood_graph_laplacian(graph = graph_bkp, alpha = model[["alpha"]], y_graph = y_graph, 
              X_cov = X_cov, maximize = FALSE, repl=repl, parameterization = "spde")
  } else if (is.character(model[["cov_function"]])) {
    if(model[["cov_function"]] %in% c("alpha1","alpha2", "GL1", "GL2")){
      model_cov <- model[["cov_function"]]
      par_names <- c("sigma", "kappa")
    } else{
      model_cov <- "isoCov"
      if(model[["cov_function"]] == "exp_covariance"){
        model[["cov_function"]] <- exp_covariance
        model[["cov_function_name"]] <- "exp_covariance"
      }
    }
    likelihood <- likelihood_graph_covariance(graph_bkp, model = model_cov, y_graph = y_graph,
                                                cov_function = model[["cov_function"]],
                                                X_cov = X_cov, repl = repl)
  } else{
    model[["cov_function_name"]] <- "other"
      likelihood <- likelihood_graph_covariance(graph_bkp, model = model_cov, y_graph = y_graph,
                                                cov_function = model[["cov_function"]],
                                                X_cov = X_cov, repl = repl)
  }

  res <- optim(start_values, 
                likelihood, method = optim_method,
                control = optim_controls,
                hessian = TRUE)

  coeff <- res$par
  coeff <- exp(c(res$par[1:3]))
  coeff <- c(coeff, res$par[-c(1:3)])

  observed_fisher <- res$hessian
  inv_fisher <- tryCatch(solve(observed_fisher), error = function(e) matrix(NA,
                                                                        nrow(observed_fisher), ncol(observed_fisher)))

  n_fixed <- ncol(X_cov)
  coeff_meas <- coeff[1]
  names(coeff_meas) <- "std. dev"
  std_err <- sqrt(diag(inv_fisher))
  std_meas <- std_err[1]
  n_random <- length(coeff) - n_fixed - 1
  coeff_fixed <- NULL
  if(n_fixed > 0){
    coeff_fixed <- coeff[(2+n_random):length(coeff)]
    std_fixed <- std_err[(2+n_random):length(coeff)]
  } else{
    std_fixed <- NULL
  }
  coeff_random <- coeff[2:(1+n_random)]
  std_random <- std_err[2:(1+n_random)]
  names(coeff_random) <- par_names



  object <- list()
  object$coeff <- list(measurement_error = coeff_meas, 
  fixed_effects = coeff_fixed, random_effects = coeff_random)
  object$std_errors <- list(std_meas = std_meas,
        std_fixed = std_fixed, std_random = std_random) 
  object$loglik <- - res$value
  object$call <- call_graph_lme
  object$terms <- list(fixed_effects = X_cov)
  object$response <- list(y = y_graph)
  object$formula <- formula
  object$estimation_method <- optim_method
  object$parameterization_latent <- parameterization_latent
  object$repl <- repl
  object$optim_controls <- optim_controls
  object$latent_model <- model
  object$loglike <- -res$value
  object$BC <- BC
  object$niter <- res$counts
  object$response <- y_term
  object$covariates <- cov_term
  object$nV_orig <- nV_orig
  if(model_matrix){
    if(ncol(X_cov)>0){
      object$model_matrix <- cbind(y_graph, X_cov)
    } else{
      object$model_matrix <- y_graph
    }
  }
  object$graph <- graph$clone()


  class(object) <- "graph_lme"
  return(object)

}


#' @name print.graph_lme
#' @title Print Method for \code{graph_lme} Objects
#' @description Provides a brief description of results related to mixed effects metric graph models.
#' @param x object of class "graph_lme" containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @return Called for its side effects.
#' @noRd
#' @method print graph_lme
#' @export 

print.graph_lme <- function(x, ...) {
  #
  model_type <- tolower(x$latent_model$type)
  call_name <- switch(model_type,
                      "whittlematern" = {paste0("Latent model - Whittle-Matern with alpha = ",x$latent_model$alpha)},
                      "graphlaplacian" = {paste0("Latent model - graph Laplacian SPDE with alpha = ",x$latent_model$alpha)},
                      "isocov" = {"Covariance-based model"}
  )

  coeff_fixed <- x$coeff$fixed_effects
  coeff_random <- x$coeff$random_effects
  
  cat("\n")
  cat(call_name)
  cat("\n\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat(paste0("Fixed effects:", "\n"))
  if(!is.null(coeff_fixed)){
    print(coeff_fixed)
  } else{
    message("No fixed effects")
  }
  cat("\n")
  cat(paste0("Random effects:", "\n"))
  print(coeff_random)
  cat("\n")
  cat(paste0("Measurement error:", "\n"))
  print(x$coeff$measurement_error)
}


#' @name summary.graph_lme
#' @title Summary Method for \code{graph_lme} Objects.
#' @description Function providing a summary of results related to metric graph mixed effects regression models.
#' @param object an object of class "graph_lme" containing results from the fitted model.
#' @param ... not used.
#' @return An object of class \code{summary_graph_lme} containing several
#' informations of a *graph_lme* object.
#' @method summary graph_lme
#' @export
summary.graph_lme <- function(object, ...) {
  ans <- list()

  nfixed <- length(object$coeff$fixed_effects)
  nrandom <- length(object$coeff$random_effects)
  model_type <- tolower(object$latent_model$type)
  call_name <- switch(model_type,
                      "whittlematern" = {paste0("Latent model - Whittle-Matern with alpha = ",object$latent_model$alpha)},
                      "graphlaplacian" = {paste0("Latent model - graph Laplacian SPDE with alpha = ",object$latent_model$alpha)},
                      "isocov" = {"Covariance-based model"}
  )

  coeff_fixed <- object$coeff$fixed_effects
  coeff_random <- object$coeff$random_effects#
  coeff_meas <- object$coeff$measurement_error

  SEr_fixed <- object$std_errors$std_fixed
  SEr_random <- object$std_errors$std_random
  SEr_meas <- object$std_errors$std_meas

  coeff <- c(coeff_fixed, coeff_random, coeff_meas)
  SEr <- c(SEr_fixed,SEr_random, SEr_meas)

  tab <- cbind(coeff, SEr, coeff / SEr, 2 * stats::pnorm(-abs(coeff / SEr)))
  colnames(tab) <- c("Estimate", "Std.error", "z-value", "Pr(>|z|)")
  rownames(tab) <- names(coeff)
  tab <- list(fixed_effects = tab[seq.int(length.out = nfixed), , drop = FALSE], random_effects = tab[seq.int(length.out = nrandom) + nfixed, , drop = FALSE], 
  meas_error = tab[seq.int(length.out = 1) + nfixed+nrandom, , drop = FALSE])

  ans$coefficients <- tab

  ans$call_name <- call_name

  ans$call <- object$call

  ans$loglike <- object$loglike

  ans$niter <- object$niter

  class(ans) <- "summary_graph_lme"
  ans
}

#' @name print.summary_graph_lme
#' @title Print Method for \code{summary_graph_lme} Objects
#' @description Provides a brief description of results related to metric graph mixed effects regression models.
#' @param x object of class "summary_graph_lme" containing results of summary method applied to a fitted model.
#' @param ... further arguments passed to or from other methods.
#' @return Called for its side effects.
#' @noRd
#' @method print summary_graph_lme
#' @export
print.summary_graph_lme <- function(x, ...) {
  tab <- x$coefficients

  #
  digits <- max(3, getOption("digits") - 3)
  #

  call_name <- x$call_name

  cat("\n")
  cat(call_name)

  cat("\n\n")
  cat("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")


  #
  #
  if (NROW(tab$fixed_effects)) {
    cat(paste0("\nFixed effects:\n"))
    stats::printCoefmat(tab[["fixed_effects"]], digits = digits, signif.legend = FALSE)
  } else {
    message("\nNo coefficients modeling the fixed effects. \n")
  }
  #
  if (NROW(tab$random_effects)) {
    cat(paste0("\nRandom effects:\n"))
    stats::printCoefmat(tab[["random_effects"]][,1:3], digits = digits, signif.legend = FALSE)
  } else {
    message("\nNo coefficients modeling the random effects. \n")
  }
  #
  cat(paste0("\nMeasurement error:\n"))
    stats::printCoefmat(tab[["meas_error"]][1,1:3,drop = FALSE], digits = digits, signif.legend = FALSE)
  #
  if (getOption("show.signif.stars")) {
    cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n\n")
  }
  #

  cat("Log-Likelihood: ", x$loglike,"\n")

  cat(paste0("Number of function calls by 'optim' = ", x$niter[1],"\n"))
}



#' @name predict.graph_lme
#' @title Prediction of a mixed effects regression model on a metric graph.
#' @param object The fitted object with the `graph_lme()` function 
#' @param data A `data.frame` or a `list` containing the covariates, the edge number and the distance on edge
#' for the locations to obtain the prediction.
#' @param mesh Obtain predictions for mesh nodes? The graph must have a mesh, and either `only_latent` is set to TRUE or the model does not have covariates.
#' @param mesh_h If the graph does not have a mesh, one will be created with this value of 'h'.
#' @param repl Which replicates to obtain the prediction. If `NULL` predictions will be obtained for all replicates. Default is `NULL`.
#' @param compute_variances Set to also TRUE to compute the kriging variances.
#' @param posterior_samples If `TRUE`, posterior samples will be returned.
#' @param n_samples Number of samples to be returned. Will only be used if `sampling` is `TRUE`.
#' @param only_latent Should the posterior samples and predictions be only given to the latent model?
#' @param edge_number Name of the variable that contains the edge number, the default is `edge_number`.
#' @param distance_on_edge Name of the variable that contains the distance on edge, the default is `distance_on_edge`.
#' @param normalized Are the distances on edges normalized?
#' @param return_as_list Should the means of the predictions and the posterior samples be returned as a list, with each replicate being an element?
#' @param return_original_order Should the results be return in the original (input) order or in the order inside the graph?
#' @param ... Not used.
#' @export
#' @method predict graph_lme

predict.graph_lme <- function(object, data = NULL, mesh = FALSE, mesh_h = 0.01, repl = NULL, compute_variances = FALSE, posterior_samples = FALSE,
                               n_samples = 100, only_latent = FALSE, edge_number = "edge_number",
                               distance_on_edge = "distance_on_edge", normalized = FALSE, return_as_list = FALSE, return_original_order = TRUE,
                               ...) {

  if(is.null(data)){
    if(!mesh){
      stop("If 'mesh' is false, you should supply data!")
    }
  }

  out <- list()

  coeff_fixed <- object$coeff$fixed_effects
  coeff_random <- object$coeff$random_effects
  coeff_meas <- object$coeff$measurement_error

  graph_bkp <- object$graph$clone()

  X_cov_initial <- stats::model.matrix(object$covariates, graph_bkp$data)
  if(ncol(X_cov_initial) > 0){
    if(mesh){
      stop("In the presence of covariates, you should provide the data, including the covariates at the prediction locations. If you only want predictions for the latent model, set 'only_latent' to TRUE.")
    }
  }


  if(sum(duplicated(cbind(data["edge_number"], data["distance_on_edge"]))) > 0){
    warning("There are duplicated locations for prediction, we will try to process the data to extract the unique locations,
    along with the corresponding covariates.")
    cov_names <- attr(object$covariates,"term.labels")
    data <- data[c(edge_number,distance_on_edge,cov_names)]
    data <- unique(data) 
    if(sum(duplicated(cbind(data["edge_number"], data["distance_on_edge"]))) > 0){
      stop("Data processing failed, please provide a data with unique locations.")
    }
  }
  
  if(!mesh){
    n_prd <- length(data[[edge_number]])
    data[[as.character(object$response)]] <- rep(NA, n_prd)
    data[["__dummy_var"]] <- rep(0, n_prd)
  } else{
    if(is.null(graph_bkp$mesh)){
      graph_bkp$build_mesh(h = mesh_h)
    }
    data <- list()
    n_prd <- nrow(graph_bkp$mesh$VtE)
    data[[as.character(object$response)]] <- rep(NA, n_prd)
    data[["__dummy_var"]] <- rep(0, n_prd)
    data[[edge_number]] <- graph_bkp$mesh$VtE[,1]
    data[[distance_on_edge]] <- graph_bkp$mesh$VtE[,2]
    normalized <- TRUE
  }

  if(return_original_order){
    ord_idx <- order(data[[edge_number]], data[[distance_on_edge]])
  }

  graph_bkp$add_observations(data = data, edge_number = edge_number, distance_on_edge = distance_on_edge, normalized = normalized)

  n <- sum(graph_bkp$data[["__group"]] == graph_bkp$data[["__group"]][1])

  ## 
  repl_vec <- graph_bkp[["data"]][["__group"]]

  if(is.null(repl)){
    u_repl <- unique(graph_bkp$data[["__group"]])
  } else{
    u_repl <- unique(repl)
  }

  ##

  X_cov_pred <- stats::model.matrix(object$covariates, graph_bkp$data)
  
  if(all(dim(X_cov_pred) == c(0,1))){
    X_cov_pred <- matrix(1, nrow = n, ncol=1)
  }
  if(ncol(X_cov_pred) > 0){
    mu <- X_cov_pred %*% coeff_fixed
  } else{
    mu <- matrix(0, nrow = n, ncol=1)
  }

  Y <- graph_bkp$data[[as.character(object$response)]] - mu

  model_type <- object$latent_model

  sigma.e <- coeff_meas[[1]]

  idx_prd <- !is.na(graph_bkp$data[["__dummy_var"]][1:n])

  n_prd <- sum(idx_prd)

  edge_nb <- graph_bkp$data[["__edge_number"]][1:n][idx_prd]
  dist_ed <- graph_bkp$data[["__distance_on_edge"]][1:n][idx_prd]

  ## construct Q

  graph_bkp$observation_to_vertex()

  if(tolower(model_type$type) == "whittlematern"){
    sigma <- object$coeff$random_effects[1]
    if(object$parameterization_latent == "spde"){
      kappa <- object$coeff$random_effects[2]
    } else{
      kappa <- sqrt(8 * 0.5) / object$coeff$random_effects[2]
    }

      if(model_type$alpha == 1){
          Q <- spde_precision(kappa = kappa, sigma = sigma,
                            alpha = 1, graph = graph_bkp)
      } else{
        PtE <- graph_bkp$get_PtE()
        n.c <- 1:length(graph_bkp$CoB$S)
        Q <- spde_precision(kappa = kappa, sigma = sigma, alpha = 2,
                            graph = graph_bkp, BC = 1)
        Qtilde <- (graph_bkp$CoB$T) %*% Q %*% t(graph_bkp$CoB$T)
        Qtilde <- Qtilde[-n.c,-n.c]
        Sigma.overdetermined  = t(graph_bkp$CoB$T[-n.c,]) %*% solve(Qtilde) %*%
          (graph_bkp$CoB$T[-n.c,])
        index.obs <- 4 * (PtE[,1] - 1) + 1.0 * (abs(PtE[, 2]) < 1e-14) +
          3.0 * (abs(PtE[, 2]) > 1e-14)
        Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
        Q <- solve(Sigma)
      }

  } else if(tolower(model_type$type) == "graphlaplacian"){
    sigma <- object$coeff$random_effects[1]
    #nV before 
    nV_temp <- object$nV_orig
    # graph_bkp$observation_to_vertex()
    if(graph_bkp$nV > nV_temp){
      warning("There are prediction locations outside of the observation locations. Refit the model with all the locations you want to obtain predictions.")
    }
    graph_bkp$compute_laplacian()
    if(object$parameterization_latent == "spde"){
      kappa <- object$coeff$random_effects[2]
    } else{
      kappa <- sqrt(8 * 0.5) / object$coeff$random_effects[2]
    }
      if(model_type$alpha == 1){
        Q <- (kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]) / sigma^2
      } else{
        Q <- kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]
        Q <- Q %*% Q / sigma^2
      }

  } else if(tolower(model_type$type) == "isocov"){
      if(is.character(model_type$cov_function)){
        sigma <- object$coeff$random_effects[1]
        kappa <- object$coeff$random_effects[2]
        if(model_type$cov_function == "alpha1"){
          # graph_bkp$observation_to_vertex()
          Q <- spde_precision(kappa = kappa, sigma = sigma,
                            alpha = 1, graph = graph_bkp)
        } else if(model_type$cov_function == "alpha2"){
          PtE <- graph_bkp$get_PtE()
          n.c <- 1:length(graph_bkp$CoB$S)
          Q <- spde_precision(kappa = kappa, sigma = sigma, alpha = 2,
                              graph = graph_bkp, BC = 1)
          Qtilde <- (graph_bkp$CoB$T) %*% Q %*% t(graph_bkp$CoB$T)
          Qtilde <- Qtilde[-n.c,-n.c]
          Sigma.overdetermined  = t(graph_bkp$CoB$T[-n.c,]) %*% solve(Qtilde) %*%
            (graph_bkp$CoB$T[-n.c,])
          index.obs <- 4 * (PtE[,1] - 1) + 1.0 * (abs(PtE[, 2]) < 1e-14) +
            3.0 * (abs(PtE[, 2]) > 1e-14)
          Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
          Q <- solve(Sigma)
        } else if(model_type$cov_function == "GL1"){
              #nV before 
              nV_temp <- object$nV_orig
              # graph_bkp$observation_to_vertex()
              if(graph_bkp$nV > nV_temp){
                warning("There are prediction locations outside of the observation locations. Refit the model with all the locations you want to obtain predictions.")
              }
              graph_bkp$compute_laplacian()        
              Q <- (kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]) / sigma^2
        } else if(model_type$cov_function == "GL2"){
              #nV before 
              nV_temp <- object$nV_orig
              # graph_bkp$observation_to_vertex()
              if(graph_bkp$nV > nV_temp){
                warning("There are prediction locations outside of the observation locations. Refit the model with all the locations you want to obtain predictions.")
              }
              graph_bkp$compute_laplacian()
              Q <- kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]
              Q <- Q %*% Q / sigma^2
        # } else if(model_type$cov_function == "exp_covariance"){
        #           print("Here 1")
        #           graph_bkp$compute_resdist(full = TRUE)
        #           Sigma <- as.matrix(exp_covariance(graph_bkp$res_dist[[1]], c(sigma,kappa)))
        #           Q <- solve(Sigma)
        # } 
        }
      } else{
        graph_bkp$compute_resdist(full = TRUE)
        cov_function <- model_type$cov_function
        Sigma <- as.matrix(cov_function(graph_bkp$res_dist[[1]], coeff_random))
        Q <- solve(Sigma)
      }
  }

  # gap <- dim(Q)[1] - n
  
  ## compute Q_x|y
  # A <- Matrix::Diagonal(dim(Q)[1])[(gap+1):dim(Q)[1], ]
  if(tolower(model_type$type) == "isocov"){
    A <- Matrix::Diagonal(dim(Q)[1])
  } else{
    A <- Matrix::Diagonal(dim(Q)[1])[graph_bkp$PtV, ]
  }

  idx_obs_full <- as.vector(!is.na(Y))
  
  # idx_obs_full <- !is.na(graph_bkp$data[[as.character(object$response)]])

  if(!return_as_list){
    out$distance_on_edge <- rep(dist_ed,length(u_repl))
    out$edge_number <- rep(edge_nb,length(u_repl))
  }

  for(repl_y in u_repl){
    if(return_as_list){
      out$distance_on_edge[[repl_y]] <- dist_ed
      out$edge_number[[repl_y]] <- edge_nb
    }
    idx_repl <- graph_bkp$data[["__group"]] == repl_y

    idx_obs <- idx_obs_full[idx_repl]

    y_repl <- Y[idx_repl]
    y_repl <- y_repl[idx_obs]

    Q_xgiveny <- t(A[idx_obs,]) %*% A[idx_obs,]/sigma_e^2 + Q

    # cov_loc <- post_Cov[idx_prd, idx_obs]
    # cov_Obs <- post_Cov[idx_obs, idx_obs]

    # mu_krig <- cov_loc %*%  solve(cov_Obs, Y[idx_obs])

    mu_krig <- solve(Q_xgiveny,as.vector(t(A[idx_obs,]) %*% y_repl / sigma_e^2))

    # mu_krig <- mu_krig[(gap+1):length(mu_krig)]
    mu_krig <- A[idx_prd,] %*% mu_krig

    if(!only_latent){
      mu_fe <- mu[idx_repl, , drop = FALSE]
      mu_krig <- mu_fe[idx_prd, , drop=FALSE] + mu_krig
    } 

    mean_tmp <- as.vector(mu_krig)

    if(return_original_order){
      mean_tmp[ord_idx] <- mean_tmp
      # var_tmp[ord_idx] <- var_tmp
    }
        
    if(!return_as_list){
      out$mean <- c(out$mean, mean_tmp)
      out$repl <- c(out$repl, rep(repl_y,n_prd))
    } else{
      out$mean[[repl_y]] <- mean_tmp
    }

    if (compute_variances) {
        post_cov <- A[idx_prd,]%*%solve(Q_xgiveny, t(A[idx_prd,]))
        var_tmp <- diag(post_cov)
      if(return_original_order){
        var_tmp[ord_idx] <- var_tmp
      }
      if(!return_as_list){
        out$variance <- rep(var_tmp, length(u_repl))
      }
      else {
          for(repl_y in u_repl){
            out$variance[[repl_y]] <- var_tmp
          }
      }
    }

    if(posterior_samples){
      mean_tmp <- as.vector(mu_krig)
      Z <- rnorm(dim(post_cov)[1] * n_samples)
      dim(Z) <- c(dim(post_cov)[1], n_samples)
      LQ <- chol(forceSymmetric(post_cov))
      X <- LQ %*% Z
      X <- X + mean_tmp
      if(!only_latent){
        X <- X + matrix(rnorm(n_samples * length(mean_tmp), sd = sigma.e), nrow = length(mean_tmp))
      } else{
        X <- X - as.vector(mu_fe[idx_prd, , drop=FALSE])
      }

      if(return_original_order){
        X[ord_idx,] <- X
      }

      if(!return_as_list){
        out$samples <- rbind(out$samples, X)
      } else{
        out$samples[[repl_y]] <- X
      }
    }
  }

  return(out)
}


