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
#' contain a parameter `cov_function`, containing the covariance function (which is `exp_covariance`, the 
#' exponential covariance, by default). Finally, for `Whittle-Matern` models, there is an additional parameter
#' `version`, which can be either 1 or 2, to tell which version of the likelihood should be used (version is 1 by default).
#' @param repl vector or list containing which replicates to consider in the model.
#' If `NULL` all replicates will be considered.
#' @param optim_method The method to be used with `optim` function.
#' @param starting_values_latent A vector containing the starting values for the latent model. If the latent model is `WhittleMatern` or `graphLaplacian`, then the starting values should be provided as a vector of the form c(sigma,kappa) or c(sigma,range) depending on the parameterization. If the model is `isoCov`, then the starting values should be provided as a vector containing the parameters of the covariance function.
#' @param parameterization_latent The parameterization for `WhittleMatern` and `graphLaplacian` models. The options are 'matern' and 'spde'. The 'matern' parameterizes as 'sigma' and 'range', whereas the 'spde' parameterization is given in terms of 'sigma' and 'kappa'.
#' @param parallel Logical. Should optimParallel be used instead of optim?
#' @param parallel_controls A list containing the parallel controls.
#' @param optim_controls Additional controls to be passed to `optim` or `optimParallel`.
#' @param BC For `WhittleMatern` models. Which boundary condition to use (0,1) 0 is no adjustment on boundary point
#'        1 is making the boundary condition stationary.
#' @return A list containing the fitted model.
#' @rdname graph_lme
#' @export
#' 

graph_lme <- function(formula, graph, 
                model = list(type = "WhittleMatern", alpha = 1, version = 1), 
                repl = NULL,
                optim_method = "L-BFGS-B", 
                starting_values_latent = NULL,
                parameterization_latent = c("matern", "spde"),
                BC = 1, 
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
  }

  if(model_type =="whittlematern"){
      if(is.null(model[["version"]])){
        model[["version"]] <- 1
      }
      if(!(model[["version"]]%in%c(1,2))){
        stop("version should be either 1 or 2.")
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

  if(model_type == "graphlaplacian"){
    graph_bkp <- graph$clone()
    graph_bkp$observation_to_vertex()
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
      if(model[["cov_function"]] == "exp_covariance"){
          model_start <- "isoExp"
      } else{
        stop("For 'isoCov' models with a non-exponential covariance, you should provide the starting values!")
      }
    }


      range_par <- ifelse(parameterization_latent == "matern",TRUE,FALSE)
      start_values <- graph_starting_values(graph = graph_bkp,
                    model = model_start, 
                    manual_data = unlist(y_graph),
                    log_scale = TRUE,
                    range_par = range_par)
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
  }

  # likelihood <- likelihood_graph_spde(graph = graph_bkp, 
  #       alpha = 1, log_scale = TRUE, maximize = FALSE,
  #           X_cov = X_cov, y = y_graph, version = version,
  #           repl = repl)

  
  res <- optim(start_values, 
                likelihood, method = optim_method,
                control = optim_controls)

  coeff <- res$par
  coeff <- exp(c(res$par[1:3]))
  coeff <- c(coeff, res$par[-c(1:3)])

  object <- list()
  object$coeff <- coeff
  object$loglik <- - res$value
  object$call <- call_graph_lme
  object$terms <- list(fixed_effects = X_cov)
  object$response <- list(y = y_graph)
  object$formula <- formula

  class(object) <- "graph_lme"
  return(object)

}
