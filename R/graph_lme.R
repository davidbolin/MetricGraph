
#' @rdname graph_lme
#' @export

graph_lme <- function(formula, graph, 
                model = list(type = "WhittleMatern", alpha = 1), 
                cov_function = NULL,
                repl = NULL,
                optim_method = "L-BFGS-B", optim_controls = list(),
                version = 2) {
  model <- model[[1]]
  if(!(model%in% c("alpha1", "alpha2", "GL1", "GL2", "isoCov"))){
    stop("The possible models are 'alpha1', 'alpha2', 'GL1', 'GL2', 'isoCov')!")
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
    stop("No data found in the graph. Please, first add observations.")
  }

  if(!(version%in%c(1,2))){
    stop("The possible choices for 'version' are 1 and 2.")
  }

  if(model %in% c("GL1", "GL2")){
    graph_bkp <- graph$clone()
    graph_bkp$observation_to_vertex()
    data <- graph_bkp$data
  } else if((model == "alpha1") && (version==2)){
    graph_bkp <- graph$clone()
    graph_bkp$observation_to_vertex()
    data <- graph_bkp$data
  } else{
    data <- graph$data
  }


  model_terms <- stats::terms(formula, data = data, rhs = 1)
  
  y_term <- stats::terms(formula)[[2]]

  y_graph <- lapply(data, function(dat){eval(y_term, envir = dat, enclos = parent.frame())})

  cov_term <- delete.response(terms(formula))

  X_cov <- lapply(data, function(dat){
    model.matrix(cov_term, dat)
    })

    # likelihood <- switch(model,
    #     alpha_1 = {
    #              likelihood_graph_spde(graph = graph, alpha = 1,
    #         X_cov = X_cov, y = y_graph, version = version)
    #     }
    # )

        likelihood <- likelihood_graph_spde(graph = graph_bkp, 
        alpha = 1, log_scale = TRUE, maximize = FALSE,
            X_cov = X_cov, y = y_graph, version = version,
            repl = repl)
  
  start_values <- graph_starting_values(graph = graph,
                    model = model, 
                    manual_data = unlist(y_graph))
  log_theta_0 <- log(start_values)
  start_theta <- c(log_theta_0, rep(0, ncol(X_cov[[1]])))
  res <- optim(start_theta, 
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
