#' Simulation of log-Gaussian Cox processes driven by Whittle-Matérn
#' fields on metric graphs
#' @param n Number of samples.
#' @param intercept Mean value of the Gaussian process.
#' @param sigma Parameter for marginal standard deviations.
#' @param range Parameter for practical correlation range.
#' @param alpha Smoothness parameter (1 or 2).
#' @param graph A `metric_graph` object.
#' @return List with Gaussian process sample and simulated points.
#' @export
graph_lgcp_sim <- function(n = 1, intercept = 0, sigma, range, alpha, graph) {

  if(is.null(graph$mesh)) {
    stop("No mesh provided")
  }

  if(n < 1 || n%%1 != 0){
    stop("n must be an integer")
  }

  if(is.null(graph$mesh$C)) {
    graph$compute_fem()
  }
  nu <- alpha - 1/2
  kappa <- sqrt(8*nu)/range
  tau <- sqrt(gamma(nu)/(gamma(alpha)*sqrt(4*pi)*kappa^(2*nu)*sigma^2))
  C <- Diagonal(dim(graph$mesh$C)[1], rowSums(graph$mesh$C))
  L <- kappa^2*C + graph$mesh$G

  if(alpha == 1) {
    Q <- tau^2*L
  } else if(alpha == 2) {
    Q <- tau^2*L%*%solve(C, L)
  } else {
    stop("not implemented yet")
  }
  R <- chol(Q)
  result <- list()
  for(i in 1:n){
    tmp <- solve(R, rnorm(dim(Q)[1]))
    u <- intercept + tmp

    lambda_max <- max(exp(u))
    domain_size <- sum(graph$edge_lengths)

    #simulate Poisson number of points
    N <- rpois(1, lambda_max*domain_size)

    #simulate locations of points from uniform distribution on edges
    p_edge <- graph$edge_lengths/domain_size
    edge_numbers <- sample(1:graph$nE,size = N, replace = TRUE, prob = p_edge)
    edge_loc <- runif(N)
    points <- cbind(edge_numbers, edge_loc)

    #Thin the sample
    lambda_loc <- exp(graph$fem_basis(points)%*%u)
    p_loc <- as.double(lambda_loc/lambda_max)
    ind_keep <- runif(N) < p_loc
    edge_numbers <- edge_loc <- NULL
    if (length(ind_keep) > 0) {
      edge_numbers <- points[ind_keep,1]
      edge_loc <- points[ind_keep,2]
    }
    result[[i]] <- list(u = u, edge_number = edge_numbers, edge_loc = edge_loc)
  }

  if(n == 1){
    return(result[[1]])
  }

  return(result)
}








#' Create a log-Gaussian Cox process model for metric graphs
#'
#' This function creates a log-Gaussian Cox process model for point pattern data on metric graphs.
#' It handles the creation of integration points and prepares the data for fitting with INLA.
#'
#' @param formula A formula object specifying the model structure
#' @param graph A metric_graph object containing the network and point pattern data
#' @param interpolate Logical; if TRUE, interpolate covariates from the graph data to integration points
#' @param manual_integration_points Data frame with columns edge_number, distance_on_edge, and E (integration weights)
#'        for manually specified integration points, or NULL to use automatic integration points
#' @param manual_covariates Named vector of covariates at integration points if interpolate is FALSE and covariates are used
#' @param use_current_mesh Logical; if TRUE, use the existing mesh in the graph as integration points
#' @param new_h Numeric; mesh size for creating a new mesh if use_current_mesh is FALSE
#' @param new_n Integer; alternative to new_h, specifies the approximate number of mesh points
#' @param repl Vector of replicates to be used in the model. For all replicates, one must use ".all".
#' @param repl_col Name of the column in the data that contains the replicates. Default is ".group".
#' @param ... Additional arguments to be passed to inla
#'
#' @return An object containing the fitted LGCP model
#' @export



lgcp_graph <- function(formula, 
                       graph, 
                       interpolate = TRUE,
                       manual_integration_points = NULL,
                       manual_covariates = NULL,
                       use_current_mesh = TRUE,
                       new_h = NULL,
                       new_n = NULL,
                       repl = ".all",
                       repl_col = ".group",
                       ...) {
          if(!is.null(manual_covariates)){
            interpolate <- FALSE
          }
          graph_bkp <- graph$clone()
          graph_bkp$.__enclos_env__$private$data[[".weights_int_points"]] <- rep(0, length(graph_bkp$.__enclos_env__$private$data[[".edge_number"]]))

          # Extract response variable name from the left-hand side of formula
          resp_variable_name <- as.character(formula[[2]])

          # Check if the response variable exists in the data
          if (!resp_variable_name %in% names(graph_bkp$.__enclos_env__$private$data)) {
            warning(paste0("Response variable '", resp_variable_name, "' not found in the data. A response variable will be created assuming counts on all available locations."))
            graph_bkp$.__enclos_env__$private$data[[resp_variable_name]] <- rep(1, length(graph_bkp$.__enclos_env__$private$data[[".edge_number"]]))
          } else{
            if(!(any(graph_bkp$.__enclos_env__$private$data[[resp_variable_name]] %in% c(0,1)))){
              stop(paste0("Response variable '", resp_variable_name, "' must be a binary variable (0 or 1)"))
            }
          }

          # Parse formula to extract covariates and their models
          formula_components <- parse_formula_components(formula)
          
          # Extract covariate names where model is a character (is_character is TRUE)
          character_covariates <- c()
          for (component in formula_components) {
            if (!is.null(component) && component$character) {
              character_covariates <- c(character_covariates, component$covariate)
            }
          }

          # Check for non-character models and validate their classes
          for (component in formula_components) {
            if (!is.null(component) && !component$character) {
              # Validate model class
              valid_classes <- c("inla_metric_graph_spde", "rspde_metric_graph")
              if (!any(valid_classes == component$model)) {
                stop(paste0("Model for '", component$covariate, 
                           "' must be one of: '", 
                           paste(valid_classes, collapse = "', '"), 
                           "', but got '", component$model, "' instead"))
              }
            }
          }

          int_points <- create_integration_points(graph = graph_bkp,
                                                 use_current_mesh = use_current_mesh,
                                                 new_h = new_h,
                                                 new_n = new_n,
                                                 interpolate = interpolate,
                                                 covariates = character_covariates,
                                                 manual_integration_points = manual_integration_points,
                                                 manual_covariates = manual_covariates,
                                                 repl = repl,
                                                 repl_col = repl_col)
          
          int_points[[resp_variable_name]] <- rep(0, nrow(int_points))

          graph_bkp$add_observations(data = int_points,
                                     edge_number = ".edge_number",
                                     distance_on_edge = ".distance_on_edge",
                                     normalized = TRUE,
                                     group = repl_col,
                                     verbose = 0)          

          # Extract model from formula components
          model_name <- NULL
          # Count non-character components (models)
          non_character_count <- sum(sapply(formula_components, function(comp) {
            !is.null(comp) && !comp$character
          }))
          
          if (non_character_count > 1) {
            stop("Only one spde-type model is allowed in the formula")
          }
          
          # Extract the model if it exists
          for (component in formula_components) {
            if (!is.null(component) && !component$character) {
              aux_spde_model <- get(component$model_name, envir = parent.frame())
              model_name <- component$covariate
              break
            }
          }

          if(inherits(aux_spde_model, "inla_metric_graph_spde")){
            aux_spde_model <- graph_spde(graph_bkp, alpha = aux_spde_model$alpha, parameterization = aux_spde_model$parameterization)
          }


          if (!is.null(aux_spde_model)) {            
            if(inherits(aux_spde_model, "inla_metric_graph_spde")){
            # Extract all covariate names from formula components
              data_spde <- graph_data_spde(aux_spde_model, name=model_name, covariates=character_covariates, repl = repl, repl_col = repl_col)
            } else{
              aux_spde_model$mesh <- graph_bkp
              data_spde <- graph_data_rspde_internal(aux_spde_model, name=model_name, covariates=character_covariates, repl = repl, repl_col = repl_col)
            }
          } else {
            aux_spde_model$mesh <- graph_bkp
            data_spde <- graph_data_linear_inla(aux_spde_model, covariates = character_covariates, repl = repl, repl_col = repl_col)
          }

          stk <- INLA::inla.stack(data = data_spde[["data"]], 
                  A = data_spde[["basis"]],
                  effects = data_spde[["index"]])

          # Update formula to use aux_spde_model instead of the original model name
          if (!is.null(model_name) && !is.null(aux_spde_model)) {
            # Extract the original formula
            formula_str <- deparse(formula)
            
            # Find and replace the model term in the formula
            for (component in formula_components) {
              if (!is.null(component) && !component$character && component$covariate == model_name) {
                # Create the pattern to match the f() term with the original model name
                pattern <- paste0("f\\(", model_name, ",\\s*model\\s*=\\s*", component$model_name)
                
                # Replace with the aux_spde_model
                replacement <- paste0("f(", model_name, ", model = aux_spde_model")
                
                # Update the formula string
                formula_str <- gsub(pattern, replacement, formula_str)
                
                # Convert back to a formula object
                formula <- as.formula(formula_str)
                break
              }
            }
          }
          
          inla_fit <- INLA::inla(formula,
                           data = INLA::inla.stack.data(stk),
                           family = "poisson",
                           control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE, link = 1),
                           E = INLA::inla.stack.data(stk)[[".weights_int_points"]],
                           ...)
          return(inla_fit)
} 
