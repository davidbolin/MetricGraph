
#' Simulate log-Gaussian Cox processes on metric graphs
#'
#' This function simulates point patterns from a log-Gaussian Cox process (LGCP) 
#' driven by Whittle-Matérn Gaussian random fields on metric graphs. The intensity 
#' function is modeled as \\eqn{\\lambda(s) = \\exp(\\beta + u(s))}{λ(s) = exp(β + u(s))}, where \\eqn{\\beta}{β} is an intercept parameter 
#' and \\eqn{u(s)}{u(s)} is a Gaussian field with Whittle-Matérn covariance.
#'
#' @param n Integer. Number of replicate point patterns to simulate. Default is 1.
#' @param intercept Numeric scalar or vector. Mean value(s) of the log-intensity 
#'   field. Can be a constant or vary spatially if provided as a vector of length 
#'   equal to the number of mesh nodes.
#' @param sigma Numeric. Marginal standard deviation parameter of the Gaussian field.
#' @param range Numeric. Practical correlation range parameter of the Gaussian field.
#' @param alpha Numeric. Smoothness parameter of the Whittle-Matérn field. Currently 
#'   supports values 1 and 2, corresponding to exponential and Matérn-3/2 covariance 
#'   functions respectively.
#' @param graph A `metric_graph` object with a built mesh. The graph must have an 
#'   existing mesh (built with `graph$build_mesh()`) and computed FEM matrices 
#'   (computed with `graph$compute_fem()`).
#'
#' @return If `n = 1`, returns a list with components:
#'   \describe{
#'     \item{u}{Numeric vector of the simulated Gaussian field values at mesh nodes}
#'     \item{edge_number}{Integer vector of edge numbers where points were simulated}
#'     \item{edge_loc}{Numeric vector of locations along edges (normalized coordinates)}
#'   }
#'   If `n > 1`, returns a list of length `n`, where each element is a list with 
#'   the above components.
#'
#' @details
#' The function implements a two-step simulation procedure:
#' \enumerate{
#'   \item Simulate the Gaussian random field u(s) using finite element methods
#'   \item Generate point locations using acceptance-rejection sampling from the 
#'         intensity \\eqn{\\lambda(s) = \\exp(\\beta + u(s))}{λ(s) = exp(β + u(s))}
#' }
#'
#' The Gaussian field is characterized by the SPDE:
#' \\eqn{(\\kappa^2 - \\Delta)^{\\alpha/2} \\tau u = \\mathcal{W}}{(κ² - Δ)^(α/2) τ u = W}
#' where \\eqn{\\kappa, \\tau}{κ, τ} are derived from the range and sigma parameters, and \\eqn{\\mathcal{W}}{W} is white noise.
#'
#' @examples
#' \dontrun{
#' # Create a metric graph and build mesh
#' graph <- metric_graph$new()
#' graph$build_mesh(h = 0.1)
#' graph$compute_fem()
#'
#' # Simulate a single point pattern
#' lgcp_data <- graph_lgcp_sim(
#'   intercept = -1, 
#'   sigma = 0.5, 
#'   range = 2, 
#'   alpha = 2, 
#'   graph = graph
#' )
#'
#' # Simulate multiple replicates
#' lgcp_replicates <- graph_lgcp_sim(
#'   n = 10,
#'   intercept = -1, 
#'   sigma = 0.5, 
#'   range = 2, 
#'   alpha = 2, 
#'   graph = graph
#' )
#'
#' # Plot the simulated intensity
#' graph$plot_function(X = exp(lgcp_data$u), vertex_size = 0)
#' }
#'
#' @seealso \code{\link{lgcp_graph}} for fitting LGCP models, 
#'   \code{\link{precompute_lgcp_graph}} for efficient model fitting
#'
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



#' Precompute expensive quantities for efficient LGCP model fitting
#'
#' This function precomputes the computationally expensive quantities needed for 
#' fitting log-Gaussian Cox process (LGCP) models on metric graphs. It enables 
#' efficient refitting of multiple models with different formulas while reusing 
#' the same spatial structure, integration points, and covariates. This is 
#' particularly beneficial for model selection, cross-validation, or sensitivity 
#' analysis scenarios.
#'
#' @param graph A `metric_graph` object containing the network structure and 
#'   point pattern data. Must have observations added via `add_observations()`.
#' @param resp_variable_name Character. Name of the response variable (typically 
#'   binary: 0 or 1) in the graph data that represents point occurrences.
#' @param model_name Character. Name to be used for the spatial field in INLA 
#'   formulas (e.g., "field" for `f(field, model = spde_model)`).
#' @param spde_model An SPDE model object of class `inla_metric_graph_spde` 
#'   (from `graph_spde()`) or `rspde_metric_graph` (from `rspde.metric_graph()`).
#' @param covariates Character vector. Names of covariates to include in the 
#'   precomputation. Only covariates listed here can be used in subsequent 
#'   model fits with the precomputed data.
#' @param interpolate Logical. If `TRUE` (default), interpolate covariate values 
#'   from graph data to integration points. If `FALSE`, use `manual_covariates`.
#' @param manual_integration_points Data frame with columns `edge_number`, 
#'   `distance_on_edge`, and `E` (integration weights). If `NULL`, automatic 
#'   integration points are created from the mesh or specified parameters.
#' @param manual_covariates Data frame or named list containing covariate values 
#'   at integration points. Required when `interpolate = FALSE`. Must include 
#'   a `.group` column for replicates if using replicated data.
#' @param use_current_mesh Logical. If `TRUE` (default), use the existing mesh 
#'   in `graph$mesh$VtE` as integration points. If `FALSE` or no mesh exists, 
#'   create a new mesh using `new_h` or `new_n`.
#' @param new_h Numeric. Mesh resolution for creating a new mesh when 
#'   `use_current_mesh = FALSE`. Smaller values create finer meshes.
#' @param new_n Integer. Alternative to `new_h`, specifies the approximate 
#'   number of mesh nodes for the new mesh.
#' @param repl Character vector or `".all"`. Specifies which replicates to 
#'   include in the model. Use `".all"` to include all available replicates.
#' @param repl_col Character. Name of the column in the graph data that contains 
#'   replicate identifiers. Default is `".group"`.
#' @param clone_graph Logical. If `TRUE` (default), clone the graph to avoid 
#'   modifying the original object. If `FALSE`, work directly on the original 
#'   graph (faster but modifies the input).
#'
#' @return A list of class `"precomputed_lgcp"` containing:
#'   \describe{
#'     \item{graph}{The modified metric_graph object with integration points}
#'     \item{covariates}{Vector of available covariate names}
#'     \item{stk}{Precomputed INLA stack object}
#'     \item{resp_variable_name}{Name of the response variable}
#'     \item{aux_spde_model}{The SPDE model object prepared for the graph}
#'     \item{type_model}{Type of SPDE model ("exact" or "rational")}
#'     \item{model_name}{Name of the spatial field for formulas}
#'   }
#'
#' @details
#' This function is most beneficial when:
#' \itemize{
#'   \item Fitting multiple models with different covariate combinations
#'   \item Performing model selection or cross-validation
#'   \item Using exact SPDE models (which have higher setup costs)
#'   \item Working with large graphs or complex spatial structures
#' }
#'
#' The speedup typically becomes apparent when fitting 3 or more models, with 
#' greater benefits for more complex spatial structures and exact SPDE models.
#'
#' @note 
#' \itemize{
#'   \item All covariates you plan to use must be specified in the initial 
#'     precomputation
#'   \item The spatial structure (mesh, integration points) is fixed during 
#'     precomputation
#'   \item Manual covariate values should correspond to mesh nodes in 
#'     `graph$mesh$VtE` when `interpolate = FALSE`
#' }
#' 
#' @examples
#' \dontrun{
#' # Setup: Create graph with mesh and add point pattern data
#' graph <- metric_graph$new()
#' graph$build_mesh(h = 0.1)
#' 
#' # Add point pattern data with covariates
#' point_data <- data.frame(
#'   y = 1,  # All observed locations have y = 1
#'   edge_number = c(1, 2, 3, 1, 2),
#'   distance_on_edge = c(0.1, 0.3, 0.8, 0.9, 0.2),
#'   elevation = c(100, 150, 200, 120, 180),
#'   temperature = c(15, 12, 8, 14, 10)
#' )
#' graph$add_observations(point_data, normalized = TRUE)
#' 
#' # Create SPDE model
#' spde_model <- graph_spde(graph, alpha = 1)
#' 
#' # Precompute for multiple model fitting scenarios
#' precomputed_data <- precompute_lgcp_graph(
#'   graph = graph,
#'   resp_variable_name = "y",
#'   model_name = "field",
#'   spde_model = spde_model,
#'   covariates = c("elevation", "temperature"),
#'   use_current_mesh = TRUE
#' )
#' 
#' # Now fit multiple models efficiently
#' # Model selection: compare different covariate combinations
#' fit1 <- lgcp_graph(y ~ elevation + f(field, model = spde_model),
#'                    graph = graph, precomputed_data = precomputed_data)
#' 
#' fit2 <- lgcp_graph(y ~ temperature + f(field, model = spde_model),
#'                    graph = graph, precomputed_data = precomputed_data)
#' 
#' fit3 <- lgcp_graph(y ~ elevation + temperature + f(field, model = spde_model),
#'                    graph = graph, precomputed_data = precomputed_data)
#' 
#' # Compare models using log marginal likelihood
#' c(fit1$mlik[1], fit2$mlik[1], fit3$mlik[1])
#' 
#' # Using manual covariates (exact values at mesh nodes)
#' manual_covs <- data.frame(
#'   elevation = graph$mesh$VtE[,1] * 50 + 100,  # Synthetic elevation
#'   temperature = 20 - graph$mesh$VtE[,1] * 10, # Synthetic temperature
#'   .group = 1
#' )
#' 
#' precomputed_manual <- precompute_lgcp_graph(
#'   graph = graph,
#'   resp_variable_name = "y", 
#'   model_name = "field",
#'   spde_model = spde_model,
#'   covariates = c("elevation", "temperature"),
#'   manual_covariates = manual_covs,
#'   interpolate = FALSE
#' )
#' 
#' # For maximum performance (modifies original graph)
#' precomputed_fast <- precompute_lgcp_graph(
#'   graph = graph,
#'   resp_variable_name = "y",
#'   model_name = "field", 
#'   spde_model = spde_model,
#'   covariates = c("elevation", "temperature"),
#'   clone_graph = FALSE
#' )
#' 
#' # Example with replicates
#' replicate_data <- data.frame(
#'   y = 1,
#'   edge_number = c(1, 2, 1, 3, 2, 3),
#'   distance_on_edge = c(0.2, 0.4, 0.7, 0.1, 0.8, 0.6),
#'   elevation = c(110, 160, 130, 190, 170, 210),
#'   replicate_id = c(1, 1, 1, 2, 2, 2)
#' )
#' 
#' graph$clear_observations()
#' graph$add_observations(replicate_data, normalized = TRUE, group = "replicate_id")
#' 
#' precomputed_reps <- precompute_lgcp_graph(
#'   graph = graph,
#'   resp_variable_name = "y",
#'   model_name = "field",
#'   spde_model = spde_model, 
#'   covariates = "elevation",
#'   repl = ".all",
#'   repl_col = "replicate_id"
#' )
#' 
#' fit_reps <- lgcp_graph(y ~ elevation + f(field, model = spde_model, 
#'                                         replicate = field.repl),
#'                        graph = graph, precomputed_data = precomputed_reps)
#' }
#'
#' @seealso 
#' \code{\link{lgcp_graph}} for fitting LGCP models with precomputed data,
#' \code{\link{graph_lgcp_sim}} for simulating LGCP data,
#' \code{\link{graph_spde}} for exact SPDE models,
#' \code{\link{spde_metric_graph_result}} for extracting spatial parameter estimates
#' 
#' @export
precompute_lgcp_graph <- function(graph, 
                           resp_variable_name,
                           model_name,
                           spde_model,
                           covariates = NULL,
                           interpolate = TRUE,
                           manual_integration_points = NULL,
                           manual_covariates = NULL,
                           use_current_mesh = TRUE,
                           new_h = NULL,
                           new_n = NULL,
                           repl = ".all",
                           repl_col = ".group",
                           clone_graph = TRUE) {


  if(inherits(spde_model, c("inla_metric_graph_spde", "inla_metric_graph_lgcp_spde"))){
    type_model <- "exact"
  } else{
    type_model <- "rational"
  }
  if(!is.null(manual_covariates)){
    interpolate <- FALSE
  }
  if(clone_graph){
    graph_bkp <- graph$clone()
  } else{
    graph_bkp <- graph
  }
  graph_bkp$.__enclos_env__$private$data[[".weights_int_points"]] <- rep(0, length(graph_bkp$.__enclos_env__$private$data[[".edge_number"]]))
  # Check if the response variable exists in the data
  if (!resp_variable_name %in% names(graph_bkp$.__enclos_env__$private$data)) {
    warning(paste0("Response variable '", resp_variable_name, "' not found in the data. A response variable will be created assuming counts on all available locations."))
    graph_bkp$.__enclos_env__$private$data[[resp_variable_name]] <- rep(1, length(graph_bkp$.__enclos_env__$private$data[[".edge_number"]]))
  } else{
    if(!(any(graph_bkp$.__enclos_env__$private$data[[resp_variable_name]] %in% c(0,1)))){
      stop(paste0("Response variable '", resp_variable_name, "' must be a binary variable (0 or 1)"))
    }
  }
  int_points <- create_integration_points(graph = graph_bkp,
                                         use_current_mesh = use_current_mesh,
                                         new_h = new_h,
                                         new_n = new_n,
                                         interpolate = interpolate,
                                         covariates = covariates,
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
  
  if(type_model == "exact"){
      spde_model <- graph_spde(graph_bkp, alpha = spde_model$alpha, parameterization = spde_model$parameterization, stationary_endpoints = spde_model$args$stationary_endpoints, directional = spde_model$directional, start_range = spde_model$args$start_range, start_kappa = spde_model$args$start_kappa, prior_kappa = spde_model$args$prior_kappa, prior_sigma = spde_model$args$prior_sigma,
      start_tau = spde_model$args$start_tau, prior_tau = spde_model$args$prior_tau, factor_start_range = spde_model$args$factor_start_range, 
      type_start_range_bbox = spde_model$args$type_start_range_bbox, shared_lib = spde_model$args$shared_lib, debug = spde_model$args$debug,
      verbose = spde_model$args$verbose)
      data_spde <- graph_data_spde(spde_model, name=model_name, covariates=covariates, repl = repl, repl_col = repl_col)
    } else{
      spde_model$mesh <- graph_bkp
      data_spde <- graph_data_rspde_internal(spde_model, name=model_name, covariates=covariates, repl = repl, repl_col = repl_col)
    }
  stk <- INLA::inla.stack(data = data_spde[["data"]], 
          A = data_spde[["basis"]],
          effects = data_spde[["index"]])

  # Create return object
  precomputed <- list(
    graph = graph_bkp,
    covariates = covariates,
    stk = stk,
    resp_variable_name = resp_variable_name,
    aux_spde_model = spde_model,
    type_model = type_model,
    model_name = model_name,
    nrow_int_points = nrow(int_points)
  )
  
  class(precomputed) <- "precomputed_lgcp"
  return(precomputed)
}






#' Fit log-Gaussian Cox process models on metric graphs
#'
#' This function fits log-Gaussian Cox process (LGCP) models for point pattern 
#' data on metric graphs using R-INLA. It handles the complex setup required for 
#' LGCP modeling, including creation of integration points, data preparation, and 
#' interface with INLA's Poisson likelihood framework. The function supports both 
#' exact and rational SPDE approximations, multiple replicates, and efficient 
#' refitting using precomputed quantities.
#'
#' @param formula A formula object specifying the model structure. Should follow 
#'   INLA syntax, e.g., `y ~ covariate + f(field, model = spde_model)` where 
#'   `field` is the spatial random effect and `spde_model` is an SPDE model object.
#' @param graph A `metric_graph` object containing the network structure and 
#'   point pattern data. Must have observations added via `add_observations()`.
#' @param interpolate Logical. If `TRUE` (default), interpolate covariate values 
#'   from graph data to integration points. If `FALSE`, use `manual_covariates`.
#' @param manual_integration_points Data frame with columns `edge_number`, 
#'   `distance_on_edge`, and `E` (integration weights). If `NULL`, automatic 
#'   integration points are created.
#' @param manual_covariates Data frame containing covariate values at integration 
#'   points when `interpolate = FALSE`. Must include a `.group` column for 
#'   replicates if using replicated data.
#' @param use_current_mesh Logical. If `TRUE` (default), use the existing mesh 
#'   in the graph as integration points. If `FALSE`, create a new mesh.
#' @param new_h Numeric. Mesh resolution for creating a new mesh when 
#'   `use_current_mesh = FALSE`. Smaller values create finer meshes.
#' @param new_n Integer. Alternative to `new_h`, specifies the approximate 
#'   number of mesh nodes for the new mesh.
#' @param repl Character vector or `".all"`. Specifies which replicates to 
#'   include in the model. Use `".all"` to include all available replicates.
#' @param repl_col Character. Name of the column in the graph data that contains 
#'   replicate identifiers. Default is `".group"`.
#' @param clone_graph Logical. If `TRUE` (default), clone the graph to avoid 
#'   modifying the original object. If `FALSE`, work directly on the original 
#'   graph (faster but modifies the input). Only used when `precomputed_data` 
#'   is `NULL`.
#' @param precomputed_data Optional object of class `"precomputed_lgcp"` from 
#'   `precompute_lgcp_graph()`. Enables efficient refitting with different 
#'   formulas using the same spatial structure and covariates.
#' @param ... Additional arguments passed to `INLA::inla()`, such as 
#'   `control.inla`, `control.predictor`, `control.compute`, etc.
#'
#' @return An object of class `"inla"` containing the fitted LGCP model. This 
#'   includes posterior marginal distributions for model parameters, fitted 
#'   values, and other standard INLA output components. Use 
#'   `spde_metric_graph_result()` to extract spatial parameter estimates in 
#'   their original scale.
#'
#' @details
#' The function implements LGCP modeling using the approach of Simpson et al. (2016), 
#' where the log-Gaussian Cox process with intensity \\eqn{\\lambda(s) = \\exp(\\eta(s))}{λ(s) = exp(η(s))} is 
#' approximated using a Poisson likelihood with carefully constructed integration 
#' points and weights.
#'
#' The key steps are:
#' \enumerate{
#'   \item Create integration points (typically mesh nodes) across the graph
#'   \item Set up data with observed points (response = 1, weights = 0) and 
#'         integration points (response = 0, weights = integration weights)
#'   \item Fit using Poisson regression with the constructed weights as exposure
#' }
#'
#' The spatial component can be modeled using:
#' \itemize{
#'   \item Exact SPDE models via `graph_spde()` (slower setup, exact likelihood)
#'   \item Rational SPDE approximations via `rspde.metric_graph()` (faster, approximate)
#' }
#'
#' @section Performance:
#' For fitting multiple models with the same spatial structure:
#' \itemize{
#'   \item Use `precompute_lgcp_graph()` first, then `lgcp_graph()` with 
#'         `precomputed_data` for substantial speedups
#'   \item Set `clone_graph = FALSE` for additional performance gains when 
#'         you don't need to preserve the original graph
#' }
#'
#' @examples
#' \dontrun{
#' # Setup: Create graph with mesh and add point pattern data
#' graph <- metric_graph$new()
#' graph$build_mesh(h = 0.1)
#' graph$add_observations(data = your_point_data, 
#'                        edge_number = "edge_id", 
#'                        distance_on_edge = "location")
#'
#' # Create SPDE model
#' spde_model <- graph_spde(graph, alpha = 1)
#' # or: rspde_model <- rspde.metric_graph(graph, nu = 1.5)
#'
#' # Fit basic LGCP model
#' fit1 <- lgcp_graph(y ~ 1 + f(field, model = spde_model), 
#'                    graph = graph)
#'
#' # Fit model with covariates
#' fit2 <- lgcp_graph(y ~ elevation + temperature + 
#'                        f(field, model = spde_model), 
#'                    graph = graph)
#'
#' # Extract spatial parameter estimates
#' spde_result <- spde_metric_graph_result(fit2, "field", spde_model)
#' summary(spde_result)
#'
#' # Efficient fitting of multiple models
#' precomputed <- precompute_lgcp_graph(
#'   graph = graph,
#'   resp_variable_name = "y",
#'   model_name = "field", 
#'   spde_model = spde_model,
#'   covariates = c("elevation", "temperature", "slope")
#' )
#'
#' # Now fit multiple models efficiently
#' fit_a <- lgcp_graph(y ~ elevation + f(field, model = spde_model), 
#'                     graph = graph, precomputed_data = precomputed)
#' fit_b <- lgcp_graph(y ~ elevation + temperature + f(field, model = spde_model), 
#'                     graph = graph, precomputed_data = precomputed)
#' fit_c <- lgcp_graph(y ~ slope + f(field, model = spde_model), 
#'                     graph = graph, precomputed_data = precomputed)
#'
#' # Model with replicates
#' fit_rep <- lgcp_graph(y ~ covariate + f(field, model = spde_model, 
#'                                        replicate = field.repl), 
#'                       graph = graph)
#' }
#'
#' @references
#' Simpson, D., Illian, J., Lindgren, F., Sørbye, S., & Rue, H. (2016). 
#' Going off grid: Computationally efficient inference for log-Gaussian Cox 
#' processes. Biometrika, 103(1), 49-70.
#'
#' @seealso 
#' \code{\link{graph_lgcp_sim}} for simulating LGCP data,
#' \code{\link{precompute_lgcp_graph}} for efficient model refitting,
#' \code{\link{graph_spde}} for exact SPDE models,
#' \code{\link{spde_metric_graph_result}} for extracting spatial parameter estimates
#'
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
                       clone_graph = TRUE,
                       precomputed_data = NULL,
                       ...) {
        # Extract response variable name from the left-hand side of formula
          resp_variable_name <- as.character(formula[[2]])

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
              valid_classes <- c("inla_metric_graph_spde", "rspde_metric_graph", "inla_metric_graph_lgcp_spde")
              if (!any(valid_classes == component$model)) {
                stop(paste0("Model for '", component$covariate, 
                           "' must be one of: '", 
                           paste(valid_classes, collapse = "', '"), 
                           "', but got '", component$model, "' instead"))
              }
            }
          }

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



          if(!is.null(precomputed_data)){
            if(!inherits(precomputed_data, "precomputed_lgcp")){
              stop("The precomputed data must be a precomputed_lgcp object")
            }
            graph_bkp <- precomputed_data$graph
            stk <- precomputed_data$stk
            nrow_int_points <- precomputed_data$nrow_int_points
            if(precomputed_data$resp_variable_name != resp_variable_name){
              warning(paste0("The response variable name in the precomputed data (", precomputed_data$resp_variable_name, ") does not match the response variable name in the formula (", resp_variable_name, "). The variable in the precomputed data will be used."))
            }
            resp_variable_name <- precomputed_data$resp_variable_name
          # Check that all covariates in the formula are available in precomputed data
          missing_covariates <- setdiff(character_covariates, precomputed_data$covariates)
          if (length(missing_covariates) > 0) {
            stop(paste0("The following covariates in the formula are not available in the precomputed data: '", 
                       paste(missing_covariates, collapse = "', '"), 
                       "'. Available covariates in precomputed data: '", 
                       paste(precomputed_data$covariates, collapse = "', '"), "'"))
          }

          # Check if model_name matches precomputed_data$model_name
          if (!is.null(model_name) && !is.null(precomputed_data$model_name)) {
            if (model_name != precomputed_data$model_name) {
              warning(paste0("The model name in the formula (", model_name, 
                           ") does not match the model name in the precomputed data (", 
                           precomputed_data$model_name, 
                           "). The model name from the precomputed data will be used."))
              model_name <- precomputed_data$model_name
            }
          } else if (is.null(model_name) && !is.null(precomputed_data$model_name)) {
            model_name <- precomputed_data$model_name
          }
            aux_spde_model <- precomputed_data$aux_spde_model

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
                formula <- formula(paste(formula_str, collapse = ""))
                break
              }
            }
          } 

          } else{
          if(!is.null(manual_covariates)){
            interpolate <- FALSE
          }
          if(clone_graph){
            graph_bkp <- graph$clone()
          } else{
            graph_bkp <- graph
          }
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
            
          nrow_int_points <- nrow(int_points)

          graph_bkp$add_observations(data = int_points,
                                     edge_number = ".edge_number",
                                     distance_on_edge = ".distance_on_edge",
                                     normalized = TRUE,
                                     group = repl_col,
                                     verbose = 0)

          if(inherits(aux_spde_model, c("inla_metric_graph_spde", "inla_metric_graph_lgcp_spde"))){
            aux_spde_model <- graph_spde(graph_bkp, alpha = aux_spde_model$alpha, parameterization = aux_spde_model$parameterization, stationary_endpoints = aux_spde_model$args$stationary_endpoints, directional = aux_spde_model$directional, start_range = aux_spde_model$args$start_range, start_kappa = aux_spde_model$args$start_kappa, prior_kappa = aux_spde_model$args$prior_kappa, prior_sigma = aux_spde_model$args$prior_sigma,
      start_tau = aux_spde_model$args$start_tau, prior_tau = aux_spde_model$args$prior_tau, factor_start_range = aux_spde_model$args$factor_start_range, 
      type_start_range_bbox = aux_spde_model$args$type_start_range_bbox, shared_lib = aux_spde_model$args$shared_lib, debug = aux_spde_model$args$debug,
      verbose = aux_spde_model$args$verbose)
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
                formula <- formula(paste(formula_str, collapse = ""))
                break
              }
            }
          }
          }
            
          inla_fit <- INLA::inla(formula,
                           data = INLA::inla.stack.data(stk),
                           family = "poisson",
                           control.predictor = list(A = INLA::inla.stack.A(stk), compute = TRUE, link = 1),
                           E = INLA::inla.stack.data(stk)[[".weights_int_points"]],
                           ...)

          if(inherits(aux_spde_model, "inla_metric_graph_spde")){
            inla_fit[["graph_lgcp_ordering"]] <- aux_spde_model$ordering
            full_inla_data <- graph_bkp$get_data()
            original_order_idx <- order(aux_spde_model$ordering)

            full_inla_data[[".predicted_field"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".predicted_field_std_dev"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".predicted_field_mode"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".predicted_field_quantile_0.025"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".predicted_field_quantile_0.975"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".predicted_field_quantile_0.5"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".linear_predictor"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".linear_predictor_std_dev"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".linear_predictor_mode"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".linear_predictor_quantile_0.025"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".linear_predictor_quantile_0.975"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".linear_predictor_quantile_0.5"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".fitted_values"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".fitted_values_std_dev"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".fitted_values_mode"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".fitted_values_quantile_0.025"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".fitted_values_quantile_0.975"]] <- numeric(nrow(full_inla_data))
            full_inla_data[[".fitted_values_quantile_0.5"]] <- numeric(nrow(full_inla_data))

            full_inla_data[[".predicted_field"]] <- inla_fit$summary.random[[model_name]]$mean[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".predicted_field_std_dev"]] <- inla_fit$summary.random[[model_name]]$sd[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".predicted_field_mode"]] <- inla_fit$summary.random[[model_name]]$mode[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".predicted_field_quantile_0.025"]] <- inla_fit$summary.random[[model_name]]$`0.025quant`[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".predicted_field_quantile_0.975"]] <- inla_fit$summary.random[[model_name]]$`0.975quant`[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".predicted_field_quantile_0.5"]] <- inla_fit$summary.random[[model_name]]$`0.5quant`[1:nrow(full_inla_data)][original_order_idx]

            full_inla_data[[".linear_predictor"]] <- inla_fit$summary.linear.predictor$mean[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".linear_predictor_std_dev"]] <- inla_fit$summary.linear.predictor$sd[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".linear_predictor_mode"]] <- inla_fit$summary.linear.predictor$mode[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".linear_predictor_quantile_0.025"]] <- inla_fit$summary.linear.predictor$`0.025quant`[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".linear_predictor_quantile_0.975"]] <- inla_fit$summary.linear.predictor$`0.975quant`[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".linear_predictor_quantile_0.5"]] <- inla_fit$summary.linear.predictor$`0.5quant`[1:nrow(full_inla_data)][original_order_idx]

            full_inla_data[[".fitted_values"]] <- inla_fit$summary.fitted.values$mean[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".fitted_values_std_dev"]] <- inla_fit$summary.fitted.values$sd[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".fitted_values_mode"]] <- inla_fit$summary.fitted.values$mode[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".fitted_values_quantile_0.025"]] <- inla_fit$summary.fitted.values$`0.025quant`[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".fitted_values_quantile_0.975"]] <- inla_fit$summary.fitted.values$`0.975quant`[1:nrow(full_inla_data)][original_order_idx]
            full_inla_data[[".fitted_values_quantile_0.5"]] <- inla_fit$summary.fitted.values$`0.5quant`[1:nrow(full_inla_data)][original_order_idx]
            
          } else{
            full_inla_data <- graph_bkp$get_data()
            idx_int_points <- which(full_inla_data[[resp_variable_name]] == 0)

            full_inla_data[[".predicted_field"]] <- rep(NA, nrow(full_inla_data))
            full_inla_data[[".predicted_field_std_dev"]] <- rep(NA, nrow(full_inla_data))
            full_inla_data[[".predicted_field_mode"]] <- rep(NA, nrow(full_inla_data))
            full_inla_data[[".predicted_field_quantile_0.025"]] <- rep(NA, nrow(full_inla_data))
            full_inla_data[[".predicted_field_quantile_0.975"]] <- rep(NA, nrow(full_inla_data))
            full_inla_data[[".predicted_field_quantile_0.5"]] <- rep(NA, nrow(full_inla_data))


            # predicted field only available at the mesh nodes
            full_inla_data[[".predicted_field"]][idx_int_points] <- inla_fit$summary.random[[model_name]]$mean[1:nrow_int_points]
            full_inla_data[[".predicted_field_std_dev"]][idx_int_points] <- inla_fit$summary.random[[model_name]]$sd[1:nrow_int_points]
            full_inla_data[[".predicted_field_mode"]][idx_int_points] <- inla_fit$summary.random[[model_name]]$mode[1:nrow_int_points]
            full_inla_data[[".predicted_field_quantile_0.025"]][idx_int_points] <- inla_fit$summary.random[[model_name]]$`0.025quant`[1:nrow_int_points]
            full_inla_data[[".predicted_field_quantile_0.975"]][idx_int_points] <- inla_fit$summary.random[[model_name]]$`0.975quant`[1:nrow_int_points]
            full_inla_data[[".predicted_field_quantile_0.5"]][idx_int_points] <- inla_fit$summary.random[[model_name]]$`0.5quant`[1:nrow_int_points]

            full_inla_data[[".linear_predictor"]] <- inla_fit$summary.linear.predictor$mean[1:nrow(full_inla_data)]
            full_inla_data[[".linear_predictor_std_dev"]] <- inla_fit$summary.linear.predictor$sd[1:nrow(full_inla_data)]
            full_inla_data[[".linear_predictor_mode"]] <- inla_fit$summary.linear.predictor$mode[1:nrow(full_inla_data)]
            full_inla_data[[".linear_predictor_quantile_0.025"]] <- inla_fit$summary.linear.predictor$`0.025quant`[1:nrow(full_inla_data)]
            full_inla_data[[".linear_predictor_quantile_0.975"]] <- inla_fit$summary.linear.predictor$`0.975quant`[1:nrow(full_inla_data)]
            full_inla_data[[".linear_predictor_quantile_0.5"]] <- inla_fit$summary.linear.predictor$`0.5quant`[1:nrow(full_inla_data)]

            full_inla_data[[".fitted_values"]] <- inla_fit$summary.fitted.values$mean[1:nrow(full_inla_data)]
            full_inla_data[[".fitted_values_std_dev"]] <- inla_fit$summary.fitted.values$sd[1:nrow(full_inla_data)]
            full_inla_data[[".fitted_values_mode"]] <- inla_fit$summary.fitted.values$mode[1:nrow(full_inla_data)]
            full_inla_data[[".fitted_values_quantile_0.025"]] <- inla_fit$summary.fitted.values$`0.025quant`[1:nrow(full_inla_data)]
            full_inla_data[[".fitted_values_quantile_0.975"]] <- inla_fit$summary.fitted.values$`0.975quant`[1:nrow(full_inla_data)]
            full_inla_data[[".fitted_values_quantile_0.5"]] <- inla_fit$summary.fitted.values$`0.5quant`[1:nrow(full_inla_data)]
          }
          
          inla_fit[["graph_fitted_data"]] <- full_inla_data
          return(inla_fit)
} 


