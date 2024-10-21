#' Metric graph linear mixed effects models
#'
#' Fitting linear mixed effects model in metric graphs. The random effects can be
#' Gaussian Whittle-Matern fields, discrete Gaussian Markov random fields based
#' on the graph Laplacian, as well as Gaussian random fields with isotropic
#' covariance functions.
#'
#' @param formula Formula object describing the relation between the response
#' variables and the fixed effects.
#' @param graph A `metric_graph` object.
#' @param model The random effects model that will be used (it also includes the
#' option of not having any random effects). It can be either a character,
#' whose options are 'lm', for linear models without random effects; 'WM1' and
#' 'WM2' for Whittle-Matern models with \eqn{\alpha}=1 and 2, with exact
#' precision matrices, respectively; 'WM' for Whittle-Matern models where one
#' also estimates the smoothness parameter via finite-element method; 'isoExp'
#' for a model with isotropic exponential covariance; 'GL1' and 'GL2' for a
#' SPDE model based on graph Laplacian with \eqn{\alpha} = 1 and 2, respectively.
#' 'WMD1' is the directed Whittle-Matern with  \eqn{\alpha}=1.
#' There is also the option to provide it as a list containing the elements
#' `type`, which can be `linearModel`, `WhittleMatern`, `graphLaplacian` or `isoCov`.
#' `linearModel` corresponds to a linear model without random effects.
#' For `WhittleMatern` models, that is, if the list contains `type = 'WhittleMatern'`,
#' one can choose between a finite element approximation of the precision matrix
#' by adding `fem = TRUE` to the list, or to use the exact precision matrix
#' (by setting `fem = FALSE`). If `fem` is `FALSE`, there is also the parameter
#' `alpha`, to determine the order of the SPDE, which is either 1 or 2. If `fem`
#' is `FALSE` and `alpha` is not specified, then the default value of `alpha=1`
#' will be used. If `fem` is `TRUE` and one does not specify `alpha`, it will be
#' estimated from the data. However, if one wants to have `alpha` fixed to some
#' value, the user can specify either `alpha` or `nu` in the list. See the
#' vignettes for examples. Finally, for type 'WhittleMatern', there is an optional
#' argument, `rspde_order`, that chooses the order of the rational approximation.
#' By default `rspde_order` is 2.
#' Finally, if one wants to fit a nonstationary model, then `fem` necessarily
#' needs to be `TRUE`, and one needs to also supply the matrices `B.tau`
#' and `B.kappa` or `B.range` and `B.sigma`.
#' For `graph-Laplacian` models, the list must also contain a parameter `alpha`
#' (which is 1 by default). For `isoCov` models, the list must
#' contain a parameter `cov_function`, containing the covariance function.
#' The function accepts a string input for the following covariance functions:
#' 'exp_covariance', 'WM1', 'WM2', 'GL1', 'GL2'. For another covariance function,
#' the function itself must be provided as the `cov_function` argument. The
#' default is 'exp_covariance', the exponential covariance. We also have
#' covariance-based versions of the Whittle-Matern and graph Laplacian models,
#' however they are much slower, they are the following (string) values for
#' 'cov_function': 'alpha1' and 'alpha2' for Whittle-Matern fields, and 'GL1'
#' and 'GL2' for graph Laplacian models. Finally, for `Whittle-Matern` models,
#' there is an additional parameter `version`, which can be either 1 or 2, to
#' tell which version of the likelihood should be used. Version is 1 by default.
#' @param which_repl Vector or list containing which replicates to consider in
#' the model. If `NULL` all replicates will be considered.
#' @param model_options A list containing additional options to be used in the model. Currently, it is possible to fix parameters during the estimation or change the starting values of the parameters. The general structure of the elements of the list is `fix_parname` and `start_parname`, where `parname` stands for the name of the parameter. If `fix_parname` is not `NULL`, then the model with be fitted with the `parname` being fixed at the value that was passed. If `start_parname` is not `NULL`, the model will be fitted using the value passed as starting value for `parname`. the For 'WM' models, the possible elements of the list are: `fix_sigma_e`, `start_sigma_e`, `fix_nu`, `start_nu`, `fix_sigma`, `start_sigma`, `fix_range`, `start_range`. Alternatively, one can use `fix_sigma_e`, `start_sigma_e`, `fix_nu`, `start_nu`, `fix_tau`, `start_tau`, `fix_kappa`, `start_kappa`. For 'WM1', 'WM2', 'isoExp', 'GL1' and 'GL2' models, the possible elements of the list are `fix_sigma_e`, `start_sigma_e`, `fix_sigma`, `start_sigma`, `fix_range`, `start_range`. Alternatively, one can use `fix_sigma_e`, `start_sigma_e`, `fix_tau`, `start_tau`, `fix_kappa`, `start_kappa`. For 'isoCov' models, the possible values are `fix_sigma_e`, `start_sigma_e`, `fix_par_vec`, `start_par_vec`. Observe that contrary to the other models, for 'isoCov' models, both `fix_par_vec` and `start_par_vec` should be given as vectors of the size of the dimension of the vector for the input of the covariance function passed to the 'isoCov' model. Furthermore, for 'isoCov' models, `fix_par_vec` is a logical vector, indicating which parameters to be fixed, and the values will be kept fixed to the values given to `start_par_vec`, one can also use `fix_sigma_e` and `start_sigma_e` for controlling the std. deviation of the measurement error.
#' @param previous_fit An object of class `graph_lme`. Use the fitted coefficients as starting values.
#' @param fix_coeff If using a previous fit, should all coefficients be fixed at the starting values?
#' @param optim_method The method to be used with `optim` function.
#' @param possible_methods Which methods to try in case the optimization fails or the hessian is not positive definite. The options are 'Nelder-Mead', 'L-BFGS-B', 'BFGS', 'CG' and 'SANN'. By default only 'L-BFGS-B' is considered.
# @param parameterization_latent The parameterization for `WhittleMatern` and `graphLaplacian` models. The options are 'matern' and 'spde'. The 'matern' parameterizes as 'sigma' and 'range', whereas the 'spde' parameterization is given in terms of 'sigma' and 'kappa'.
#' @param BC For `WhittleMatern` models, decides which boundary condition to use
#' (0,1). Here, 0 is Neumann boundary conditions and 1 specifies stationary boundary
#' conditions.
# @param model_matrix logical indicating whether the model matrix should be returned as component of the returned value.
#' @param parallel logical. Indicating whether to use `optimParallel()` or not.
#' @param n_cores Number of cores to be used if parallel is true.
#' @param optim_controls Additional controls to be passed to `optim()` or `optimParallel()`.
#' @param improve_hessian Should a more precise estimate of the hessian be obtained?
#' Turning on might increase the overall time.
#' @param hessian_args List of controls to be used if `improve_hessian` is `TRUE`.
#' The list can contain the arguments to be passed to the `method.args` argument
#' in the `hessian` function. See the help of the `hessian` function in 'numDeriv'
#' package for details. Observet that it only accepts the "Richardson" method for
#' now, the method "complex" is not supported.
#' @param check_euclidean Check if the graph used to compute the resistance distance has Euclidean edges? The graph used to compute the resistance distance has the observation locations as vertices.
#' @return A list containing the fitted model.
#' @rdname graph_lme
#' @export
#'
graph_lme <- function(formula, graph,
                model = list(type = "linearModel"),
                which_repl = NULL,
                optim_method = "L-BFGS-B",
                possible_methods = "L-BFGS-B",
                model_options = list(),
                BC = 1,
                previous_fit = NULL,
                fix_coeff = FALSE,
                parallel = FALSE,
                n_cores = parallel::detectCores()-1,
                optim_controls = list(),
                improve_hessian = FALSE,
                hessian_args = list(),
                check_euclidean = TRUE) {

  if(!is.list(model)){
    if(!is.character(model)){
      stop("The 'model' argument must be either a list or a character (string).")
    }
    model <- model[[1]]
    model <- tolower(model)
    if(!(model%in% c("lm", "wm", "wm1", "wm2", "isoexp", "gl1", "gl2","wmd1"))){
      stop("If model is a character (string), the options are 'lm', 'WM','WMD1', 'WM1', 'WM2', 'isoExp', 'GL1' or 'GL2'.")
    }
    model <- switch(model,
            "lm" = list(type = "linearModel"),
            "wm1" = list(type = "WhittleMatern", fem = FALSE, alpha = 1, version = 1, directional=0),
            "wm2" = list(type = "WhittleMatern", fem = FALSE, alpha = 2),
            "isoexp" = list(type = "isoCov"),
            "gl1" = list(type = "graphLaplacian", alpha = 1),
            "gl2" = list(type = "graphLaplacian", alpha = 2),
            "wm" = list(type = "WhittleMatern", fem = TRUE),
            'wmd1' = list(type = "WhittleMatern", fem = FALSE, alpha = 1, directional=1)
            )
  }

  possible_methods <- unique(c(optim_method,possible_methods))
  possible_methods <- intersect(possible_methods, c("Nelder-Mead", "L-BFGS-B", "BFGS", "SANN", "CG"))

  if(length(possible_methods) == 0){
    possible_methods <- optim_method[[1]]
  }

  model_type <- model[["type"]]
  model_type <- tolower(model_type)

  model_type <- switch(model_type,
            "lm" = "linearModel",
            "wm1" = "WhittleMatern",
            "wm2" = "WhittleMatern",
            "wmd1" = "WhittleMatern",
            "isoexp" = "isoCov",
            "gl1" = "graphLaplacian",
            "gl2" = "graphLaplacian",
            "wm" = "WhittleMatern",
            model_type
            )

  start_previous <- NULL     
  par_vec <- FALSE

  if(!is.null(previous_fit)){
    if(!inherits(previous_fit, "graph_lme")){
      warning("previous_fit is not a 'graph_lme' object, thus will be ignored.")
    }
    if(tolower(previous_fit$latent_model$type) %in% c("whittlematern", "graphlaplacian")){
      par_vec <- FALSE
      if(tolower(model_type) == "isocov"){
        warning("previous_fit is of type 'isoCov' and thus will be ignored.")
      } else {
        if(fix_coeff){
          model_options <- list(fix_kappa = previous_fit$coeff$random_effects[2],
                              fix_tau = previous_fit$coeff$random_effects[1],
                              fix_sigma_e = previous_fit$coeff$measurement_error)
        } else{
          model_options <- list(start_kappa = previous_fit$coeff$random_effects[2],
                              start_tau = previous_fit$coeff$random_effects[1],
                              start_sigma_e = previous_fit$coeff$measurement_error)
        }
      }
    } else if(tolower(previous_fit$latent_model$type) == "isocov"){
      par_vec <- TRUE
      if(tolower(model_type) %in% c("whittlematern", "graphlaplacian")){
        warning("previous_fit is not of type 'isoCov' and thus will be ignored")
      } else{ 
        if(fix_coeff){
          model_options <- list(start_par_vec = previous_fit$coeff$random_effects,
                                fix_par_vec = rep(TRUE, length(previous_fit$coeff$random_effects)),
                              fix_sigma_e = previous_fit$coeff$measurement_error)
        } else{
          model_options <- list(start_par_vec = previous_fit$coeff$random_effects,
                              start_sigma_e = previous_fit$coeff$measurement_error)
        }
      }
    } 
  } 

  check_model_options(model_options)

  if(model_type == "isocov"){
    if(is.null(graph$characteristics)){
      warning("No check for Euclidean edges have been perfomed on this graph. The isotropic covariance models are only known to work for graphs with Euclidean edges. You can check if the graph has Euclidean edges by running the `check_euclidean()` method. See the vignette https://davidbolin.github.io/MetricGraph/articles/isotropic_noneuclidean.html for further details.")
    } else if(is.null(graph$characteristics$euclidean)){
            warning("No check for Euclidean edges have been perfomed on this graph. The isotropic covariance models are only known to work for graphs with Euclidean edges. You can check if the graph has Euclidean edges by running the `check_euclidean()` method. See the vignette https://davidbolin.github.io/MetricGraph/articles/isotropic_noneuclidean.html for further details.")
    } else if(!graph$characteristics$euclidean){
                  warning("This graph DOES NOT have Euclidean edges. The isotropic covariance models are NOT guaranteed to work for this graph! See the vignette https://davidbolin.github.io/MetricGraph/articles/isotropic_noneuclidean.html for further details.")
    }
  }

  # parameterization_latent <- parameterization_latent[[1]]

  # if(!(parameterization_latent%in%c("matern", "spde"))){
  #   stop("The possible values for 'parameterization_latent' are 'matern' and 'spde'!")
  # }

  if(!(BC%in%c(0,1))){
    stop("The possible values for 'BC' are 0 and 1!")
  }

  if(!(model_type%in% c("whittlematern", "graphlaplacian", "isocov", "linearmodel"))){
    stop("The possible models are 'linearModel', 'WhittleMatern', 'graphLaplacian', 'isoCov')!")
  }

  if(is.null(which_repl)){
    which_repl <- unique(graph$.__enclos_env__$private$data[[".group"]])
  }

  if(model_type%in% c("whittlematern", "graphlaplacian")){
    # if(parameterization_latent == "spde"){
      par_names <- c("tau", "kappa")
    # } else{
    #   par_names <- c("sigma", "range")
    # }
  }

  if(model_type == "graphlaplacian"){
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

  if(model_type == "whittlematern"){
    if(is.null(model[["fem"]])){
      fem <- FALSE
      if(is.null(model[["alpha"]])){
        alpha <- 1
      }
    }
  }

  if(model_type == "whittlematern"){
    if(is.null(model[["fem"]])){
      model[["fem"]] <- FALSE
    }
    if(is.null(model[["directional"]])){
      model[["directional"]] <- 0
    }
    if(model[["fem"]]){
      if(is.null(model[["rspde_order"]])){
        rspde_order <- 2
      } else{
        rspde_order <- model[["rspde_order"]]
        if(!is.numeric(rspde_order)){
          stop("rspde_order must be numeric.")
        }
        if(rspde_order<=0){
          stop("rspde_order must be positive.")
        }
        if(!(rspde_order%%1 == 0) || (rspde_order > 8)){
          stop("rspde_order must be an integer between 1 and 8.")
        }
        if(length(rspde_order)>1){
          stop("rspde_order must be a number, not a vector.")
        }
      }


      if(!is.null(model[["B.tau"]]) && !is.null(model[["B.kappa"]])){
              rspde_object <- rSPDE::spde.matern.operators(graph = graph,
                                                m = rspde_order,
                                                parameterization = "spde",
                                                B.tau = model[["B.tau"]],
                                                B.kappa = model[["B.kappa"]])

      } else if(!is.null(model[["B.sigma"]]) && !is.null(model[["B.range"]])){
              rspde_object <- rSPDE::spde.matern.operators(graph = graph,
                                                m = rspde_order,
                                                parameterization = "matern",
                                                B.sigma = model[["B.sigma"]],
                                                B.range = model[["B.range"]])
      } else if ( (!is.null(model[["B.tau"]]) && is.null(model[["B.kappa"]])) ||
       (is.null(model[["B.tau"]]) && !is.null(model[["B.kappa"]])) ||
       (!is.null(model[["B.sigma"]]) && is.null(model[["B.range"]])) ||
       (is.null(model[["B.sigma"]]) && !is.null(model[["B.range"]]))){
        stop("You must either define both B.tau and B.kappa or both B.sigma and B.range.")
      } else{
        rspde_object <- rSPDE::matern.operators(graph = graph,
                                                m = rspde_order,
                                                parameterization = "spde")
      }

      if(!is.null(model[["alpha"]])){
        if(!is.numeric(model[["alpha"]])){
          stop("alpha must be numeric.")
        }
        nu <- model[["alpha"]] - 0.5

        if(nu <= 0){
          stop("nu = alpha - 0.5 must be positive.")
        }
        if(length(nu)>1){
          stop("nu must be a number, not a vector.")
        }
      } else if(!is.null(model[["nu"]])){
        nu <- model[["nu"]]
        if(!is.numeric(nu)){
          stop("alpha must be numeric.")
        }

        if(nu <= 0){
          stop("nu = alpha - 0.5 must be positive.")
        }
        if(length(nu)>1){
          stop("nu must be a number, not a vector.")
        }

      } else{
        nu <- NULL
      }

      df_data <- graph$.__enclos_env__$private$data

      y_term <- stats::terms(formula)[[2]]

      cov_term <- stats::delete.response(terms(formula))

      X_cov <- stats::model.matrix(cov_term, df_data)

      cov_names <- NULL

      if(!is.null(X_cov)){
        cov_names <- attr(cov_term, "term.labels")
      }

      names_temp <- c(as.character(y_term), cov_names, c(".edge_number", ".distance_on_edge", ".group", ".coord_x", ".coord_y"))

      df_data <- lapply(names_temp, function(i){df_data[[i]]})
      names(df_data) <- names_temp

      idx_notanyNA <- idx_not_any_NA(df_data)

      df_data <- lapply(df_data, function(dat){dat[idx_notanyNA]})

      fit <- rSPDE::rspde_lme(formula = formula,
                            loc = cbind(df_data[[".edge_number"]],
                            df_data[[".distance_on_edge"]]),
                            model = rspde_object,
                            repl = df_data[[".group"]],
                            nu = nu, which_repl = which_repl,
                            optim_method = optim_method,
                            data = df_data,
                            use_data_from_graph = FALSE,
                            parallel = parallel,
                            n_cores = n_cores,
                            starting_values_latent = NULL,
                            start_sigma_e = NULL,
                            optim_controls = optim_controls,
                            improve_hessian = improve_hessian,
                            hessian_args = hessian_args)
      fit$call <- call_graph_lme

      graph_bkp <- graph$clone()
      graph_bkp$.__enclos_env__$private$data <- lapply(names_temp, function(i){graph_bkp$.__enclos_env__$private$data[[i]]})
      names(graph_bkp$.__enclos_env__$private$data) <- names_temp
      fit$graph <- graph_bkp

      # if(fit$estimate_nu){
      #   names(fit$coeff$random_effects)[1] <- "alpha"
      #   fit$coeff$random_effects[1] <- fit$coeff$random_effects[1] + 0.5
      # }
      class(fit) <- c("graph_lme", class(fit))
      return(fit)
    } else{
      if(!(model[["alpha"]] %in% c(1,2))){
        stop("For WhittleMatern models, alpha must be either 1 or 2. For different values, set 'fem' to 'TRUE' instead.")
      }
    }
  }


  nV_orig <- NULL

  graph_bkp <- graph$clone()

  likelihood_new <- NULL
  new_likelihood <- NULL

  if(model_type %in% c("graphlaplacian", "isocov")){
    graph_bkp$observation_to_vertex(mesh_warning=FALSE)
    nV_orig <- graph_bkp$nV
    data <- graph_bkp$.__enclos_env__$private$data
  } else if((model_type == "whittlematern") && (model[["alpha"]] == 1) && (model[["version"]]==2)){
    graph_bkp$observation_to_vertex(mesh_warning=FALSE)
    data <- graph_bkp$.__enclos_env__$private$data
  } else{
    data <- graph$.__enclos_env__$private$data
  }

  y_term <- stats::terms(formula)[[2]]

  y_graph <- eval(y_term, envir = data, enclos = parent.frame())
  y_graph <- as.numeric(y_graph)

  cov_term <- stats::delete.response(terms(formula))

  X_cov <- stats::model.matrix(cov_term, data)

  cov_names <- NULL

  if(!is.null(X_cov)){
    n_cov <- ncol(X_cov)
    cov_names <- attr(cov_term, "term.labels")
  } else{
    n_cov <- 0
  }

  names_temp <- NULL

  if(all(dim(X_cov) == c(0,1))){
    names_temp <- colnames(X_cov)
    X_cov <- matrix(1, nrow = length(y_graph))
    colnames(X_cov) <- names_temp
  }


  names_temp <- c(as.character(y_term), cov_names, c(".edge_number", ".distance_on_edge", ".group", ".coord_x", ".coord_y"))

  graph_bkp$.__enclos_env__$private$data <- lapply(names_temp, function(i){graph_bkp$.__enclos_env__$private$data[[i]]})
  names(graph_bkp$.__enclos_env__$private$data) <- names_temp

  time_build_likelihood_start <- Sys.time()

  if( (is.null(model_options$start_par_vec) && !par_vec) || tolower(model_type) != "isocov"){
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
    } else if(model_type == "isocov"){
      graph_bkp$res_dist <- NULL
      if(is.character(model[["cov_function"]])){
        if(model[["cov_function"]] == "exp_covariance"){
          model_start <- "isoExp"
          par_names <- c("tau", "kappa")
        } else if(model[["cov_function"]] %in% c("WM1","WM2", "GL1", "GL2")){
          model_start <- model[["cov_function"]]
        }
      } else{
        stop("For 'isoCov' models with a non-exponential covariance, that are not 'WM1', 'WM2', 'GL1' or 'GL2', you should provide the starting values!")
      }
    }


    if(model_type != "linearmodel"){
      # range_par <- FALSE
      # if(model_type == "whittlematern"){
      #   range_par <- ifelse(parameterization_latent == "matern",TRUE,FALSE)
      # }
      if(model_type %in% c("whittlematern", "graphlaplacian")){
        rec_tau <- TRUE
      } else{
        rec_tau <- FALSE
      }
      start_fixed_values <- graph_starting_values(graph = graph_bkp,
                    model = model_start,
                    manual_data = unlist(y_graph),
                    log_scale = TRUE,
                    # range_par = range_par)
                    range_par = FALSE,
                    model_options = model_options,
                    rec_tau = rec_tau)
                
      start_values <- start_fixed_values$start_values
      start_values_orig <- start_values
      fixed_values <- start_fixed_values$fixed_values
    }
  } else if(model_type != "linearmodel") {
    if(is.null(model_options$start_sigma_e) && is.null(model_options$fix_sigma_e)){
      start_values <- c(log(0.1*sd(y_graph)),log(model_options$start_par_vec))
    } else if(!is.null(model_options$fix_sigma_e)) {
      start_values <- c(log(model_options$fix_sigma_e),log(model_options$start_par_vec))
    } else{
      start_values <- c(log(model_options$start_sigma_e),log(model_options$start_par_vec))
    }

    start_values_orig <- start_values
    
    par_names <- names(model_options$start_par_vec)
    if(is.null(model_options$fix_sigma_e)){
        fixed_values <- FALSE
    } else{
        fixed_values <- TRUE
    }
    if(is.null(model_options$fix_par_vec)){
      fixed_values <- c(fixed_values, rep(FALSE, length(model_options$start_par_vec)))
    } else{
      fixed_values <- c(fixed_values, model_options$fix_par_vec)
    }
  }

  if(ncol(X_cov)>0 && model_type != "linearmodel"){
    names_tmp <- colnames(X_cov)
    data_tmp <- cbind(y_graph, X_cov)
    data_tmp <- na.omit(data_tmp)
    temp_coeff <- lm(data_tmp[,1] ~ data_tmp[,-1] - 1)$coeff
    names(temp_coeff) <- names_tmp
    start_values <- c(start_values, temp_coeff)
    rm(data_tmp)
  }


  if(model_type == "whittlematern"){

    fix_v <- create_fix_vec_val(fixed_values)

    fix_vec <- fix_v$fix_vec
    fix_v_val <- fix_v$fix_v_val

    if(model[["directional"]] == 1){
      if(model[["alpha"]] == 1){
        likelihood <- function(theta){
            if(!is.null(X_cov)){
                  n_cov <- ncol(X_cov)
            } else{
                  n_cov <- 0
            }
          fix_v_val_full <- c(fix_v_val, rep(NA, n_cov))
          fix_vec_full <- c(fix_vec, rep(FALSE, n_cov))
          new_theta <- fix_v_val_full
          new_theta[!fix_vec_full] <- theta
          return(-likelihood_alpha1_directional(theta = new_theta, graph = graph_bkp,
                                    data_name = NULL, manual_y = y_graph,
                                    X_cov = X_cov, repl = which_repl,
                                    parameterization = "spde")) # , parameterization = parameterization_latent))
        }
      }

    }else{
      if(model[["alpha"]] == 1){
        if(model[["version"]] == 2){
          likelihood <- function(theta){
            if(!is.null(X_cov)){
                  n_cov <- ncol(X_cov)
            } else{
                  n_cov <- 0
            }
            fix_v_val_full <- c(fix_v_val, rep(NA, n_cov))
            fix_vec_full <- c(fix_vec, rep(FALSE, n_cov))
            new_theta <- fix_v_val_full
            new_theta[!fix_vec_full] <- theta
            return(-likelihood_alpha1_v2(theta = new_theta, graph = graph_bkp,
                X_cov = X_cov, y = y_graph, repl = which_repl, BC = BC,
                parameterization = "spde")) # parameterization = parameterization_latent))
          }
        } else {
          likelihood <- function(theta){
            if(!is.null(X_cov)){
                  n_cov <- ncol(X_cov)
            } else{
                  n_cov <- 0
            }
            fix_v_val_full <- c(fix_v_val, rep(NA, n_cov))
            fix_vec_full <- c(fix_vec, rep(FALSE, n_cov))
            new_theta <- fix_v_val_full
            new_theta[!fix_vec_full] <- theta

            return(-likelihood_alpha1(theta = new_theta, graph = graph_bkp,
                                      data_name = NULL, manual_y = y_graph,
                               X_cov = X_cov, repl = which_repl, BC = BC,
                               parameterization = "spde")) # , parameterization = parameterization_latent))
          }
        }
      } else{
        likelihood <- function(theta){
            if(!is.null(X_cov)){
                  n_cov <- ncol(X_cov)
            } else{
                  n_cov <- 0
            }
            fix_v_val_full <- c(fix_v_val, rep(NA, n_cov))
            fix_vec_full <- c(fix_vec, rep(FALSE, n_cov))
            new_theta <- fix_v_val_full
            new_theta[!fix_vec_full] <- theta       
            return(-likelihood_alpha2(theta = new_theta, graph = graph_bkp,
                                      data_name = NULL, manual_y = y_graph,
                               X_cov = X_cov, repl = which_repl, BC = BC,
                               parameterization = "spde")) # , parameterization = parameterization_latent))
          }
      }
    }
  } else if (model_type == "graphlaplacian"){

    fix_v <- create_fix_vec_val(fixed_values)

    fix_vec <- fix_v$fix_vec
    fix_v_val <- fix_v$fix_v_val

    likelihood <- likelihood_graph_laplacian(graph = graph_bkp,
                                               alpha = model[["alpha"]],
                                               y_graph = y_graph,
                                               X_cov = X_cov, maximize = FALSE,
                                               repl=which_repl,
                                               parameterization = "spde",
                                               fix_vec = fix_vec,
                                               fix_v_val = fix_v_val)
  } else if(model_type == "isocov") {
  if (is.character(model[["cov_function"]]) && !par_vec) {
    if(model[["cov_function"]] %in% c("WM1","WM2", "GL1", "GL2")){
      model_cov <- model[["cov_function"]]
      par_names <- c("tau", "kappa")
    } else if (is.character(model[["cov_function"]])){
      model_cov <- "isoCov"
      if(model[["cov_function"]] == "exp_covariance"){
        model[["cov_function"]] <- exp_covariance
        model[["cov_function_name"]] <- "exp_covariance"
      }
    }

    fix_v <- create_fix_vec_val(fixed_values)

    fix_vec <- fix_v$fix_vec
    fix_v_val <- fix_v$fix_v_val

    likelihood <- likelihood_graph_covariance(graph_bkp, model = model_cov,
                                              y_graph = y_graph,
                                              cov_function = model[["cov_function"]],
                                              X_cov = X_cov, repl = which_repl,
                                              fix_vec = fix_vec,
                                              fix_v_val = fix_v_val,
                                               check_euclidean = check_euclidean)
    } else{
    model[["cov_function_name"]] <- "other"
    model_cov <- "isoCov"
    if(is.character(model[["cov_function"]])){
      if(model[["cov_function"]] == "exp_covariance"){
          model[["cov_function"]] <- exp_covariance
          model[["cov_function_name"]] <- "exp_covariance"
      }    
    }

    # fix_vec <- model_options$fix_par_vec
    # fix_v_val <- model_options$start_par_values
    fix_vec <- fixed_values
    fix_vec_full <- c(fix_vec, rep(FALSE, ncol(X_cov)))
    fix_v_val <- start_values_orig
    start_values <- start_values[!fix_vec_full]

    likelihood <- likelihood_graph_covariance(graph_bkp, model = model_cov,
                                                y_graph = y_graph,
                                                cov_function = model[["cov_function"]],
                                                X_cov = X_cov, repl = which_repl,
                                                fix_vec = fix_vec,
                                                fix_v_val = fix_v_val,
                                               check_euclidean = check_euclidean)
    }
    }

  if(model_type != "linearmodel"){
      time_build_likelihood_end <- Sys.time()

      time_build_likelihood <- time_build_likelihood_end - time_build_likelihood_start

      hessian <- TRUE

      if(improve_hessian){
        hessian <- FALSE
      }

      time_par <- NULL

      likelihood_new <- function(theta){
        l_tmp <- tryCatch(likelihood(theta),
                            error = function(e){return(NULL)})
          if(is.null(l_tmp)){
              return(10^100)
          }
          return(l_tmp)
        }

      # if(!is.null(fix_par_vec)){
      #   likelihood_new2 <- function_factory_fix_var(likelihood_new, fix_par_vec, length(start_values), get_fixed_values(start_values, fix_par_vec, n_cov), n_cov)
      # } else{
      #   fix_par_vec <- rep(FALSE, length(starting_values_latent))
      # }


      if(parallel){
        start_par <- Sys.time()
        cl <- parallel::makeCluster(n_cores)
        parallel::setDefaultCluster(cl = cl)
        parallel::clusterExport(cl, "y_graph", envir = environment())
        parallel::clusterExport(cl, "graph_bkp", envir = environment())
        parallel::clusterExport(cl, "X_cov", envir = environment())
        parallel::clusterExport(cl, "which_repl", envir = environment())
        parallel::clusterExport(cl, "model", envir = environment())
        parallel::clusterExport(cl, "fix_vec", envir = environment())        
        parallel::clusterExport(cl, "fix_v_val", envir = environment())      
        # parallel::clusterExport(cl, "y_list", envir = environment())
        parallel::clusterExport(cl, "likelihood_graph_covariance",
                       envir = as.environment(asNamespace("MetricGraph")))
        parallel::clusterExport(cl, "likelihood_graph_laplacian",
                       envir = as.environment(asNamespace("MetricGraph")))
        parallel::clusterExport(cl, "likelihood_alpha2",
                       envir = as.environment(asNamespace("MetricGraph")))
        parallel::clusterExport(cl, "likelihood_alpha1",
                       envir = as.environment(asNamespace("MetricGraph")))
        parallel::clusterExport(cl, "likelihood_alpha1_v2",
                       envir = as.environment(asNamespace("MetricGraph")))

        end_par <- Sys.time()
        time_par <- end_par - start_par

          start_fit <- Sys.time()
          # res <- optimParallel::optimParallel(start_values_fix(start_values, fix_par_vec, n_cov),
          #               likelihood_new2, method = optim_method,
          #               control = optim_controls,
          #               hessian = hessian,
          #               parallel = list(forward = FALSE, cl = cl,
          #                   loginfo = FALSE))

          res <- optimParallel::optimParallel(start_values,
                        likelihood_new,
                        control = optim_controls,
                        hessian = hessian,
                        parallel = list(forward = FALSE, cl = cl,
                            loginfo = FALSE))

        end_fit <- Sys.time()
        time_fit <- end_fit-start_fit
        parallel::stopCluster(cl)


        time_hessian <- NULL

        if(!is.na(res[1])){
          if(!improve_hessian){
            observed_fisher <- res$hessian
          } else{
            if(!is.list(hessian_args)){
              stop("hessian_controls must be a list")
            }

            start_hessian <- Sys.time()
            observed_fisher <- numDeriv::hessian(likelihood_new, res$par,
                                                 method.args = hessian_args)
            end_hessian <- Sys.time()
            time_hessian <- end_hessian-start_hessian
          }
          if(nrow(observed_fisher)>0){
            eig_hes <- eigen(observed_fisher)$value
            cond_pos_hes <- (min(eig_hes) > 1e-15)
          }
        } else{
          stop("Could not fit the model. Please, try another method with 'parallel' set to FALSE.")
        }

         if(min(eig_hes) < 1e-15){
          warning("The optimization failed to provide a numerically positive-definite Hessian. You can try to obtain a positive-definite Hessian by setting 'improve_hessian' to TRUE or by setting 'parallel' to FALSE, which allows other optimization methods to be used.")
        }

      } else{
        # possible_methods <- c("Nelder-Mead", "L-BFGS-B", "BFGS", "CG")  
        start_fit <- Sys.time()
        # res <- withCallingHandlers(tryCatch(optim(start_values_fix(start_values, fix_par_vec, n_cov),
        #           likelihood_new2, method = optim_method,
        #           control = optim_controls,
        #           hessian = hessian), error = function(e){return(NA)}),
        #           warning = function(w){invokeRestart("muffleWarning")})

        res <- withCallingHandlers(tryCatch(optim(start_values,
                  likelihood_new, method = optim_method,
                  control = optim_controls,
                  hessian = hessian), error = function(e){return(NA)}),
                  warning = function(w){invokeRestart("muffleWarning")})

        end_fit <- Sys.time()
        time_fit <- end_fit-start_fit

        cond_pos_hes <- FALSE
        time_hessian <- NULL

        if(!is.na(res[1])){
          if(!improve_hessian){
            observed_fisher <- res$hessian
          } else{
            if(!is.list(hessian_args)){
              stop("hessian_controls must be a list")
            }

            start_hessian <- Sys.time()
            observed_fisher <- numDeriv::hessian(likelihood_new, res$par,
                                                 method.args = hessian_args)
            end_hessian <- Sys.time()
            time_hessian <- end_hessian-start_hessian
          }
          if(nrow(observed_fisher) > 0){
            eig_hes <- eigen(observed_fisher)$value
            cond_pos_hes <- (min(eig_hes) > 1e-15)
          }
        }



        problem_optim <- list()

        if(is.na(res[1]) || !cond_pos_hes){
          problem_optim[[optim_method]] <- list()
          if(is.na(res[1])){
            problem_optim[[optim_method]][["lik"]] <- NA
            problem_optim[[optim_method]][["res"]] <- res
            problem_optim[[optim_method]][["hess"]] <- NA
            problem_optim[[optim_method]][["time_hessian"]] <- NA
            problem_optim[[optim_method]][["time_fit"]] <- NA
          } else{
            problem_optim[[optim_method]][["lik"]] <- -res$value
            problem_optim[[optim_method]][["res"]] <- res
            problem_optim[[optim_method]][["hess"]] <- observed_fisher
            problem_optim[[optim_method]][["time_hessian"]] <- time_hessian
            problem_optim[[optim_method]][["time_fit"]] <- time_fit
          }
        }

        ok_optim <- FALSE
        orig_optim_method <- optim_method
        test_other_optim <- (is.na(res[1]) || !cond_pos_hes)

         if(test_other_optim){
              tmp_method <- optim_method
              while(length(possible_methods)>1){
                possible_methods <- setdiff(possible_methods, tmp_method)
                new_method <- possible_methods[1]
                time_fit <- NULL
                start_fit <- Sys.time()
                  res <- withCallingHandlers(tryCatch(optim(start_values,
                            likelihood_new, method = new_method,
                            control = optim_controls,
                            hessian = hessian), error = function(e){return(NA)}),
                            warning = function(w){invokeRestart("muffleWarning")})
                end_fit <- Sys.time()
                time_fit <- end_fit-start_fit
                tmp_method <- new_method
                cond_pos_hes <- FALSE
                if(is.na(res[1])){
                  problem_optim[[tmp_method]][["lik"]] <- NA
                  problem_optim[[tmp_method]][["res"]] <- res
                  problem_optim[[tmp_method]][["hess"]] <- NA
                  problem_optim[[tmp_method]][["time_hessian"]] <- NA
                  problem_optim[[tmp_method]][["time_fit"]] <- NA
                } else{
                  if(!improve_hessian){
                    observed_fisher <- res$hessian
                  } else{
                    if(!is.list(hessian_args)){
                      stop("hessian_controls must be a list")
                    }

                    start_hessian <- Sys.time()
                    observed_fisher <- numDeriv::hessian(likelihood, res$par,
                                                         method.args = hessian_args)
                    end_hessian <- Sys.time()
                    time_hessian <- end_hessian-start_hessian
                  }
                  if(nrow(observed_fisher)>0){
                    eig_hes <- eigen(observed_fisher)$value
                    cond_pos_hes <- (min(eig_hes) > 1e-15)
                  }

                  problem_optim[[tmp_method]][["lik"]] <- -res$value
                  problem_optim[[tmp_method]][["res"]] <- res
                  problem_optim[[tmp_method]][["hess"]] <- observed_fisher
                  problem_optim[[tmp_method]][["time_hessian"]] <- time_hessian
                  problem_optim[[tmp_method]][["time_fit"]] <- time_fit
                }

                cond_ok <- ((!is.na(res[1])) && cond_pos_hes)
                if(cond_ok){
                  optim_method <- new_method
                  ok_optim <- TRUE
                  break
                }
              }
          }

          if(test_other_optim){
              lik_val <- lapply(problem_optim, function(dat){dat[["lik"]]})
                if(all(is.na(lik_val))){
                  stop("The model could not be fitted. All optimizations method failed.")
                } else if(ok_optim){
                  warning(paste("optim method",orig_optim_method,"failed to provide a positive-definite Hessian. Another optimization method was used."))
                  rm(problem_optim)
                } else{
                  max_method <- which.max(lik_val)
                  res <- problem_optim[[max_method]][["res"]]
                  observed_fisher <- problem_optim[[max_method]][["hess"]]
                  time_hessian <- problem_optim[[max_method]][["time_hessian"]]
                  time_fit <- problem_optim[[max_method]][["time_fit"]]
                  if(length(fix_vec) != sum(fix_vec)){
                    warning("All optimization methods failed to provide a positive-definite Hessian. The optimization method with largest likelihood was chosen. You can try to obtain a positive-definite Hessian by setting 'improve_hessian' to TRUE.")
                  }
                }
          }
      }


  # res <- optim(start_values,
  #               likelihood, method = optim_method,
  #               control = optim_controls,
  #               hessian = TRUE)

  if(model_type %in% c("graphlaplacian", "whittlematern")){
    n_par_coeff <- 3
  } else if(model[["cov_function_name"]] == "exp_covariance"){
    n_par_coeff <- 3
  } else{
    n_par_coeff <- length(model_options$start_par_vec) + 1
  }

  coeff <- res$par

  if(all(res$par == start_values)){
    likelihood(start_values)
    if(length(fix_vec) != sum(fix_vec)){
      stop("The optimizer is stuck at the starting values! This probably indicates that the precision matrix is not positive-definite.")
    }
  }
  
  if(sum(!fix_vec)>0){
    coeff <- exp(c(res$par[1:(sum(!fix_vec))]))
    coeff <- c(coeff, res$par[-c(1:(sum(!fix_vec)))])
  } else{
    coeff <- res$par
  }

  loglik <- -res$value

  if(!is.null(X_cov)){
        n_fixed <- ncol(X_cov)
  } else{
        n_fixed <- 0
  }

  # n_random <- length(coeff) - n_fixed - 1
  n_random <- length(fix_vec) - 1

  if(sum(!fix_vec)>0){
    tmp_vec <- c(exp(-c(res$par[1:(sum(!fix_vec))])), rep(1,n_fixed))
  } else{
    tmp_vec <- rep(1,n_fixed)
  }

  if(length(tmp_vec)>1){
      par_change <- diag(tmp_vec)
  } else{
    par_change <- tmp_vec
  }

  observed_fisher <- par_change %*% observed_fisher %*% par_change

  new_likelihood <- NULL

  if(model_type %in% c("graphlaplacian", "whittlematern")){


      tmp_coeff <- rep(NA,3)
      tmp_coeff[!fix_vec] <- coeff[1:(sum(!fix_vec))]
      tmp_coeff[fix_vec] <- exp(fix_v_val[!is.na(fix_v_val)])

      tmp_coeff[2] <- 1/tmp_coeff[2]

      if(sum(!fix_vec)>0){
        grad_tmp <- diag(c(c(1,-1/(tmp_coeff[2]^2), 1)[!fix_vec], rep(1,n_fixed)))
        observed_fisher <- grad_tmp %*% observed_fisher %*% grad_tmp
      }



      time_matern_par_start <- Sys.time()
      # new_likelihood <- function(theta){
      #   new_par <- res$par
      #   new_par[2:3] <- theta
      #   new_par[2] <- -new_par[2]
      #   return(likelihood(new_par))
      # }
      
      fix_vec_full <- c(fix_vec, rep(FALSE, n_fixed))

      observed_tmp <- matrix(nrow = length(fix_vec_full), ncol = length(fix_vec_full))

      observed_tmp[!fix_vec_full, !fix_vec_full] <- observed_fisher

      # coeff_tmp <- coeff[2:3]
      # new_observed_fisher <- observed_fisher[2:3,2:3]

      coeff_tmp <- tmp_coeff[2:3]

      new_observed_fisher <- observed_tmp[2:3,2:3]

      change_par <- change_parameterization_graphlme( #new_likelihood,
                                                     model[["alpha"]]-0.5,
                                              coeff_tmp,
                                              hessian = new_observed_fisher, fix_vec[2:3]
                                              )
      
      matern_coeff <- list()
      matern_coeff$random_effects <- change_par$coeff
      names(matern_coeff$random_effects) <- c("sigma", "range")
      matern_coeff$std_random <- change_par$std_random
      time_matern_par_end <- Sys.time()
      time_matern_par <- time_matern_par_end - time_matern_par_start
    } else{

      tmp_coeff <- rep(NA,length(fix_vec))
      tmp_coeff[!fix_vec] <- coeff[1:(sum(!fix_vec))]
      tmp_coeff[fix_vec] <- exp(fix_v_val[!is.na(fix_v_val)])

      matern_coeff <- NULL
      time_matern_par <- NULL
    }

  inv_fisher <- tryCatch(solve(observed_fisher),
                         error = function(e) matrix(NA,
                                                    nrow(observed_fisher),
                                                    ncol(observed_fisher)))
  std_err <- sqrt(diag(inv_fisher))

  fix_vec_full <- c(fix_vec, rep(FALSE,n_fixed))

  std_err_tmp <- rep(NA, length(fix_vec_full))
  std_err_tmp[!fix_vec_full] <- std_err

  std_err <- std_err_tmp

  coeff_random <- tmp_coeff[2:(1+n_random)]
  std_random <- std_err[2:(1+n_random)]

  names(coeff_random) <- par_names

  coeff_meas <- tmp_coeff[1]
  names(coeff_meas) <- "std. dev"

  std_meas <- std_err[1]

  coeff_fixed <- NULL
  if(n_fixed > 0){
    if(sum(!fix_vec)>0){
      coeff_fixed <- res$par[-c(1:(sum(!fix_vec)))]
    } else{
      coeff_fixed <- res$par
    }
    std_fixed <- std_err[(2+n_random):length(fix_vec_full)]
  } else{
    std_fixed <- NULL
  }

  } else{
    time_build_likelihood <- NULL

    coeff_random <- NULL
    std_random <- NULL

    matern_coeff <- NULL
    time_matern_par <- NULL

    if(ncol(X_cov) == 0){
      stop("The model does not have either random nor fixed effects.")
    }

    names_tmp <- colnames(X_cov)
    data_tmp <- cbind(y_graph, X_cov)
    data_tmp <- na.omit(data_tmp)
    start_time <- Sys.time()
    res <- lm(data_tmp[,1] ~ data_tmp[,-1] - 1)
    end_time <- Sys.time()
    time_fit <- end_time - start_time
    time_hessian <- NULL
    time_par <- NULL
    coeff_fixed <- res$coeff
    names(coeff_fixed) <- names_tmp
    sm_temp <- summary(res)
    std_fixed <- sm_temp$coefficients
    rownames(std_fixed) <- names_tmp
    coeff_meas <- sm_temp$sigma
    names(coeff_meas) <- "std. dev"
    std_meas <- NULL
    loglik <- logLik(res)[[1]]
    start_values <- NULL
    fixed_values <- NULL

  }

  if(is.null(coeff_fixed) && is.null(coeff_random)){
    stop("The model does not have either random nor fixed effects.")
  }



  object <- list()
  object$coeff <- list(measurement_error = coeff_meas,
  fixed_effects = coeff_fixed, random_effects = coeff_random)
  object$std_errors <- list(std_meas = std_meas,
        std_fixed = std_fixed, std_random = std_random)
  object$call <- call_graph_lme
  object$terms <- list(fixed_effects = X_cov)
  object$response_data <- list(y = y_graph)
  object$formula <- formula
  object$estimation_method <- optim_method
  # object$parameterization_latent <- parameterization_latent
  object$which_repl <- which_repl
  object$nobs <- sum(graph$.__enclos_env__$private$data[[".group"]] %in% which_repl)
  object$optim_controls <- optim_controls
  object$latent_model <- model
  object$loglik <- loglik
  object$BC <- BC
  object$niter <- res$counts
  object$response_var <- y_term
  object$matern_coeff <- matern_coeff
  object$time_matern_par <- time_matern_par
  object$optim_method <- optim_method
  object$covariates <- cov_term
  object$nV_orig <- nV_orig
  object$fitting_time <- time_fit
  object$time_likelihood <- time_build_likelihood
  object$improve_hessian <- improve_hessian
  object$time_hessian <- time_hessian
  object$parallel <- parallel
  object$time_par <- time_par
  object$start_values <- start_values
  object$fixed_values <- fixed_values
  if(!is.null(graph_bkp$res_dist)){
    object$euclidean <- attr(graph_bkp$res_dist[[".complete"]], "euclidean")
  }
  # if(model_matrix){
    if(ncol(X_cov)>0){
      object$model_matrix <- cbind(y_graph, X_cov)
    } else{
      object$model_matrix <- y_graph
    }
  # }
  rm(graph_bkp)
  graph_bkp <- graph$clone()
  graph_bkp$.__enclos_env__$private$data <- lapply(names_temp, function(i){graph_bkp$.__enclos_env__$private$data[[i]]})
  names(graph_bkp$.__enclos_env__$private$data) <- names_temp
  object$graph <- graph_bkp
  object$df.residual <- object$nobs -(1 + length(object$coeff$fixed_effects) + length(object$coeff$random_effects))
  object$lik_fun <- likelihood_new
  object$mle_par_orig <- res$par


  class(object) <- "graph_lme"
  return(object)

}

#' @name logLik.graph_lme
#' @title Log-likelihood for \code{graph_lme} objects
#' @description computes the log-likelihood for a fitted mixed effects model on
#' metric graphs.
#' @param x Object of class `graph_lme` containing results from the fitted model.
#' @param ... Currently not used.
#' @return Log-likelihood value.
#' @noRd
#' @method logLik graph_lme
#' @export
logLik.graph_lme <- function(object, ...){
  ll <- object$loglik
  attr(ll, "df") <- 1 + length(object$coeff$fixed_effects) + length(object$coeff$random_effects)
  return(ll)
}

#' @name nobs.graph_lme
#' @title Number of observations for \code{graph_lme} objects
#' @description Gets the number of observations of the fitted object.
#' @param x Object of class `graph_lme` containing results from the fitted model.
#' @param ... Currently not used.
#' @return The number of observations.
#' @noRd
#' @method nobs graph_lme
#' @export
nobs.graph_lme <- function(object, ...){
  return(object$nobs)
}


#' @name deviance.graph_lme
#' @title Deviance for \code{graph_lme} objects
#' @description Gets the deviance of `graph_lme` fitted object.
#' @param x Object of class `graph_lme` containing results from the fitted model.
#' @param ... Currently not used.
#' @return Deviance
#' @noRd
#' @method deviance graph_lme
#' @export
deviance.graph_lme <- function(object, ...){
  if(length(object$coeff$random_effects) > 0){
    return(-2*object$loglik)
  } else{
    df_temp <- as.data.frame(object$model_matrix)
    colnames(df_temp)[1] <- as.character(object$response_var)

    fit_lm <- stats::lm(object$formula, data = df_temp)
    return(stats::deviance(fit_lm))
  }
}


#' @name glance.graph_lme
#' @title Glance at a \code{graph_lme} object
#' @aliases glance glance.graph_lme
#' @description Glance accepts a \code{graph_lme} object and returns a
#' [tidyr::tibble()] with exactly one row of model summaries.
#' The summaries are the square root of the estimated variance of the measurement error, residual
#' degrees of freedom, AIC, BIC, log-likelihood,
#' the type of latent model used in the fit and the total number of observations.
#' @param x A \code{graph_lme} object.
#' @param ... Additional arguments. Currently not used.
#' @return A [tidyr::tibble()] with exactly one row and columns:
#' \itemize{
#'   \item `nobs` Number of observations used.
#'   \item `sigma` the square root of the estimated residual variance
#'   \item `logLik` The log-likelihood of the model.
#'   \item `AIC` Akaike's Information Criterion for the model.
#'   \item `BIC` Bayesian Information Criterion for the model.
#'   \item `deviance` Deviance of the model.
#'   \item `df.residual` Residual degrees of freedom.
#'   \item `model.type` Type of latent model fitted.
#'   }
#' @seealso [augment.graph_lme]
#' @method glance graph_lme
#' @export

glance.graph_lme <- function(x, ...){
  alpha <- NULL
  if(x$latent_model$type == "Covariance-Based Matern SPDE Approximation"){
    alpha <- x$coeff$random_effects[[1]]
  } else if(!is.null(x$latent_model$alpha)){
    alpha <- x$latent_model$alpha
  }
  tidyr::tibble(nobs = stats::nobs(x),
                  sigma = as.numeric(x$coeff$measurement_error[[1]]),
                   logLik = as.numeric(stats::logLik(x)), AIC = stats::AIC(x),
                   BIC = stats::BIC(x), deviance = stats::deviance(x),
                   df.residual = stats::df.residual(x),
                   model = x$latent_model$type,
                   alpha = alpha,
                   cov_function = x$latent_model$cov_function_name)
}


#' @name augment.graph_lme
#' @title Augment data with information from a `graph_lme` object
#' @aliases augment augment.graph_lme
#' @description Augment accepts a model object and a dataset and adds information about each observation in the dataset. It includes
#' predicted values in the `.fitted` column, residuals in the `.resid` column, and standard errors for the fitted values in a `.se.fit` column.
#' It also contains the New columns always begin with a . prefix to avoid overwriting columns in the original dataset.
#' @param x A `graph_lme` object.
#' @param newdata A `data.frame` or a `list` containing the covariates, the edge
#' number and the distance on edge for the locations to obtain the prediction. If `NULL`, the fitted values will be given for the original locations where the model was fitted.
#' @param which_repl  Which replicates to obtain the prediction. If `NULL` predictions
#' will be obtained for all replicates. Default is `NULL`.
#' @param edge_number Name of the variable that contains the edge number, the
#' default is `edge_number`.
#' @param distance_on_edge Name of the variable that contains the distance on
#' edge, the default is `distance_on_edge`.
#' @param coord_x Column (or entry on the list) of the `data` that contains
#' the x coordinate. If not supplied, the column with name "coord_x" will be
#' chosen. Will not be used if `Spoints` is not `NULL` or if `data_coords` is
#' `PtE`.
#' @param coord_y Column (or entry on the list) of the `data` that contains
#' the y coordinate. If not supplied, the column with name "coord_x" will be
#' chosen. Will not be used if `Spoints` is not `NULL` or if `data_coords` is
#' `PtE`.
#' @param data_coords To be used only if `Spoints` is `NULL`. It decides which
#' coordinate system to use. If `PtE`, the user must provide `edge_number` and
#' `distance_on_edge`, otherwise if `spatial`, the user must provide
#' `coord_x` and `coord_y`.
#' @param normalized Are the distances on edges normalized?
#' @param sd_post_re Logical indicating whether or not a .sd_post_re column should be added to the augmented output containing the posterior standard deviations of the random effects. 
#' @param se_fit Logical indicating whether or not a .se_fit column should be added to the augmented output containing the standard errors of the fitted values. If `TRUE`, the posterior standard deviations of the random effects will also be returned.
#' @param conf_int Logical indicating whether or not confidence intervals for the posterior mean of the random effects should be built.
#' @param pred_int Logical indicating whether or not prediction intervals for the fitted values should be built. If `TRUE`, the confidence intervals for the posterior random effects will also be built.
#' @param level Level of confidence and prediction intervals if they are constructed.
#' @param no_nugget Should the prediction be done without nugget?
#' @param check_euclidean Check if the graph used to compute the resistance distance has Euclidean edges? The graph used to compute the resistance distance has the observation locations as vertices.
#' @param ... Additional arguments.
#'
#' @return A [tidyr::tibble()] with columns:
#' \itemize{
#'   \item `.fitted` Fitted or predicted value.
#'   \item `.relwrconf` Lower bound of the confidence interval of the random effects, if conf_int = TRUE
#'   \item `.reuprconf` Upper bound of the confidence interval of the random effects, if conf_int = TRUE
#'   \item `.fittedlwrpred` Lower bound of the prediction interval, if conf_int = TRUE
#'   \item `.fitteduprpred` Upper bound of the prediction interval, if conf_int = TRUE
#'   \item `.fixed` Prediction of the fixed effects.
#'   \item `.random` Prediction of the random effects.
#'   \item `.resid` The ordinary residuals, that is, the difference between observed and fitted values.
#'   \item `.std_resid` The standardized residuals, that is, the ordinary residuals divided by the standard error of the fitted values (by the prediction standard error), if se_fit = TRUE or pred_int = TRUE.
#'   \item `.se_fit` Standard errors of fitted values, if se_fit = TRUE.
#'    \item `.sd_post_re` Standard deviation of the posterior mean of the random effects, if se_fit = TRUE.
#'   }
#'
#' @seealso [glance.graph_lme]
#' @method augment graph_lme
#' @export
augment.graph_lme <- function(x, newdata = NULL, which_repl = NULL, sd_post_re = FALSE, se_fit = FALSE, conf_int = FALSE, pred_int = FALSE, level = 0.95, edge_number = "edge_number", distance_on_edge = "distance_on_edge", coord_x = "coord_x", coord_y = "coord_y", data_coords = c("PtE", "spatial"),  normalized = FALSE, no_nugget = FALSE, check_euclidean = FALSE,...) {

  .resid <- FALSE
  if(is.null(newdata)){
    .resid <-  TRUE
  }


  level <- level[[1]]
  if(!is.numeric(level)){
    stop("'level' must be numeric!")
  }
  if(level > 1 || level < 0){
    stop("'level' must be between 0 and 1!")
  }

  if(is.null(newdata)){
    newdata = x$graph$get_data(group = x$graph$get_groups()[1], format="tibble")
  } else{
    newdata <- x$graph$process_data(data = newdata,
    edge_number = edge_number, distance_on_edge = distance_on_edge,
    coord_x = coord_x, coord_y = coord_y, data_coords = data_coords, group = NULL,
    format="tibble", normalized = normalized)
  }

 if(pred_int || se_fit){
    pred <- stats::predict(x, newdata = newdata, which_repl = which_repl, compute_pred_variances = TRUE,
                  posterior_samples = FALSE, edge_number = ".edge_number",
                  distance_on_edge = ".distance_on_edge", normalized = TRUE, return_original_order = FALSE, return_as_list = FALSE, no_nugget = no_nugget, check_euclidean = check_euclidean)
 } else if(conf_int || sd_post_re){
    pred <- stats::predict(x, newdata = newdata, which_repl = which_repl, compute_variances = TRUE,
                  posterior_samples = FALSE, edge_number = ".edge_number",
                  distance_on_edge = ".distance_on_edge", normalized = TRUE, return_original_order = FALSE, return_as_list = FALSE, no_nugget = no_nugget, check_euclidean = check_euclidean)
  } else{
      pred <- stats::predict(x, newdata = newdata, which_repl = which_repl, compute_variances = FALSE, posterior_samples = FALSE, edge_number = ".edge_number",
                  distance_on_edge = ".distance_on_edge", normalized = TRUE, return_original_order = FALSE, return_as_list = FALSE, no_nugget = no_nugget, check_euclidean = check_euclidean)
  }

  newdata[[".fitted"]] <- pred$mean
  newdata[[".fixed"]] <- pred$fe_mean
  newdata[[".random"]] <- pred$re_mean  
  if(.resid){
      newdata[[".resid"]] <- pred$mean - newdata[[as.character(x$response_var)]]
      if(se_fit || pred_int){
        newdata[[".std_resid"]] <- (pred$mean - newdata[[as.character(x$response_var)]])/sqrt(pred$pred_variance)
      }
  }  

  if(se_fit || pred_int){
    newdata[[".se_fit"]] <- sqrt(pred$pred_variance)
  } 

  if(conf_int || sd_post_re || se_fit || pred_int){
    newdata[[".sd_post_re"]] <- sqrt(pred$variance)
  }   


  if(conf_int || pred_int){
    newdata$.relwrconf <- newdata[[".random"]] + stats::qnorm( (1-level)/2 )*newdata[[".sd_post_re"]]
    newdata$.reuprconf <- newdata[[".random"]] + stats::qnorm( (1+level)/2 )*newdata[[".sd_post_re"]]
  }

  if(pred_int){
    newdata$.fittedlwrpred <- newdata[[".fitted"]] + stats::qnorm( (1-level)/2 )*newdata[[".se_fit"]]
    newdata$.fitteduprpred <- newdata[[".fitted"]] + stats::qnorm( (1+level)/2 )*newdata[[".se_fit"]]
  }

  attr(newdata, "euclidean") <- pred$euclidean
  newdata
}




#' @name print.graph_lme
#' @title Print Method for \code{graph_lme} Objects
#' @description Provides a brief description of results related to mixed effects
#' metric graph models.
#' @param x object of class `graph_lme` containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @return No return value. Called for its side effects.
#' @noRd
#' @method print graph_lme
#' @export
print.graph_lme <- function(x, ...) {
  #
  if(inherits(x, "rspde_lme")){
    class(x) <- "rspde_lme"
    return(print(x))
  }
  model_type <- tolower(x$latent_model$type)
  call_name <- switch(model_type,
                      "whittlematern" = {paste0("Latent model - Whittle-Matern with alpha = ",x$latent_model$alpha)},
                      "graphlaplacian" = {paste0("Latent model - graph Laplacian SPDE with alpha = ",x$latent_model$alpha)},
                      "isocov" = {"Latent model - Covariance-based model"},
                      "linearmodel" = {"Linear regression model"}
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
  if(!is.null(coeff_random)){
    print(coeff_random)
    if(!is.null(x$matern_coeff$random_effects)){
      cat(paste0("\n", "Random effects (Matern parameterization):", "\n"))
      print(x$matern_coeff$random_effects)
    }
  } else{
    message("No random effects")
  }
  cat("\n")
  cat(paste0("Measurement error:", "\n"))
  print(x$coeff$measurement_error)
}


#' @name summary.graph_lme
#' @title Summary Method for \code{graph_lme} Objects
#' @description Function providing a summary of results related to metric graph
#' mixed effects regression models.
#' @param object an object of class `graph_lme` containing results from the
#' fitted model.
#' @param all_times Show all computed times.
#' @param ... not used.
#' @return An object of class \code{summary_graph_lme} containing information
#' about a *graph_lme* object.
#' @method summary graph_lme
#' @export
summary.graph_lme <- function(object, all_times = FALSE, ...) {
  if(inherits(object, "rspde_lme")){
    class(object) <- "rspde_lme"
    return(summary(object))
  }
  ans <- list()

  nfixed <- length(object$coeff$fixed_effects)
  nrandom <- length(object$coeff$random_effects)
  model_type <- tolower(object$latent_model$type)
  call_name <- switch(model_type,
                      "whittlematern" = {paste0("Latent model - Whittle-Matern with alpha = ",object$latent_model$alpha)},
                      "graphlaplacian" = {paste0("Latent model - graph Laplacian SPDE with alpha = ",object$latent_model$alpha)},
                      "isocov" = {"Latent model - Covariance-based model"},
                      "linearmodel" = {"Linear regression model"}
  )

  coeff_fixed <- object$coeff$fixed_effects
  coeff_random <- object$coeff$random_effects#
  coeff_meas <- object$coeff$measurement_error

  SEr_fixed <- object$std_errors$std_fixed
  SEr_random <- object$std_errors$std_random
  SEr_meas <- object$std_errors$std_meas

  if(model_type %in% c("whittlematern", "graphlaplacian")){
    coeff <- c(coeff_fixed, coeff_random, object$matern_coeff$random_effects, coeff_meas)
    SEr <- c(SEr_fixed,SEr_random, object$matern_coeff$std_random, SEr_meas)
  } else{
    coeff <- c(coeff_fixed, coeff_random, coeff_meas)
    SEr <- c(SEr_fixed,SEr_random, SEr_meas)
  }

  if(model_type != "linearmodel"){
    tab <- cbind(coeff, SEr, coeff / SEr, 2 * stats::pnorm(-abs(coeff / SEr)))
    colnames(tab) <- c("Estimate", "Std.error", "z-value", "Pr(>|z|)")
    rownames(tab) <- names(coeff)
    if(model_type %in% c("whittlematern", "graphlaplacian")){
      tab <- list(fixed_effects = tab[seq.int(length.out = nfixed), , drop = FALSE], random_effects = tab[seq.int(length.out = nrandom) + nfixed, , drop = FALSE],
      random_effects_matern = tab[seq.int(length.out = nrandom) + nrandom + nfixed, , drop = FALSE],
      meas_error = tab[seq.int(length.out = 1) + nfixed+2*nrandom, , drop = FALSE])
    } else{
      tab <- list(fixed_effects = tab[seq.int(length.out = nfixed), , drop = FALSE], random_effects = tab[seq.int(length.out = nrandom) + nfixed, , drop = FALSE],
      meas_error = tab[seq.int(length.out = 1) + nfixed+nrandom, , drop = FALSE])
    }
  } else{
    tab <- list(fixed_effects = SEr_fixed, coeff_meas = coeff_meas)
  }




  ans$coefficients <- tab

  ans$all_times <- all_times

  ans$model_type <- model_type

  ans$call_name <- call_name

  ans$call <- object$call

  ans$loglik <- object$loglik

  ans$optim_method <- object$optim_method

  ans$niter <- object$niter

  ans$fitting_time <- object$fitting_time

  ans$improve_hessian <- object$improve_hessian

  ans$time_hessian <- object$time_hessian

  ans$parallel <- object$parallel

  ans$time_par <- object$time_par

  ans$time_matern_par <- object$time_matern_par

  ans$time_likelihood <- object$time_likelihood


  class(ans) <- "summary_graph_lme"
  ans
}

#' @name print.summary_graph_lme
#' @title Print method for \code{summary_graph_lme} objects
#' @description Provides a brief description of results related to metric graph
#' mixed effects regression models.
#' @param x object of class "summary_graph_lme" containing results of summary
#' method applied to a fitted model.
#' @param ... further arguments passed to or from other methods.
#' @return No return value. Called for its side effects.
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
  model_type <- tolower(x$model_type)
  #
  if(model_type != "linearmodel"){
      if (NROW(tab$fixed_effects)) {
        cat(paste0("\nFixed effects:\n"))
        stats::printCoefmat(tab[["fixed_effects"]], digits = digits,
                            signif.legend = FALSE)
      } else {
        message("\nNo fixed effects. \n")
      }
      #
      if (NROW(tab$random_effects)) {
        cat(paste0("\nRandom effects:\n"))
        stats::printCoefmat(tab[["random_effects"]][,1:3], digits = digits,
                            signif.legend = FALSE)
      } else {
        message("\nNo random effects. \n")
      }
      if (NROW(tab$random_effects_matern)) {
        cat(paste0("\nRandom effects (Matern parameterization):\n"))
        stats::printCoefmat(tab[["random_effects_matern"]][,1:3], digits = digits,
                            signif.legend = FALSE)
      }
      #
      cat(paste0("\nMeasurement error:\n"))
        stats::printCoefmat(tab[["meas_error"]][1,1:3,drop = FALSE], digits = digits,
                            signif.legend = FALSE)
  } else{
        cat(paste0("\nFixed effects:\n"))
        stats::printCoefmat(tab[["fixed_effects"]], digits = digits,
                            signif.legend = FALSE)

        message("\nNo random effects. \n")
        cat(paste0("\nMeasurement error:\n"))
        print(tab$coeff_meas)

  }
  #
  if (getOption("show.signif.stars")) {
    cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n\n")
  }
  #

  cat("Log-Likelihood: ", x$loglik,"\n")
  if(model_type != "linearmodel"){
    cat(paste0("Number of function calls by 'optim' = ", x$niter[1],"\n"))
    cat(paste0("Optimization method used in 'optim' = ", x$optim_method,"\n"))
    cat(paste0("\nTime used to:"))
    if(x$all_times){
      cat("\t build the likelihood = ",
          paste(trunc(x$time_likelihood[[1]] * 10^5)/10^5, attr(x$time_likelihood, "units"),"\n"))
      cat("\t compute Matern parameterization = ",
          paste(trunc(x$time_matern_par[[1]] * 10^5)/10^5,attr(x$time_likelihood, "units"),"\n"))
    }
    cat("\t fit the model = ",
        paste(trunc(x$fitting_time[[1]] * 10^5)/10^5,attr(x$fitting_time, "units"),"\n"))
    if(x$improve_hessian){
    cat(paste0("\t compute the Hessian = ",
               paste(trunc(x$time_hessian[[1]] * 10^5)/10^5,attr(x$time_hessian, "units"),"\n")))
    }
    if(x$parallel){
    cat(paste0("\t set up the parallelization = ",
               paste(trunc(x$time_par[[1]] * 10^5)/10^5,attr(x$time_par, "units"),"\n")))
    }
  }
}



#' @name predict.graph_lme
#' @title Prediction for a mixed effects regression model on a metric graph
#' @param object The fitted object with the `graph_lme()` function.
#' @param newdata A `data.frame` or a `list` containing the covariates, the edge
#' number and the distance on edge for the locations to obtain the prediction. Observe that you should not provide the locations for each replicate. Only a single set of locations and covariates, and the predictions for the different replicates will be obtained for this same set of locations.
#' @param mesh Obtain predictions for mesh nodes? The graph must have a mesh and should not have covariates.
#' @param mesh_h If the graph does not have a mesh, one will be created with this
#' value of 'h'.
#' @param which_repl Which replicates to obtain the prediction. If `NULL` predictions
#' will be obtained for all replicates. Default is `NULL`.
#' @param compute_variances Set to TRUE to compute the kriging variances.
#' @param compute_pred_variances Set to TRUE to compute the prediction variances. Will only be computed if newdata is `NULL`.
#' @param posterior_samples If `TRUE`, posterior samples for the random effect will be returned.
#' @param pred_samples If `TRUE`, prediction samples for the response variable will be returned. Will only be computed if newdata is `NULL`.
#' @param n_samples Number of samples to be returned. Will only be used if
#' `sampling` is `TRUE`.
#' @param edge_number Name of the variable that contains the edge number, the
#' default is `edge_number`.
#' @param distance_on_edge Name of the variable that contains the distance on
#' edge, the default is `distance_on_edge`.
#' @param normalized Are the distances on edges normalized?
#' @param no_nugget Should the prediction be carried out without the nugget?
#' @param return_as_list Should the means of the predictions and the posterior
#' samples be returned as a list, with each replicate being an element?
#' @param return_original_order Should the results be return in the original
#' (input) order or in the order inside the graph?
#' @param check_euclidean Check if the graph used to compute the resistance distance has Euclidean edges? The graph used to compute the resistance distance has the observation locations as vertices.
#' @param ... Not used.
#' @param data `r lifecycle::badge("deprecated")` Use `newdata` instead.
#' @return A list with elements `mean`, which contains the means of the
#' predictions, `fe_mean`, which is the prediction for the fixed effects, `re_mean`, which is the prediction for the random effects, `variance` (if `compute_variance` is `TRUE`), which contains the
#' posterior variances of the random effects, `samples` (if `posterior_samples` is `TRUE`),
#' which contains the posterior samples.
#' @export
#' @method predict graph_lme

predict.graph_lme <- function(object,
                              newdata = NULL,
                              mesh = FALSE,
                              mesh_h = 0.01,
                              which_repl = NULL,
                              compute_variances = FALSE,
                              compute_pred_variances = FALSE,
                              posterior_samples = FALSE,
                              pred_samples = FALSE,
                              n_samples = 100,
                              edge_number = "edge_number",
                              distance_on_edge = "distance_on_edge",
                              normalized = FALSE,
                              no_nugget = FALSE,
                              return_as_list = FALSE,
                              return_original_order = TRUE,
                              check_euclidean = TRUE,
                               ...,
                               data = deprecated()) {

                                
  repl <- which_repl
  if (lifecycle::is_present(data)) {
    if (is.null(newdata)) {
      lifecycle::deprecate_warn("1.2.0", "predict(data)", "predict(newdata)",
        details = c("`data` was provided but not `newdata`. Setting `newdata <- data`.")
      )
      newdata <- data
    } else {
      lifecycle::deprecate_warn("1.2.0", "predict(data)", "predict(newdata)",
        details = c("Both `newdata` and `data` were provided. Only `newdata` will be considered.")
      )
    }
    data <- NULL
  }

  data <- newdata
  if(is.null(data)){
    if(!mesh){
      stop("If 'mesh' is false, you should supply data!")
    }
  }

  if(inherits(object, "rspde_lme")){
    class(object) <- "rspde_lme"
    return(stats::predict(object = object,
              newdata = newdata,
              mesh = mesh, which_repl = which_repl,
              compute_variances = compute_variances, posterior_samples = posterior_samples,
                               n_samples = n_samples, edge_number = edge_number,
                               distance_on_edge = distance_on_edge, normalized = normalized, return_as_list = return_as_list, return_original_order = return_original_order)
              )
  }

  out <- list()

  coeff_fixed <- object$coeff$fixed_effects
  coeff_random <- object$coeff$random_effects
  coeff_meas <- object$coeff$measurement_error

  BC <- object$BC

  graph_bkp <- object$graph$clone()


  X_cov_initial <- stats::model.matrix(object$covariates, graph_bkp$.__enclos_env__$private$data)
  if(ncol(X_cov_initial) > 0){
    if(mesh){
      stop("In the presence of covariates, you should provide the data, including the covariates at the prediction locations.")
    }
  }


  if(sum(duplicated(cbind(data[[edge_number]], data[[distance_on_edge]]))) > 0){
    warning("There are duplicated locations for prediction, we will try to process the data to extract the unique locations,
    along with the corresponding covariates.")
    cov_names <- attr(object$covariates,"term.labels")
    data <- data[c(edge_number,distance_on_edge,cov_names)]
    data <- unique(data)
    if(sum(duplicated(cbind(data[[edge_number]], data[[distance_on_edge]]))) > 0){
      stop("Data processing failed, please provide a data with unique locations.")
    }
  }

  if(!mesh){
    n_prd <- length(data[[edge_number]])
    data[["__dummy_var"]] <- rep(0, n_prd)
    # Convert data to normalized
    if(!normalized){
      data[[distance_on_edge]] <- data[[distance_on_edge]] / graph_bkp$edge_lengths[data[[edge_number]]]
    }
  } else{
    if(is.null(graph_bkp$mesh)){
      graph_bkp$build_mesh(h = mesh_h)
    }
    data <- list()
    n_prd <- nrow(graph_bkp$mesh$VtE)
    data[["__dummy_var"]] <- rep(0, n_prd)
    data[[edge_number]] <- graph_bkp$mesh$VtE[,1]
    data[[distance_on_edge]] <- graph_bkp$mesh$VtE[,2]
    normalized <- TRUE
  }


    ord_idx <- order(data[[edge_number]], data[[distance_on_edge]])


  if(!is.null(data[[as.character(object$response_var)]])){
    data[[as.character(object$response_var)]] <- NULL
  }

  data_graph_temp <- list()
  idx_group1 <-  graph_bkp$.__enclos_env__$private$data[[".group"]] == graph_bkp$.__enclos_env__$private$data[[".group"]][1]
  data_graph_temp[[edge_number]] <- graph_bkp$.__enclos_env__$private$data[[".edge_number"]][idx_group1]
  data_graph_temp[[distance_on_edge]] <- graph_bkp$.__enclos_env__$private$data[[".distance_on_edge"]][idx_group1]
  data_graph_temp[[as.character(object$response_var)]] <- graph_bkp$.__enclos_env__$private$data[[as.character(object$response_var)]][idx_group1]
  # data_graph_temp <- as.data.frame(data_graph_temp)
  # colnames(data_graph_temp)[1:2] <- c(edge_number, distance_on_edge)

  data_prd_temp <- list()
  data_prd_temp[[edge_number]] <- data[[edge_number]]
  data_prd_temp[[distance_on_edge]] <- data[[distance_on_edge]]
  data_prd_temp[["included"]] <- TRUE

  temp_merge <- merge(data_prd_temp, data_graph_temp, all = TRUE)

  temp_merge <- temp_merge[!is.na(temp_merge[["included"]]),]

  temp_merge[["included"]] <- NULL

  data <- merge(temp_merge, data)

  rm(temp_merge)
  # rm(data_prd_temp)
  rm(data_graph_temp)

  old_data <- graph_bkp$.__enclos_env__$private$data

  data[[".group"]] <- old_data[[".group"]][1]

  graph_bkp$clear_observations()

  graph_bkp$add_observations(data = data, edge_number = edge_number,
                             distance_on_edge = distance_on_edge,
                             normalized = TRUE, group = ".group", verbose = 0,
                  suppress_warnings = TRUE)

  graph_bkp$add_observations(data = old_data, edge_number = ".edge_number",
                             distance_on_edge = ".distance_on_edge",
                             group = ".group", normalized = TRUE, verbose = 0,
                  suppress_warnings = TRUE)

  graph_bkp$.__enclos_env__$private$data[["__dummy_ord_var"]] <- 1:length(graph_bkp$.__enclos_env__$private$data[[".edge_number"]])


  model_type <- object$latent_model

  sigma.e <- coeff_meas[[1]]
  sigma_e <- sigma.e

  # graph_bkp$observation_to_vertex(mesh_warning=FALSE)

  n <- sum(graph_bkp$.__enclos_env__$private$data[[".group"]] == graph_bkp$.__enclos_env__$private$data[[".group"]][1])

  ##
  repl_vec <- graph_bkp$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    u_repl <- unique(graph_bkp$.__enclos_env__$private$data[[".group"]])
  } else{
    u_repl <- unique(repl)
  }

  ##

  X_cov_pred <- stats::model.matrix(object$covariates, graph_bkp$.__enclos_env__$private$data)

  if(all(dim(X_cov_pred) == c(0,1))){
    X_cov_pred <- matrix(1, nrow = length(graph_bkp$.__enclos_env__$private$data[[".group"]]), ncol=1)
  }
  if(ncol(X_cov_pred) > 0){
    mu <- X_cov_pred %*% coeff_fixed
  } else{
    mu <- matrix(0, nrow = length(graph_bkp$.__enclos_env__$private$data[[".group"]]), ncol=1)
  }

  Y <- graph_bkp$.__enclos_env__$private$data[[as.character(object$response_var)]] - mu


  if(!is.null(graph_bkp$.__enclos_env__$private$data[["__dummy_var"]])){
      idx_prd <- !is.na(graph_bkp$.__enclos_env__$private$data[["__dummy_var"]][1:n])
  } else {
      idx_prd <- !is.na(graph_bkp$.__enclos_env__$private$data[["X__dummy_var"]][1:n])
  }

  n_prd <- sum(idx_prd)

  edge_nb <- graph_bkp$.__enclos_env__$private$data[[".edge_number"]][1:n][idx_prd]
  dist_ed <- graph_bkp$.__enclos_env__$private$data[[".distance_on_edge"]][1:n][idx_prd]

  ## construct Q


  if(tolower(model_type$type) == "whittlematern"){
    tau <- object$coeff$random_effects[1]
    # if(object$parameterization_latent == "spde"){
      kappa <- object$coeff$random_effects[2]

      if (compute_variances || posterior_samples || no_nugget || compute_pred_variances){
              graph_bkp$observation_to_vertex(mesh_warning=FALSE)
      }
    # } else{
    #   kappa <- sqrt(8 * 0.5) / object$coeff$random_effects[2]
    # }

      # if(model_type$alpha == 1){
      #     Q <- spde_precision(kappa = kappa, sigma = sigma,
      #                       alpha = 1, graph = graph_bkp)
      # }
      # else{
      #   PtE <- graph_bkp$get_PtE()
      #   n.c <- 1:length(graph_bkp$CoB$S)
      #   Q <- spde_precision(kappa = kappa, sigma = sigma, alpha = 2,
      #                       graph = graph_bkp, BC = BC)
      #   Qtilde <- (graph_bkp$CoB$T) %*% Q %*% t(graph_bkp$CoB$T)
      #   Qtilde <- Qtilde[-n.c,-n.c]
      #   Sigma.overdetermined  = t(graph_bkp$CoB$T[-n.c,]) %*% solve(Qtilde) %*%
      #     (graph_bkp$CoB$T[-n.c,])
      #   index.obs <- 4 * (PtE[,1] - 1) + 1.0 * (abs(PtE[, 2]) < 1e-16) +
      #     3.0 * (abs(PtE[, 2]) > 1e-16)
      #   Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
      #   Q <- solve(Sigma)
      # }

  } else if(tolower(model_type$type) == "graphlaplacian"){
    tau <- object$coeff$random_effects[1]
    #nV before
    nV_temp <- object$nV_orig
    graph_bkp$observation_to_vertex(mesh_warning=FALSE)
    if(graph_bkp$nV > nV_temp){
      warning("There are prediction locations outside of the observation locations. Refit the model with all the locations you want to obtain predictions.")
    }
    graph_bkp$compute_laplacian(full=TRUE)
    # if(object$parameterization_latent == "spde"){
      kappa <- object$coeff$random_effects[2]
    # } else{
    #   kappa <- sqrt(8 * 0.5) / object$coeff$random_effects[2]
    # }

      if(model_type$alpha == 1){
        # Q <- (kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]) * tau^2
        Q <- (kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]) * tau^2
      } else{
        # Q <- kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]
        Q <- kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]
        Q <- Q %*% Q * tau^2
      }

  } else if(tolower(model_type$type) == "isocov"){
      graph_bkp$observation_to_vertex(mesh_warning=FALSE)
      if(is.character(model_type$cov_function)){
        sigma <- object$coeff$random_effects[1]
        kappa <- object$coeff$random_effects[2]
        # if(model_type$cov_function == "alpha1"){
        #   Q <- spde_precision(kappa = kappa, sigma = sigma,
        #                     alpha = 1, graph = graph_bkp)
        # } else if(model_type$cov_function == "alpha2"){
        #   PtE <- graph_bkp$get_PtE()
        #   n.c <- 1:length(graph_bkp$CoB$S)
        #   Q <- spde_precision(kappa = kappa, sigma = sigma, alpha = 2,
        #                       graph = graph_bkp, BC = BC)
        #   Qtilde <- (graph_bkp$CoB$T) %*% Q %*% t(graph_bkp$CoB$T)
        #   Qtilde <- Qtilde[-n.c,-n.c]
        #   Sigma.overdetermined  = t(graph_bkp$CoB$T[-n.c,]) %*% solve(Qtilde) %*%
        #     (graph_bkp$CoB$T[-n.c,])
        #   index.obs <- 4 * (PtE[,1] - 1) + 1.0 * (abs(PtE[, 2]) < 1e-14) +
        #     3.0 * (abs(PtE[, 2]) > 1e-14)
        #   Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
        #   Q <- solve(Sigma)
        # } else
        if(model_type$cov_function == "GL1"){
              #nV before
              tau <- object$coeff$random_effects[1]
              nV_temp <- object$nV_orig
              # graph_bkp$observation_to_vertex(mesh_warning=FALSE)
              if(graph_bkp$nV > nV_temp){
                warning("There are prediction locations outside of the observation locations. Refit the model with all the locations you want to obtain predictions.")
              }
              graph_bkp$compute_laplacian(full=TRUE)
              Q <- (kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]) * tau^2
        } else if(model_type$cov_function == "GL2"){
              #nV before
              tau <- object$coeff$random_effects[1]
              nV_temp <- object$nV_orig
              # graph_bkp$observation_to_vertex(mesh_warning=FALSE)
              if(graph_bkp$nV > nV_temp){
                warning("There are prediction locations outside of the observation locations. Refit the model with all the locations you want to obtain predictions.")
              }
              graph_bkp$compute_laplacian(full=TRUE)
              Q <- kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]
              Q <- Q %*% Q * tau^2
        # } else if(model_type$cov_function == "exp_covariance"){
        #           graph_bkp$compute_resdist(full = TRUE)
        #           Sigma <- as.matrix(exp_covariance(graph_bkp$res_dist[[1]], c(sigma,kappa)))
        #           Q <- solve(Sigma)
        # }
        }
      } else{
        graph_bkp$compute_resdist(full = TRUE, check_euclidean = check_euclidean)
        cov_function <- model_type$cov_function
        Sigma <- as.matrix(cov_function(graph_bkp$res_dist[[".complete"]], coeff_random))
      }
  }



  ord_var <- graph_bkp$.__enclos_env__$private$data[["__dummy_ord_var"]]

  # gap <- dim(Q)[1] - n

  ## compute Q_x|y
  # A <- Matrix::Diagonal(dim(Q)[1])[(gap+1):dim(Q)[1], ]
  # if(tolower(model_type$type) == "isocov"){
  #   # A <- Matrix::Diagonal(dim(Q)[1])
  #   # A[graph_bkp$data[["__dummy_ord_var"]],] <- A
  #   A <- Matrix::Diagonal(dim(Q)[1])[graph_bkp$data[["__dummy_ord_var"]], ]
  #   A <- Matrix::Diagonal(dim(Q)[1])[graph_bkp$PtV, ]

  #   print(graph_bkp$PtV)
  # } else{
  #   A <- Matrix::Diagonal(dim(Q)[1])[graph_bkp$PtV, ]
  # }

  cond_aux1 <- (tolower(model_type$type) == "whittlematern")
  cond_aux2 <- (tolower(model_type$type) == "isocov" && is.character(model_type$cov_function))
  if(cond_aux2){
    cond_aux2 <- (model_type$cov_function == "WM2" || model_type$cov_function == "WM1")
  }
  cond_wm <- cond_aux1 || cond_aux2

  cond_isocov <- (tolower(model_type$type) == "isocov" && !is.character(model_type$cov_function))

  if(!cond_wm && !cond_isocov){
    A <- Matrix::Diagonal(dim(Q)[1])[graph_bkp$PtV, ]
  }

  idx_obs_full <- as.vector(!is.na(Y))

  # idx_obs_full <- !is.na(graph_bkp$data[[as.character(object$response)]])

  if(return_original_order){
          dist_ed[ord_idx] <- dist_ed
          edge_nb[ord_idx] <- edge_nb
  }

  if(!return_as_list){
    out$distance_on_edge <- rep(dist_ed,length(u_repl))
    out$edge_number <- rep(edge_nb,length(u_repl))
  }

  if(cond_wm){
    PtE_full <- graph_bkp$get_PtE()
    PtE_pred <- PtE_full[idx_prd,]
  }

  cond_alpha2 <- FALSE
  cond_alpha1 <- FALSE
  if(cond_aux1){
    if(model_type$alpha == 2){
      cond_alpha2 <- TRUE
    } else {
      cond_alpha1 <- TRUE
    }
  }
  if(cond_aux2){
    if(model_type$cov_function == "WM2"){
      cond_alpha2 <- TRUE
    } else{
      cond_alpha1 <- TRUE
    }
  }

  if(compute_variances || posterior_samples || no_nugget || compute_pred_variances){
    if(cond_wm){
      tau <- object$coeff$random_effects[1]
      if(cond_alpha1){
        Q <- spde_precision(kappa = kappa, tau = tau,
                          alpha = 1, graph = graph_bkp, BC = BC)
        A <- Matrix::Diagonal(dim(Q)[1])[graph_bkp$PtV, ]
      } else{
        if(is.null(graph_bkp$CoB)){
          graph_bkp$buildC(2)
        }
        PtE <- graph_bkp$get_PtE()
        n.c <- 1:length(graph_bkp$CoB$S)
        Q <- spde_precision(kappa = kappa, tau = tau, alpha = 2,
                            graph = graph_bkp, BC = BC)
        Qtilde <- (graph_bkp$CoB$T) %*% Q %*% t(graph_bkp$CoB$T)
        Qtilde <- Qtilde[-n.c,-n.c]
        Sigma.overdetermined  = t(graph_bkp$CoB$T[-n.c,]) %*% solve(Qtilde) %*%
          (graph_bkp$CoB$T[-n.c,])
        index.obs <- 4 * (PtE[,1] - 1) + 1.0 * (abs(PtE[, 2]) < 1e-14) +
          3.0 * (abs(PtE[, 2]) > 1e-14)
        Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
        # Q <- solve(Sigma)
        # A <- Matrix::Diagonal(dim(Q)[1]) #[graph_bkp2$PtV, ]
      }
    }
  }



  for(repl_y in u_repl){
    if(return_as_list){
      out$distance_on_edge[[repl_y]] <- dist_ed
      out$edge_number[[repl_y]] <- edge_nb
    }
    idx_repl <- graph_bkp$.__enclos_env__$private$data[[".group"]] == repl_y

    idx_obs <- idx_obs_full[idx_repl]

    y_repl <- Y[idx_repl]
    y_repl <- y_repl[idx_obs]

    if(!cond_wm && !cond_isocov){

        if(!no_nugget){
          Q_xgiveny <- t(A[idx_obs,]) %*% A[idx_obs,]/sigma_e^2 + Q
          mu_krig <- solve(Q_xgiveny,as.vector(t(A[idx_obs,]) %*% y_repl / sigma_e^2))
          mu_krig <- A[idx_prd,] %*% mu_krig
        } else{
          QiAt <- solve(Q, t(A[idx_obs,]))
          AQiA <- A[idx_obs,] %*% QiAt
          mu_krig <- solve(Q, t(A[idx_obs,]) %*% solve(AQiA, y_repl))
          mu_krig <- as.vector(A[idx_prd,] %*% mu_krig)
        }


        # mu_krig <- mu_krig[(gap+1):length(mu_krig)]

        mu_fe <- mu[idx_repl, , drop = FALSE]
        mu_fe <- mu_fe[idx_prd, , drop=FALSE]
        mu_re <- mu_krig[ord_idx]

        mu_krig <- mu_fe + mu_re
    } else if (cond_wm){

      PtE_obs <- PtE_full[idx_obs,]

      if(cond_alpha2){
        if(is.null(graph_bkp$CoB)){
          graph_bkp$buildC(2)
        }
          if(!no_nugget){
                      mu_krig <- posterior_mean_obs_alpha2(c(sigma.e,tau,kappa),
                        graph = graph_bkp, PtE_resp = PtE_obs, resp = y_repl,
                        PtE_pred = PtE_pred, no_nugget = no_nugget)
                      mu_re <- mu_krig #[ord_idx]
          } else{
                cov_loc <- Sigma[idx_prd, idx_obs]
                cov_Obs <- Sigma[idx_obs, idx_obs]    
                mu_krig <- cov_loc %*%  solve(cov_Obs, y_repl)
                mu_re <- mu_krig
          }

        mu_fe <- mu[idx_repl, , drop = FALSE]
        mu_fe <- mu_fe[idx_prd, , drop=FALSE]


        mu_krig <- mu_fe + mu_re

      } else{          

        if(!no_nugget){
          mu_krig <- posterior_mean_obs_alpha1(c(sigma.e,tau,kappa),
                        graph = graph_bkp, PtE_resp = PtE_obs, resp = y_repl,
                        PtE_pred = PtE_pred, no_nugget = no_nugget)

                        mu_re <- mu_krig[ord_idx]      
               
        } else{
          QiAt <- solve(Q, t(A[idx_obs,]))
          AQiA <- A[idx_obs,] %*% QiAt
          mu_krig <- solve(Q, t(A[idx_obs,]) %*% solve(AQiA, y_repl))
          mu_krig <- as.vector(A[idx_prd,] %*% mu_krig)
          mu_re <- mu_krig
        }

          # Q <- spde_precision(kappa = kappa, tau = tau,
          #                 alpha = 1, graph = graph_bkp)
          # A <- Matrix::Diagonal(dim(Q)[1])[graph_bkp$PtV, ]

          # Q_xgiveny <- t(A[idx_obs,]) %*% A[idx_obs,]/sigma_e^2 + Q
          # mu_krig <- solve(Q_xgiveny,as.vector(t(A[idx_obs,]) %*% y_repl / sigma_e^2))
          # mu_krig <- A[idx_prd,] %*% mu_krig

          # print(mu_krig)
      

          mu_re <- mu_krig[ord_idx]      
          mu_re <- mu_krig   
          mu_fe <- mu[idx_repl, , drop = FALSE]
          mu_fe <- mu_fe[idx_prd, , drop=FALSE]

          mu_krig <- mu_fe + mu_re
      }
    } else {

        graph_bkp$compute_resdist(full=TRUE, check_euclidean=TRUE)

        Sigma <- as.matrix(cov_function(graph_bkp$res_dist[[".complete"]], coeff_random))

        # A <- Matrix::Diagonal(dim(Sigma)[1])[graph_bkp$PtV, ]

        if(!no_nugget){
          # Q_tmp <- solve(Sigma)
          # Q_tmp <- Q_tmp + t(A[idx_obs,]) %*% A[idx_obs,] /sigma_e^2   
          # mu_krig <- A[idx_prd,]%*%solve(Q_tmp, t(A[idx_prd,])%*%y_repl/sigma_e^2)

          # nV <- graph_bkp$nV - nrow(graph_bkp$get_PtE())          
          # idx_obs_tmp <- c(rep(FALSE,nV), idx_obs)
          # idx_prd_tmp <- c(rep(FALSE,nV), idx_prd)
          idx_obs_tmp <- idx_obs
          idx_prd_tmp <- idx_prd          
          cov_loc <- Sigma[idx_prd_tmp, idx_obs_tmp]
          cov_Obs <- Sigma[idx_obs_tmp, idx_obs_tmp]    
          diag(cov_Obs) <- diag(cov_Obs) + sigma_e^2      
          mu_krig <- cov_loc %*%  solve(cov_Obs, y_repl)
        } else{
          nV <- graph_bkp$nV - nrow(graph_bkp$get_PtE())          
          # idx_obs_tmp <- c(rep(FALSE,nV), idx_obs)
          # idx_prd_tmp <- c(rep(FALSE,nV), idx_prd)
          idx_obs_tmp <- idx_obs
          idx_prd_tmp <- idx_prd                      
          cov_loc <- Sigma[idx_prd_tmp, idx_obs_tmp]
          cov_Obs <- Sigma[idx_obs_tmp, idx_obs_tmp]          
          mu_krig <- cov_loc %*%  solve(cov_Obs, y_repl)
        }

        # Observe that the "fixed-effects" mean has been subtracted from y_repl


        mu_re <- mu_krig
        mu_fe <- mu[idx_repl, , drop = FALSE]
        mu_fe <- mu_fe[idx_prd, , drop=FALSE]

        mu_krig <- mu_fe + mu_re

    }


    mean_tmp <- as.vector(mu_krig)
    mean_fe_tmp <- as.vector(mu_fe)
    mean_re_tmp <- as.vector(mu_re)

    if(return_original_order){
        mean_tmp[ord_idx] <- mean_tmp
        mean_fe_tmp[ord_idx] <- mean_fe_tmp
        mean_re_tmp[ord_idx] <- mean_re_tmp
      # var_tmp[ord_idx] <- var_tmp
    }

    if(!return_as_list){
      out$fe_mean <- c(out$fe_mean, mean_fe_tmp)
      out$re_mean <- c(out$re_mean, mean_re_tmp)
      out$mean <- c(out$mean, mean_tmp)
      out$repl <- c(out$repl, rep(repl_y,n_prd))

    } else{
      out$mean[[repl_y]] <- mean_tmp
      out$fe_mean[[repl_y]] <- mean_fe_tmp
      out$re_mean[[repl_y]] <- mean_re_tmp
    }

    if(compute_variances || posterior_samples || compute_pred_variances){
      if(cond_wm){
        if(!cond_alpha2){
            Q_xgiveny <- t(A[idx_obs,]) %*% A[idx_obs,]/sigma_e^2 + Q
        }
      }
    }

    if (compute_variances || compute_pred_variances) {
        if(!cond_isocov){
          if(cond_alpha2){
            if(!no_nugget){
                cov_loc <- Sigma[idx_prd, idx_obs]
                cov_Obs <- Sigma[idx_obs, idx_obs]    
                diag(cov_Obs) <- diag(cov_Obs) + sigma_e^2    
                post_cov_tmp <- cov_loc %*%  solve(cov_Obs, t(cov_loc))
                var_tmp <- diag(Sigma[idx_prd, idx_prd] - post_cov_tmp)
                if(compute_pred_variances  || pred_samples) {
                  pred_cov <- Matrix::Diagonal(dim(cov_loc)[1], x=sigma_e^2) -cov_loc + post_cov_tmp
                  pred_var_tmp <- diag(pred_cov)
                }
            } else{
              cov_tmp <- Sigma[idx_prd, idx_prd] - Sigma[idx_prd, idx_obs] %*% solve(Sigma[idx_obs, idx_obs],t(Sigma[idx_prd, idx_obs]))
              var_tmp <- diag(cov_tmp)
              var_tmp <- ifelse(var_tmp < 0, 0, var_tmp) # possible numerical errors
              if(compute_pred_variances  || pred_samples) {
                  pred_cov <- cov_tmp
                  pred_var_tmp <- diag(pred_cov)
                }
            }
          } else{
            if(!no_nugget){
              post_cov <- A[idx_prd,]%*%solve(Q_xgiveny, t(A[idx_prd,]))
              var_tmp <- diag(post_cov)
              # var_tmp <- var_tmp * (var_tmp > 0)
                if(compute_pred_variances  || pred_samples) {
                  Q_x_tmp <- A[idx_prd,]%*%solve(Q, t(A[idx_prd,]))
                  pred_cov <- sigma_e^2 - Q_x_tmp + Q_x_tmp %*% solve(Q_x_tmp + Matrix::Diagonal(dim(Q_x_tmp)[1], x=sigma_e^2), t(Q_x_tmp))
                  pred_var_tmp <- diag(pred_cov)
                }           
            } else{
                QiAt <- solve(Q, t(A[idx_obs,]))
                AQiA <- A[idx_obs,] %*% QiAt            
                M <- Q - QiAt %*% solve(AQiA, t(QiAt))
                cov_tmp <- A[idx_prd,] %*% M %*% t(A[idx_prd,])
                var_tmp <- diag(cov_tmp)
              if(compute_pred_variances  || pred_samples) {
                  pred_cov <- cov_tmp
                  pred_var_tmp <- diag(pred_cov)
                }                
            }
          }
        } else{
          # nV <- graph_bkp$nV - nrow(graph_bkp$get_PtE())          
          # idx_obs_tmp <- c(rep(FALSE,nV), idx_obs)
          # idx_prd_tmp <- c(rep(FALSE,nV), idx_prd)   
          # idx_obs_tmp <- idx_obs
          # idx_prd_tmp <- idx_prd

          if(!no_nugget){
            # Sigma_prd <- solve(Q_tmp)
            # var_tmp <- diag(Sigma_prd[idx_prd_tmp,idx_prd_tmp])
            # print(var_tmp)
            # var_tmp <- diag(Sigma[idx_prd_tmp, idx_prd_tmp] - cov_loc %*%  solve(cov_Obs, t(cov_loc)))
            post_cov_tmp <- cov_loc %*%  solve(cov_Obs, t(cov_loc))
            var_tmp <- diag(Sigma[idx_prd, idx_prd] - post_cov_tmp)
            if(compute_pred_variances  || pred_samples) {
                  pred_cov <- Matrix::Diagonal(dim(cov_loc)[1], x= sigma_e^2) -cov_loc + post_cov_tmp
                  pred_var_tmp <- diag(pred_cov) 
            }

            
          } else{
            cov_tmp <- Sigma[idx_prd_tmp, idx_prd_tmp] - Sigma[idx_prd_tmp, idx_obs_tmp] %*% solve(Sigma[idx_obs_tmp, idx_obs_tmp],t(Sigma[idx_prd_tmp, idx_obs_tmp]))
            var_tmp <- diag(cov_tmp)
            var_tmp <- ifelse(var_tmp < 0, 0, var_tmp) # possible numerical errors
            if(compute_pred_variances  || pred_samples) {
                  pred_cov <- cov_tmp
                  pred_var_tmp <- diag(cov_tmp)
              }                  
          }
        }

        # var_tmp[graph_bkp$data[["__dummy_ord_var"]]] <- var_tmp

      if(return_original_order){
        var_tmp[ord_idx] <- var_tmp
        if(compute_pred_variances  || pred_samples) {
          pred_var_tmp[ord_idx] <- pred_var_tmp
        }        
      }
      if(!return_as_list){
        out$variance <- rep(var_tmp, length(u_repl))
        if(compute_pred_variances  || pred_samples) {
            out$pred_variance <- rep(pred_var_tmp, length(u_repl))
        }
      }
      else {
          for(repl_y in u_repl){
            out$variance[[repl_y]] <- var_tmp
            if(compute_pred_variances  || pred_samples) {
              out$pred_variance[[repl_y]] <- pred_var_tmp     
            }
          }
      }
    }

    if(posterior_samples){
      mean_tmp <- as.vector(mu_re[idx_prd, , drop=FALSE])
        if(cond_isocov){
          if(!no_nugget){
              post_cov <- Sigma[idx_prd_tmp, idx_prd_tmp] - cov_loc %*%  solve(cov_Obs, t(cov_loc))
          } else{
              post_cov <- Sigma[idx_prd_tmp, idx_prd_tmp] - Sigma[idx_prd_tmp, idx_obs_tmp] %*% solve(Sigma[idx_obs_tmp, idx_obs_tmp],t(Sigma[idx_prd_tmp, idx_obs_tmp]))
          }
        } else if(cond_alpha2){
              if(!no_nugget){
                cov_loc <- Sigma[idx_prd, idx_obs]
                cov_Obs <- Sigma[idx_obs, idx_obs]    
                diag(cov_Obs) <- diag(cov_Obs) + sigma_e^2    
                post_cov <- Sigma[idx_prd, idx_prd] - cov_loc %*%  solve(cov_Obs, t(cov_loc))
            } else{
              post_cov <- Sigma[idx_prd, idx_prd] - Sigma[idx_prd, idx_obs] %*% solve(Sigma[idx_obs, idx_obs],t(Sigma[idx_prd, idx_obs]))
            }
        } else{
          if(!no_nugget){
            post_cov <- A[idx_prd,]%*%solve(Q_xgiveny, t(A[idx_prd,]))
            } else{
          post_cov <- A[idx_prd,]  %*% M %*% t(A[idx_prd,])
          }
        }
      Z <- rnorm(dim(post_cov)[1] * n_samples)
      dim(Z) <- c(dim(post_cov)[1], n_samples)
      LQ <- chol(forceSymmetric(post_cov))
      X <- t(LQ) %*% Z
      X <- X + mean_tmp

      if(return_original_order){
        X[ord_idx,] <- X
      }

      if(!return_as_list){
        out$samples <- rbind(out$samples, X)
      } else{
        out$samples[[repl_y]] <- X
      }
    }
    
    if(pred_samples){
      mean_tmp <- as.vector(mu_krig[idx_prd, , drop=FALSE])
      Z <- rnorm(dim(pred_cov)[1] * n_samples)
      dim(Z) <- c(dim(pred_cov)[1], n_samples)
      LQ <- chol(forceSymmetric(pred_cov))
      X <- t(LQ) %*% Z
      X <- X + mean_tmp

      if(return_original_order){
        X[ord_idx,] <- X
      }

      if(!return_as_list){
        out$pred_samples <- rbind(out$pred_samples, X)
      } else{
        out$pred_samples[[repl_y]] <- X
      }
    }


  }


  if(check_euclidean){
    if(!is.null(graph_bkp$res_dist)){
      out$euclidean <- attr(graph_bkp$res_dist[[".complete"]], "euclidean")
    }
  }

  return(out)
}

#' @noRd 
# theta depends on the model
# object is a graph_lme object

get_covariance_precision <- function(object){
  graph_bkp <- object$graph$clone()
  graph_bkp$.__enclos_env__$private$data[["__dummy_ord_var"]] <- 1:length(graph_bkp$.__enclos_env__$private$data[[".edge_number"]])
  graph_bkp$observation_to_vertex(mesh_warning=FALSE)
  model_type <- object$latent_model
  coeff_fixed <- object$coeff$fixed_effects
  coeff_random <- object$coeff$random_effects
  coeff_meas <- object$coeff$measurement_error
  sigma_e <- coeff_meas[[1]]
  if(tolower(model_type$type) == "whittlematern"){
    tau <- object$coeff$random_effects[1]
    kappa <- object$coeff$random_effects[2] 
    if(model_type$alpha == 1){
        BC = object$BC
        Q <- spde_precision(kappa = kappa, tau = tau,
                          alpha = 1, graph = graph_bkp,
                          BC = BC)
        prec_cov <- Q
        attr(prec_cov, "prec_cov") <- "prec"
      } else{
        if(is.null(graph_bkp$CoB)){
          graph_bkp$buildC(2)
        }
        BC = object$BC
        PtE <- graph_bkp$get_PtE()
        n.c <- 1:length(graph_bkp$CoB$S)
        Q <- spde_precision(kappa = kappa, tau = tau, alpha = 2,
                            graph = graph_bkp, BC = BC)
        Qtilde <- (graph_bkp$CoB$T) %*% Q %*% t(graph_bkp$CoB$T)
        Qtilde <- Qtilde[-n.c,-n.c]
        Sigma.overdetermined  = t(graph_bkp$CoB$T[-n.c,]) %*% solve(Qtilde) %*%
          (graph_bkp$CoB$T[-n.c,])
        index.obs <- 4 * (PtE[,1] - 1) + 1.0 * (abs(PtE[, 2]) < 1e-14) +
          3.0 * (abs(PtE[, 2]) > 1e-14)
        Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
        prec_cov <- Sigma
        attr(prec_cov, "prec_cov") <- "cov" 
      }    
  } else if(tolower(model_type$type) == "graphlaplacian"){
    tau <- object$coeff$random_effects[1]
    graph_bkp$compute_laplacian(full=TRUE)
    kappa <- object$coeff$random_effects[2]
    if(model_type$alpha == 1){
        Q <- (kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]) * tau^2
      } else{
        Q <- kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]
        Q <- Q %*% Q * tau^2
      }
    prec_cov <- Q
    attr(prec_cov, "prec_cov") <- "prec"
  }  else if(tolower(model_type$type) == "isocov"){
      if(is.character(model_type$cov_function)){
        sigma <- object$coeff$random_effects[1]
        kappa <- object$coeff$random_effects[2]
        if(model_type$cov_function == "GL1"){
              tau <- object$coeff$random_effects[1]
              nV_temp <- object$nV_orig
              graph_bkp$compute_laplacian(full=TRUE)
              Q <- (kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]) * tau^2
        prec_cov <- Q
        attr(prec_cov, "prec_cov") <- "prec"              
        } else if(model_type$cov_function == "GL2"){
              tau <- object$coeff$random_effects[1]
              graph_bkp$compute_laplacian(full=TRUE)
              Q <- kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]
              Q <- Q %*% Q * tau^2
              prec_cov <- Q
              attr(prec_cov, "prec_cov") <- "prec" 
        }
      } else{
        graph_bkp$compute_resdist(full = TRUE, check_euclidean = FALSE)
        cov_function <- model_type$cov_function
        Sigma <- as.matrix(cov_function(graph_bkp$res_dist[[".complete"]], coeff_random))
        prec_cov <- Sigma
        attr(prec_cov, "prec_cov") <- "cov" 
      }
  }

  if(attr(prec_cov,"prec_cov")  == "prec"){
    A <- Matrix::Diagonal(dim(prec_cov)[1])[graph_bkp$PtV, ]
  } else{
    A <- NULL
  }
  return(list(order_vec = graph_bkp$.__enclos_env__$private$data[["__dummy_ord_var"]], A = A, prec_cov = prec_cov))
}


#' @title Simulation of models on metric graphs
#'
#' @description The function samples a Gaussian random field based on a
#' fitted model using `graph_lme()`.
#'
#' @param object A `graph_lme` object
#' @param nsim The number of simulations. 
#' @param seed an object specifying if and how the random number generator should be initialized (‘seeded’).
#' @param sample_latent If `FALSE`, samples for the response variable will be generated. If `TRUE`, samples for the latent model will be generated. The default is `FALSE`.
#' @param posterior Should posterior samples be generated? If `FALSE`, samples will be computed based on the estimated prior distribution. The default is `FALSE`.
#' @param which_repl Which replicates to generate the samples. If `NULL` samples will
#' be generated for all replicates. Default is `NULL`.
#' @param ... Currently not used.
#'
#' @return A list containing elements `samples`, `edge_number` and `distance_on_edge`. Each of them is a list, whose indexes are the replicates, and in `samples` a matrix is given with `nsim` columns, each one being a sample. `edge_number` and `distance_on_edges` contain the respective edge numbers and distances on edge for each sampled element. The locations of the samples are the location of the data in which the model was fitted.
#' @export
#' @method simulate graph_lme
#'
simulate.graph_lme <- function(object,
                              nsim = 1,
                              seed = NULL,                              
                              sample_latent = FALSE,
                              posterior = FALSE,
                              which_repl = NULL,
                              ...) {

  repl <- which_repl                                
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  if (!inherits(object, "graph_lme")) {
    stop("input object is not of class graph_lme")
  }

  coeff_fixed <- object$coeff$fixed_effects
  coeff_random <- object$coeff$random_effects
  coeff_meas <- object$coeff$measurement_error

  model_type <- object$latent_model

  sigma_e <- coeff_meas[[1]]

  # We dont need to clone it since we will not modify it
  graph_bkp <- object$graph  

  n <- sum(graph_bkp$.__enclos_env__$private$data[[".group"]] == graph_bkp$.__enclos_env__$private$data[[".group"]][1])

  ##
  repl_vec <- graph_bkp$.__enclos_env__$private$data[[".group"]]

  if(is.null(repl)){
    u_repl <- unique(graph_bkp$.__enclos_env__$private$data[[".group"]])
  } else{
    u_repl <- unique(repl)
  }

  X_cov_pred <- stats::model.matrix(object$covariates, graph_bkp$.__enclos_env__$private$data)

  if(all(dim(X_cov_pred) == c(0,1))){
    X_cov_pred <- matrix(1, nrow = length(graph_bkp$.__enclos_env__$private$data[[".group"]]), ncol=1)
  }
  if(ncol(X_cov_pred) > 0){
    mu <- X_cov_pred %*% coeff_fixed
  } else{
    mu <- matrix(0, nrow = length(graph_bkp$.__enclos_env__$private$data[[".group"]]), ncol=1)
  }

  edge_nb <- graph_bkp$.__enclos_env__$private$data[[".edge_number"]][1:n]
  dist_ed <- graph_bkp$.__enclos_env__$private$data[[".distance_on_edge"]][1:n]  

  Y <- graph_bkp$.__enclos_env__$private$data[[as.character(object$response_var)]] - mu

  out <- list()


  for(repl_y in u_repl){
    out$distance_on_edge[[repl_y]] <- dist_ed
    out$edge_number[[repl_y]] <- edge_nb    

    idx_repl <- graph_bkp$.__enclos_env__$private$data[[".group"]] == repl_y
    y_repl <- Y[idx_repl]
    y_repl <- y_repl
    
    mu_fe <- mu[idx_repl, , drop = FALSE]    

    prec_cov_prior <- get_covariance_precision(object)
    A <- prec_cov_prior$A
    order_idx <- prec_cov_prior$order_vec
    prec_cov_prior <- prec_cov_prior$prec_cov
    type_prec_cov <- attr(prec_cov_prior, "prec_cov")
    
    if(!posterior){
      mu_re <- 0
      if(type_prec_cov == "prec"){
        krig_cov <- A%*%solve(prec_cov_prior, t(A))   
      } else{
        krig_cov <- prec_cov_prior
      }
      pred_cov <- Matrix::Diagonal(dim(krig_cov)[1], x=sigma_e^2)
    } else{
      if(type_prec_cov == "prec"){
          krig_Q <- t(A) %*% A/sigma_e^2 + prec_cov_prior
          mu_re <- solve(krig_Q,as.vector(t(A) %*% y_repl / sigma_e^2))
          mu_re <- A %*% mu_re
          krig_cov <- A%*%solve(krig_Q, t(A))          
          pred_cov <- sigma_e^2 - krig_cov + krig_cov %*% solve(krig_cov + Matrix::Diagonal(dim(krig_cov)[1], x=sigma_e^2), t(krig_cov))          
      } else{
          cov_loc <- prec_cov_prior
          cov_Obs <- prec_cov_prior
          diag(cov_Obs) <- diag(cov_Obs) + sigma_e^2      
          mu_re <- cov_loc %*%  solve(cov_Obs, y_repl)
          post_cov_tmp <- cov_loc %*%  solve(cov_Obs, t(cov_loc))
          krig_cov <- cov_loc - post_cov_tmp
          pred_cov <- Matrix::Diagonal(dim(cov_loc)[1], x=sigma_e^2) -cov_loc + post_cov_tmp
      }
    }

    mu_krig <- mu_fe + mu_re

    mean_tmp <- as.vector(mu_krig)

    if(sample_latent){
      Z <- rnorm(dim(krig_cov)[1] * nsim)
      dim(Z) <- c(dim(krig_cov)[1], nsim)
      LQ <- chol(forceSymmetric(krig_cov))
      X <- t(LQ) %*% Z
      X[order_idx,] <- X
      out$samples[[repl_y]] <- X
    } else{
      if(posterior){
        Z <- rnorm(dim(pred_cov)[1] * nsim)
        dim(Z) <- c(dim(pred_cov)[1], nsim)
        LQ <- chol(forceSymmetric(pred_cov))
        X <- t(LQ) %*% Z
        X <- X + mean_tmp
        X[order_idx,] <- X      
        out$samples[[repl_y]] <- X      
      } else{
        Z <- rnorm(dim(krig_cov)[1] * nsim)
        dim(Z) <- c(dim(krig_cov)[1], nsim)
        LQ <- chol(forceSymmetric(krig_cov))
        X <- t(LQ) %*% Z
        X <- X + mean_tmp
        Z <- rnorm(dim(pred_cov)[1] * nsim) * sigma_e
        dim(Z) <- c(dim(pred_cov)[1], nsim)
        X <- X + Z
        X[order_idx,] <- X
        out$samples[[repl_y]] <- X
      }
    }
  }
  return(out)
}
