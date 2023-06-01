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
#' There is also the option to provide it as a list containing the elements
#' `type`, which can be `linearModel`, `WhittleMatern`, `graphLaplacian` or `isoCov`.
#' `linearModel` corresponds to a linear model without random effects.
#' For `WhittleMatern` models, that is, if the list contains `type = 'WhittleMatern'`,
#' one can choose between a finite element approximation of the precision matrix
#' by adding `fem = TRUE` to the list, or to use the exact precision matrix
#' (by setting `fem = FALSE`). If `fem` is `FALSE`, there is also the parameter
#' `alpha`, to determine the order of the SPDE, which is either 1 or 2. If `fem`
#' is `TRUE` and `alpha` is not specified, then the default value of `alpha=1`
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
#' @param optim_method The method to be used with `optim` function.
#' @param starting_values_latent A vector containing the starting values for the
#' latent model. If the latent model is `WhittleMatern` or `graphLaplacian`, then
#' the starting values should be provided as a vector of the form c(sigma,kappa)
#' or c(sigma,range) depending on the parameterization. If the model is `isoCov`,
#' then the starting values should be provided as a vector containing the parameters
#' of the covariance function.
#' @param start_sigma_e Starting value for the standard deviation of the measurament
#' error.
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
#' @return A list containing the fitted model.
#' @rdname graph_lme
#' @export
#'
graph_lme <- function(formula, graph,
                model = list(type = "linearModel"),
                which_repl = NULL,
                optim_method = "L-BFGS-B",
                starting_values_latent = NULL,
                start_sigma_e = NULL,
                # parameterization_latent = c("matern", "spde"),
                BC = 1,
                # model_matrix = TRUE,
                parallel = FALSE,
                n_cores = parallel::detectCores()-1,
                optim_controls = list(),
                improve_hessian = FALSE,
                hessian_args = list()) {

  if(!is.list(model)){
    if(!is.character(model)){
      stop("The 'model' argument must be either a list or a character (string).")
    }
    model <- model[[1]]
    model <- tolower(model)
    if(!(model%in% c("lm", "wm", "wm1", "wm2", "isoexp", "gl1", "gl2"))){
      stop("If model is a character (string), the options are 'lm', 'WM', 'WM1', 'WM2', 'isoExp', 'GL1' or 'GL2'.")
    }
    model <- switch(model,
            "lm" = list(type = "linearModel"),
            "wm1" = list(type = "WhittleMatern", fem = FALSE, alpha = 1, version = 1),
            "wm2" = list(type = "WhittleMatern", fem = FALSE, alpha = 2),
            "isoexp" = list(type = "isoCov"),
            "gl1" = list(type = "graphLaplacian", alpha = 1),
            "gl2" = list(type = "graphLaplacian", alpha = 2),
            "wm" = list(type = "WhittleMatern", fem = TRUE)
            )
  }

  model_type <- model[["type"]]
  model_type <- tolower(model_type)

  model_type <- switch(model_type,
            "lm" = "linearModel",
            "wm1" = "WhittleMatern",
            "wm2" = "WhittleMatern",
            "isoexp" = "isoCov",
            "gl1" = "graphLaplacian",
            "gl2" = "graphLaplacian",
            "wm" = "WhittleMatern",
            model_type
            )

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

  if(is.null(graph$data)){
    stop("No data found in the graph. Please, add observations before fitting the model.")
  }

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
      } else{ rspde_object <- rSPDE::matern.operators(graph = graph,
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

      fit <- rSPDE::rspde_lme(formula = formula, model = rspde_object,
                            nu = nu, which_repl = which_repl,
                            optim_method = optim_method,
                            use_data_from_graph = TRUE,
                            parallel = parallel,
                            n_cores = n_cores,
                            starting_values_latent = starting_values_latent,
                            start_sigma_e = start_sigma_e,
                            optim_controls = optim_controls,
                            improve_hessian = improve_hessian,
                            hessian_args = hessian_args)
      fit$call <- call_graph_lme
      # if(fit$estimate_nu){
      #   names(fit$coeff$random_effects)[1] <- "alpha"
      #   fit$coeff$random_effects[1] <- fit$coeff$random_effects[1] + 0.5
      # }
      class(fit) <- c(class(fit), "graph_lme")
      return(fit)
    } else{
      if(!(model[["alpha"]] %in% c(1,2))){
        stop("For WhittleMatern models, alpha must be either 1 or 2. For different values, set 'fem' to 'TRUE' instead.")
      }
    }
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

  time_build_likelihood_start <- Sys.time()

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
      range_par <- FALSE
      # if(model_type == "whittlematern"){
      #   range_par <- ifelse(parameterization_latent == "matern",TRUE,FALSE)
      # }
      start_values <- graph_starting_values(graph = graph_bkp,
                    model = model_start,
                    manual_data = unlist(y_graph),
                    log_scale = TRUE,
                    # range_par = range_par)
                    range_par = FALSE)
    }
  } else if(model_type != "linearmodel") {
    start_values <- c(log(0.1*sd(y_graph),log(starting_values_latent)))
    par_names <- names(starting_values_latent)

    if(!is.null(start_sigma_e)){
        start_values[1] <- log(start_sigma_e)
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
    if(model[["alpha"]] == 1){
      if(model[["version"]] == 2){
        likelihood <- function(theta){
          return(-likelihood_alpha1_v2(theta = theta, graph = graph_bkp,
              X_cov = X_cov, y = y_graph, repl = which_repl, BC = BC,
              parameterization = "spde")) # parameterization = parameterization_latent))
        }
      } else {
        likelihood <- function(theta){
          return(-likelihood_alpha1(theta = theta, graph = graph_bkp,
                                    data_name = NULL, manual_y = y_graph,
                             X_cov = X_cov, repl = which_repl, BC = BC,
                             parameterization = "spde")) # , parameterization = parameterization_latent))
        }
      }
    } else{
      likelihood <- function(theta){
          return(-likelihood_alpha2(theta = theta, graph = graph_bkp,
                                    data_name = NULL, manual_y = y_graph,
                             X_cov = X_cov, repl = which_repl, BC = BC,
                             parameterization = "spde")) # , parameterization = parameterization_latent))
        }
    }
  } else if (model_type == "graphlaplacian"){
      likelihood <- likelihood_graph_laplacian(graph = graph_bkp,
                                               alpha = model[["alpha"]],
                                               y_graph = y_graph,
                                               X_cov = X_cov, maximize = FALSE,
                                               repl=which_repl,
                                               parameterization = "spde")
  } else if(model_type == "isocov") {

  if (is.character(model[["cov_function"]])) {
    if(model[["cov_function"]] %in% c("WM1","WM2", "GL1", "GL2")){
      model_cov <- model[["cov_function"]]
      par_names <- c("tau", "kappa")
    } else{
      model_cov <- "isoCov"
      if(model[["cov_function"]] == "exp_covariance"){
        model[["cov_function"]] <- exp_covariance
        model[["cov_function_name"]] <- "exp_covariance"
      }
    }
    likelihood <- likelihood_graph_covariance(graph_bkp, model = model_cov,
                                              y_graph = y_graph,
                                              cov_function = model[["cov_function"]],
                                              X_cov = X_cov, repl = which_repl)
    } else{
    model[["cov_function_name"]] <- "other"
      likelihood <- likelihood_graph_covariance(graph_bkp, model = model_cov,
                                                y_graph = y_graph,
                                                cov_function = model[["cov_function"]],
                                                X_cov = X_cov, repl = which_repl)
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

      if(parallel){
        start_par <- Sys.time()
        cl <- parallel::makeCluster(n_cores)
        parallel::setDefaultCluster(cl = cl)
        parallel::clusterExport(cl, "y_graph", envir = environment())
        parallel::clusterExport(cl, "graph_bkp", envir = environment())
        parallel::clusterExport(cl, "X_cov", envir = environment())
        parallel::clusterExport(cl, "which_repl", envir = environment())
        parallel::clusterExport(cl, "model", envir = environment())
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
          res <- optimParallel::optimParallel(start_values,
                        likelihood, method = optim_method,
                        control = optim_controls,
                        hessian = hessian,
                        parallel = list(forward = FALSE, cl = cl,
                            loginfo = FALSE))
        end_fit <- Sys.time()
        time_fit <- end_fit-start_fit
        parallel::stopCluster(cl)
      } else{
        start_fit <- Sys.time()
            res <- optim(start_values,
                        likelihood, method = optim_method,
                        control = optim_controls,
                        hessian = hessian)
        end_fit <- Sys.time()
        time_fit <- end_fit-start_fit
      }


  time_hessian <- NULL

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


  # res <- optim(start_values,
  #               likelihood, method = optim_method,
  #               control = optim_controls,
  #               hessian = TRUE)

  if(model_type %in% c("graphlaplacian", "whittlematern")){
    n_par_coeff <- 3
  } else if(model[["cov_function_name"]] == "exp_covariance"){
    n_par_coeff <- 3
  } else{
    n_par_coeff <- length(starting_values_latent) + 1
  }

  coeff <- res$par
  coeff <- exp(c(res$par[1:n_par_coeff]))
  coeff <- c(coeff, res$par[-c(1:n_par_coeff)])

  loglik <- -res$value

  n_fixed <- ncol(X_cov)
  n_random <- length(coeff) - n_fixed - 1

  n_coeff_nonfixed <- length(coeff) - n_fixed

  par_change <- diag(c(exp(-c(res$par[1:n_coeff_nonfixed])), rep(1,n_fixed)))
  observed_fisher <- par_change %*% observed_fisher %*% par_change

  if(model_type %in% c("graphlaplacian", "whittlematern")){

      coeff[2] <- 1/coeff[2]

      grad_tmp <- diag(c(1,-1/(coeff[2]^2), 1, rep(1,n_fixed)))

      observed_fisher <- grad_tmp %*% observed_fisher %*% grad_tmp

      time_matern_par_start <- Sys.time()
      new_likelihood <- function(theta){
        new_par <- res$par
        new_par[2:3] <- theta
        new_par[2] <- -new_par[2]
        return(likelihood(new_par))
      }

      coeff_tmp <- coeff[2:3]
      new_observed_fisher <- observed_fisher[2:3,2:3]
      change_par <- change_parameterization_graphlme(new_likelihood,
                                                     model[["alpha"]]-0.5,
                                              coeff_tmp,
                                              hessian = new_observed_fisher
                                              )
      matern_coeff <- list()
      matern_coeff$random_effects <- change_par$coeff
      names(matern_coeff$random_effects) <- c("sigma", "range")
      matern_coeff$std_random <- change_par$std_random
      time_matern_par_end <- Sys.time()
      time_matern_par <- time_matern_par_end - time_matern_par_start
    } else{
      matern_coeff <- NULL
      time_matern_par <- NULL
    }


  inv_fisher <- tryCatch(solve(observed_fisher),
                         error = function(e) matrix(NA,
                                                    nrow(observed_fisher),
                                                    ncol(observed_fisher)))
  std_err <- sqrt(diag(inv_fisher))

  coeff_random <- coeff[2:(1+n_random)]
  std_random <- std_err[2:(1+n_random)]
  names(coeff_random) <- par_names

  coeff_meas <- coeff[1]
  names(coeff_meas) <- "std. dev"

  std_meas <- std_err[1]

  coeff_fixed <- NULL
  if(n_fixed > 0){
    coeff_fixed <- coeff[(2+n_random):length(coeff)]
    std_fixed <- std_err[(2+n_random):length(coeff)]
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
  object$response <- list(y = y_graph)
  object$formula <- formula
  object$estimation_method <- optim_method
  # object$parameterization_latent <- parameterization_latent
  object$which_repl <- which_repl
  object$optim_controls <- optim_controls
  object$latent_model <- model
  object$loglik <- loglik
  object$BC <- BC
  object$niter <- res$counts
  object$response <- y_term
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
  # if(model_matrix){
    if(ncol(X_cov)>0){
      object$model_matrix <- cbind(y_graph, X_cov)
    } else{
      object$model_matrix <- y_graph
    }
  # }
  object$graph <- graph$clone()


  class(object) <- "graph_lme"
  return(object)

}

#' @name logLik.graph_lme
#' @title Log-likelihood for \code{graph_lme} objects
#' @description computes the log-likelihood for a fitted mixed effects model on
#' metric graphs.
#' @param x Object of class `graph_lme` containing results from the fitted model.
#' @param ... further arguments passed to or from other methods.
#' @return Log-likelihood value.
#' @noRd
#' @method logLik graph_lme
#' @export
logLik.graph_lme <- function(object, ...){
  return(object$loglik)
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
    cat(paste0("\n", "Random effects (Matern parameterization):", "\n"))
    print(x$matern_coeff$random_effects)
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
#' @param data A `data.frame` or a `list` containing the covariates, the edge
#' number and the distance on edge for the locations to obtain the prediction.
#' @param mesh Obtain predictions for mesh nodes? The graph must have a mesh,
#' and either `only_latent` is set to TRUE or the model does not have covariates.
#' @param mesh_h If the graph does not have a mesh, one will be created with this
#' value of 'h'.
#' @param repl Which replicates to obtain the prediction. If `NULL` predictions
#' will be obtained for all replicates. Default is `NULL`.
#' @param compute_variances Set to also TRUE to compute the kriging variances.
#' @param posterior_samples If `TRUE`, posterior samples will be returned.
#' @param n_samples Number of samples to be returned. Will only be used if
#' `sampling` is `TRUE`.
#' @param only_latent Should the posterior samples and predictions be only given
#' to the latent model?
#' @param edge_number Name of the variable that contains the edge number, the
#' default is `edge_number`.
#' @param distance_on_edge Name of the variable that contains the distance on
#' edge, the default is `distance_on_edge`.
#' @param normalized Are the distances on edges normalized?
#' @param return_as_list Should the means of the predictions and the posterior
#' samples be returned as a list, with each replicate being an element?
#' @param return_original_order Should the results be return in the original
#' (input) order or in the order inside the graph?
#' @param ... Not used.
#' @return A list with elements `mean`, which contains the means of the
#' predictions, `variance` (if `compute_variance` is `TRUE`), which contains the
#' variances of the predictions, `samples` (if `posterior_samples` is `TRUE`),
#' which contains the posterior samples.
#' @export
#' @method predict graph_lme

predict.graph_lme <- function(object,
                              data = NULL,
                              mesh = FALSE,
                              mesh_h = 0.01,
                              repl = NULL,
                              compute_variances = FALSE,
                              posterior_samples = FALSE,
                              n_samples = 100,
                              only_latent = FALSE,
                              edge_number = "edge_number",
                              distance_on_edge = "distance_on_edge",
                              normalized = FALSE,
                              return_as_list = FALSE,
                              return_original_order = TRUE,
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

  BC <- object$BC

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


  if(!is.null(data[[as.character(object$response)]])){
    data[[as.character(object$response)]] <- NULL
  }

  data_graph_temp <- list()
  idx_group1 <-  graph_bkp$data[["__group"]] == graph_bkp$data[["__group"]][1]
  data_graph_temp[[edge_number]] <- graph_bkp$data[["__edge_number"]][idx_group1]
  data_graph_temp[[distance_on_edge]] <- graph_bkp$data[["__distance_on_edge"]][idx_group1]
  data_graph_temp[[as.character(object$response)]] <- graph_bkp$data[[as.character(object$response)]][idx_group1]
  data_graph_temp <- as.data.frame(data_graph_temp)

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

  old_data <- graph_bkp$data

  data[["__group"]] <- old_data[["__group"]][1]

  graph_bkp$clear_observations()

  graph_bkp$add_observations(data = data, edge_number = edge_number,
                             distance_on_edge = distance_on_edge,
                             normalized = TRUE, group = "__group")

  graph_bkp$add_observations(data = old_data, edge_number = "__edge_number",
                             distance_on_edge = "__distance_on_edge",
                             group = "__group", normalized = TRUE)

  graph_bkp$data[["__dummy_ord_var"]] <- 1:length(graph_bkp$data[["__edge_number"]])

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
    X_cov_pred <- matrix(1, nrow = length(graph_bkp$data[["__group"]]), ncol=1)
  }
  if(ncol(X_cov_pred) > 0){
    mu <- X_cov_pred %*% coeff_fixed
  } else{
    mu <- matrix(0, nrow = length(graph_bkp$data[["__group"]]), ncol=1)
  }

  Y <- graph_bkp$data[[as.character(object$response)]] - mu

  model_type <- object$latent_model

  sigma.e <- coeff_meas[[1]]
  sigma_e <- sigma.e

  if(!is.null(graph_bkp$data[["__dummy_var"]])){
      idx_prd <- !is.na(graph_bkp$data[["__dummy_var"]][1:n])
  } else {
      idx_prd <- !is.na(graph_bkp$data[["X__dummy_var"]][1:n])
  }

  n_prd <- sum(idx_prd)

  edge_nb <- graph_bkp$data[["__edge_number"]][1:n][idx_prd]
  dist_ed <- graph_bkp$data[["__distance_on_edge"]][1:n][idx_prd]

  ## construct Q

  # graph_bkp$data <- lapply(graph_bkp$data, function(dat){dat_temp <- dat
  #                       dat_temp[graph_bkp$data[["__dummy_ord_var"]]] <- dat
  #                                       return(dat_temp)})


  if(tolower(model_type$type) == "whittlematern"){
    tau <- object$coeff$random_effects[1]
    # if(object$parameterization_latent == "spde"){
      kappa <- object$coeff$random_effects[2]
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
    graph_bkp$observation_to_vertex()
    tau <- object$coeff$random_effects[1]
    #nV before
    nV_temp <- object$nV_orig
    # graph_bkp$observation_to_vertex()
    if(graph_bkp$nV > nV_temp){
      warning("There are prediction locations outside of the observation locations. Refit the model with all the locations you want to obtain predictions.")
    }
    graph_bkp$compute_laplacian()
    # if(object$parameterization_latent == "spde"){
      kappa <- object$coeff$random_effects[2]
    # } else{
    #   kappa <- sqrt(8 * 0.5) / object$coeff$random_effects[2]
    # }
      if(model_type$alpha == 1){
        Q <- (kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]) * tau^2
      } else{
        Q <- kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]
        Q <- Q %*% Q * tau^2
      }

  } else if(tolower(model_type$type) == "isocov"){
      if(is.character(model_type$cov_function)){
        sigma <- object$coeff$random_effects[1]
        kappa <- object$coeff$random_effects[2]
        # if(model_type$cov_function == "alpha1"){
        #   # graph_bkp$observation_to_vertex()
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
              graph_bkp$observation_to_vertex()
              if(graph_bkp$nV > nV_temp){
                warning("There are prediction locations outside of the observation locations. Refit the model with all the locations you want to obtain predictions.")
              }
              graph_bkp$compute_laplacian()
              Q <- (kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]) * tau^2
        } else if(model_type$cov_function == "GL2"){
              #nV before
              tau <- object$coeff$random_effects[1]
              nV_temp <- object$nV_orig
              graph_bkp$observation_to_vertex()
              if(graph_bkp$nV > nV_temp){
                warning("There are prediction locations outside of the observation locations. Refit the model with all the locations you want to obtain predictions.")
              }
              graph_bkp$compute_laplacian()
              Q <- kappa^2 * Matrix::Diagonal(graph_bkp$nV, 1) + graph_bkp$Laplacian[[1]]
              Q <- Q %*% Q * tau^2
        # } else if(model_type$cov_function == "exp_covariance"){
        #           graph_bkp$compute_resdist(full = TRUE)
        #           Sigma <- as.matrix(exp_covariance(graph_bkp$res_dist[[1]], c(sigma,kappa)))
        #           Q <- solve(Sigma)
        # }
        }
      } else{
        graph_bkp$compute_resdist(full = TRUE)
        cov_function <- model_type$cov_function
        Sigma <- as.matrix(cov_function(graph_bkp$res_dist[[1]], coeff_random))
      }
  }

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

  if(compute_variances || posterior_samples){
    if(cond_wm){
      tau <- object$coeff$random_effects[1]
      graph_bkp2 <- graph_bkp$clone()
      graph_bkp2$observation_to_vertex()
      if(cond_alpha1){
        Q <- spde_precision(kappa = kappa, tau = tau,
                          alpha = 1, graph = graph_bkp2)
      } else{
        PtE <- graph_bkp2$get_PtE()
        n.c <- 1:length(graph_bkp2$CoB$S)
        Q <- spde_precision(kappa = kappa, tau = tau, alpha = 2,
                            graph = graph_bkp2, BC = BC)
        Qtilde <- (graph_bkp2$CoB$T) %*% Q %*% t(graph_bkp2$CoB$T)
        Qtilde <- Qtilde[-n.c,-n.c]
        Sigma.overdetermined  = t(graph_bkp2$CoB$T[-n.c,]) %*% solve(Qtilde) %*%
          (graph_bkp2$CoB$T[-n.c,])
        index.obs <- 4 * (PtE[,1] - 1) + 1.0 * (abs(PtE[, 2]) < 1e-14) +
          3.0 * (abs(PtE[, 2]) > 1e-14)
        Sigma <-  as.matrix(Sigma.overdetermined[index.obs, index.obs])
        Q <- solve(Sigma)
      }
      A <- Matrix::Diagonal(dim(Q)[1])[graph_bkp2$PtV, ]
      rm(graph_bkp2)
    }
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

    if(!cond_wm && !cond_isocov){
        Q_xgiveny <- t(A[idx_obs,]) %*% A[idx_obs,]/sigma_e^2 + Q

        mu_krig <- solve(Q_xgiveny,as.vector(t(A[idx_obs,]) %*% y_repl / sigma_e^2))

        # mu_krig <- mu_krig[(gap+1):length(mu_krig)]
        mu_krig <- A[idx_prd,] %*% mu_krig

        if(!only_latent){
          mu_fe <- mu[idx_repl, , drop = FALSE]
          mu_krig <- mu_fe[idx_prd, , drop=FALSE] + mu_krig
        }
    } else if (cond_wm){

      PtE_obs <- PtE_full[idx_obs,]

      if(cond_alpha2){
          mu_krig <- posterior_mean_obs_alpha2(c(sigma.e,tau,kappa),
                        graph = graph_bkp, PtE_resp = PtE_obs, resp = y_repl,
                        PtE_pred = cbind(data_prd_temp[[edge_number]],
                                         data_prd_temp[[distance_on_edge]]))
        if(!only_latent){
          mu_fe <- mu[idx_repl, , drop = FALSE]
          mu_krig <- mu_fe[idx_prd, , drop=FALSE] + mu_krig[ord_idx]
        } else{
            mu_krig <- mu_krig[ord_idx]
          }
      } else{
          mu_krig <- posterior_mean_obs_alpha1(c(sigma.e,tau,kappa),
                        graph = graph_bkp, PtE_resp = PtE_obs, resp = y_repl,
                        PtE_pred = cbind(data_prd_temp[[edge_number]],
                                         data_prd_temp[[distance_on_edge]]))
                        # PtE_pred = cbind(edge_nb, dist_ed))
          if(!only_latent){
            mu_fe <- mu[idx_repl, , drop = FALSE]
            mu_krig <- mu_fe[idx_prd, , drop=FALSE] + mu_krig[ord_idx]
          } else{
            mu_krig <- mu_krig[ord_idx]
          }
      }
    } else {
        Sigma <- as.matrix(cov_function(graph_bkp$res_dist[[1]], coeff_random))

        cov_loc <- Sigma[idx_prd, idx_obs]
        cov_Obs <- Sigma[idx_obs, idx_obs]

        # Observe that the "fixed-effects" mean has been subtracted from y_repl

        mu_krig <- cov_loc %*%  solve(cov_Obs, y_repl)

        if(!only_latent){
            mu_fe <- mu[idx_repl, , drop = FALSE]
            mu_krig <- mu_fe[idx_prd, , drop=FALSE] + mu_krig
          } else{
            mu_krig <- mu_krig
          }

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

    if(compute_variances || posterior_samples){
      if(cond_wm){
            Q_xgiveny <- t(A[idx_obs,]) %*% A[idx_obs,]/sigma_e^2 + Q
      }
    }

    if (compute_variances) {
      if(!cond_isocov){
        post_cov <- A[idx_prd,]%*%solve(Q_xgiveny, t(A[idx_prd,]))
        var_tmp <- max(diag(post_cov),0)
      } else{
        var_tmp <- diag(Sigma[idx_prd, idx_prd] - Sigma[idx_prd, idx_obs] %*% solve(Sigma[idx_obs, idx_obs],t(Sigma[idx_prd, idx_obs])))
        var_tmp <- ifelse(var_tmp < 0, 0, var_tmp) # possible numerical errors
      }

        # var_tmp[graph_bkp$data[["__dummy_ord_var"]]] <- var_tmp

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
      if(cond_isocov){
        post_cov <- Sigma[idx_prd, idx_prd] - Sigma[idx_prd, idx_obs] %*% solve(Sigma[idx_obs, idx_obs], t(Sigma[idx_prd,  idx_obs]))
      } else{
        post_cov <- A[idx_prd,]%*%solve(Q_xgiveny, t(A[idx_prd,]))
      }
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


