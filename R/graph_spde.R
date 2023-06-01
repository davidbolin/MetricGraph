
#' 'INLA' implementation of Whittle-Matérn fields for metric graphs
#'
#' This function creates an 'INLA' object that can be used
#' in 'INLA' or 'inlabru' to fit Whittle-Matérn fields on metric graphs.
#'
#' @param graph_object A `metric_graph` object.
#' @param alpha The order of the SPDE.
#' @param stationary_endpoints Which vertices of degree 1 should contain
#' stationary boundary conditions?
#' @param parameterization Which parameterization to be used? The options are
#' 'matern' (sigma and range) and 'spde' (sigma and kappa).
#' @param start_sigma Starting value for sigma.
#' @param prior_sigma a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation of sigma on the log scale.
#' @param start_range Starting value for range parameter.
#' @param prior_range a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation of the range parameter on
#' the log scale. Will not be used if prior.kappa is non-null.
#' @param start_kappa Starting value for kappa.
#' @param prior_kappa a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation of kappa on the log scale.
#' @param shared_lib Which shared lib to use for the cgeneric implementation?
#' If "detect", it will check if the shared lib exists locally, in which case it will
#' use it. Otherwise it will use 'INLA's shared library.
#' If 'INLA', it will use the shared lib from 'INLA's installation. If 'MetricGraph', then
#' it will use the local installation (does not work if your installation is from CRAN).
#' Otherwise, you can directly supply the path of the .so (or .dll) file.
#' @param debug Should debug be displayed?
#'
#' @return An 'INLA' object.
#' @details
#' This function is used to construct a Matern SPDE model on a metric graph.
#' The latent field \eqn{u} is the solution of the SPDE
#' \deqn{(\kappa^2 - \Delta)^\alpha u = \sigma W,} where \eqn{W} is Gaussian
#' white noise on the metric graph. This model implements exactly
#' the cases in which \eqn{\alpha = 1} or \eqn{\alpha = 2}. For a finite
#' element approximation for general \eqn{\alpha} we refer the reader to the
#' 'rSPDE' package and to the Whittle--Matérn fields with general smoothness vignette.
#'
#' We also have the alternative parameterization \eqn{\rho = \frac{\sqrt{8(\alpha-0.5)}}{\kappa}},
#' which can be interpreted as a range parameter.
#'
#' Let \eqn{\kappa_0} and \eqn{\sigma_0} be the starting values for \eqn{\kappa} and
#' \eqn{\sigma}, we write \eqn{\sigma = \exp\{\theta_1\}} and \eqn{\kappa = \exp\{\theta_2\}}.
#' We assume priors on \eqn{\theta_1} and \eqn{\theta_2} to be normally distributed
#' with mean, respectively, \eqn{\log(\sigma_0)} and \eqn{\log(\kappa_0)}, and variance 10.
#' Similarly, if we let \eqn{\rho_0} be the starting value for \eqn{\rho}, then
#' we write \eqn{\rho = \exp\{\theta_2\}} and assume a normal prior for \eqn{\theta_2},
#' with mean \eqn{\log(\rho_0)} and variance 10.
#'
#' @export
graph_spde <- function(graph_object,
                       alpha = 1,
                       stationary_endpoints = "all",
                       parameterization = c("matern", "spde"),
                       start_range = NULL,
                       prior_range = NULL,
                       start_kappa = NULL,
                       start_sigma = NULL,
                       prior_kappa = NULL,
                       prior_sigma = NULL,
                       shared_lib = "detect",
                       debug = FALSE){

  graph_spde <- graph_object$clone()

  graph_spde$observation_to_vertex()

  parameterization <- parameterization[[1]]
  if(!(alpha%in%c(1,2))){
    stop("alpha must be either 1 or 2!")
  }
  if(alpha == 2){
    stop("Only alpha=1 implemented.")
  }
  nu <- alpha - 0.5
  V <- graph_spde$V
  EtV <- graph_spde$E
  El <- graph_spde$edge_lengths

  i_ <- j_ <- rep(0, dim(V)[1]*4)
  nE <- dim(EtV)[1]
  count <- 0
  for(i in 1:nE){
    l_e <- El[i]

    if(EtV[i,1]!=EtV[i,2]){

      i_[count + 1] <- EtV[i,1] - 1
      j_[count + 1] <- EtV[i,1] - 1

      i_[count + 2] <- EtV[i,2] - 1
      j_[count + 2] <- EtV[i,2] - 1


      i_[count + 3] <- EtV[i,1] - 1
      j_[count + 3] <- EtV[i,2] - 1

      i_[count + 4] <- EtV[i,2] - 1
      j_[count + 4] <- EtV[i,1] - 1
      count <- count + 4
    }else{
      i_[count + 1] <- EtV[i,1] - 1
      j_[count + 1] <- EtV[i,1] - 1
      count <- count + 1
    }
  }
  n.v <- dim(V)[1]

  if(stationary_endpoints == "all"){
    i.table <- table(i_[1:count])
    index <- as.integer(names(which(i.table<3)))
    index <- index
  } else if(stationary_endpoints == "none"){
    index <- NULL
  } else{
    index <- stationary_endpoints - 1
  }
    if(!is.null(index)){
    #does this work for circle?
        i_ <- c(i_[1:count], index)
        j_ <- c(j_[1:count], index)
        count <- count + length(index)
    }

    if(is.null(index)){
        index <- -1
    }

    EtV2 <- EtV[,1]
    EtV3 <- EtV[,2]
    El <- as.vector(El)


  idx_ij <- order(i_, j_)
  j_ <- j_[idx_ij]
  i_ <- i_[idx_ij]

  idx_sub <- which(i_ <= j_)
  j_ <- j_[idx_sub]
  i_ <- i_[idx_sub]

  graph_matrix_1 <- cbind(i_, j_)
  graph_matrix <- unique(graph_matrix_1)

  count_idx <- numeric(nrow(graph_matrix))

  row_tmp <- 1
  for(i in 1:nrow(graph_matrix)){
    count_tmp <- 0
    j <- row_tmp
    while(j <= nrow(graph_matrix_1) && all(graph_matrix[i,]==graph_matrix_1[j,])){
      count_tmp <- count_tmp + 1
      j <- j + 1
    }
    row_tmp <- j
    count_idx[i] <- count_tmp
  }

  i_ <- graph_matrix[,1]
  j_ <- graph_matrix[,2]

  idx_ij <- idx_ij[idx_sub]
  idx_ij <- sort(idx_ij, index.return=TRUE)
  idx_ij <- idx_ij$ix

  idx_ij <- idx_ij - 1


    if(is.null(prior_kappa$meanlog) && is.null(prior_range$meanlog)){
      model_start <- ifelse(alpha==1,"alpha1", "alpha2")
      start_values_vector <- graph_starting_values(graph_spde,
                      model = model_start, data=FALSE)

      prior_kappa$meanlog <- log(start_values_vector[3])
      prior_range$meanlog <- log(sqrt(8 * nu)) - prior_kappa$meanlog
    } else if(is.null(prior_range$meanlog)){
      prior_range$meanlog <- log(sqrt(8 * nu)) -
      prior_kappa$meanlog
    } else{
      prior_kappa$meanlog <- log(sqrt(8 * nu)) -
      prior_range$meanlog
    }

  if (is.null(prior_kappa$sdlog)) {
    prior_kappa$sdlog <- sqrt(10)
  }

  if(is.null(prior_range$sdlog)){
    prior_range$sdlog <- sqrt(10)
  }

  if(is.null(prior_sigma$meanlog)){
    prior_sigma$meanlog <- 0
  }

  if(is.null(prior_sigma$sdlog)){
    prior_sigma$sdlog <- sqrt(10)
  }

  if (is.null(start_kappa)) {
    start_lkappa <- prior_kappa$meanlog
  } else{
    start_lkappa <- log(start_kappa)
  }
  if (is.null(start_sigma)) {
    start_lsigma <- prior_sigma$meanlog
  } else{
    start_lsigma <- log(start_sigma)
  }
  if(is.null(start_range)){
    start_lrange <- prior_range$meanlog
  } else{
    start_lrange <- log(start_range)
  }

  if(parameterization == "matern"){
    start_theta <- start_lrange
    prior_theta <- prior_range
  } else{
    start_theta <- start_lkappa
    prior_theta <- prior_kappa
  }

  ### Location of object files

  gpgraph_lib <- shared_lib

  if(shared_lib == "INLA"){
    gpgraph_lib <- INLA::inla.external.lib('rSPDE')
  } else if(shared_lib == "MetricGraph"){
    gpgraph_lib <- system.file('shared', package='MetricGraph')
    if(Sys.info()['sysname']=='Windows') {
		gpgraph_lib <- paste0(gpgraph_lib, "/gpgraph_cgeneric_models.dll")
            } else {
		gpgraph_lib <- paste0(gpgraph_lib, "/gpgraph_cgeneric_models.so")
            }
  } else if(shared_lib == "detect"){
    gpgraph_lib_local <- system.file('shared', package='MetricGraph')
    if(Sys.info()['sysname']=='Windows') {
		gpgraph_lib_local <- paste0(gpgraph_lib_local, "/gpgraph_cgeneric_models.dll")
            } else {
		gpgraph_lib_local <- paste0(gpgraph_lib_local, "/gpgraph_cgeneric_models.so")
            }
    if(file.exists(gpgraph_lib_local)){
      gpgraph_lib <- gpgraph_lib_local
    } else{
      gpgraph_lib <- INLA::inla.external.lib('rSPDE')
    }
  }


  model <-
        do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_gpgraph_alpha1_model",
            shlib=gpgraph_lib,
            n=as.integer(n.v), debug=debug,
            prec_graph_i = as.integer(i_),
            prec_graph_j = as.integer(j_),
            index_graph = as.integer(idx_ij),
            count_idx = as.integer(count_idx),
            EtV2 = EtV2,
            EtV3 = EtV3,
            El = El,
            stationary_endpoints = as.integer(index),
            start_theta = start_theta,
            start_lsigma = start_lsigma,
            prior_theta_meanlog = prior_theta$meanlog,
            prior_theta_sdlog = prior_theta$sdlog,
            prior_sigma_meanlog = prior_sigma$meanlog,
            prior_sigma_sdlog = prior_sigma$sdlog,
            parameterization = parameterization))
model$graph_spde <- graph_spde
model$data_PtE <- graph_object$get_PtE()
model$parameterization <- parameterization
class(model) <- c("inla_metric_graph_spde", class(model))
return(model)
}


#'  Model index vector generation for metric graph models
#'
#' Generates a list of named index vectors for
#' 'INLA'-based metric graph models.
#'
#' @param name A character string with the base name of the effect.
#' @param graph_spde An `inla_metric_graph_spde` object built with the
#' `graph_spde()` function.
#' @param n.group Number of groups.
#' @param n.repl Number of replicates.
#' @param ... Currently not being used.
#'
#' @return A list of indexes.
#' @export
graph_spde_make_index <- function (name,
                                   graph_spde,
                                   n.group = 1,
                                   n.repl = 1,
                                   ...) {
    n.spde <- dim(graph_spde$graph_spde$V)[1]
    name.group <- paste(name, ".group", sep = "")
    name.repl <- paste(name, ".repl", sep = "")
    out <- list()
    out[[name]] <- rep(rep(1:n.spde, times = n.group), times = n.repl)
    out[[name.group]] <- rep(rep(1:n.group, each = n.spde), times = n.repl)
    out[[name.repl]] <- rep(1:n.repl, each = n.spde * n.group)
    return(out)
}


#' Observation/prediction matrices for 'rSPDE' models
#'
#' Constructs observation/prediction weight matrices
#' for metric graph models.
#'
#' @param graph_spde An `inla_metric_graph_spde` object built with the
#' `graph_spde()` function.
#' @param repl Which replicates? If there is no replicates, or to
#' use all replicates, one can set to `NULL`.
#' @return The observation matrix.
#' @export

graph_spde_make_A <- function (graph_spde, repl = NULL) {
   return(graph_spde$graph_spde$A(group = repl))
}


#' Data extraction for 'rSPDE' models
#'
#' Extracts data from metric graphs to be used by 'INLA' and 'inlabru'.
#'
#' @param graph_spde An `inla_metric_graph_spde` object built with the
#' `graph_spde()` function or an `rspde_metric_graph` object built with the
#' `rspde.metric_graph()` function from the 'rSPDE' package.
#' @param repl Which replicates? If there is no replicates, one
#' can set `repl` to `NULL`. If one wants all replicates,
#' then one sets to `repl` to `__all`.
#' @param only_pred Should only return the `data.frame` to the prediction data?
#' @param loc Character with the name of the location variable to be used in
#' 'inlabru' prediction.
#' @return An 'INLA' and 'inlabru' friendly list with the data.
#' @export

graph_data_spde <- function (graph_spde, repl = NULL,
                                only_pred = FALSE,
                                loc = NULL){
  graph_tmp <- graph_spde$graph_spde$clone()
  if(only_pred){
    idx_allNA <- !idx_not_all_NA(graph_tmp$data)
    graph_tmp$data <- lapply(graph_tmp$data, function(dat){return(dat[idx_allNA])})
  }

  if(is.null(repl)){
    groups <- graph_tmp$data[["__group"]]
    repl <- groups[1]
    ret <- select_group(graph_tmp$data, repl)
  } else if(repl[1] == "__all") {
    ret <- graph_tmp$data
  } else {
    ret <- select_group(graph_tmp$data, repl)
  }
  if(!is.null(loc)){
    ret[[loc]] <- cbind(ret[["__edge_number"]],
                          ret[["__distance_on_edge"]])
  }
  return(ret)
}

#' Extraction of vector of replicates for 'inlabru'
#'
#' Extracts the vector of replicates from an 'rSPDE'
#' model object for 'inlabru'
#'
#' @param graph_spde An `rspde_metric_graph` object built with the
#' `rspde.metric_graph()` function from the 'rSPDE' package.
#' @param repl Which replicates? If there is no replicates, one
#' can set `repl` to `NULL`. If one wants all replicates,
#' then one sets to `repl` to `__all`.
#' @return The vector of replicates.
#' @export

graph_repl_spde <- function (graph_spde, repl = NULL){
  graph_tmp <- graph_spde$graph_spde$clone()
  if(is.null(repl)){
    groups <- graph_tmp$data[["__group"]]
    repl <- groups[1]
    ret <- select_group(graph_tmp$data, repl)
  } else if(repl[1] == "__all") {
    ret <- graph_tmp$data
  } else {
    ret <- select_group(graph_tmp$data, repl)
  }
  return(ret[["__group"]])
}



#' @name spde_metric_graph_result
#' @title Metric graph SPDE result extraction from 'INLA' estimation results
#' @description Extract field and parameter values and distributions
#' for a metric graph spde effect from an 'INLA' result object.
#' @param inla An 'INLA' object obtained from a call to `inla()`.
#' @param name A character string with the name of the 'rSPDE' effect
#' in the model.
#' @param metric_graph_spde The `inla_metric_graph_spde` object used for the
#' random effect in the model.
#' @param compute.summary Should the summary be computed?
#' @param n_samples The number of samples to be used if parameterization is `matern`.
#' @param n_density The number of equally spaced points to estimate the density.
#' @return If the model was fitted with `matern` parameterization (the default),
#' it returns a list containing:
#' \item{marginals.range}{Marginal densities for the range parameter.}
#' \item{marginals.log.range}{Marginal densities for log(range).}
#' \item{marginals.sigma}{Marginal densities for std. deviation.}
#' \item{marginals.log.sigma}{Marginal densities for log(std. deviation).}
#' \item{marginals.values}{Marginal densities for the field values.}
#' \item{summary.log.range}{Summary statistics for log(range).}
#' \item{summary.log.sigma}{Summary statistics for log(std. deviation).}
#' \item{summary.values}{Summary statistics for the field values.}
#' If `compute.summary` is `TRUE`, then the list will also contain
#' \item{summary.kappa}{Summary statistics for kappa.}
#' \item{summary.tau}{Summary statistics for tau.}
#' If the model was fitted with the `spde` parameterization, it returns a list containing:
#' \item{marginals.kappa}{Marginal densities for kappa.}
#' \item{marginals.log.kappa}{Marginal densities for log(kappa).}
#' \item{marginals.log.tau}{Marginal densities for log(tau).}
#' \item{marginals.tau}{Marginal densities for tau.}
#' \item{marginals.values}{Marginal densities for the field values.}
#' \item{summary.log.kappa}{Summary statistics for log(kappa).}
#' \item{summary.log.tau}{Summary statistics for log(tau).}
#' \item{summary.values}{Summary statistics for the field values.}
#' If `compute.summary` is `TRUE`, then the list will also contain
#' \item{summary.kappa}{Summary statistics for kappa.}
#' \item{summary.tau}{Summary statistics for tau.}
#' @export

spde_metric_graph_result <- function(inla, name,
                                     metric_graph_spde,
                                     compute.summary = TRUE,
                                     n_samples = 5000,
                                     n_density = 1024) {
  if(!inherits(metric_graph_spde, "inla_metric_graph_spde")){
    stop("You should provide an inla_metric_graph_spde object!")
  }

  result <- list()

  parameterization <- metric_graph_spde$parameterization

    if(parameterization == "spde"){
      row_names <- c("tau", "kappa")
    } else{
      row_names <- c("sigma", "range")
    }

  result$summary.values <- inla$summary.random[[name]]

  if (!is.null(inla$marginals.random[[name]])) {
    result$marginals.values <- inla$marginals.random[[name]]
  }

  if(parameterization == "spde"){
    name_theta1 <- "reciprocal_tau"
    name_theta1_t <- "tau"
    name_theta2 <- "kappa"
  } else{
    name_theta1 <- "reciprocal_tau"
    name_theta1_t <- "sigma"
    name_theta2 <- "range"
  }


  result[[paste0("summary.log.",name_theta1)]] <- INLA::inla.extract.el(
    inla$summary.hyperpar,
    paste("Theta1 for ", name, "$", sep = "")
  )
  rownames(  result[[paste0("summary.log.",name_theta1)]]) <- paste0("log(",name_theta1,")")

  result[[paste0("summary.log.",name_theta2)]] <- INLA::inla.extract.el(
    inla$summary.hyperpar,
    paste("Theta2 for ", name, "$", sep = "")
  )
  rownames(result[[paste0("summary.log.",name_theta2)]]) <- paste0("log(", name_theta2,")")

  if (!is.null(inla$marginals.hyperpar[[paste0("Theta1 for ", name)]])) {
    result[[paste0("marginals.log.",name_theta1)]] <- INLA::inla.extract.el(
      inla$marginals.hyperpar,
      paste("Theta1 for ", name, "$", sep = "")
    )
    names(result[[paste0("marginals.log.",name_theta1)]]) <- name_theta1
    result[[paste0("marginals.log.",name_theta2)]] <- INLA::inla.extract.el(
      inla$marginals.hyperpar,
      paste("Theta2 for ", name, "$", sep = "")
    )
    names(result[[paste0("marginals.log.",name_theta2)]]) <- name_theta2

    if(parameterization == "spde"){
            result[[paste0("marginals.",name_theta1_t)]] <- lapply(
              result[[paste0("marginals.log.",name_theta1_t)]],
              function(x) {
                INLA::inla.tmarginal(
                  function(y) exp(-y),
                  x
                )
              }
            )
            result[[paste0("marginals.",name_theta2)]] <- lapply(
              result[[paste0("marginals.log.",name_theta2)]],
              function(x) {
                INLA::inla.tmarginal(
                  function(y) exp(y),
                  x
                )
              }
            )
    } else{
            hyperpar_sample <- INLA::inla.hyperpar.sample(n_samples, inla)
            reciprocal_tau_est <- exp(hyperpar_sample[, paste0('Theta1 for ',name)])
            tau_est <- 1/reciprocal_tau_est
            range_est <- exp(hyperpar_sample[, paste0('Theta2 for ',name)])
            kappa_est <- sqrt(4)/range_est
            sigma_est <- sqrt(gamma(0.5) / (tau_est^2 * kappa_est^(2 * 0.5) *
                    (4 * pi)^(1 / 2) * gamma(0.5 + 1 / 2)))

            density_sigma <- stats::density(sigma_est, n = n_density)

            result[[paste0("marginals.",name_theta1_t)]] <- list()
            result[[paste0("marginals.",name_theta1_t)]][[name_theta1_t]] <- cbind(density_sigma$x, density_sigma$y)
            colnames(result[[paste0("marginals.",name_theta1_t)]][[name_theta1_t]]) <- c("x","y")

            result[[paste0("marginals.",name_theta2)]] <- lapply(
              result[[paste0("marginals.log.",name_theta2)]],
              function(x) {
                INLA::inla.tmarginal(
                  function(y) exp(y),
                  x
                )
              }
            )
    }
  }

  if (compute.summary) {
    norm_const <- function(density_df) {
      min_x <- min(density_df[, "x"])
      max_x <- max(density_df[, "x"])
      denstemp <- function(x) {
        dens <- sapply(x, function(z) {
          if (z < min_x) {
            return(0)
          } else if (z > max_x) {
            return(0)
          } else {
            return(approx(x = density_df[, "x"],
            y = density_df[, "y"], xout = z)$y)
          }
        })
        return(dens)
      }
      norm_const <- stats::integrate(
        f = function(z) {
          denstemp(z)
        }, lower = min_x, upper = max_x,
        subdivisions = nrow(density_df)
      )$value
      return(norm_const)
    }

    norm_const_theta1 <- norm_const(result[[paste0("marginals.",name_theta1_t)]][[name_theta1_t]])
    result[[paste0("marginals.",name_theta1_t)]][[name_theta1_t]][, "y"] <-
    result[[paste0("marginals.",name_theta1_t)]][[name_theta1_t]][, "y"] / norm_const_theta1

    norm_const_theta2 <- norm_const(result[[paste0("marginals.",name_theta2)]][[name_theta2]])
    result[[paste0("marginals.",name_theta2)]][[name_theta2]][, "y"] <-
    result[[paste0("marginals.",name_theta2)]][[name_theta2]][, "y"] / norm_const_theta2




    result[[paste0("summary.",name_theta1_t)]] <- create_summary_from_density(result[[paste0("marginals.",name_theta1_t)]][[name_theta1_t]],
    name = name_theta1_t)
    result[[paste0("summary.",name_theta2)]] <-
    create_summary_from_density(result[[paste0("marginals.",name_theta2)]][[name_theta2]], name = name_theta2)
  }

  class(result) <- "metric_graph_spde_result"
  result$params <- c(name_theta1_t,name_theta2)
  return(result)
}


#' Data frame for metric_graph_spde_result objects to be used in 'ggplot2'
#'
#' Returns a 'ggplot2'-friendly data-frame with the marginal posterior densities.
#' @aliases gg_df gg_df.metric_graph_spde_result
#' @param result A metric_graph_spde_result object.
#' @param parameter Vector. Which parameters to get the posterior density in the
#' data.frame? The options are `sigma`, `range` or `kappa`.
#' @param transform Should the posterior density be given in the original scale?
#' @param restrict_x_axis Variables to restrict the range of x axis based on quantiles.
#' @param restrict_quantiles List of quantiles to restrict x axis.
#' @param ... Not being used.
#'
#' @return A `data.frame` containing the posterior densities.
#' @export
gg_df.metric_graph_spde_result <- function(result,
                          parameter = result$params,
                          transform = TRUE,
                          restrict_x_axis = parameter,
                          restrict_quantiles = list(sigma = c(0,1),
                          range = c(0,1),
                          kappa = c(0,1),
                          sigma = c(0,1)),...) {
      parameter <- intersect(parameter, c("kappa", "range", "sigma"))
      if(length(parameter) == 0){
        stop("You should choose at least one of the parameters 'kappa', 'range' or 'sigma'!")
      }

  spde_result <- result
  param <- parameter[[1]]
  if(transform){
    param <- paste0("marginals.", param)
  } else{
      param <- paste0("marginals.log.", param)
  }
  ret_df <- data.frame(x = spde_result[[param]][[parameter[1]]][,1],
  y = spde_result[[param]][[parameter[1]]][,2],
  parameter = parameter[[1]])

  if(parameter[[1]] %in% restrict_x_axis){
    if(is.null( restrict_quantiles[[parameter[[1]]]])){
      warning("If you want to restrict x axis you should provide a quantile for the parameter!")
       restrict_quantiles[[parameter[[1]]]] <- c(0,1)
    }
    d_t <- c(0,diff(ret_df$x))
    emp_cdf <- cumsum(d_t*ret_df$y)
    lower_quant <- restrict_quantiles[[parameter[[1]]]][1]
    upper_quant <- restrict_quantiles[[parameter[[1]]]][2]
    filter_coord <- (emp_cdf >= lower_quant) * (emp_cdf <= upper_quant)
    filter_coord <- as.logical(filter_coord)
    ret_df <- ret_df[filter_coord, ]
  }

  if(length(parameter) > 1){
  for(i in 2:length(parameter)){
  param <- parameter[[i]]
  if(transform){
    param <- paste0("marginals.", param)
  } else{
      param <- paste0("marginals.log.", param)
  }
    tmp <- data.frame(x = spde_result[[param]][[parameter[i]]][,1],
      y = spde_result[[param]][[parameter[i]]][,2],
      parameter = parameter[[i]])

    if(parameter[[i]] %in% restrict_x_axis){
    if(is.null( restrict_quantiles[[parameter[[i]]]])){
      warning(paste("No quantile for", parameter[[i]]))
      warning("If you want to restrict x axis you should provide a quantile for the parameter!")
       restrict_quantiles[[parameter[[i]]]] <- c(0,1)
    }
    d_t <- c(0,diff(tmp$x))
    emp_cdf <- cumsum(d_t*tmp$y)
    lower_quant <- restrict_quantiles[[parameter[[i]]]][1]
    upper_quant <- restrict_quantiles[[parameter[[i]]]][2]
    filter_coord <- (emp_cdf >= lower_quant) * (emp_cdf <= upper_quant)
    filter_coord <- as.logical(filter_coord)
    tmp <- tmp[filter_coord, ]
  }

    ret_df <- rbind(ret_df, tmp)
  }
  }
  return(ret_df)
}


#' Summary for posteriors of field parameters for an `inla_rspde`
#' model from a `rspde.result` object
#'
#' Summary for posteriors of 'rSPDE' field parameters in
#' their original scales.
#'
#' @param object A `rspde.result` object.
#' @param digits Integer, used for number formatting with signif()
#' @param ... Currently not used.
#'
#' @return A `data.frame` containing the summary.
#' @export
#' @method summary metric_graph_spde_result

summary.metric_graph_spde_result <- function(object,
                                 digits = 6,
                                 ...) {

  if (is.null(object[[paste0("summary.",object$params[1])]])) {
    warning("The summary was not computed, rerun spde_metric_graph_result with
    compute.summary set to TRUE.")
  } else {
    out <- object[[paste0("summary.",object$params[1])]]
    out <- rbind(out, object[[paste0("summary.",object$params[2])]])
    return(signif(out, digits = digits))
  }
}


#'
#' @title Metric graph 'inlabru' mapper
#' @name bru_mapper.inla_metric_graph_spde
#' @param model An `inla_metric_graph_spde` for which to construct or extract a mapper
#' @param \dots Arguments passed on to other methods
#' @rdname bru_mapper.inla_metric_graph_spde
#' @rawNamespace if (getRversion() >= "3.6.0") {
#'   S3method(inlabru::bru_get_mapper, inla_metric_graph_spde)
#'   S3method(inlabru::ibm_n, bru_mapper_inla_metric_graph_spde)
#'   S3method(inlabru::ibm_values, bru_mapper_inla_metric_graph_spde)
#'   S3method(inlabru::ibm_jacobian, bru_mapper_inla_metric_graph_spde)
#' }
#'

bru_get_mapper.inla_metric_graph_spde <- function(model, ...){
  mapper <- list(model = model)
  inlabru::bru_mapper_define(mapper, new_class = "bru_mapper_inla_metric_graph_spde")
}

#' @param mapper A `bru_mapper.inla_metric_graph_spde` object
#' @rdname bru_mapper.inla_metric_graph_spde
ibm_n.bru_mapper_inla_metric_graph_spde <- function(mapper, ...) {
  model <- mapper[["model"]]
  n_groups <- length(unique(model$graph_spde$data[["__group"]]))
  return(model$f$n)
}
#' @rdname bru_mapper.inla_metric_graph_spde
ibm_values.bru_mapper_inla_metric_graph_spde <- function(mapper, ...) {
  seq_len(inlabru::ibm_n(mapper))
}
#' @param input The values for which to produce a mapping matrix
#' @rdname bru_mapper.inla_metric_graph_spde
ibm_jacobian.bru_mapper_inla_metric_graph_spde <- function(mapper, input, ...) {
  model <- mapper[["model"]]
  pte_tmp <- model$graph_spde$get_PtE()
  input_list <- lapply(1:nrow(input), function(i){input[i,]})
  pte_tmp_list <- lapply(1:nrow(pte_tmp), function(i){pte_tmp[i,]})
  idx_tmp <- match(input_list, pte_tmp_list)
  A_tmp <- model$graph_spde$A()
  return(A_tmp[idx_tmp,])
}





#' @name create_summary_from_density
#' @title Creates a summary from a density data frame
#' @description Auxiliary function to create summaries from density data frames.
#' @param density_df A density data frame.
#' @param name Name of the parameter.
#' @return A data frame containing a basic summary.
#' @noRd

create_summary_from_density <- function(density_df, name) {
  min_x <- min(density_df[, "x"])
  max_x <- max(density_df[, "x"])
  denstemp <- function(x) {
    dens <- sapply(x, function(z) {
      if (z < min_x) {
        return(0)
      } else if (z > max_x) {
        return(0)
      } else {
        return(approx(x = density_df[, "x"], y = density_df[, "y"], xout = z)$y)
      }
    })
    return(dens)
  }

  ptemp <- function(q) {
    prob_temp <- sapply(q, function(v) {
      if (v <= min_x) {
        return(0)
      } else if (v >= max_x) {
        return(1)
      } else {
        stats::integrate(
          f = denstemp, lower = min_x, upper = v,
          subdivisions = min(nrow(density_df), 500)
        )$value
      }
    })
    return(prob_temp)
  }

  mean_temp <- stats::integrate(
    f = function(z) {
      denstemp(z) * z
    }, lower = min_x, upper = max_x,
    subdivisions = nrow(density_df)
  )$value

  sd_temp <- sqrt(stats::integrate(
    f = function(z) {
      denstemp(z) * (z - mean_temp)^2
    }, lower = min_x, upper = max_x,
    subdivisions = nrow(density_df)
  )$value)

  mode_temp <- density_df[which.max(density_df[, "y"]), "x"]

  qtemp <- function(p) {
    quant_temp <- sapply(p, function(x) {
      if (x < 0 | x > 1) {
        return(NaN)
      } else {
        return(stats::uniroot(function(y) {
          ptemp(y) - x
        }, lower = min_x, upper = max_x)$root)
      }
    })
    return(quant_temp)
  }

  out <- data.frame(
    mean = mean_temp, sd = sd_temp, `0.025quant` = qtemp(0.025),
    `0.5quant` = qtemp(0.5), `0.975quant` = qtemp(0.975), mode = mode_temp
  )
  rownames(out) <- name
  colnames(out) <- c("mean", "sd", "0.025quant",
  "0.5quant", "0.975quant", "mode")
  return(out)
}


#' @name bru_graph_rep
#' @title Creates a vector of replicates to be used with 'inlabru'
#' @description Auxiliary function to create a vector of replicates to be used
#' with 'inlabru'.
#' @param repl A vector of replicates. If set to `__all`, a vector
#' for all replicates will be generated.
#' @param graph_spde Name of the field.
#' @return A vector of replicates to be used with 'inlabru'.
#' @export

bru_graph_rep <- function(repl, graph_spde){
  groups <- unique(graph_spde$graph_spde$data[["__group"]])
  if(repl[1] == "__all"){
    repl <- groups
  }
  n_groups <- length(groups)
  length_resp <- sum(graph_spde$graph_spde$data[["__group"]] == groups[1])
  return(rep(repl, each = length_resp ))
}

#' @name predict.inla_metric_graph_spde
#' @title Predict method for 'inlabru' fits on Metric Graphs
#' @description Auxiliar function to obtain predictions of the field
#' using 'inlabru'.
#' @param object An `inla_metric_graph_spde` object built with the `graph_spde()`
#' function.
#' @param cmp The 'inlabru' component used to fit the model.
#' @param bru_fit A fitted model using 'inlabru' or 'INLA'.
#' @param data A data.frame of covariates needed for the prediction. The
#' locations must be normalized PtE.
#' @param formula A formula where the right hand side defines an R expression to
#' evaluate for each generated sample. If NULL, the latent and hyperparameter
#' states are returned as named list elements. See Details for more information.
#' @param data_coords It decides which coordinate system to use. If `PtE`, the
#' user must provide the locations as a data frame with the first column being
#' the edge number and the second column as the distance on edge, otherwise if
#' `euclidean`, the user must provide a data frame with the first column being
#' the `x` Euclidean coordinates and the second column being the `y` Euclidean
#' coordinates.
#' @param  normalized if `TRUE`, then the distances in distance on edge are
#' assumed to be normalized to (0,1). Default TRUE. Will not be
#' used if `data_coords` is `euclidean`.
#' @param n.samples Integer setting the number of samples to draw in order to
#' calculate the posterior statistics. The default is rather low but provides a
#' quick approximate result.
#' @param seed Random number generator seed passed on to `inla.posterior.sample()`
#' @param probs	A numeric vector of probabilities with values in the standard
#' unit interval to be passed to stats::quantile
#' @param return_original_order Should the predictions be returned in the
#' original order?
#' @param num.threads	Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to 'INLA'. When seed != 0, overridden to "1:1"
#' @param include	Character vector of component labels that are needed by the
#' predictor expression; Default: NULL (include all components that are not
#' explicitly excluded)
#' @param exclude	Character vector of component labels that are not used by the
#' predictor expression. The exclusion list is applied to the list as determined
#' by the include parameter; Default: NULL (do not remove any components from
#' the inclusion list)
#' @param drop logical; If keep=FALSE, data is a SpatialDataFrame, and the
#' prediciton summary has the same number of rows as data, then the output is a
#' SpatialDataFrame object. Default FALSE.
#' @param... Additional arguments passed on to `inla.posterior.sample()`.
#' @return A list with predictions.
#' @export

predict.inla_metric_graph_spde <- function(object,
                                           cmp,
                                           bru_fit,
                                           data = NULL,
                                           formula = NULL,
                                           data_coords = c("PtE", "euclidean"),
                                           normalized = TRUE,
                                           n.samples = 100,
                                           seed = 0L,
                                           probs = c(0.025, 0.5, 0.975),
                                           return_original_order = TRUE,
                                           num.threads = NULL,
                                           include = NULL,
                                           exclude = NULL,
                                           drop = FALSE,
                                           ...){
  data_coords <- data_coords[[1]]
  if(!(data_coords %in% c("PtE", "euclidean"))){
    stop("data_coords must be either 'PtE' or 'euclidean'!")
  }
  graph_tmp <- object$graph_spde$get_initial_graph()
  # graph_tmp <- object$graph_spde$clone()
  name_locations <- bru_fit$bru_info$model$effects$field$main$input$input
  original_data <- object$graph_spde$data
  original_data[["__edge_number"]] <- object$data_PtE[,1]
  original_data[["__distance_on_edge"]] <- object$data_PtE[,2]

  graph_tmp$add_observations(data = original_data,
                  edge_number = "__edge_number",
                  distance_on_edge = "__distance_on_edge",
                  data_coords = "PtE",
                  normalized = TRUE)

  new_data <- data
  new_data[[name_locations]] <- NULL
  n_locations <- nrow(data[[name_locations]])
  names_columns <- names(original_data)
  names_columns <- setdiff(names_columns, c("__group", "__coord_x",
                                            "__coord_y", "__edge_number",
                                            "__distance_on_edge"))
  # for(name_column in names_columns){
  #   new_data[[name_column]] <- rep(NA, n_locations)
  # }
  if(data_coords == "PtE"){
    new_data[["__edge_number"]] <- data[[name_locations]][,1]
    new_data[["__distance_on_edge"]] <- data[[name_locations]][,2]
  } else{
    new_data[["__coord_x"]] <- data[[name_locations]][,1]
    new_data[["__coord_y"]] <- data[[name_locations]][,2]
  }

  new_data[["__dummy_var"]] <- 1:length(new_data[["__edge_number"]])

  graph_tmp$add_observations(data = new_data,
                  edge_number = "__edge_number",
                  distance_on_edge = "__distance_on_edge",
                  coord_x = "__coord_x",
                  coord_y = "__coord_y",
                  data_coords = data_coords,
                  normalized = normalized)

  dummy1 <- graph_tmp$data[["__dummy_var"]]

  graph_tmp$data[["__dummy_var2"]] <- 1:length(graph_tmp$data[["__dummy_var"]])

  pred_PtE <- cbind(graph_tmp$data[["__edge_number"]],
                          graph_tmp$data[["__distance_on_edge"]])

  # pred_PtE <- pred_PtE[!is.na(dummy1),]

  # Adding the original data

  # graph_tmp$add_observations(data = original_data,
  #                   coord_x = "__coord_x",
  #                   coord_y = "__coord_y",
  #                   data_coords = "euclidean")

  graph_tmp$observation_to_vertex()

  # tmp_list2 <- cbind(graph_tmp$data[["__coord_x"]],
  #                                       graph_tmp$data[["__coord_y"]])
  # tmp_list2 <- lapply(1:nrow(tmp_list2), function(i){tmp_list2[i,]})
  # idx_list <- match(tmp_list, tmp_list2)

  new_data_list <- graph_tmp$data

  idx_list <- !is.na(new_data_list[["__dummy_var"]])

  new_data_list <- lapply(new_data_list, function(dat){dat[idx_list]})

  pred_PtE <- pred_PtE[graph_tmp$data[["__dummy_var2"]],][idx_list,]

  # new_data_list[[name_locations]] <- cbind(graph_tmp$data[["__edge_number"]][idx_list],
  #                                             graph_tmp$data[["__distance_on_edge"]][idx_list])

  new_data_list[[name_locations]] <- cbind(new_data_list[["__edge_number"]],
                                              new_data_list[["__distance_on_edge"]])

  spde____model <- graph_spde(graph_tmp)
  cmp_c <- as.character(cmp)
  name_model <- deparse(substitute(object))
  cmp_c[3] <- sub(name_model, "spde____model", cmp_c[3])
  cmp <- as.formula(paste(cmp_c[2], cmp_c[1], cmp_c[3]))
  bru_fit_new <- inlabru::bru(cmp,
          data = graph_data_spde(spde____model, loc = name_locations))
  pred <- predict(object = bru_fit_new,
                    data = new_data_list,
                    formula = formula,
                    n.samples = n.samples,
                    seed = seed,
                    probs = probs,
                    num.threads = num.threads,
                    include = include,
                    exclude = exclude,
                    drop = drop,
                    ...)
  pred_list <- list()
  pred_list[["pred"]] <- pred
  pred_list[["PtE_pred"]] <- pred_PtE
  if(return_original_order){
    ord <- graph_tmp$data[["__dummy_var"]][idx_list]
    pred_list[["pred"]][ord,] <- pred
    pred_list[["PtE_pred"]][ord,] <- pred_PtE
  }
  pred_list[["initial_graph"]] <- graph_tmp$get_initial_graph()
  # pred_list[["new_model"]] <- spde____model
  # pred_list[["new_fit"]] <- bru_fit_new
  # pred_list[["new_graph"]] <- graph_tmp

  class(pred_list) <- "graph_bru_pred"
  return(pred_list)
}


#' @name plot.graph_bru_pred
#' @title Plot of predicted values with 'inlabru'
#' @description Auxiliary function to obtain plots of the predictions of the field
#' using 'inlabru'.
#' @param x A predicted object obtained with the `predict` method.
#' @param y Not used.
#' @param vertex_size Size of the vertices.
#' @param ... Additional parameters to be passed to the plot function.
#' @return A 'ggplot2' object.
#' @export

plot.graph_bru_pred <- function(x, y = NULL, vertex_size = 0, ...){
  m_prd_bru <- x$pred$mean
  PtE_prd <- x$PtE_pred
  p <- x$initial_graph$plot_function(X = cbind(PtE_prd, m_prd_bru),
                                      vertex_size = vertex_size, ...)
  p
}

#' @name predict.rspde_metric_graph
#' @title Predict method for 'inlabru' fits on Metric Graphs for 'rSPDE' models
#' @description Auxiliar function to obtain predictions of the field
#' using 'inlabru' and 'rSPDE'.
#' @param object An `rspde_metric_graph` object built with the
#' `rspde.metric_graph()` function.
#' @param cmp The 'inlabru' component used to fit the model.
#' @param bru_fit A fitted model using 'inlabru' or 'INLA'.
#' @param data A data.frame of covariates needed for the prediction. The locations
#' must be normalized PtE.
#' @param formula A formula where the right hand side defines an R expression to
#' evaluate for each generated sample. If NULL, the latent and hyperparameter
#' states are returned as named list elements. See Details for more information.
#' @param data_coords It decides which coordinate system to use. If `PtE`, the
#' user must provide the locations as a data frame with the first column being
#' the edge number and the second column as the distance on edge, otherwise if
#' `euclidean`, the user must provide a data frame with the first column being
#' the `x` Euclidean coordinates and the second column being the `y` Euclidean
#' coordinates.
#' @param  normalized if `TRUE`, then the distances in distance on edge are
#' assumed to be normalized to (0,1). Default TRUE. Will not be used if
#' `data_coords` is `euclidean`.
#' @param n.samples Integer setting the number of samples to draw in order to
#' calculate the posterior statistics. The default is rather low but provides a
#' quick approximate result.
#' @param seed Random number generator seed passed on to inla.posterior.sample
#' @param probs	A numeric vector of probabilities with values in the standard
#' unit interval to be passed to stats::quantile.
#' @param num.threads	Specification of desired number of threads for parallel
#' computations. Default NULL, leaves it up to 'INLA'. When seed != 0, overridden
#' to "1:1"
#' @param include	Character vector of component labels that are needed by the
#' predictor expression; Default: NULL (include all components that are not
#' explicitly excluded)
#' @param exclude	Character vector of component labels that are not used by the
#' predictor expression. The exclusion list is applied to the list as determined
#' by the include parameter; Default: NULL (do not remove any components from the
#' inclusion list)
#' @param drop logical; If keep=FALSE, data is a SpatialDataFrame, and the
#' prediciton summary has the same number of rows as data, then the output is a
#' SpatialDataFrame object. Default FALSE.
#' @param... Additional arguments passed on to inla.posterior.sample.
#' @return A list with predictions.
#' @export

predict.rspde_metric_graph <- function(object,
                                           cmp,
                                           bru_fit,
                                           data = NULL,
                                           formula = NULL,
                                           data_coords = c("PtE", "euclidean"),
                                           normalized = TRUE,
                                           n.samples = 100,
                                           seed = 0L,
                                           probs = c(0.025, 0.5, 0.975),
                                           num.threads = NULL,
                                           include = NULL,
                                           exclude = NULL,
                                           drop = FALSE,
                                           ...){
  data_coords <- data_coords[[1]]
  if(!(data_coords %in% c("PtE", "euclidean"))){
    stop("data_coords must be either 'PtE' or 'euclidean'!")
  }
  graph_tmp <- object$graph_spde$get_initial_graph()
  name_locations <- bru_fit$bru_info$model$effects$field$main$input$input

  if(data_coords == "PtE"){
    if(normalized){
      pred_PtE <- data[[name_locations]]
    } else{
      pred_PtE <- data[[name_locations]]
      pred_PtE[,2] <- pred_PtE[,2]/graph_tmp$edge_lengths[pred_PtE[, 1]]
    }
  } else{
    pred_PtE <- graph_tmp$coordinates(XY = data[[name_locations]])
  }

  pred <- predict(object = bru_fit,
                    data = data,
                    formula = formula,
                    n.samples = n.samples,
                    seed = seed,
                    probs = probs,
                    num.threads = num.threads,
                    include = include,
                    exclude = exclude,
                    drop = drop,
                    ...)
  pred_list <- list()
  pred_list[["pred"]] <- pred
  pred_list[["PtE_pred"]] <- pred_PtE
  pred_list[["initial_graph"]] <- graph_tmp$get_initial_graph()

  class(pred_list) <- "graph_bru_pred"
  return(pred_list)
}


#' @name graph_bru_process_data
#' @title Prepare data frames or data lists to be used with 'inlabru' in metric
#' graphs
#' @param data A `data.frame` or a `list` containing the covariates, the edge
#' number and the distance on edge for the locations to obtain the prediction.
#' @param loc character. Name of the locations to be used in 'inlabru' component.
#' @param edge_number Name of the variable that contains the edge number, the
#' default is `edge_number`.
#' @param distance_on_edge Name of the variable that contains the distance on
#' edge, the default is `distance_on_edge`.
#' @return A list containing the processed data to be used in a user-friendly
#' manner by 'inlabru'.
#' @export

graph_bru_process_data <- function(data, edge_number = "edge_number",
                                        distance_on_edge = "distance_on_edge",
                                        loc = "loc"){

                                        data[[loc]] <- cbind(data[[edge_number]], data[[distance_on_edge]])
                                        data[[edge_number]] <- NULL
                                        data[[distance_on_edge]] <- NULL
                                        return(data)
                                        }
