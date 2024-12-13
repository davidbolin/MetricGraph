
#' 'INLA' implementation of Whittle-Matérn fields for metric graphs
#'
#' This function creates an 'INLA' object that can be used
#' in 'INLA' or 'inlabru' to fit Whittle-Matérn fields on metric graphs.
#'
#' @param graph_object A `metric_graph` object.
#' @param alpha The order of the SPDE.
#' @param directional Should a directional model be used? Currently only implemented for `alpha=1`. 
#' @param stationary_endpoints Which vertices of degree 1 should contain
#' stationary boundary conditions? Set to "all" for all vertices of degree 1, "none" for none of the vertices of degree 1, or pass the indices of the vertices of degree 1 for which stationary conditions are desired.
#' @param parameterization Which parameterization to be used? The options are
#' 'matern' (sigma and range) and 'spde' (sigma and kappa).
#' @param start_sigma Starting value for sigma.
#' @param prior_sigma a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation of sigma on the log scale.
#' @param start_tau Starting value for tau.
#' @param prior_tau a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation of tau on the log scale.
#' @param start_range Starting value for range parameter.
#' @param prior_range a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation of the range parameter on
#' the log scale. Will not be used if prior.kappa is non-null.
#' @param start_kappa Starting value for kappa.
#' @param prior_kappa a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation of kappa on the log scale.
#' @param factor_start_range Factor to multiply the max/min dimension of the bounding box to obtain a starting value for range. Default is 0.3.
#' @param type_start_range_bbox Which dimension from the bounding box should be used? The options are 'diag', the default, 'max' and 'min'.
#' @param shared_lib Which shared lib to use for the cgeneric implementation?
#' If "detect", it will check if the shared lib exists locally, in which case it will
#' use it. Otherwise it will use 'INLA's shared library.
#' If 'INLA', it will use the shared lib from 'INLA's installation. If 'rSPDE', then
#' it will use the local installation of the rSPDE package (does not work if your installation is from CRAN).
#' Otherwise, you can directly supply the path of the .so (or .dll) file.
#' @param debug Should debug be displayed?
#' @param verbose Level of verbosity. 0 is silent, 1 prints basic information, 2 prints more.
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
                       directional = FALSE,
                       stationary_endpoints = "all",
                       parameterization = c("matern", "spde"),
                       start_range = NULL,
                       prior_range = NULL,
                       start_kappa = NULL,
                       prior_kappa = NULL,
                       start_sigma = NULL,
                       prior_sigma = NULL,
                       start_tau = NULL,
                       prior_tau = NULL,
                       factor_start_range = 0.3,
                       type_start_range_bbox = "diag",
                       shared_lib = "detect",
                       debug = FALSE,
                       verbose = 0){

  if(!(alpha%in%c(1,2))){
    stop("alpha must be either 1 or 2!")
  }

  graph_spde <- graph_object$clone()

  if(verbose>0){
    message("Turning observations into vertices...")
  }

  if(!is.null(graph_spde$.__enclos_env__$private$data)){
    graph_spde$observation_to_vertex(mesh_warning=FALSE, verbose = verbose)
  }

  parameterization <- parameterization[[1]]
  if(!(alpha%in%c(1,2))){
    stop("alpha must be either 1 or 2!")
  }

  A_tmp <- NULL
  Tc <- NULL
  nu <- alpha - 0.5
  V <- graph_spde$V
  EtV <- graph_spde$E
  El <- graph_spde$edge_lengths

  if(verbose>0){
    message("Constructing the sparsity graph...")
  }

  i_ <- j_ <- rep(0, dim(V)[1]*4)
  nE <- dim(EtV)[1]
  count <- 0
  if(alpha == 1){
    if(!directional){
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
    } else{
    weights_directional <- 0
    Q_tmp <- Qalpha1_edges(theta = c(1,1), graph = graph_spde, BC=BC, stationary_points=stationary_endpoints, w = weights_directional)
    
    if(is.character(stationary_endpoints)){
      stationary_endpoints <- stationary_endpoints[[1]]
      if(!(stationary_endpoints %in% c("all", "none"))){
        stop("If stationary_endpoints is a string, it must be either 'all' or 'none', otherwise it must be a numeric vector.")
      }
      stat_indices <- which(graph_spde$get_degrees("indegree")==0)
    } else{
      stat_indices <- stationary_endpoints
      if(!is.numeric(stat_indices)){
        stop("stationary_endpoints must be either numeric or a string.")
      }
    }
    if(stationary_endpoints == "none"){
      BC <- 0
    } else{
      BC <- 1
    }

    ind_stat_indices <- NULL

    for (v in stat_indices) {
      edge <- which(graph_spde$E[,1]==v)[1] #only put stationary of one of indices
      ind_stat_indices <- c(ind_stat_indices, 2 *  (edge-1))
    }
    
    if(is.null(graph_spde$CoB)){
      graph_spde$buildDirectionalConstraints(alpha = 1)
    } else if(graph_spde$CoB$alpha == 2){
      graph_spde$buildDirectionalConstraints(alpha = 1)
    }
    n_const <- length(graph_spde$CoB$S)
    ind.const <- c(1:n_const)
    Tc <- graph_spde$CoB$T[-ind.const, ]                      
    Q_tmp <- Tc%*%Q_tmp%*%t(Tc)

    Q_tmp <- INLA::inla.as.sparse(Q_tmp)
    ii <- Q_tmp@i
    Q_tmp@i <- Q_tmp@j
    Q_tmp@j <- ii
    idx <- which(Q_tmp@i <= Q_tmp@j)
    Q_tmp@i <- Q_tmp@i[idx]
    Q_tmp@j <- Q_tmp@j[idx]
    Q_tmp@x <- Q_tmp@x[idx]
    
    i_ <- Q_tmp@i
    j_ <- Q_tmp@j

    Tc <- as(Tc, "TsparseMatrix")
    i_Tc <- Tc@i
    j_Tc <- Tc@j
    x_Tc <- Tc@x
    }
  } else if(alpha == 2){
    if(directional){
      stop("Directional models are currently not implemented for 'alpha=2'.")
    }
    if(stationary_endpoints == "all"){
        i.table <- table(c(graph_spde$E))
        index <- as.integer(names(which(i.table == 1)))
        BC = 1
    } else if(stationary_endpoints == "none"){
      index <- NULL
      BC = 0
    } else{
        index <- stationary_endpoints - 1
        BC = 1
    }

    Q_tmp <- Qalpha2(theta = c(1,1), graph = graph_spde, BC=BC, stationary_points=index)
    if(verbose>0){
      message("Checking/Computing constraint matrix...")
    }
    if(is.null(graph_spde$CoB)){
      graph_spde$buildC(2, edge_constraint = BC)
    } else if(graph_spde$CoB$alpha == 1){
      graph_spde$buildC(2, edge_constraint = BC)
    }
    n_const <- length(graph_spde$CoB$S)
    ind.const <- c(1:n_const)
    Tc <- graph_spde$CoB$T[-ind.const, ]                      
    Q_tmp <- Tc%*%Q_tmp%*%t(Tc)

    Q_tmp <- INLA::inla.as.sparse(Q_tmp)
    ii <- Q_tmp@i
    Q_tmp@i <- Q_tmp@j
    Q_tmp@j <- ii
    idx <- which(Q_tmp@i <= Q_tmp@j)
    Q_tmp@i <- Q_tmp@i[idx]
    Q_tmp@j <- Q_tmp@j[idx]
    Q_tmp@x <- Q_tmp@x[idx]
    
    i_ <- Q_tmp@i
    j_ <- Q_tmp@j

    if(is.null(index)){
      index <- -1
    }
    Tc <- as(Tc, "TsparseMatrix")
    i_Tc <- Tc@i
    j_Tc <- Tc@j
    x_Tc <- Tc@x
    lower.edges <- NULL
    upper.edges <- NULL
    if(!is.null(index)){
      lower.edges <- which(graph_spde$E[, 1] %in% index)
      upper.edges <- which(graph_spde$E[, 2] %in% index)
    }

    lower_edges_len <- length(lower.edges)
    upper_edges_len <- length(upper.edges)
    if(length(lower.edges) == 0){
      lower.edges <- -1
    }
    if(length(upper.edges)==0){
      upper.edges <- -1
    }
  }


  ## End of alpha = 2 
  
    if(is.null(prior_kappa$meanlog) && is.null(prior_range$meanlog)){
      model_start <- ifelse(alpha==1,"alpha1", "alpha2")
      if(verbose>0){
          message("Computing starting values...")
      }
      start_values_vector <- graph_starting_values(graph_spde,
                      model = model_start, data=FALSE,
                      factor_start_range = factor_start_range,
                      type_start_range_bbox = type_start_range_bbox)$start_values

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
    # the prior for sigma
    prior_sigma$meanlog <- 0
  }
  # converting to reciprocal tau
  const_tmp <-  sqrt(gamma(nu) / (exp(prior_kappa$meanlog)^(2 * nu) * (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))   
  prior_sigma$meanlog <- log(exp(prior_sigma$meanlog)/const_tmp)

  if(!is.null(prior_tau$meanlog)){
    prior_sigma$meanlog <- -prior_tau$meanlog
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
  # converting to reciprocal tau
  const_tmp <-  sqrt(gamma(nu) / (exp(start_lkappa)^(2 * nu) * (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))   
  start_lsigma <- log(exp(start_lsigma)/const_tmp)
  if(!is.null(start_tau)){
    start_lsigma <- -log(start_tau)
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
  } else if(shared_lib == "rSPDE"){
    gpgraph_lib <- system.file('shared', package='rSPDE')
    if(Sys.info()['sysname']=='Windows') {
		gpgraph_lib <- paste0(gpgraph_lib, "/rspde_cgeneric_models.dll")
            } else {
		gpgraph_lib <- paste0(gpgraph_lib, "/rspde_cgeneric_models.so")
            }
  } else if(shared_lib == "detect"){
    gpgraph_lib_local <- system.file('shared', package='rSPDE')
    if(Sys.info()['sysname']=='Windows') {
		gpgraph_lib_local <- paste0(gpgraph_lib_local, "/rspde_cgeneric_models.dll")
            } else {
		gpgraph_lib_local <- paste0(gpgraph_lib_local, "/rspde_cgeneric_models.so")
            }
    if(file.exists(gpgraph_lib_local)){
      gpgraph_lib <- gpgraph_lib_local
    } else{
      gpgraph_lib <- INLA::inla.external.lib('rSPDE')
    }
  }

  if(verbose > 0){
    message("Creating INLA model...")
  }

if(alpha == 1){
  if(!directional){
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
  } else{
      model <-
        do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_gpgraph_alpha1_directional_model",
            shlib=gpgraph_lib,
            n=dim(Q_tmp)[1], debug=debug,
            prec_graph_i = as.integer(i_),
            prec_graph_j = as.integer(j_),
            Tc = Tc,
            El = El,
            start_theta = start_theta,
            start_lsigma = start_lsigma,
            prior_theta_meanlog = prior_theta$meanlog,
            prior_theta_sdlog = prior_theta$sdlog,
            prior_sigma_meanlog = prior_sigma$meanlog,
            prior_sigma_sdlog = prior_sigma$sdlog,
            parameterization = parameterization,
            w = weights_directional,
            BC = as.integer(BC),
            ind_stat_indices = as.integer(ind_stat_indices)))
  }

} else{
    model <-
        do.call(eval(parse(text='INLA::inla.cgeneric.define')),
        list(model="inla_cgeneric_gpgraph_alpha2_model",
            shlib=gpgraph_lib,
            n=dim(Q_tmp)[1], debug=debug,
            prec_graph_i = as.integer(i_),
            prec_graph_j = as.integer(j_),
            Tc = Tc,
            El = El,
            upper_edges = as.integer(upper.edges),
            lower_edges = as.integer(lower.edges),
            upper_edges_len = upper_edges_len,
            lower_edges_len = lower_edges_len,
            start_theta = start_theta,
            start_lsigma = start_lsigma,
            prior_theta_meanlog = prior_theta$meanlog,
            prior_theta_sdlog = prior_theta$sdlog,
            prior_sigma_meanlog = prior_sigma$meanlog,
            prior_sigma_sdlog = prior_sigma$sdlog,
            parameterization = parameterization))
}
model$graph_spde <- graph_spde
model$directional <- directional
model$data_PtE <- suppressWarnings(graph_object$get_PtE())
model$original_data <- graph_object$.__enclos_env__$private$data
model$parameterization <- parameterization
model$Tc <- Tc
model$alpha <- alpha
if(alpha == 2){
    A_tmp <- t(Tc)
    index.obs1 <- sapply(graph_spde$PtV, function(i){idx_temp <- i == graph_spde$E[,1]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
    index.obs1 <- (index.obs1-1)*4+1
    index.obs2 <- NULL
    na_obs1 <- is.na(index.obs1)
    if(any(na_obs1)){
          idx_na <- which(na_obs1)
          PtV_NA <- graph_spde$PtV[idx_na]
          index.obs2 <- sapply(PtV_NA, function(i){idx_temp <- i == graph_spde$E[,2]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
          index.obs1[na_obs1] <- (index.obs2-1)*4 + 3                                                                      
          }
    A_tmp <- A_tmp[index.obs1,] #A matrix for alpha=    
} else if(directional){
    A_tmp <- t(Tc)
    index.obs1 <- sapply(graph_spde$PtV, function(i){idx_temp <- i == graph_spde$E[,1]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
    index.obs1 <- (index.obs1-1)*2+1
    index.obs2 <- NULL
    na_obs1 <- is.na(index.obs1)
    if(any(na_obs1)){
          idx_na <- which(na_obs1)
          PtV_NA <- graph_spde$PtV[idx_na]
          index.obs2 <- sapply(PtV_NA, function(i){idx_temp <- i == graph_spde$E[,2]
                                                                      idx_temp <- which(idx_temp)
                                                                      return(idx_temp[1])})
          index.obs1[na_obs1] <- (index.obs2-1)*2 + 2                                                                      
          }
    A_tmp <- A_tmp[index.obs1,] #A matrix for alpha=   
}
model$A <- A_tmp
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
#' @noRd
graph_spde_make_index <- function (name,
                                   graph_spde,
                                   n.group = 1,
                                   n.repl = 1,
                                   ...) {
    n.spde <- graph_spde$f$n
    name.group <- paste(name, ".group", sep = "")
    name.repl <- paste(name, ".repl", sep = "")
    out <- list()
    out[[name]] <- rep(rep(1:n.spde, times = n.group), times = n.repl)
    out[[name.group]] <- rep(rep(1:n.group, each = n.spde), times = n.repl)
    out[[name.repl]] <- rep(1:n.repl, each = n.spde * n.group)
    return(out)
}


#' Deprecated - Observation/prediction matrices for 'SPDE' models
#'
#' Constructs observation/prediction weight matrices
#' for metric graph models.
#'
#' @param graph_spde An `inla_metric_graph_spde` object built with the
#' `graph_spde()` function.
#' @param repl Which replicates? If there is no replicates, or to
#' use all replicates, one can set to `NULL`.
#' @param drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`.
#' @param drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`.
#' @return The observation matrix.
#' @export

graph_spde_basis <- function (graph_spde, repl = NULL, drop_na = FALSE, drop_all_na = TRUE) {
    lifecycle::deprecate_warn("1.3.0", "graph_spde_basis()", "graph_data_spde()")
    warning("This function only works for alpha = 1. Use graph_data_spde() function for alpha = 1 and alpha = 2.")
   return(graph_spde$graph_spde$.__enclos_env__$private$A(group = repl, drop_na = drop_na, drop_all_na = drop_all_na))
}

#' Deprecated - Observation/prediction matrices for 'SPDE' models
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

graph_spde_make_A <- function(graph_spde, repl = NULL){
  lifecycle::deprecate_warn("1.2.0", "graph_spde_make_A()", "graph_spde_basis()")
    warning("This function only works for alpha = 1. Use graph_data_spde() function for alpha = 1 and alpha = 2.")
  return(graph_spde_basis(graph_spde, repl = repl, drop_na = FALSE, drop_all_na = FALSE))
}


#' Data extraction for 'spde' models
#'
#' Extracts data from metric graphs to be used by 'INLA' and 'inlabru'.
#'
#' @param graph_spde An `inla_metric_graph_spde` object built with the
#' `graph_spde()` function.
#' @param name A character string with the base name of the effect.
#' @param repl Which replicates? If there is no replicates, one
#' can set `repl` to `NULL`. If one wants all replicates,
#' then one sets to `repl` to `.all`.
#' @param repl_col Column containing the replicates. If the replicate is the internal group variable, set the replicates
#' to ".group". If not replicates, set to `NULL`.
#' @param group Which groups? If there is no groups, one
#' can set `group` to `NULL`. If one wants all groups,
#' then one sets to `group` to `.all`.
#' @param group_col Which "column" of the data contains the group variable?
#' @param likelihood_col If only a single likelihood, this variable should be `NULL`. In case of multiple likelihoods, which column contains the variable indicating the number of the likelihood to be considered?
#' @param resp_col If only a single likelihood, this variable should be `NULL`. In case of multiple likelihoods, column containing the response variable.
#' @param covariates Vector containing the column names of the covariates. If no covariates, then it should be `NULL`.
#' @param only_pred Should only return the `data.frame` to the prediction data?
#' @param loc `r lifecycle::badge("deprecated")` Use `loc_name` instead.
#' @param loc_name Character with the name of the location variable to be used in
#' 'inlabru' prediction.
#' @param tibble Should the data be returned as a `tidyr::tibble`?
#' @param drop_na Should the rows with at least one NA for one of the columns be removed? DEFAULT is `FALSE`. This option is turned to `FALSE` if `only_pred` is `TRUE`.
#' @param drop_all_na Should the rows with all variables being NA be removed? DEFAULT is `TRUE`. This option is turned to `FALSE` if `only_pred` is `TRUE`.
#' @return An 'INLA' and 'inlabru' friendly list with the data.
#' @export

graph_data_spde <- function (graph_spde, name = "field", repl = NULL, repl_col = NULL, group = NULL, 
                                group_col = NULL,
                                likelihood_col = NULL,
                                resp_col = NULL,
                                covariates = NULL,
                                only_pred = FALSE,
                                loc_name = NULL,
                                tibble = FALSE,
                                drop_na = FALSE, drop_all_na = TRUE,
                                loc = deprecated()){

        if (lifecycle::is_present(loc)) {
         if (is.null(loc_name)) {
           lifecycle::deprecate_warn("1.2.0", "graph_data_spde(loc)", "graph_data_spde(loc_name)",
             details = c("`loc` was provided but not `loc_name`. Setting `loc_name <- loc`.")
           )
           loc_name <- loc
         } else {
           lifecycle::deprecate_warn("1.2.0", "graph_data_spde(loc)", "graph_data_spde(loc_name)",
             details = c("Both `loc_name` and `loc` were provided. Only `loc_name` will be considered.")
           )
         }
         loc <- NULL
       }  

  
  if(!is.null(likelihood_col)){
    # Processing likelihoods
    if(any(is.na(graph_spde$graph_spde$.__enclos_env__$private$data[[likelihood_col]]))){
      tmp_group_col_val <- graph_spde$graph_spde$.__enclos_env__$private$data[[".group"]]
      tmp_unique_group_val <- unique(tmp_group_col_val)
      for(tmp_i in tmp_unique_group_val){
        idx_tmp <- (tmp_group_col_val == tmp_i)
        idx_like_val_temp <- !is.na(graph_spde$graph_spde$.__enclos_env__$private$data[[likelihood_col]][idx_tmp])
        like_val_temp <- unique(graph_spde$graph_spde$.__enclos_env__$private$data[[likelihood_col]][idx_tmp][idx_like_val_temp])
        if(length(like_val_temp)>1){
          stop("Likelihood processing error. There was something wrong when grouping the likelihood data.")
        }
        graph_spde$graph_spde$.__enclos_env__$private$data[[likelihood_col]][idx_tmp] <- like_val_temp
      }
    }

    like_val <- unique(graph_spde$graph_spde$.__enclos_env__$private$data[[likelihood_col]])
    if(is.null(resp_col)){
      stop("If likelihood_col is non-NULL, then resp_col should be non-NULL!")
    }
  } else{
    like_val <- 1
  }

  if(is.null(repl_col)){
    graph_spde$graph_spde$.__enclos_env__$private$data[[".dummy_repl_col"]] <- rep(1,length(graph_spde$graph_spde$.__enclos_env__$private$data[[".group"]]))
    repl_col <- ".dummy_repl_col"
  }

    if(is.null(group_col)){
    graph_spde$graph_spde$.__enclos_env__$private$data[[".dummy_group_col"]] <- rep(1,length(graph_spde$graph_spde$.__enclos_env__$private$data[[".group"]]))
    group_col <- ".dummy_group_col"
  }

  alpha <- graph_spde$alpha

  ret_list <- list()

  for(lik in like_val){
    ret <- list()
    
    graph_tmp <- graph_spde$graph_spde$clone()

    if(is.null(repl) || repl[1] == ".all") {
      groups <- graph_tmp$.__enclos_env__$private$data[[repl_col]]
      repl <- unique(groups)
    }     

    lik_tmp <- lik

    if(length(like_val) == 1){
      lik_tmp <- NULL
    } 

    graph_tmp$.__enclos_env__$private$data <- select_repl_group(graph_tmp$.__enclos_env__$private$data, repl = repl, repl_col = repl_col, group = lik_tmp, group_col = likelihood_col)   

    if(is.null((graph_tmp$.__enclos_env__$private$data))){
      stop("The graph has no data!")
    }
    if(only_pred){
      idx_anyNA <- !idx_not_any_NA(graph_tmp$.__enclos_env__$private$data)
      graph_tmp$.__enclos_env__$private$data <- lapply(graph_tmp$.__enclos_env__$private$data, function(dat){return(dat[idx_anyNA])})
      drop_na <- FALSE
      drop_all_na <- FALSE
    }
  
     ret[["data"]] <- select_repl_group(graph_tmp$.__enclos_env__$private$data, repl = repl, repl_col = repl_col,  group = group, group_col = group_col)   

    n.repl <- length(unique(repl))
  
    if(is.null(group)){
      if(!is.null(group_col)){
        n.group <- length(unique(graph_tmp$.__enclos_env__$private$data[[group_col]]))
        group <- unique(graph_tmp$.__enclos_env__$private$data[[group_col]])
      } else{
        n.group <- 1
      }
    } else if (group[1] == ".all"){
      n.group <- length(unique(graph_tmp$.__enclos_env__$private$data[[group_col]]))
      group <- unique(graph_tmp$.__enclos_env__$private$data[[group_col]])
    } else{
      n.group <- length(unique(group))
    }
  
    A <- Matrix::Diagonal(0)  

    for(i in 1:n.repl){
     for(j in 1:n.group){
         data_group_repl <- select_repl_group(ret[["data"]], repl = repl[i], repl_col = repl_col, group = group[j], group_col = group_col)
         if(drop_na){
           idx_notna <- idx_not_any_NA(data_group_repl)
         } else if(drop_all_na){
           idx_notna <- idx_not_all_NA(data_group_repl)
         } else{
           idx_notna <- rep(TRUE, length(data_group_repl[[repl_col]]))
         }
         # nV_tmp <- sum(idx_notna)        
         if(alpha == 1){
          if(!graph_spde$directional){
            A_tmp <- Matrix::Diagonal(graph_tmp$nV)[graph_tmp$PtV[idx_notna], ]
          } else{
           A_tmp <- graph_spde$A 
           A_tmp <- A_tmp[idx_notna,]
          }
         } else{
           A_tmp <- graph_spde$A 
           A_tmp <- A_tmp[idx_notna,]
         }
         A <- Matrix::bdiag(A, A_tmp)
     }
    }
   
    if(tibble){
      ret[["data"]] <-tidyr::as_tibble(ret[["data"]])
    }
  
    if(drop_all_na){
      is_tbl <- inherits(ret, "tbl_df")
        idx_temp <- idx_not_all_NA(ret[["data"]])
        ret[["data"]] <- lapply(ret[["data"]], function(dat){dat[idx_temp]}) 
        if(is_tbl){
          ret[["data"]] <- tidyr::as_tibble(ret[["data"]])
        }
    }    
    if(drop_na){
      if(!inherits(ret[["data"]], "tbl_df")){
        idx_temp <- idx_not_any_NA(ret[["data"]])
        ret[["data"]] <- lapply(ret[["data"]], function(dat){dat[idx_temp]})
      } else{
        ret[["data"]] <- tidyr::drop_na(ret[["data"]])
      }
    }
    
    if(!is.null(loc_name)){
        ret[["data"]][[loc_name]] <- cbind(ret[["data"]][[".edge_number"]],
                            ret[["data"]][[".distance_on_edge"]])
    }
  
    if(!inherits(ret[["data"]], "metric_graph_data")){
      class(ret[["data"]]) <- c("metric_graph_data", class(ret))
    }
  
    ret[["repl"]] <- bru_graph_rep(repl = repl, graph_spde = graph_spde, repl_col = repl_col)
  
     ret[["group"]] <- bru_graph_rep(repl = group, graph_spde = graph_spde, repl_col = group_col)
  
     ret[["index"]] <- graph_spde_make_index(name = name, graph_spde = graph_spde,
                                     n.group = n.group,
                                     n.repl = n.repl)
    
    if(!is.null(covariates)){
      cov_tmp <- list()
      for(cov_var in covariates){
        cov_tmp[[cov_var]] <- ret[["data"]][[cov_var]]
        if(is.null(loc_name)){
          ret[["data"]][[cov_var]] <- NULL
        }
      }
      ret[["index"]] <- list(ret[["index"]], cov_tmp)
      ret[["basis"]] <- list(A, 1)
    } else{
      ret[["basis"]] <- A
    }
        

    ret_list[[lik]] <- ret
  }
  
  if(is.null(likelihood_col)){
    return(ret)
  } else if(is.null(loc_name)){
    count <- 1
    for(lik in like_val){
      tmp_data <- matrix(nrow = length(ret_list[[lik]][["data"]][[resp_col]]), ncol = length(like_val))
      tmp_data[,count] <- ret_list[[lik]][["data"]][[resp_col]]
      ret_list[[lik]][["data"]][[resp_col]] <- tmp_data
      count <- count + 1
    }
  }
  
  return(ret_list)
}


#' Select replicate and group
#' @noRd
#'
select_repl_group <- function(data_list, repl, repl_col, group, group_col){
    if(!is.null(group) && is.null(group_col)){
      stop("If you specify group, you need to specify group_col!")
    }
    if(!is.null(group)){
      grp <- data_list[[group_col]]
      grp <- which(grp %in% group)
      data_result <- lapply(data_list, function(dat){dat[grp]})
      replicates <- data_result[[repl_col]]
      replicates <- which(replicates %in% repl)
      data_result <- lapply(data_result, function(dat){dat[replicates]})
      return(data_result)
    } else{
      replicates <- data_list[[repl_col]]
      replicates <- which(replicates %in% repl)
      data_result <- lapply(data_list, function(dat){dat[replicates]})
      return(data_result)
    }
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

  alpha <- metric_graph_spde$alpha
  nu <- alpha - 0.5

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
              result[[paste0("marginals.log.",name_theta1)]],
              function(x) {
                INLA::inla.tmarginal(
                  function(y) exp(-y),
                  x
                )
              }
            )
            names(result[[paste0("marginals.",name_theta1_t)]]) <- name_theta1_t
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
            kappa_est <- sqrt(8*nu)/range_est
            sigma_est <- sqrt(gamma(nu) / (tau_est^2 * kappa_est^(2 * nu) *
                    (4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))

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
                  stop.on.error = FALSE,
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
  if(model$alpha == 1 && !(model$directional)){
    pte_tmp <- model$graph_spde$get_PtE()
    input_list <- lapply(1:nrow(input), function(i){input[i,]})
    pte_tmp_list <- lapply(1:nrow(pte_tmp), function(i){pte_tmp[i,]})
    idx_tmp <- match(input_list, pte_tmp_list)
    A_tmp <- model$graph_spde$.__enclos_env__$private$A()
    return(A_tmp[idx_tmp, , drop=FALSE])
  } else{
    pte_tmp <- model$graph_spde$get_PtE()
    input_list <- lapply(1:nrow(input), function(i){input[i,]})
    pte_tmp_list <- lapply(1:nrow(pte_tmp), function(i){pte_tmp[i,]})
    idx_tmp <- match(input_list, pte_tmp_list)
    A_tmp <- model$A
    return(A_tmp[idx_tmp, , drop=FALSE])    
  }
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
                  stop.on.error = FALSE,
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
                  stop.on.error = FALSE,
    subdivisions = nrow(density_df)
  )$value

  sd_temp <- sqrt(stats::integrate(
    f = function(z) {
      denstemp(z) * (z - mean_temp)^2
    }, lower = min_x, upper = max_x,
                  stop.on.error = FALSE,
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
#' @param repl A vector of replicates. If set to `.all`, a vector
#' for all replicates will be generated.
#' @param graph_spde Name of the field.
#' @return A vector of replicates to be used with 'inlabru'.
#' @noRd

bru_graph_rep <- function(repl, graph_spde, repl_col){
  groups <- unique(graph_spde$graph_spde$.__enclos_env__$private$data[[repl_col]])
  if(repl[1] == ".all"){
    repl <- groups
  }
  n_groups <- length(groups)
  length_resp <- sum(graph_spde$graph_spde$.__enclos_env__$private$data[[repl_col]] == groups[1])
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
#' @param newdata A data.frame of covariates needed for the prediction. The
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
#' @param repl Which replicates? If there is no replicates, one
#' can set `repl` to `NULL`. If one wants all replicates,
#' then one sets to `repl` to `.all`.
#' @param repl_col Column containing the replicates. If the replicate is the internal group variable, set the replicates
#' to ".group". If not replicates, set to `NULL`.
#' @param group Which groups? If there is no groups, one
#' can set `group` to `NULL`. If one wants all groups,
#' then one sets to `group` to `.all`.
#' @param group_col Which "column" of the data contains the group variable?
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
#' @param tolerance_merge Tolerance for merging prediction points into original points to increase stability.
#' @param... Additional arguments passed on to `inla.posterior.sample()`.
#' @param data `r lifecycle::badge("deprecated")` Use `newdata` instead.
#' @return A list with predictions.
#' @export

predict.inla_metric_graph_spde <- function(object,
                                           cmp,
                                           bru_fit,
                                           newdata = NULL,
                                           formula = NULL,
                                           data_coords = c("PtE", "euclidean"),
                                           normalized = TRUE,
                                           repl = NULL,
                                           repl_col = NULL,
                                           group = NULL,
                                           group_col = NULL,
                                           n.samples = 100,
                                           seed = 0L,
                                           probs = c(0.025, 0.5, 0.975),
                                           return_original_order = TRUE,
                                           num.threads = NULL,
                                           include = NULL,
                                           exclude = NULL,
                                           drop = FALSE,
                                           tolerance_merge = 1e-5,
                                           ...,
                                           data = deprecated()){
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
  data_coords <- data_coords[[1]]
  if(!(data_coords %in% c("PtE", "euclidean"))){
    stop("data_coords must be either 'PtE' or 'euclidean'!")
  }
  graph_tmp <- object$graph_spde$get_initial_graph()
  graph_tmp$clear_observations()
  # graph_tmp <- object$graph_spde$clone()
  name_locations <- bru_fit$bru_info$model$effects$field$main$input$input
  original_data <- object$original_data

  group_variables <- attr(object$graph_spde$.__enclos_env__$private$data, "group_variable")

  if(group_variables == ".none"){
  graph_tmp$add_observations(data = original_data,
                  edge_number = ".edge_number",
                  distance_on_edge = ".distance_on_edge",
                  data_coords = "PtE",
                  normalized = TRUE,
                  verbose=0,
                  suppress_warnings = TRUE)
  } else{
      graph_tmp$add_observations(data = original_data,
                  edge_number = ".edge_number",
                  distance_on_edge = ".distance_on_edge",
                  data_coords = "PtE",
                  normalized = TRUE,
                  verbose=0,
                  group = group_variables,
                  suppress_warnings = TRUE)
  }

  new_data <- data
  new_data[[name_locations]] <- NULL
  n_locations <- nrow(data[[name_locations]])
  names_columns <- names(original_data)
  names_columns <- setdiff(names_columns, c(".group", ".coord_x",
                                            ".coord_y", ".edge_number",
                                            ".distance_on_edge"))

  # for(name_column in names_columns){
  #   new_data[[name_column]] <- rep(NA, n_locations)
  # }
  if(data_coords == "PtE"){
    new_data[[".edge_number"]] <- data[[name_locations]][,1]
    new_data[[".distance_on_edge"]] <- data[[name_locations]][,2]
  } else{
    new_data[[".coord_x"]] <- data[[name_locations]][,1]
    new_data[[".coord_y"]] <- data[[name_locations]][,2]
  }

  new_data[["__dummy_var"]] <- 1:length(new_data[[".edge_number"]])

  if(group_variables == ".none"){
    group_variables <- NULL
  } else if(group_variables == ".group"){
    if(!(".group" %in% names(new_data))){
      if(length(unique(original_data[[".group"]])) == 1){
        new_data[[".group"]] <- rep(unique(original_data[[".group"]]), length(new_data[[".edge_number"]]))
      }
    }
  } 

  if(!is.null(group_variables)){
    if(!(group_variables %in% names(new_data))){
      warning("The replicate variable was not found in newdata. Predictions were only given for the first replicate.")
      new_data[[".group"]] <- rep(original_data[[group_variables]][1], length(new_data[[".edge_number"]]))
    }
  }


  graph_tmp$add_observations(data = new_data,
                  edge_number = ".edge_number",
                  distance_on_edge = ".distance_on_edge",
                  coord_x = ".coord_x",
                  coord_y = ".coord_y",
                  data_coords = data_coords,
                  normalized = normalized,
                  group = group_variables,
                  verbose=0,
                  suppress_warnings = TRUE,
                  tolerance_merge = tolerance_merge)

  dummy1 <- graph_tmp$.__enclos_env__$private$data[["__dummy_var"]]

  graph_tmp$.__enclos_env__$private$data[["__dummy_var2"]] <- 1:length(graph_tmp$.__enclos_env__$private$data[["__dummy_var"]])

  pred_PtE <- cbind(graph_tmp$.__enclos_env__$private$data[[".edge_number"]],
                          graph_tmp$.__enclos_env__$private$data[[".distance_on_edge"]])

  # pred_PtE <- pred_PtE[!is.na(dummy1),]

  # Adding the original data

  # graph_tmp$add_observations(data = original_data,
  #                   coord_x = ".coord_x",
  #                   coord_y = ".coord_y",
  #                   data_coords = "euclidean", verbose=0)

  graph_tmp$observation_to_vertex(mesh_warning=FALSE)

  # tmp_list2 <- cbind(graph_tmp$data[[".coord_x"]],
  #                                       graph_tmp$data[[".coord_y"]])
  # tmp_list2 <- lapply(1:nrow(tmp_list2), function(i){tmp_list2[i,]})
  # idx_list <- match(tmp_list, tmp_list2)

  new_data_list <- graph_tmp$.__enclos_env__$private$data

  idx_list <- !is.na(new_data_list[["__dummy_var"]])

  new_data_list <- lapply(new_data_list, function(dat){dat[idx_list]})

  pred_PtE <- pred_PtE[graph_tmp$.__enclos_env__$private$data[["__dummy_var2"]],][idx_list,]

  # new_data_list[[name_locations]] <- cbind(graph_tmp$data[[".edge_number"]][idx_list],
  #                                             graph_tmp$data[[".distance_on_edge"]][idx_list])

  new_data_list[[name_locations]] <- cbind(new_data_list[[".edge_number"]],
                                              new_data_list[[".distance_on_edge"]])                                      

  spde____model <- graph_spde(graph_tmp, alpha = object$alpha, directional = object$directional)

  cmp_c <- as.character(cmp)
  name_model <- deparse(substitute(object))
  cmp_c[3] <- sub(name_model, "spde____model", cmp_c[3])
  cmp <- as.formula(paste(cmp_c[2], cmp_c[1], cmp_c[3]))

  new_data_tmp <- graph_data_spde(spde____model, loc_name = name_locations, 
                        drop_all_na = FALSE, drop_na = FALSE, group = group, group_col = group_col,
                        repl_col = repl_col, repl = repl)[["data"]]


  info <- bru_fit[["bru_info"]]
  info[["options"]] <- inlabru::bru_call_options(inlabru::bru_options(info[["options"]]))

  bru_fit_new <- inlabru::bru(cmp,
          data = new_data_tmp, options = info[["options"]])
  
  pred <- predict(object = bru_fit_new,
                    newdata = new_data_list,
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
    ord <- graph_tmp$.__enclos_env__$private$data[["__dummy_var"]][idx_list]
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
#' @param ... Additional parameters to be passed to plot_function.
#' @return A 'ggplot2' object.
#' @export

plot.graph_bru_pred <- function(x, y = NULL, vertex_size = 0, ...){
  m_prd_bru <- x$pred$mean
  PtE_prd <- x$PtE_pred
  newdata <- data.frame("edge_number" = PtE_prd[,1],
                        "distance_on_edge" = PtE_prd[,2],
                        "pred_y" = m_prd_bru)
  newdata <- x$initial_graph$process_data(data = newdata, normalized = TRUE)
  
  p <- x$initial_graph$plot_function(data = "pred_y", newdata=newdata, vertex_size = vertex_size,...)
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
#' @param newdata A data.frame of covariates needed for the prediction. The locations
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
#' @param data `r lifecycle::badge("deprecated")` Use `newdata` instead.
#' @return A list with predictions.
#' @export

predict.rspde_metric_graph <- function(object,
                                           cmp,
                                           bru_fit,
                                           newdata = NULL,
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
                                           ...,
                                           data = deprecated()){
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
  data_coords <- data_coords[[1]]
  if(!(data_coords %in% c("PtE", "euclidean"))){
    stop("data_coords must be either 'PtE' or 'euclidean'!")
  }
  graph_tmp <- object$mesh$get_initial_graph()
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
                    newdata = newdata,
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


#' @name process_rspde_predictions
#' @title Process predictions of `rspde_metric_graph` objects obtained by using `inlabru`
#' @description Auxiliar function to transform the predictions of the field into a plot friendly object.
#' @param pred The predictions of the field obtained by using `inlabru`
#' @param graph The original `metric_graph` object in which the predictions were obtained.
#' @param PtE Normalized locations of the points on the edge.
#' @return A list with predictions.
#' @export

process_rspde_predictions <- function(pred,
                                        graph,
                                        PtE = NULL){
  pred_list <- list()
  pred_list[["pred"]] <- pred
  pred_list[["PtE_pred"]] <- PtE
  pred_list[["graph"]] <- graph

  class(pred_list) <- "graph_bru_proc_pred"
  return(pred_list)
}


#' @name plot.graph_bru_proc_pred
#' @title Plot of processed predicted values with 'inlabru'
#' @description Auxiliary function to obtain plots of the processed predictions of the field
#' using 'inlabru'.
#' @param x A processed predicted object obtained with the `process_rspde_predictions` function.
#' @param y Not used.
#' @param vertex_size Size of the vertices.
#' @param ... Additional parameters to be passed to plot_function.
#' @return A 'ggplot2' object.
#' @export

plot.graph_bru_proc_pred <- function(x, y = NULL, vertex_size = 0, ...){
  m_prd_bru <- x$pred$mean
  PtE_prd <- x$PtE_pred
  graph <- x$graph
  newdata <- data.frame("edge_number" = PtE_prd[,1],
                        "distance_on_edge" = PtE_prd[,2],
                        "pred_y" = m_prd_bru)
  newdata <- graph$process_data(data = newdata, normalized = TRUE)
  
  p <- graph$plot_function(data = "pred_y", newdata=newdata, vertex_size = vertex_size,...)
  p
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
                                        if(inherits(data, "metric_graph_data")){
                                          edge_number <- ".edge_number"
                                          distance_on_edge = ".distance_on_edge"
                                        }
                                        if(is.null(data[[edge_number]])){
                                          stop(paste("No column",edge_number,"was found in data."))
                                        }
                                        if(is.null(data[[distance_on_edge]])){
                                          stop(paste("No column",distance_on_edge,"was found in data."))
                                        }                                        
                                        data[[loc]] <- cbind(data[[edge_number]], data[[distance_on_edge]])
                                        data[[edge_number]] <- NULL
                                        data[[distance_on_edge]] <- NULL
                                        return(data)
                                        }
