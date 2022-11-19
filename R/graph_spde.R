

gpgraph_spde <- function(graph_object, alpha = 1, stationary_endpoints = "all",
 start_kappa = NULL, start_sigma = NULL, prior_kappa = NULL,
 prior_sigma = NULL, debug = FALSE){

  nu <- alpha - 0.5
  V <- graph_object$V
  EtV <- graph_object$E
  El <- graph_object$edge_lengths 

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

  gpgraph_lib <- system.file('shared', package='GPGraph')

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

  if(is.null(prior_kappa$meanlog)){
      if(is.null(graph_object$geo.dist)){
        graph_object$compute_geodist()
      }
      prior.range.nominal <- max(graph_object$geo.dist) * 0.2
      prior_kappa$meanlog <- log(sqrt(8 *
      exp(nu) / prior.range.nominal))
    }
  
  if (is.null(prior_kappa$sdlog)) {
    prior_kappa$sdlog <- sqrt(10)
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

  model <- do.call(
        'inla.cgeneric.define',
        list(model="inla_cgeneric_gpgraph_alpha1_model",
            shlib=paste0(gpgraph_lib, '/gpgraph_cgeneric_models.so'),
            n=as.integer(n.v), debug=debug,
            prec_graph_i = as.integer(i_),
            prec_graph_j = as.integer(j_),
            index_graph = as.integer(idx_ij),
            count_idx = as.integer(count_idx),
            EtV2 = EtV2,
            EtV3 = EtV3,
            El = El,
            stationary_endpoints = as.integer(index),
            start_lkappa = start_lkappa,
            start_lsigma = start_lsigma,
            prior_kappa_meanlog = prior_kappa$meanlog,
            prior_kappa_sdlog = prior_kappa$sdlog,
            prior_sigma_meanlog = prior_sigma$meanlog,
            prior_sigma_sdlog = prior_sigma$sdlog))
model$graph_obj <- graph_object
class(model) <- c("inla_metric_graph_spde", class(model))
return(model)
}


graph_spde_make_index <- function (name, graph, n.group = 1, n.repl = 1, ...) {
   
    n.spde <- dim(graph$V)[1]
    name.group <- paste(name, ".group", sep = "")
    name.repl <- paste(name, ".repl", sep = "")
    out <- list()
    out[[name]] <- rep(rep(1:n.spde, times = n.group), times = n.repl)
    out[[name.group]] <- rep(rep(1:n.group, each = n.spde), times = n.repl)
    out[[name.repl]] <- rep(1:n.repl, each = n.spde * n.group)
    return(out)
}


#'
#' @title metric graph inlabru mapper
#' @name bru_mapper.inla_metric_graph_spde
#' @param model An `inla_metric_graph_spde` for which to construct or extract a mapper
#' @param \dots Arguments passed on to other methods
#' @rdname bru_mapper.inla_metric_graph_spde
#' @rawNamespace if (getRversion() >= "3.6.0") {
#'   S3method(inlabru::bru_mapper, inla_metric_graph_spde)
#'   # S3method(inlabru::bru_get_mapper, inla_metric_graph_spde)
#'   S3method(inlabru::ibm_n, bru_mapper_inla_metric_graph_spde)
#'   S3method(inlabru::ibm_values, bru_mapper_inla_metric_graph_spde)
#'   S3method(inlabru::ibm_jacobian, bru_mapper_inla_metric_graph_spde)
#' }
#' 
bru_mapper.inla_metric_graph_spde <- function(model,...) {
  mapper <- list(model = model)
  # Note 1: From inlabru > 2.5.3, use bru_mapper_define instead.
  # Note 2: bru_mapper.default is not exported from inlabru, so
  # must call the generic bru_mapper()
  inlabru::bru_mapper(mapper, new_class = "bru_mapper_inla_metric_graph_spde")
}

#' @param mapper A `bru_mapper.inla_metric_graph_spde` object
#' @rdname bru_mapper.inla_metric_graph_spde
ibm_n.bru_mapper_inla_metric_graph_spde <- function(mapper, ...) {
  model <- mapper[["model"]]
  model$f$n
}
#' @rdname bru_mapper.inla_rspde
ibm_values.bru_mapper_inla_metric_graph_spde <- function(mapper, ...) {
  seq_len(inlabru::ibm_n(mapper))
}
#' @param input The values for which to produce a mapping matrix
#' @rdname bru_mapper.inla_rspde
ibm_jacobian.bru_mapper_inla_metric_graph_spde <- function(mapper, input, ...) {
  if (is.null(input)) {
    return(Matrix::Matrix(0, 0, inlabru::ibm_n(mapper)))
  }
  if (!is.matrix(input) && !inherits(input, "Spatial")) {
    input <- as.matrix(input)
  }
  model <- mapper[["model"]]
  if(is.null(model$graph_obj$A)){
    model$graph_obj$observation_to_vertex()
  }
  model$graph_obj$A
}

#' @rdname bru_mapper.inla_metric_graph_spde
bru_get_mapper.inla_metric_graph_spde <- function(model, ...){
 inlabru::bru_mapper(model)
}
