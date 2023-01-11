
#' INLA implementation of Whittle-Matérn fields for metric graphs.
#'
#' This function creates an inla object that can be used
#' in `INLA` or `inlabru` to fit Whittle-Matérn fields
#' on metric graphs.
#'
#' @param graph_object A `metric_graph` object.
#' @param alpha The order of the SPDE.
#' @param stationary_endpoints Which vertices of degree 1 should contain stationary boundary conditions?
#' @param parameterization Which parameterization to be used? The options are 'matern' (sigma and range) and 'spde' (sigma and kappa).
#' @param start_range Starting value for range parameter.
#' @param prior_range a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation on the log scale. Will not be used if prior.kappa is non-null.
#' @param start_kappa Starting value for kappa.
#' @param start_sigma Starting value for sigma.
#' @param prior_kappa a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation on the log scale.
#' @param prior_sigma a `list` containing the elements `meanlog` and
#' `sdlog`, that is, the mean and standard deviation on the log scale.
#' @param debug Should debug be displayed?
#'
#' @return An inla object.
#' @export
graph_spde <- function(graph_object, alpha = 1, stationary_endpoints = "all",
 parameterization = c("matern", "spde"),
 start_range = NULL, prior_range = NULL,
 start_kappa = NULL, start_sigma = NULL,
 prior_kappa = NULL,
 prior_sigma = NULL, debug = FALSE){

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

  gpgraph_lib <- system.file('shared', package='MetricGraph')
  
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
            start_theta = start_theta,
            start_lsigma = start_lsigma,
            prior_theta_meanlog = prior_theta$meanlog,
            prior_theta_sdlog = prior_theta$sdlog,
            prior_sigma_meanlog = prior_sigma$meanlog,
            prior_sigma_sdlog = prior_sigma$sdlog,
            parameterization = parameterization))
model$graph_spde <- graph_spde
model$parameterization <- parameterization
class(model) <- c("inla_metric_graph_spde", class(model))
return(model)
}



#'  Model index vector generation for metric graph models
#'
#' Generates a list of named index vectors for
#' `INLA`-based metric graph models.
#'
#' @param name A character string with the base name of the effect.
#' @param graph_spde An `inla_metric_graph_spde` object built with the `graph_spde()` function.
#' @param n.group Number of groups.
#' @param n.repl Number of replicates.
#' @param ... Currently not being used.
#'
#' @return A list of indexes.
#' @export
graph_spde_make_index <- function (name, graph_spde, n.group = 1, n.repl = 1, ...) {
    n.spde <- dim(graph_spde$graph_spde$V)[1]
    name.group <- paste(name, ".group", sep = "")
    name.repl <- paste(name, ".repl", sep = "")
    out <- list()
    out[[name]] <- rep(rep(1:n.spde, times = n.group), times = n.repl)
    out[[name.group]] <- rep(rep(1:n.group, each = n.spde), times = n.repl)
    out[[name.repl]] <- rep(1:n.repl, each = n.spde * n.group)
    return(out)
}


#' Observation/prediction matrices for rSPDE models.
#'
#' Constructs observation/prediction weight matrices
#' for metric graph models.
#'
#' @param graph_spde An `inla_metric_graph_spde` object built with the `graph_spde()` function.
#' @param repl Which replicates? If there is no replicates, or to
#' use all replicates, one can set to `NULL`.
#' @return The observation matrix
#' @export

graph_spde_make_A <- function (graph_spde, repl = NULL) {
   return(graph_spde$graph_spde$A(group = repl))
}


#' Observation/prediction matrices for rSPDE models.
#'
#' Constructs observation/prediction weight matrices
#' for metric graph models.
#'
#' @param graph_spde An `inla_metric_graph_spde` object built with the `graph_spde()` function or
#' an `rspde_metric_graph` object built with the `rspde.metric_graph()` function from the `rSPDE` package.
#' @param repl Which replicates? If there is no replicates, one
#' can set `repl` to `NULL`. If one wants all replicates,
#' then one sets to `repl` to `__all`.
#' @param only_pred Should only return the `data.frame` to the prediction data?
#' @param loc Character with the name of the location variable to be used in `inlabru`'s prediction.
#' @return The observation matrix
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
    ret[[loc]] <- cbind(graph_tmp$data[["__edge_number"]],
                          graph_tmp$data[["__distance_on_edge"]])
  }
  return(ret)
}



#' @name spde_metric_graph_result
#' @title metric graph SPDE result extraction from INLA estimation results
#' @description Extract field and parameter values and distributions
#' for a metric graph spde effect from an inla result object.
#' @param inla An `inla` object obtained from a call to
#' `inla()`.
#' @param name A character string with the name of the rSPDE effect
#' in the inla formula.
#' @param metric_graph_spde The `inla_metric_graph_spde` object used for the effect in
#' the inla formula.
#' @param compute.summary Should the summary be computed?
#' @return If the model was fitted with `matern` parameterization (the default), it returns a list containing:
#' \item{marginals.range}{Marginal densities for the range parameter}
#' \item{marginals.log.range}{Marginal densities for log(range)}
#' \item{marginals.sigma}{Marginal densities for std. deviation}
#' \item{marginals.log.sigma}{Marginal densities for log(std. deviation)}
#' \item{marginals.values}{Marginal densities for the field values}
#' \item{summary.log.range}{Summary statistics for log(range)}
#' \item{summary.log.sigma}{Summary statistics for log(std. deviation)}
#' \item{summary.values}{Summary statistics for the field values}
#' If `compute.summary` is `TRUE`, then the list will also contain
#' \item{summary.kappa}{Summary statistics for kappa}
#' \item{summary.tau}{Summary statistics for tau}
#' If the model was fitted with the `spde` parameterization, it returns a list containing:
#' \item{marginals.kappa}{Marginal densities for kappa}
#' \item{marginals.log.kappa}{Marginal densities for log(kappa)}
#' \item{marginals.log.tau}{Marginal densities for log(tau)}
#' \item{marginals.tau}{Marginal densities for tau}
#' \item{marginals.values}{Marginal densities for the field values}
#' \item{summary.log.kappa}{Summary statistics for log(kappa)}
#' \item{summary.log.tau}{Summary statistics for log(tau)}
#' \item{summary.values}{Summary statistics for the field values}
#' If `compute.summary` is `TRUE`, then the list will also contain
#' \item{summary.kappa}{Summary statistics for kappa}
#' \item{summary.tau}{Summary statistics for tau}
#' @export

spde_metric_graph_result <- function(inla, name, metric_graph_spde, compute.summary = TRUE) {
  if(!inherits(metric_graph_spde, "inla_metric_graph_spde")){
    stop("You should provide an inla_metric_graph_spde object!")
  }

  result <- list()

  parameterization <- metric_graph_spde$parameterization

    if(parameterization == "spde"){
      row_names <- c("sigma", "kappa")
    } else{
      row_names <- c("sigma", "range")
    }

  result$summary.values <- inla$summary.random[[name]]

  if (!is.null(inla$marginals.random[[name]])) {
    result$marginals.values <- inla$marginals.random[[name]]
  }

  if(parameterization == "spde"){
    name_theta1 <- "sigma"
    name_theta2 <- "kappa"
  } else{
    name_theta1 <- "sigma"
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

    result[[paste0("marginals.",name_theta1)]] <- lapply(
      result[[paste0("marginals.log.",name_theta1)]],
      function(x) {
        INLA::inla.tmarginal(
          function(y) exp(y),
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

    norm_const_theta1 <- norm_const(result[[paste0("marginals.",name_theta1)]][[name_theta1]])
    result[[paste0("marginals.",name_theta1)]][[name_theta1]][, "y"] <-
    result[[paste0("marginals.",name_theta1)]][[name_theta1]][, "y"] / norm_const_theta1

    norm_const_theta2 <- norm_const(result[[paste0("marginals.",name_theta2)]][[name_theta2]])
    result[[paste0("marginals.",name_theta2)]][[name_theta2]][, "y"] <-
    result[[paste0("marginals.",name_theta2)]][[name_theta2]][, "y"] / norm_const_theta2




    result[[paste0("summary.",name_theta1)]] <- create_summary_from_density(result[[paste0("marginals.",name_theta1)]][[name_theta1]],
    name = name_theta1)
    result[[paste0("summary.",name_theta2)]] <-
    create_summary_from_density(result[[paste0("marginals.",name_theta2)]][[name_theta2]], name = name_theta2)
  }

  class(result) <- "metric_graph_spde_result"
  result$params <- c(name_theta1,name_theta2)
  return(result)
}


#' Data frame for metric_graph_spde_result objects to be used in ggplot2
#'
#' Returns a ggplot-friendly data-frame with the marginal posterior densities.
#' @aliases gg_df gg_df.metric_graph_spde_result
#' @param result An metric_graph_spde_result object.
#' @param parameter Vector. Which parameters to get the posterior density in the data.frame? The options are `sigma`, `range` or `kappa`.
#' @param transform Should the posterior density be given in the original scale?
#' @param restrict_x_axis Variables to restrict the range of x axis based on quantiles.
#' @param restrict_quantiles List of quantiles to restrict x axis.
#' @param ... Not being used.
#'
#' @return A data frame containing the posterior densities.
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
#' Summary for posteriors of rSPDE field parameters in
#' their original scales.
#'
#' @param object A `rspde.result` object.
#' @param digits integer, used for number formatting with signif()
#' @param ... Currently not used.
#'
#' @return Returns a `data.frame`
#' containing the summary.
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
#' @title metric graph inlabru mapper
#' @name bru_mapper.inla_metric_graph_spde
#' @param model An `inla_metric_graph_spde` for which to construct or extract a mapper
#' @param \dots Arguments passed on to other methods
#' @rdname bru_mapper.inla_metric_graph_spde
#' @rawNamespace if (getRversion() >= "3.6.0") {
#'   S3method(inlabru::bru_mapper, inla_metric_graph_spde)
#'   S3method(inlabru::bru_get_mapper, inla_metric_graph_spde)
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
  if(is.null(input)){
    return(model$graph_spde$A())
  } else if(input[1] == "__all"){
    return(model$graph_spde$A(group="__all"))
  } else{
    A_list <- list()
    for(repl_ in unique(input)){
      graph_tmp <- model$graph_spde$get_initial_graph()
      data_tmp <- graph_data_spde(model, 
            repl=repl_)
      graph_tmp$add_observations(data = data_tmp,
                    coord_x = "__coord_x",
                    coord_y = "__coord_y",
                    data_coords = "euclidean")
      graph_tmp$observation_to_vertex()
      A_list <- c(A_list, graph_tmp$A(group=repl_))
    }
    A <- do.call(rbind, A_list)
    return(A)
  }
}

#' @rdname bru_mapper.inla_metric_graph_spde
bru_get_mapper.inla_metric_graph_spde <- function(model, ...){
 inlabru::bru_mapper(model)
}



#' @name create_summary_from_density
#' @title Creates a summary from a density data frame
#' @description Auxiliar function to create summaries from density data drames
#' @param density_df A density data frame
#' @param name Name of the parameter
#' @return A data frame containing a basic summary
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
#' @title Creates a vector of replicates to be used with inlabru
#' @description Auxiliar function to create a vector of replicates to be used with inlabru
#' @param repl A vector of replicates. If set to `__all`, a vector
#' for all replicates will be generated.
#' @param graph_spde Name of the field
#' @return A vector of replicates to be used with inlabru
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


#' @name inlabru_predict
#' @title Kriging with inlabru
#' @description Auxiliar function to obtain predictions of the field
#' using inlabru.
#' @param bru_model  An `inla_metric_graph_spde` object built with the `graph_spde()` function.
#' @param bru_fit A fitted model using `inlabru` or `inla`.
#' @param cmp An inlabru component.
#' @param XY Euclidean coordinates of the prediction locations.
#' @param PtE Relative positions on the edge to obtain predictions.
#' @param data_pred A data.frame containing the prediction locations along with the response variables (with NA if they are missing).
#' @return A list with predictions.
#' @export

inlabru_predict <- function(bru_model, bru_fit, cmp, XY = NULL, PtE = NULL,
                            data_pred = NULL){
  if(is.null(XY) && is.null(PtE) && is.null(data_pred)){
    stop("No location to predict was provided!")
  }
  graph_tmp <- bru_model$graph_spde$get_initial_graph()
  repl <- unique(bru_model$graph_spde$data[["__group"]])

  data_tmp <- graph_data_spde(bru_model, 
            repl=repl[1])
  graph_tmp$add_observations(data = data_tmp,
                    coord_x = "__coord_x",
                    coord_y = "__coord_y",
                    data_coords = "euclidean")

  if(is.null(data_pred)){
    if(!is.null(XY)){
      PtE <- graph_tmp$coordinates(XY = XY)
    }
    resp <- as.character(cmp[2])
    data_pred <- list()
    data_pred[[resp]] <- rep(NA, nrow(PtE))
    data_pred[["edge_number"]] <- PtE[,1]
    data_pred[["distance_on_edge"]] <- PtE[,2]
  } 


  graph_tmp$add_observations(data = data_pred)
  graph_tmp$observation_to_vertex()
  spde____model <- graph_spde(graph_tmp)
  cmp_c <- as.character(cmp)
  name_model <- deparse(substitute(bru_model))
  cmp_c[3] <- sub(name_model, "spde____model", cmp_c[3])
  cmp <- as.formula(paste(cmp_c[2], cmp_c[1], cmp_c[3]))
  bru_fit_new <- inlabru::bru(cmp, 
          data = graph_data_spde(spde____model))#,
          # options = list(
          #   control.mode = list(
          #     theta = bru_fit$mode$theta,
          #     fixed = TRUE
          #   )
          # ))
  fitted_values <- bru_fit_new$summary.fitted.values
  
idx_prd <- which(is.na(graph_data_spde(spde____model)[[resp]]))
return(fitted_values[idx_prd,])
}


