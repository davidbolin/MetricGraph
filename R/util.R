#' internal function for checking metric_graph inputs
#' @noRd
check_graph <- function(graph)
{
  if (!inherits(graph, "metric_graph")) {
    stop("The graph object is not a metric graph")
  }
  out <- list(has.mesh = FALSE,
              has.obs = FALSE)
  if(!is.null(graph$mesh)){
    out$has.mesh = TRUE
  }
  if(!is.null(graph$data))
    out$has.data = TRUE
  return(out)
}

#'
#' computes the covariance of free space (r'=0) neumann boundary
#' @param s (n x 1) location
#' @param t (n x 1) location
#' @param kappa a range-like parameter (from the SPDE)
#' @param sigma the standard deviation
#' @param nu the smoothness parameter
#' @param L interval length
#' @param deriv a vector containing the order of the derivatives
#' @noRd
matern_neumann_free <- function(s, t, kappa, sigma, nu=3/2, L = 1, deriv = c(0,0)){

    if(nu==3/2){

    }else if(nu==1/2){

    }else{
        stop('nu not yet implimented')
    }
    D <- outer(s, t, "-")
    phi_t <- phi1(t, kappa, sigma, nu, L, deriv = deriv[2])
    phi_s <- phi1(s, kappa, sigma, nu, L, deriv = deriv[1])
    A <- corrector_neumann_free(kappa, sigma, nu, L)
    if(sum(deriv)==0){
        r0 <- matern.covariance(D,kappa=kappa,nu=nu,sigma=sigma)
    }else{
        r0 <- (-1)^(deriv[2]+2)*matern_derivative(D,kappa=kappa,nu=nu,sigma=sigma, deriv = sum(deriv))
    }
    r <- r0 - t(phi_s)%*%A%*%phi_t
    return(r)
}



#' computes the covariance of free space (r'=0) neumann boundary with arbitrary
#' covarians on X'
#' @param s (n x 1) location
#' @param t (n x 1) location
#' @param C (2 x 2) covarians of X'
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
#' @param deriv  derivative order s, t
#' @noRd
matern_neumann_free2 <- function(s, t, C, kappa, sigma=1, nu=3/2, L = 1, deriv = c(0,0)){

    if(nu==3/2){

    }else{
        stop('nu not yet implimented')
    }
    D <- outer(s, t, "-")
    phi_t <- phi(t, kappa, sigma, nu, L, deriv = deriv[2])
    phi_s <- phi(s, kappa, sigma, nu, L, deriv = deriv[1])
    A <- corrector(kappa, sigma, nu, L)
    if(sum(deriv)==0){
        r0 <- rSPDE::matern.covariance(D,kappa=kappa,nu=nu,sigma=sigma)
    }else{
        r0 <- (-1)^(deriv[2]+2)*matern_derivative(D,kappa=kappa,nu=nu,sigma=sigma, deriv = sum(deriv))
    }
    r <- r0 - t(phi_s)%*%A%*%phi_t

    # extra term
    phi2_t <- phi_t[3:4,]
    phi2_s <- phi_s[3:4,]
    phi1_t <- phi_t[1:2,]
    phi1_s <- phi_s[1:2,]
    Ainv.e <- corrector_inverse_e(kappa, sigma, nu, L)
    B11.inv.B12 <- solve(Ainv.e$B11,t(Ainv.e$B12))
    Sigma_X1_tilde <- Ainv.e$B22 - t(Ainv.e$B12)%*%solve(Ainv.e$B11,Ainv.e$B12)
    Sigma_Xt_X1_tilde <- -t(phi2_t) + t(phi1_t)%*%B11.inv.B12
    Sigma_Xs_X1_tilde <- -t(phi2_s) + t(phi1_s)%*%B11.inv.B12
    Sigma_X1_tilde.inv <- solve(Sigma_X1_tilde)
    A2 <- Sigma_X1_tilde.inv%*%C%*%Sigma_X1_tilde.inv
    r <- r + Sigma_Xs_X1_tilde%*%A2%*%t(Sigma_Xt_X1_tilde)
    return(r)
}

#' Derivatives of the Matern covariance
#'
#' @param h distances where the covariance should be evaluated
#' @param kappa range parameter
#' @param nu smoothness parameter
#' @param sigma standard deviation
#' @param deriv order of derivative
#'
#' @return covariance function
#' @noRd
matern_derivative <- function(h, kappa, nu, sigma,deriv=1)
{
  if(deriv==1){
    C = h*rSPDE::matern.covariance(h,kappa=kappa,nu=nu-1,sigma=sigma)
    C[h==0] = 0
  } else if (deriv == 2){
    C = rSPDE::matern.covariance(h,kappa=kappa,nu=nu-1,sigma=sigma)+
      h*matern_derivative(h,kappa=kappa,nu=nu-1,sigma=sigma,deriv=1)

  } else {
    C = (deriv-1)*matern_derivative(h,kappa=kappa,nu=nu-1,sigma=sigma,deriv=deriv-2) +
      h*matern_derivative(h,kappa=kappa,nu=nu-1,sigma=sigma,deriv=deriv-1)
  }
  return(-(kappa^2/(2*(nu-1)))*as.matrix(C))
}

#' The corrector matrix, A, such that free space neumann is:
#' r(t,s)  =r_0(t-s) - phi(t) A phi(s)
#' where r_0 is stationary matern
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
#' @noRd
corrector_neumann_free <- function(kappa, sigma, nu=3/2, L = 1){
    if(nu==3/2){
        D <- matrix(c(0,L,L,0),2,2)
        B11 <- -matern.covariance(D,kappa=kappa,nu=3/2,sigma=sigma)
        B11[1,2] <- B11[2,1] <-  -B11[1,2]
    }else{
        stop('nu not yet implimented')
    }
    return(solve(B11))
}

#' The corrector matrix, A, such that neumann is:
#' r(t,s)  =r_0(t-s) - phi(t) A phi(s)
#' where r_0 is stationary matern
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
#' @return The corrector matrix
#' @noRd
corrector <- function(kappa, sigma, nu=3/2, L = 1){

    B <- corrector_inverse_e(kappa, sigma, nu, L)
    if(nu==1/2){
        B <- B$B11
    }else if(nu==3/2){
        B <- cbind( rbind(B$B11, B$B12) , rbind(t(B$B12), B$B22) )
    }else{
        stop('nu not yet implimented')
    }
    A <- solve(B)
    return(A)
}

#'
#' simple dim corrector function
#' @param t (n x 1) where to evaluate phi (in `[0,L]`)
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
#' @param deriv  how many derivatives
#' @return phi - (n x m)
#' @noRd
phi1 <- function(t, kappa, sigma, nu=3/2, L=1, deriv=0){
    n <- length(t)
    r <- matrix(0, nrow=2, ncol=n)
    if(deriv==0){
        r[1,] <- matern.covariance(t, kappa=kappa, nu = nu, sigma = sigma )
        r[2,] <- matern.covariance(t-L, kappa=kappa, nu = nu, sigma = sigma )
    }else{
        r[1,] <- matern_derivative(t, kappa=kappa, nu = nu, sigma = sigma , deriv = deriv)
        r[2,] <-  matern_derivative(t-L, kappa=kappa, nu = nu, sigma = sigma, deriv = deriv )

    }
    return(r)
}




#'
#' simple dim corrector function v2
#' @param t (n x 1) where to evaluate phi (in `[0,L]`)
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
#' @param deriv  how many derivatives
#' @return phi - (n x m)
#' @noRd
phi2 <- function(t, kappa, sigma, nu=3/2, L=1, deriv=0){
    n <- length(t)
    if(nu==3/2){
        r <- matrix(0, nrow=2, ncol=n)
        r[1,] <- -matern_derivative(t, kappa=kappa, nu = 3/2, sigma = sigma , deriv=deriv+1)
        r[2,] <- -matern_derivative(t-L, kappa=kappa, nu = 3/2, sigma = sigma, deriv=deriv+1 )

        return(r)
    }else{
        stop('nu not yet implimented')
    }
}



#'
#' one dim corrector function
#' @param t (n x 1) where to evaluate phi (in `[0,L]`)
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
#' @param deriv  how many derivatives
#' @return phi - (n x m)
#' @noRd
phi <- function(t, kappa, sigma, nu=3/2, L=1, deriv=0)
{
    n <- length(t)
    if(nu==1/2){
        r <- matrix(0, nrow=2, ncol=n)
        r[1:2,] <- phi1(t, kappa=kappa, nu = 1/2, sigma = sigma, L=L,deriv)

    }else if(nu==3/2){
        r <- matrix(0, nrow=4, ncol=n)
        r[1:2,] <- phi1(t, kappa=kappa, nu = 3/2, sigma = sigma, L=L,deriv)
        r[3:4,] <- -phi1(t, kappa=kappa, nu = 3/2, sigma = sigma,L=L ,deriv + 1)
    }else if(nu==5/2){
        alpha <- nu + 1/2
        r <- matrix(0, nrow=alpha*2, ncol=n)
        for(i in 1:alpha){
            r[2*(i-1)+(1:2),] <- (-1)^(i+1) * phi1(t, kappa=kappa, nu = nu, sigma = sigma, L=L,deriv + (i-1))
        }
    }else{

        stop('nu not yet implimented')
    }
    return(r)
}


#' inverse corrector elements
#' builds the elements of the inverse of the corrector matrix
#' r(t,s)  =r_0(t-s) - phi(t) A phi(s)
#' where r_0 is stationary matern
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
#' @noRd
corrector_inverse_e <- function(kappa, sigma, nu=3/2, L = 1){
    B.element <- list()
    if(nu ==1/2){
        D <- matrix(c(0,L,L,0),2,2)
        B11 <- -matern.covariance(D,kappa=kappa,nu=1/2,sigma=sigma)
        B11[1,2] <- B11[2,1] <-  -B11[1,2]
        B.element$B11 <- B11
    }else if(nu==3/2){
        D <- matrix(c(0,L,L,0),2,2)
        B11 <- -matern.covariance(D,kappa=kappa,nu=3/2,sigma=sigma)
        B11[1,2] <- B11[2,1] <-  -B11[1,2]
        B.element$B11 <- B11
        B12 <- matern_derivative(D,kappa=kappa,nu=3/2,sigma=sigma,1)
        B12[1,2] <-  -B12[1,2]
        B.element$B12 <- B12
        B22 <- -matern_derivative(D,kappa=kappa,nu=3/2,sigma=sigma,2)
        B.element$B22 <- B22
    }else{
        stop('nu not yet implimented')
    }
    return(B.element)
}



#' Starting values for random field models on metric graphs
#'
#' Computes appropriate starting values for optimization of Gaussian random
#' field models on metric graphs.
#'
#' @param graph A `metric_graph` object.
#' @param model Type of model, "alpha1", "alpha2", "isoExp", "GL1", and "GL2"
#' are supported.
#' @param data Should the data be used to obtain improved starting values?
#' @param data_name The name of the response variable in `graph$data`.
#' @param manual_data A vector (or matrix) of response variables.
#' @param range_par Should an initial value for range parameter be returned
#' instead of for kappa?
#' @param nu Should an initial value for nu be returned?
#' @param like_format Should the starting values be returned with sigma.e as the
#' last element? This is the format for the likelihood constructor from the
#' 'rSPDE' package.
#' @param log_scale Should the initial values be returned in log scale?
#'
#' @return A vector, `c(start_sigma_e, start_sigma, start_kappa)`
#' @export
graph_starting_values <- function(graph,
                                  model = c("alpha1", "alpha2", "isoExp",
                                            "GL1", "GL2"),
                                  data = TRUE,
                                  data_name = NULL,
                                  range_par = FALSE,
                                  nu = FALSE,
                                  manual_data = NULL,
                                  like_format = FALSE,
                                  log_scale = FALSE){

  check_graph(graph)

  model <- model[[1]]
  if((!model%in%c("alpha1", "alpha2", "isoExp", "GL1", "GL2"))){
    stop("The model should be one of 'alpha1', 'alpha2', 'isoExp',
      'GL1' or 'GL2'!")
  }
  if(data){
    if(is.null(graph$data) && is.null(manual_data)) {
      stop("No data provided, if you want the version without data set the 'data' argument to FALSE!")
    }
    if(is.null(data_name) && is.null(manual_data)){
      stop("If data is true, you must either supply the column data or manual data.")
    }
    if(!is.null(data_name)){
      y <- graph$data[[data_name]]
      y <- na.omit(y)
    }
    if(!is.null(manual_data)){
      y <- manual_data
      y <- na.omit(y)
    }
    data_std <- sqrt(var(as.vector(y)))
  } else{
    data_std <- NA
  }

  if(is.null(graph$geo_dist)){
        graph$compute_geodist(obs=FALSE)
  }
  finite_geodist <- is.finite(graph$geo_dist[[".vertices"]])
  finite_geodist <- graph$geo_dist[[".vertices"]][finite_geodist]
  prior.range.nominal <- max(finite_geodist) * 0.2

  if (model == "alpha1") {
    start_kappa <- sqrt(8 * 0.5) / prior.range.nominal
    #variance is sigma^2/2 kappa
    if(data){
      start_sigma <- sqrt(2*start_kappa) * data_std
    } else{
      start_sigma <- 1
    }
  } else if (model == "alpha2") {
    start_kappa <- sqrt(8 * 1.5) / prior.range.nominal
    if(data){
      #variance is sigma^2/(4 * kappa^3)
      start_sigma <- sqrt(4*start_kappa^3) * data_std
    } else{
      start_sigma <- 1
    }
  } else if (model == "isoExp") {
    start_kappa <- sqrt(8 * 0.5) / prior.range.nominal
    if(data){
      start_sigma <- data_std
    } else{
      start_sigma <- 1
    }
  } else if (model == "GL1") {
    if(is.null(graph$Laplacian)) {
      graph$compute_laplacian()
    }
    h <- mean(graph$edge_lengths)
    k <- sqrt(8 * 0.5) / prior.range.nominal
    start_kappa <- exp(-k*h)/(1-exp(-2*k*h)) + 2*exp(-k*h) - 2

    if(data){
      Q <- start_kappa^2*Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]]
      v <- rep(0,graph$nV)
      v[1] <- 1
      s2 <- solve(Q,v)[1]
      start_sigma <- data_std / sqrt(s2)
    } else{
      start_sigma <- 1
    }

  } else if (model == "GL2") {
    if(is.null(graph$Laplacian)) {
      graph$compute_laplacian()
    }
    h <- mean(graph$edge_lengths)
    k <- sqrt(8 * 0.5) / prior.range.nominal
    start_kappa <- exp(-k*h)/(1-exp(-2*k*h)) + 2*exp(-k*h) - 2
    if(data){
      Q <- start_kappa^2*Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]]
      v <- rep(0,graph$nV)
      v[1] <- 1
      s2 <- solve(Q %*% Q,v)[1]
      start_sigma <- data_std / sqrt(s2)
    } else{
      start_sigma <- 1
    }

  } else {
    stop("wrong model choice")
  }
  if(like_format){
      if(!nu){
        out_vec <- start_sigma
      } else{
        out_vec <- 1
      }

      if(range_par){
        out_vec <- c(out_vec, prior.range.nominal)
      } else{
        out_vec <- c(out_vec, start_kappa)
      }
      if(nu){
        out_vec <- c(out_vec,1)
      }
      out_vec <- c(out_vec, 0.1 * data_std)
  } else{
      if(!nu){
        out_vec <- c(0.1 * data_std,start_sigma)
      } else{
        out_vec <- c(0.1 * data_std,1)
      }

      if(range_par){
        out_vec <- c(out_vec, prior.range.nominal)
      } else{
        out_vec <- c(out_vec, start_kappa)
      }

      if(nu){
        out_vec <- c(out_vec,1)
      }
  }

  if(log_scale){
    out_vec <- log(out_vec)
  }

  return(out_vec)
}




#' Exponential covariance function
#'
#' Evaluates the exponential covariance function
#' \deqn{C(h) = \sigma^2 \exp\{-kappa h\}}
#'
#' @param h Distances to evaluate the covariance function at.
#' @param theta A vector `c(sigma, kappa)`, where `sigma` is the standard
#' deviation and `kappa` is a range-like parameter.
#'
#' @return A vector with the values of the covariance function.
#' @export
exp_covariance <- function(h, theta){
  sigma <- theta[1]
  kappa <- theta[2]
  return(sigma^2 * exp(-kappa * h))
}


#' Processing data to be used in add_observations
#' @noRd
process_data_add_obs <- function(PtE, new_data, old_data, group_vector){
  new_data[[".edge_number"]] <- PtE[,1]
  new_data[[".distance_on_edge"]] <- PtE[,2]

  if (is.null(group_vector)) {
    if(!is.null(old_data)){
      group_vector <- rep(old_data[[".group"]][1], length(PtE[,1]))
    } else{
      group_vector <- rep(1, length(PtE[,1]))
    }
  }

  if(is.null(old_data)){
      group_val <- unique(group_vector)
  } else {
      group_val <- unique(union(old_data[[".group"]], group_vector))
  }

  if (is.null(old_data)) {
    full_colnames <- names(new_data)
    data_coords_new <- data.frame(PtE1 = PtE[,1], PtE2 = PtE[,2])
    data_coords <- unique(data_coords_new)
    data_coords <- data_coords[order(data_coords$PtE1, data_coords$PtE2), ]

    data_coords_tmp <- data_coords
    # group_val <- unique(group_vector)
    n_group <- length(group_val)
    # data_coords[["group"]] <- group_val[[1]]
    # if (n_group>1) {
    #   for (i in 2:n_group) {
    #     tmp_coords <- data_coords_tmp
    #     tmp_coords[["group"]] <- group_val[[i]]
    #     data_coords <- rbind(data_coords, tmp_coords)
    #   }
    # }

    data_coords_1 <- rep(data_coords_tmp[,1], times = n_group)
    data_coords_2 <-  rep(data_coords_tmp[,2], times = n_group)
    data_coords_3 <- rep(group_val, each = length(data_coords_tmp[,1]))

    data_coords <- data.frame(PtE1 = data_coords_1, PtE2 = data_coords_2, group = data_coords_3)
    rm(data_coords_1)
    rm(data_coords_2)
    rm(data_coords_3)

    data_coords_new[["group"]] <- group_vector
    data_coords[["idx"]] <- 1:nrow(data_coords)
    idx_new_entries <- merge(data_coords_new, data_coords, all=FALSE,
                             sort = FALSE)

    idx_new_entries <- idx_new_entries[["idx"]]
    list_result <- vector(mode = "list", length(full_colnames))
    names(list_result) <- full_colnames
    list_result[1:length(list_result)] <- full_colnames
    list_result[[".edge_number"]] <- NULL
    list_result[[".distance_on_edge"]] <- NULL
    list_result[[".group"]] <- NULL
    new_data <- lapply(list_result, function(col_name){
          mode_vector <- typeof(new_data[[col_name]])
          tmp <- vector(mode=mode_vector, length = nrow(data_coords))
          is.na(tmp) <- 1:length(tmp)
          #  for(i in 1:length(idx_new_entries)){
            for(i in 1:length(new_data[[col_name]])){
               tmp[[idx_new_entries[i]]] <- new_data[[col_name]][[i]]
            }
            return(tmp)
          })

    new_data[[".edge_number"]] <- data_coords[["PtE1"]]
    new_data[[".distance_on_edge"]] <- data_coords[["PtE2"]]
    new_data[[".group"]] <- data_coords[["group"]]
    return(new_data)
  } else {
    old_colnames <- names(old_data)
    new_colnames <- names(new_data)
    full_colnames <- union(old_colnames, new_colnames)

    new_df <- data.frame(PtE1 = PtE[,1], PtE2 = PtE[,2])

    old_df <- data.frame(PtE1 = old_data[[".edge_number"]],
                         PtE2 = old_data[[".distance_on_edge"]])

    data_coords <- unique(rbind(old_df, new_df))
    data_coords <- data_coords[order(data_coords$PtE1, data_coords$PtE2), ]

    data_coords_tmp <- data_coords
    # group_val <- unique(group_vector)
    n_group <- length(group_val)
    data_coords[["group"]] <- group_val[[1]]
    if (n_group>1) {
      for (i in 2:n_group) {
          tmp_coords <- data_coords_tmp
          tmp_coords[["group"]] <- group_val[[i]]
          data_coords <- rbind(data_coords, tmp_coords)
      }
    }
    data_coords <- as.data.frame(data_coords)
    new_df[["group"]] <- group_vector
    old_df[["group"]] <- old_data[[".group"]]
    data_coords[["idx"]] <- 1:nrow(data_coords)

    idx_new_entries <- merge(new_df, data_coords, all = FALSE, sort = FALSE)
    idx_new_entries <- idx_new_entries[["idx"]]
    idx_old_entries <- merge(old_df, data_coords, all = FALSE, sort = FALSE)
    idx_old_entries <- idx_old_entries[["idx"]]
    list_result <- vector(mode = "list", length(full_colnames))
    names(list_result) <- full_colnames
    list_result[1:length(list_result)] <- full_colnames
    list_result[[".edge_number"]] <- NULL
    list_result[[".distance_on_edge"]] <- NULL
    list_result[[".group"]] <- NULL
    list_result <- lapply(list_result, function(col_name){
        if(!is.null(new_data[[col_name]])){
          mode_vector <- typeof(new_data[[col_name]])
        } else{
          mode_vector <- typeof(old_data[[col_name]])
        }
        tmp <- vector(mode=mode_vector, length = nrow(data_coords))
        is.na(tmp) <- 1:length(tmp)

        if(length(idx_new_entries)>0){
          for(i in 1:length(idx_new_entries)){
            if(!is.null(new_data[[col_name]][[i]])){
              tmp[[idx_new_entries[i]]] <- new_data[[col_name]][[i]]
            } else{
              tmp[[idx_new_entries[i]]] <- rep(NA, length(idx_new_entries[i]))
            }
          }
        }
        if(length(idx_old_entries)>0){
          for(i in 1:length(idx_old_entries)){
            if(!is.null(old_data[[col_name]][[i]])){
              tmp[[idx_old_entries[i]]] <- old_data[[col_name]][[i]]
            }
          }
        }
        return(tmp)
      })
  list_result[[".edge_number"]] <- data_coords[["PtE1"]]
  list_result[[".distance_on_edge"]] <- data_coords[["PtE2"]]
  list_result[[".group"]] <- data_coords[["group"]]
  return(list_result)
  }
}

#' find indices of the rows with all NA's in lists
#' @noRd
#'
idx_not_all_NA <- function(data_list){
     data_list[[".edge_number"]] <- NULL
     data_list[[".distance_on_edge"]] <- NULL
     data_list[[".coord_x"]] <- NULL
     data_list[[".coord_y"]] <- NULL
     data_list[[".group"]] <- NULL
     data_names <- names(data_list)
     n_data <- length(data_list[[data_names[1]]])
     idx_non_na <- logical(n_data)
     for(i in 1:n_data){
        na_idx <- lapply(data_list, function(dat){
          return(is.na(dat[i]))
        })
        idx_non_na[i] <- !all(unlist(na_idx))
     }
     return(idx_non_na)
}

#' find indices of the rows with at least one NA's in lists
#' @noRd
#'
idx_not_any_NA <- function(data_list){
     data_list[[".edge_number"]] <- NULL
     data_list[[".distance_on_edge"]] <- NULL
     data_list[[".coord_x"]] <- NULL
     data_list[[".coord_y"]] <- NULL
     data_list[[".group"]] <- NULL
     data_names <- names(data_list)
     n_data <- length(data_list[[data_names[1]]])
     idx_non_na <- logical(n_data)
     for(i in 1:n_data){
        na_idx <- lapply(data_list, function(dat){
          return(is.na(dat[i]))
        })
        idx_non_na[i] <- !any(unlist(na_idx))
     }
     return(idx_non_na)
}


#' Select replicate
#' @noRd
#'
select_group <- function(data_list, group){
    grp <- data_list[[".group"]]
    grp <- which(grp %in% group)
    data_result <- lapply(data_list, function(dat){dat[grp]})
    return(data_result)
}

#' Create lines for package name
#'
#' @return `SpatialLines` object with package name.
#' @export
logo_lines <- function(){
  n <- 100
  #G
  theta <- seq(from=pi,to=3*pi/2,length.out = n)
  line1 <- cbind(1+sin(theta),2+2*cos(theta))

  theta <- seq(from=pi/2,to=pi,length.out = n)
  line2 <- cbind(1+sin(theta),1.5+1.5*cos(theta))

  theta <- seq(from=3*pi/2,to=2*pi,length.out = n)
  line3 <- cbind(2+2*sin(theta),2+2*cos(theta))

  line4 <- rbind(c(1,1.5),c(2,1.5))

  #R
  line5 <- rbind(c(2,0),c(2,4))
  line6 <- rbind(c(2,4),c(3,4))
  theta <- seq(from=0,to=pi,length.out = n)
  line7 <- cbind(3+sin(theta),3+cos(theta))
  line8 <- rbind(c(3,2),c(2,2))
  line9 <- rbind(c(2,2),c(4,0))

  #A
  line10 <- rbind(c(4,0),c(5,4))
  line11 <- rbind(c(5,4),c(6,0))
  line12 <- rbind(c(4.5,2),c(5.5,2))

  #P
  line13 <- rbind(c(6,0),c(6,4))
  line14 <- rbind(c(6,4),c(7,4))
  theta <- seq(from=0,to=pi,length.out = n)
  line15 <- cbind(7+sin(theta),3+cos(theta))
  line16 <- rbind(c(7,2),c(6,2))

  #H
  line17 <- rbind(c(8,0),c(8,4))
  line18 <- rbind(c(10,0),c(10,4))
  line19 <- rbind(c(8,2),c(10,2))

  #M
  line20 <- rbind(c(0,4),c(0.75,8))
  line21 <- rbind(c(0.75,8),c(1.5,5))
  line22 <- rbind(c(1.5,5),c(2.25,8))
  line23 <- rbind(c(2.25,8),c(3,4))

  # E
  line24 <- rbind(c(3,4),c(3,8))
  line25 <- rbind(c(3,8),c(4,8))
  line26 <- rbind(c(3,6),c(4,6))
  line27 <- rbind(c(3,4),c(5,4))

  # T
  line28 <- rbind(c(5,4),c(5,8))
  line29 <- rbind(c(4,8),c(6,8))


  # R
  line30 <- rbind(c(6,4),c(6,8))
  line31 <- rbind(c(6,8),c(7,8))
  theta <- seq(from=0,to=pi,length.out = n)
  line32 <- cbind(7+sin(theta),7+cos(theta))
  line33 <- rbind(c(7,6),c(6,6))
  line34 <- rbind(c(6,6),c(8,4))

  # I
  line35 <- rbind(c(8,4),c(8,8))

  # C
  theta <- seq(from=pi,to=3*pi/2,length.out = n)
  line36 <- cbind(10+2*sin(theta),6+2*cos(theta))
  theta <- seq(from=3*pi/2,to=2*pi,length.out = n)
  line37 <- cbind(10+2*sin(theta),6+2*cos(theta))

  return(list(line1,
                               line2,
                               line3,
                               line4,
                               line5,
                               line6,
                               line7,
                               line8,
                               line9,
                               line10,
                               line11,
                               line12,
                               line13,
                               line14,
                               line15,
                               line16,
                               line17,
                               line18,
                               line19,
                               line20,
                               line21,
                               line22,
                               line23,
                               line24,
                               line25,
                               line26,
                               line27,
                               line28,
                               line29,
                               line30,
                               line31,
                               line32,
                               line33,
                               line34,
                               line35,
                               line36,
                               line37))
}



#' @noRd

projectVecLine2 <- function(lines, points, normalized = FALSE){
  return(projectVecLine(lines, points, normalized))
}

#' @noRd

distance2 <- function(points, lines, byid=FALSE, longlat, crs){
  if(!is.null(points)){
    class(points) <- setdiff(class(points), "metric_graph_edge")
  }
  if(!is.null(lines)){
    class(lines) <- setdiff(class(lines), "metric_graph_edge")
  }

  if(!longlat){
    points_sf <- sf::st_as_sf(as.data.frame(points), coords = 1:2)
  } else{
    points_sf <- sf::st_as_sf(as.data.frame(points), coords = 1:2, crs = crs)
  }
  if(!is.list(lines)){
    lines <- list(lines)
  }

  if(!longlat){
    lines_sf <- sf::st_sfc(lapply(lines, function(i){sf::st_linestring(i)}))
  } else{
    lines_sf <- sf::st_sfc(lapply(lines, function(i){sf::st_linestring(i)}), crs = crs)
  }

  dist_result <- sf::st_distance(points_sf, lines_sf)
  if(byid){
    ID_names <- 1:length(lines)
    dist_result <- t(dist_result)
    row.names(dist_result) <- ID_names
    return(dist_result)
  } else{
    return(min(dist_result))
  }
}

#' @noRd

intersection2 <- function(lines1, lines2){
  lines1_sf <- sf::st_as_sf(lines1)
  lines2_sf <- sf::st_as_sf(lines2)
  inter_lines <- sf::st_intersection(lines1_sf, lines2_sf)
  inter_lines <- unique(inter_lines)
  return(inter_lines)
}

#' @noRd

intersection3 <- function(lines1_sf, lines2_sf){
  inter_lines <- sf::st_intersection(lines1_sf, lines2_sf)
  inter_lines <- unique(inter_lines)
  return(sf::st_as_sfc(inter_lines))
}

#' @noRd

interpolate2 <- function(lines, pos, normalized = FALSE, get_idx = FALSE){
  if(!get_idx){
    return(interpolate2_aux(lines, pos, normalized)[["coords"]])
  } else{
    return(interpolate2_aux(lines, pos, normalized))
  }
}

#' @noRd

make_Aprd <- function(graph, edge_number, distance_on_edge){
  X_loc <- cbind(edge_number, distance_on_edge)
  order_X <- order(X_loc[,1], X_loc[,2])
  X_loc <- X_loc[order_X,]
  edge_n <- sort(unique(edge_number))
  idx_tmp <- (graph$data[[".edge_number"]] == edge_n[1])
  mesh_tmp <- graph$data[[".distance_on_edge"]][idx_tmp]
  loc <- (X_loc[X_loc[,1] == edge_n[1], 2])
  A_prd <- rSPDE::rSPDE.A1d(mesh_tmp, )
}

#' @noRd

change_parameterization_graphlme <- function(likelihood, nu, par, hessian
){
  tau <- par[1]
  kappa <- par[2]

  C1 <- sqrt(8*nu)
  C2 <- sqrt(gamma(nu) / ((4 * pi)^(1 / 2) * gamma(nu + 1 / 2)))

  sigma <- C2 /(tau * kappa^nu)
  range <- C1/kappa

  grad_par <- matrix(c(-C2/(kappa^nu * sigma^2),0,
                    nu * range^(nu-1) * C2/(sigma * C1^nu),
                    -C1/range^2), nrow = 2, ncol=2)


  new_observed_fisher <- t(grad_par) %*% hessian %*% (grad_par)

  # No need to include the additional term as the gradient is approximately zero.
  # from some numerical experiments, the approximation without the additional term
  # seems to be better in general.

  inv_fisher <- tryCatch(solve(new_observed_fisher),
                         error = function(e) matrix(NA, nrow(new_observed_fisher),
                                                    ncol(new_observed_fisher)))

  std_err <- sqrt(diag(inv_fisher))

  return(list(coeff = c(sigma, range), std_random = std_err))
}

#' @noRd 

process_factor_unit <- function(vertex_unit, length_unit){
  if(is.null(vertex_unit) && is.null(length_unit)){
    return(1)
  } else if(is.null(vertex_unit) || is.null(length_unit)){
          stop("If one of 'vertex_unit' or 'length_unit' is NULL, then the other must also be NULL.")
  }
  if(vertex_unit == length_unit){
    return(1)
  } else if(vertex_unit == "degrees"){
    fact <- switch(length_unit, "km" = 1,
                        "m" = 1000,
                        "miles" = 0.621371192)
    return(fact)
  } else if(vertex_unit == "km"){
    fact <- switch(length_unit, "km" = 1,
                        "m" = 1000,
                        "miles" = 0.621371192) 
    return(fact) 
  } else if(vertex_unit == "m"){
    fact <- switch(length_unit, "km" = 1e-3,
                        "m" = 1,
                        "miles" = 0.621371192*1e-3) 
    return(fact)
  } else if(vertex_unit == "miles"){
    fact <- switch(length_unit, "km" = 1.609344,
                        "m" = 1.609344*1e3,
                        "miles" = 1) 
    return(fact)
  }
}


#' code from https://gist.github.com/MansMeg/1ec56b54e1d9d238b4fd
#' 
#' Message progress bar
#' 
#' @description 
#' A simple progress bar to use in R packages where messages are prefered to console output.
#' 
#' @field iter Total number of iterations
#' @field i Current iteration
#' @field width Width of the R console
#' @field width_bar Width of the progress bar
#' @field progress The number of character printed (continous)
#' @field progress_step Addition to progress per iteration
#' 
#' @examples
#' test_bar <- function(i = 10){
#'  bar <- msg_progress_bar(i)
#'  for(j in 1:i){
#'    bar$increment()
#'    Sys.sleep(4/i)
#'    }
#'  }
#'  test_bar(100)
#'   
#' @author Mans Magnusson (MansMeg @ github)
#'   
#' @noRd 
msg_progress_bar <- 
  setRefClass(
    Class = "msg_progress_bar", 
    fields = list(iter = "numeric",
                  i = "numeric",
                  progress = "numeric",
                  progress_step = "numeric",
                  width = "numeric",
                  width_bar = "numeric"),
    
    methods = list(
      initialize = function(iter){
        'Initialize a messagebar object'
        .self$width <- getOption("width")
        .self$iter <- iter
        .self$i <- 0
        .self$progress <- 0
        white_part <- paste(rep(" ", (.self$width - 11) %/% 4), collapse="")
        init_length <- .self$width - ((.self$width - 11) %/% 4) * 4 - 11
        white_init <- paste(rep(" ", init_length), collapse="")
        .self$width_bar <- .self$width - init_length - 2 + 0.1
        .self$progress_step <- .self$width_bar / .self$iter
        message(paste(white_init, "|", white_part, "25%", white_part, "50%", white_part, "75%", white_part, "|","\n", white_init, "|", sep=""), appendLF = FALSE)
      },
      
      increment = function(){
        'A messagebar object.'
        if(.self$i > .self$iter) return(invisible(NULL))
        new_progress <- .self$progress + .self$progress_step
        diff_in_char <- floor(new_progress) - floor(.self$progress)
        if(diff_in_char > 0) {
          message(paste(rep("=", diff_in_char),collapse=""), appendLF = FALSE)
        }
        
        .self$progress <- new_progress
        .self$i <- .self$i + 1
        if(.self$i == .self$iter) message("|\n", appendLF = FALSE)
        
      }
    )
  )

  #' @noRd
  #' 
  
  get_rel_pos_prune <- function(which_line_starts, Line_1, Line_2, start_1, end_1, start_2, end_2, length_line1, length_line2){

      if(Line_1 != Line_2){
        total_length <- length_line2 + length_line1
        if(which_line_starts == 1){
          # if(start_2 == 0) {
          #   end <- (end_2*length_line2+length_line1)/total_length
          # } else {
          #   end <- ((end_2 - start_2)*length_line2+length_line1)/total_length
          # }
          # if(end_1 == 1){
          #   start <- (start_1 * length_line1)/total_length
          # } else{
          #   start <- ((end_1 - start_1)*length_line1+length_line2)/total_length
          # }
          start <- start_1 * length_line1 / total_length
          end <- (end_2 * length_line2 + length_line1)/ total_length
        } else{
          # if(start_1 == 0) {
          #   end <- (end_1*length_line1+length_line2)/total_length
          # } else {
          #   end <- ((end_1 - start_1)*length_line1+length_line2)/total_length
          # }
          # if(end_2 == 1){
          #   start <- (start_2 * length_line2)/total_length
          # } else{
          #   start <- ((end_2 - start_2)*length_line2+length_line1)/total_length
          # }
          start <- start_2 * length_line2 / total_length
          end <- (end_1 * length_line1 + length_line2)/ total_length
        }
      } else {
          end <- max(end_1,end_2)
          start <- min(start_1,start_2)
      }

      return(list(start = start, end = end))
  }



#' @noRd 
#' 

get_vertex_pos_in_line <- function(V, coords_line){
    return(which.min(sapply(1:nrow(coords_line), function(i){norm(as.matrix(V - coords_line[i,]))})))
}


#' @noRd 
#' 

check_lines_input <- function(lines){
  is_matrix <- sapply(lines, function(i){is.matrix(i)})
  is_data_frame <- sapply(lines, function(i){is.data.frame(i)})
  if(!all(is_matrix | is_data_frame)) {
    stop("The list must contain either matrices of data.frames!")
  }
    n_cols <- sapply(lines, ncol)
  if(any(n_cols != 2)){
    stop("The elements in the list must have two columns!")
  }
  lines <- lapply(lines, function(i){as.matrix(i)})
  return(lines)
}

#' @noRd 
#' 

compute_line_lengths <- function(edge, longlat, unit, crs, proj4string, which_longlat, vertex_unit, project_data){
  if(!is.null(edge)){
      class(edge) <- setdiff(class(edge), "metric_graph_edge")
  }
    if(!longlat || project_data){
      fact <- process_factor_unit(vertex_unit, unit)
      return(compute_length(edge) * fact)
    } else if(which_longlat == "sf"){
      if(!is.null(edge)){
        linestring <- sf::st_sfc(sf::st_linestring(edge), crs = crs)
        length <- sf::st_length(linestring)
        units(length) <- unit
        units(length) <- NULL
        return(length)
      } else{
        return(0)
      }
    } else{
      Line <- sp::Line(edge)
      length <- sp::LineLength(Line, longlat = longlat)
      fact <- process_factor_unit(vertex_unit, unit)
      return(length * fact)
    }
}


#' @noRd 
#' 

compute_aux_distances <- function(lines, crs, longlat, proj4string, points = NULL, fact, which_longlat, length_unit){
  if(!is.null(points)){
    class(points) <- setdiff(class(points), "metric_graph_edge")
  }
  if(!is.null(lines)){
    class(lines) <- setdiff(class(lines), "metric_graph_edge")
  }
    if(!longlat){
      if(is.null(points)){
        dists <- dist(lines) * fact
      } else{
        if(nrow(lines)>nrow(points)){
          if(nrow(points)>1){
            stop("Points must have either the same number of rows as lines, or only 1 row!")
          }
        }
        dists <- sqrt((lines[,1] - points[,1])^2 + (lines[,2]-points[,2])^2) * fact
      }
    } else if (which_longlat == "sf") {
        sf_points <- sf::st_as_sf(as.data.frame(lines), coords = 1:2, crs = crs)
        if(is.null(points)){
          dists <- sf::st_distance(sf_points, which = "Great Circle")
        } else{
          sf_p_points <- sf::st_as_sf(as.data.frame(points), coords = 1:2, crs = crs)
          dists <- sf::st_distance(x = sf_points, y = sf_p_points, which = "Great Circle", by_element = TRUE)
        }
        units(dists) <- length_unit
        units(dists) <- NULL
    } else{
        sp_points <- sp::SpatialPoints(coords = lines, proj4string = proj4string) 
        if(is.null(points)){
          dists <- sp::spDists(sp_points, longlat = TRUE) * fact
        } else{
          sp_p_points <- sp::SpatialPoints(coords = points, proj4string = proj4string) 
          dists <- sp::spDists(x = sp_points, y=sp_p_points, longlat = TRUE, diagonal = TRUE) * fact
        }
    }
    return(dists)
}





#' A version of `dplyr::select()` function for datasets on metric graphs
#'
#' Selects columns on metric graphs, while keeps the spatial positions.
#'
#' @aliases select select.metric_graph_data
#' @param .data The data list or `tidyr::tibble` obtained from a metric graph object.
#' @param ... Additional parameters to be passed to `dplyr::select()`.
#' @return A `tidyr::tibble` with the resulting selected columns.
#' @method select metric_graph_data
#' @export
#' 
select.metric_graph_data <- function(.data, ...){
    bkp <- list()
    bkp[[".group"]] <- .data[[".group"]] 
    bkp[[".edge_number"]] <- .data[[".edge_number"]]
    bkp[[".distance_on_edge"]] <- .data[[".distance_on_edge"]]
    bkp[[".coord_x"]] <- .data[[".coord_x"]]
    bkp[[".coord_y"]] <- .data[[".coord_y"]]

    data_res <- dplyr::select(.data = tidyr::as_tibble(.data), ...)
    data_res[[".group"]] <- bkp[[".group"]] 
    data_res[[".edge_number"]] <- bkp[[".edge_number"]]
    data_res[[".distance_on_edge"]] <- bkp[[".distance_on_edge"]]
    data_res[[".coord_x"]] <- bkp[[".coord_x"]]
    data_res[[".coord_y"]] <- bkp[[".coord_y"]]
    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }    
    return(data_res)
}

#' A version of `dplyr::mutate()` function for datasets on metric graphs
#'
#' Applies `dplyr::mutate()` function for datasets obtained from a metric graph object.
#'
#' @aliases mutate mutate.metric_graph_data
#' @param .data The data list or `tidyr::tibble` obtained from a metric graph object.
#' @param ... Additional parameters to be passed to `dplyr::mutate()`.
#' @return A `tidyr::tibble` with the resulting selected columns.
#' @method mutate metric_graph_data
#' @export
#' 
mutate.metric_graph_data <- function(.data, ...){
    data_res <- dplyr::mutate(.data = tidyr::as_tibble(.data), ...)
    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }    
    return(data_res)
}


#' A version of `tidyr::drop_na()` function for datasets on metric graphs
#'
#' Applies `tidyr::drop_na()` function for datasets obtained from a metric graph object.
#'
#' @aliases drop_na drop_na.metric_graph_data
#' @param data The data list or `tidyr::tibble` obtained from a metric graph object.
#' @param ... Additional parameters to be passed to `tidyr::drop_na()`.
#' @return A `tidyr::tibble` with the resulting selected columns.
#' @method drop_na metric_graph_data
#' @export
#' 
drop_na.metric_graph_data <- function(data, ...){
    data_res <- tidyr::drop_na(data = tidyr::as_tibble(data), ...)
    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }    
    return(data_res)
}


#' A version of `dplyr::filter()` function for datasets on metric graphs
#'
#' Applies `dplyr::filter()` function for datasets obtained from a metric graph object.
#'
#' @aliases filter filter.metric_graph_data
#' @param .data The data list or `tidyr::tibble` obtained from a metric graph object.
#' @param ... Additional parameters to be passed to `dplyr::filter()`.
#' @return A `tidyr::tibble` with the resulting selected columns.
#' @method filter metric_graph_data
#' @export
#' 
filter.metric_graph_data <- function(.data, ...){
    data_res <- dplyr::filter(.data = tidyr::as_tibble(.data), ...)
    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }    
    return(data_res)   
}


#' A version of `dplyr::summarise()` function for datasets on metric graphs
#'
#' Creates summaries, while keeps the spatial positions.
#'
#' @aliases summarise summarise.metric_graph_data
#' @param .data The data list or `tidyr::tibble` obtained from a metric graph object.
#' @param ... Additional parameters to be passed to `dplyr::summarise()`.
#' @param .include_graph_groups Should the internal graph groups be included in the grouping variables? The default is `FALSE`. This means that, when summarising, the data will be grouped by the internal group variable together with the spatial locations.
#' @param .groups A vector of strings containing the names of the columns to be additionally grouped, when computing the summaries. The default is `NULL`.
#' @return A `tidyr::tibble` with the resulting selected columns.
#' @method summarise metric_graph_data
#' @export
#' 
summarise.metric_graph_data <- function(.data, ..., .include_graph_groups = FALSE, .groups = NULL){
    group_vars <- c(".edge_number", ".distance_on_edge", ".coord_x", ".coord_y")
    if(.include_graph_groups){
      group_vars <- c(".group", group_vars)
    }
    group_vars <- c(.groups, group_vars)
    previous_groups <- as.character(dplyr::groups(.data))
    group_vars <- c(previous_groups, group_vars)

    data_res <- dplyr::group_by_at(.tbl = tidyr::as_tibble(.data), .vars = group_vars)
    data_res <- dplyr::summarise(.data = data_res, ...)
    data_res <- dplyr::ungroup(data_res)
    if(is.null(data_res[[".group"]])){
      data_res[[".group"]] <- 1
    }

    ord_data <- order(data_res[[".group"]], data_res[[".edge_number"]], data_res[[".distance_on_edge"]])

    data_res <- data_res[ord_data,]

    if(!inherits(data_res, "metric_graph_data")){
      class(data_res) <- c("metric_graph_data", class(data_res))
    }        
    return(data_res)
}


#' Pipe operator
#'
#' See \code{\link[magrittr]{%>%}} for more details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @usage lhs \%>\% rhs
NULL


#' @name summary.metric_graph
#' @title Summary Method for \code{metric_graph} Objects
#' @description Function providing a summary of several informations/characteristics of a metric graph object.
#' @param object an object of class `metric_graph`.
#' @param messages Should message explaining how to build the results be given for missing quantities?
#' @param compute_characteristics Should the characteristics of the graph be computed?
#' @param check_euclidean Check if the graph has Euclidean edges?
#' @param check_distance_consistency Check the distance consistency assumption?#' 
#' @param ... not used.
#' @return An object of class \code{summary_graph_lme} containing information
#' about a *metric_graph* object.
#' @method summary metric_graph
#' @export
summary.metric_graph <- function(object, messages = FALSE, compute_characteristics = TRUE, check_euclidean = TRUE, check_distance_consistency = TRUE, ...){
  object$summary(messages = messages, compute_characteristics = compute_characteristics, check_euclidean = check_euclidean, check_distance_consistency = check_distance_consistency)
}





#' @name print.metric_graph_vertices
#' @title Print Method for \code{metric_graph_vertices} Objects
#' @description Provides a brief description of the vertices of a metric graph
#' @param x object of class `metric_graph_vertices`.
#' @param n number of rows to show
#' @param ... Currently not used.
#' @return No return value. Called for its side effects.
#' @noRd
#' @method print metric_graph_vertices
#' @export
print.metric_graph_vertices <- function(x, n = 10, ...) {
  cat("Vertices of the metric graph\n\n")
  cat("Longitude and Latitude coordinates:", attr(x[[1]], "longlat"), "\n")
  if(attr(x[[1]], "longlat")){
    cat("Coordinate reference system:",attr(x[[1]], "crs"), "\n")
  }
  if(attr(x[[1]], "longlat")){
    lab_x = "Longitude"
    lab_y = "Latitude"
  } else{
    lab_x <- "x"
    lab_y <- "y"
  }
  cat("\nSummary: \n")
  coord_tmp <- matrix(nrow = length(x), ncol = 6)
  coord_tmp <- as.data.frame(coord_tmp)
  for(i in 1:length(x)){
    coord_tmp[i,1:5] <- c(x[[i]], attr(x[[i]], "degree"),attr(x[[i]], "indegree"), attr(x[[i]], "outdegree"))
    coord_tmp[i,6] <- attr(x[[i]], "problematic")
  }
  rownames(coord_tmp) <- 1:length(x)
  colnames(coord_tmp) <- c(lab_x, lab_y, "Degree", "Indegree", "Outdegree", "Problematic")
  print(coord_tmp[1:min(n, nrow(coord_tmp)),])
  if(n < nrow(coord_tmp)){
    message(paste("#", nrow(coord_tmp)-n,"more rows"))
    message("# Use `print(n=...)` to see more rows")
  }
}




#' @name print.metric_graph_edges
#' @title Print Method for \code{metric_graph_edges} Objects
#' @description Provides a brief description of the edges of a metric graph
#' @param x object of class `metric_graph_edges`.
#' @param n number of edges to show
#' @param ... Currently not used.
#' @return No return value. Called for its side effects.
#' @noRd
#' @method print metric_graph_edges
#' @export
print.metric_graph_edges <- function(x, n = 4, ...) {
  cat("Edges of the metric graph\n\n")
  cat("Longitude and Latitude coordinates:", attr(x[[1]], "longlat"), "\n")
  if(attr(x[[1]], "longlat")){
    cat("Coordinate reference system:",attr(x[[1]], "crs"), "\n")
  }
  if(attr(x[[1]], "longlat")){
    lab_x = "Longitude"
    lab_y = "Latitude"
  } else{
    lab_x <- "x"
    lab_y <- "y"
  }
  edge_lengths <- 
  cat("\nSummary: \n\n")
  for(i in 1:min(n,length(x))){
    edge <- x[[i]]
    edge_df <- data.frame(x = edge[,1], y = edge[,2])
    n_edge_df <- nrow(edge_df)
    edge_df <- edge_df[c(1,n_edge_df),]
    colnames(edge_df) <- c(lab_x,lab_y)
    cat(paste0("Edge ",i," (first and last coordinates):"),"\n")
    print(edge_df, row.names=FALSE)
    cat("Total number of coordinates:",n_edge_df,"\n")
    if(!is.null(attr(attr(x[[i]],"length"),"units"))){
      cat("Edge length:", attr(x[[i]], "length"),units(attr(x[[i]], "length"))$numerator,"\n")
    } else{
      cat("Edge length:", attr(x[[i]], "length"),"\n")
    }
    if(is.data.frame(attr(x[[i]], "weight"))){
      cat("Weights: \n")
      print(attr(x[[i]], "weight"), row.names=FALSE)
      cat("\n")
    } else{
      cat("Weight:", attr(x[[i]], "weight"),"\n\n")
    }
    
  }
  if(n < length(x)){
    message(paste("#", length(x)-n,"more edges"))
    message("# Use `print(n=...)` to see more edges")
  }
}

#' @name print.metric_graph_edge
#' @title Print Method for \code{metric_graph_edge} Objects
#' @description Provides a brief description of the chosen edge of a metric graph
#' @param x object of class `metric_graph_edge`.
#' @param n number of coordinates to show
#' @param ... Currently not used.
#' @return No return value. Called for its side effects.
#' @noRd
#' @method print metric_graph_edge
#' @export
print.metric_graph_edge <- function(x, n = 4, ...) {
  cat("Edge",attr(x,"id"),"of the metric graph\n\n")
  cat("Longitude and Latitude coordinates:", attr(x, "longlat"), "\n")
  if(attr(x, "longlat")){
    cat("Coordinate reference system:",attr(x, "crs"), "\n")
  }
  if(attr(x, "longlat")){
    lab_x = "Longitude"
    lab_y = "Latitude"
  } else{
    lab_x <- "x"
    lab_y <- "y"
  }
  edge_lengths <- 
  cat("\nCoordinates of the vertices of the edge: \n")
  edge_df <- data.frame(a = x[,1], b = x[,2]) 
  n_edge_df <- nrow(edge_df)
  edge_df <- edge_df[c(1,n_edge_df),]
  colnames(edge_df) <- c(lab_x,lab_y)
  print(edge_df, row.names=FALSE)

  cat("\n")

  cat("Coordinates of the edge:\n")
  edge_df <- data.frame(a = x[,1], b = x[,2]) 
  colnames(edge_df) <- c(lab_x,lab_y)
  print(edge_df[1:min(n,nrow(edge_df)),], row.names=FALSE)
  if(n < nrow(edge_df)){
    message(paste("#", nrow(x)-n,"more coordinates"))
    message("# Use `print(n=...)` to see more coordinates")
  }
  
  cat("\n")

  if(is.null(attr(x, "PtE"))){
    message("Relative positions of the coordinates on the graph edges were not computed.")
    message("To compute them, run the `compute_PtE_edges()` method.")
  } else{
  cat("Relative positions of the edge:\n")
  PtE <- attr(x, "PtE")
  PtE_df <- data.frame(a = PtE[,1], b = PtE[,2]) 
  colnames(PtE_df) <- c("Edge number","Distance on edge")
  print(PtE_df[1:min(n,nrow(edge_df)),], row.names=FALSE)
  if(n < nrow(PtE_df)){
    message(paste("#", nrow(PtE_df)-n,"more relative positions"))
    message("# Use `print(n=...)` to see more relative positions")
  }
  }
  
  cat("\n")
  cat("Total number of coordinates:",nrow(edge_df),"\n")
    if(!is.null(attr(attr(x,"length"),"units"))){
      cat("Edge length:", attr(x, "length"),units(attr(x, "length"))$numerator,"\n")
    } else{
      cat("Edge length:", attr(x, "length"),"\n")
    }
    if(is.data.frame(attr(x, "weight"))){
      cat("Weights: \n")
      print(attr(x, "weight"), row.names=FALSE)
      cat("\n")
    } else{
      cat("Weight:", attr(x, "weight"),"\n")
    }

}




#' @name print.metric_graph_vertex
#' @title Print Method for \code{metric_graph_vertice} Objects
#' @description Provides a brief description of the chosen vertex of a metric graph
#' @param x object of class `metric_graph_vertex`.
#' @param n number of rows to show
#' @param ... Currently not used.
#' @return No return value. Called for its side effects.
#' @noRd
#' @method print metric_graph_vertex
#' @export
print.metric_graph_vertex <- function(x, n = 10, ...) {
  cat("Vertex", attr(x, "id"),"of the metric graph\n\n")
  cat("Longitude and Latitude coordinates:", attr(x, "longlat"), "\n")
  if(attr(x, "longlat")){
    cat("Coordinate reference system:",attr(x, "crs"), "\n")
  }
  if(attr(x, "longlat")){
    lab_x = "Longitude"
    lab_y = "Latitude"
  } else{
    lab_x <- "x"
    lab_y <- "y"
  }
  cat("\nSummary: \n")
  coord_tmp <- matrix(nrow = 1, ncol = 6)
  coord_tmp <- as.data.frame(coord_tmp)
  coord_tmp[1,1:5] <- c(x, attr(x, "degree"),attr(x, "indegree"), attr(x, "outdegree"))
  coord_tmp[1,6] <- attr(x, "problematic")
  colnames(coord_tmp) <- c(lab_x, lab_y, "Degree", "Indegree", "Outdegree", "Problematic")
  print(coord_tmp, row.names = FALSE)
}



#' @noRd 

na.const <- function(x){
  if(!any(is.na(x))){
    return(x)
  }
  not_na <- which(!is.na(x))
  min_nonna <- min(not_na)
  max_nonna <- max(not_na)
  if(min_nonna > 1){
    x[1:(min_nonna-1)] <- x[min_nonna]
  }
  if(max_nonna < length(x)){
    x[(max_nonna+1):length(x)] <- x[max_nonna]
  }
  return(x)
}
