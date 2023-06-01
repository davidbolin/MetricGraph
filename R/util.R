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
  finite_geodist <- is.finite(graph$geo_dist[["__vertices"]])
  finite_geodist <- graph$geo_dist[["__vertices"]][finite_geodist]
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
  new_data[["__edge_number"]] <- PtE[,1]
  new_data[["__distance_on_edge"]] <- PtE[,2]

  if (is.null(group_vector)) {
    if(!is.null(old_data)){
      group_vector <- rep(old_data[["__group"]][1], length(PtE[,1]))
    } else{
      group_vector <- rep(1, length(PtE[,1]))
    }
  }

  if(is.null(old_data)){
      group_val <- unique(group_vector)
  } else {
      group_val <- unique(union(old_data[["__group"]], group_vector))
  }

  if (is.null(old_data)) {
    full_colnames <- names(new_data)
    data_coords_new <- data.frame(PtE1 = PtE[,1], PtE2 = PtE[,2])
    data_coords <- unique(data_coords_new)
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
    data_coords_new[["group"]] <- group_vector
    data_coords[["idx"]] <- 1:nrow(data_coords)
    idx_new_entries <- merge(data_coords_new, data_coords, all=FALSE,
                             sort = FALSE)
    idx_new_entries <- idx_new_entries[["idx"]]
    list_result <- vector(mode = "list", length(full_colnames))
    names(list_result) <- full_colnames
    list_result[1:length(list_result)] <- full_colnames
    list_result[["__edge_number"]] <- NULL
    list_result[["__distance_on_edge"]] <- NULL
    list_result[["__group"]] <- NULL
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
    new_data[["__edge_number"]] <- data_coords[["PtE1"]]
    new_data[["__distance_on_edge"]] <- data_coords[["PtE2"]]
    new_data[["__group"]] <- data_coords[["group"]]
    return(new_data)
  } else {
    old_colnames <- names(old_data)
    new_colnames <- names(new_data)
    full_colnames <- union(old_colnames, new_colnames)

    new_df <- data.frame(PtE1 = PtE[,1], PtE2 = PtE[,2])

    old_df <- data.frame(PtE1 = old_data[["__edge_number"]],
                         PtE2 = old_data[["__distance_on_edge"]])

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
    old_df[["group"]] <- old_data[["__group"]]
    data_coords[["idx"]] <- 1:nrow(data_coords)

    idx_new_entries <- merge(new_df, data_coords, all = FALSE, sort = FALSE)
    idx_new_entries <- idx_new_entries[["idx"]]
    idx_old_entries <- merge(old_df, data_coords, all = FALSE, sort = FALSE)
    idx_old_entries <- idx_old_entries[["idx"]]
    list_result <- vector(mode = "list", length(full_colnames))
    names(list_result) <- full_colnames
    list_result[1:length(list_result)] <- full_colnames
    list_result[["__edge_number"]] <- NULL
    list_result[["__distance_on_edge"]] <- NULL
    list_result[["__group"]] <- NULL
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
  list_result[["__edge_number"]] <- data_coords[["PtE1"]]
  list_result[["__distance_on_edge"]] <- data_coords[["PtE2"]]
  list_result[["__group"]] <- data_coords[["group"]]
  return(list_result)
  }
}

#' find indices of the rows with all NA's in lists
#' @noRd
#'
idx_not_all_NA <- function(data_list){
     data_list[["__edge_number"]] <- NULL
     data_list[["__distance_on_edge"]] <- NULL
     data_list[["__coord_x"]] <- NULL
     data_list[["__coord_y"]] <- NULL
     data_list[["__group"]] <- NULL
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

#' Select replicate
#' @noRd
#'
select_group <- function(data_list, group){
    grp <- data_list[["__group"]]
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
  line1 <- Line(cbind(1+sin(theta),2+2*cos(theta)))

  theta <- seq(from=pi/2,to=pi,length.out = n)
  line2 <- Line(cbind(1+sin(theta),1.5+1.5*cos(theta)))

  theta <- seq(from=3*pi/2,to=2*pi,length.out = n)
  line3 <- Line(cbind(2+2*sin(theta),2+2*cos(theta)))

  line4 <- Line(rbind(c(1,1.5),c(2,1.5)))

  #R
  line5 <- Line(rbind(c(2,0),c(2,4)))
  line6 <- Line(rbind(c(2,4),c(3,4)))
  theta <- seq(from=0,to=pi,length.out = n)
  line7 <- Line(cbind(3+sin(theta),3+cos(theta)))
  line8 <- Line(rbind(c(3,2),c(2,2)))
  line9 <- Line(rbind(c(2,2),c(4,0)))

  #A
  line10 <- Line(rbind(c(4,0),c(5,4)))
  line11 <- Line(rbind(c(5,4),c(6,0)))
  line12 <- Line(rbind(c(4.5,2),c(5.5,2)))

  #P
  line13 <- Line(rbind(c(6,0),c(6,4)))
  line14 <- Line(rbind(c(6,4),c(7,4)))
  theta <- seq(from=0,to=pi,length.out = n)
  line15 <- Line(cbind(7+sin(theta),3+cos(theta)))
  line16 <- Line(rbind(c(7,2),c(6,2)))

  #H
  line17 <- Line(rbind(c(8,0),c(8,4)))
  line18 <- Line(rbind(c(10,0),c(10,4)))
  line19 <- Line(rbind(c(8,2),c(10,2)))

  #M
  line20 <- Line(rbind(c(0,4),c(0.75,8)))
  line21 <- Line(rbind(c(0.75,8),c(1.5,5)))
  line22 <- Line(rbind(c(1.5,5),c(2.25,8)))
  line23 <- Line(rbind(c(2.25,8),c(3,4)))

  # E
  line24 <- Line(rbind(c(3,4),c(3,8)))
  line25 <- Line(rbind(c(3,8),c(4,8)))
  line26 <- Line(rbind(c(3,6),c(4,6)))
  line27 <- Line(rbind(c(3,4),c(5,4)))

  # T
  line28 <- Line(rbind(c(5,4),c(5,8)))
  line29 <- Line(rbind(c(4,8),c(6,8)))


  # R
  line30 <- Line(rbind(c(6,4),c(6,8)))
  line31 <- Line(rbind(c(6,8),c(7,8)))
  theta <- seq(from=0,to=pi,length.out = n)
  line32 <- Line(cbind(7+sin(theta),7+cos(theta)))
  line33 <- Line(rbind(c(7,6),c(6,6)))
  line34 <- Line(rbind(c(6,6),c(8,4)))

  # I
  line35 <- Line(rbind(c(8,4),c(8,8)))

  # C
  theta <- seq(from=pi,to=3*pi/2,length.out = n)
  line36 <- Line(cbind(10+2*sin(theta),6+2*cos(theta)))
  theta <- seq(from=3*pi/2,to=2*pi,length.out = n)
  line37 <- Line(cbind(10+2*sin(theta),6+2*cos(theta)))

  return(sp::SpatialLines(list(Lines(list(line1),ID="1"),
                               Lines(list(line2),ID="2"),
                               Lines(list(line3),ID="3"),
                               Lines(list(line4),ID="4"),
                               Lines(list(line5),ID="5"),
                               Lines(list(line6),ID="6"),
                               Lines(list(line7),ID="7"),
                               Lines(list(line8),ID="8"),
                               Lines(list(line9),ID="9"),
                               Lines(list(line10),ID="10"),
                               Lines(list(line11),ID="11"),
                               Lines(list(line12),ID="12"),
                               Lines(list(line13),ID="13"),
                               Lines(list(line14),ID="14"),
                               Lines(list(line15),ID="15"),
                               Lines(list(line16),ID="16"),
                               Lines(list(line17),ID="17"),
                               Lines(list(line18),ID="18"),
                               Lines(list(line19),ID="19"),
                               Lines(list(line20),ID="20"),
                               Lines(list(line21),ID="21"),
                               Lines(list(line22),ID="22"),
                               Lines(list(line23),ID="23"),
                               Lines(list(line24),ID="24"),
                               Lines(list(line25),ID="25"),
                               Lines(list(line26),ID="26"),
                               Lines(list(line27),ID="27"),
                               Lines(list(line28),ID="28"),
                               Lines(list(line29),ID="29"),
                               Lines(list(line30),ID="30"),
                               Lines(list(line31),ID="31"),
                               Lines(list(line32),ID="32"),
                               Lines(list(line33),ID="33"),
                               Lines(list(line34),ID="34"),
                               Lines(list(line35),ID="35"),
                               Lines(list(line36),ID="36"),
                               Lines(list(line37),ID="37"))))
}



#' @noRd

projectVecLine2 <- function(lines, points, normalized = FALSE){
  lines <- lines@lines[[1]]@Lines[[1]]@coords
  points <- points@coords
  return(projectVecLine(lines, points, normalized))
}

#' @noRd

distance2 <- function(points, lines, byid=FALSE){
  rownames(points@coords) <- 1:nrow(points@coords)
  points_sf <- sf::st_as_sf(points)
  lines_sf <- sf::st_as_sf(lines)
  dist_result <- sf::st_distance(points_sf, lines_sf)
  if(byid){
    return(matrix(apply(dist_result,1,min),nrow=1))
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
  return(inter_lines)
}

#' @noRd

interpolate2 <- function(lines, pos, normalized = FALSE){
    lines <- lines@lines[[1]]@Lines[[1]]@coords
    return(interpolate2_aux(lines, pos, normalized))
}

#' @noRd

make_Aprd <- function(graph, edge_number, distance_on_edge){
  X_loc <- cbind(edge_number, distance_on_edge)
  order_X <- order(X_loc[,1], X_loc[,2])
  X_loc <- X_loc[order_X,]
  edge_n <- sort(unique(edge_number))
  idx_tmp <- (graph$data[["__edge_number"]] == edge_n[1])
  mesh_tmp <- graph$data[["__distance_on_edge"]][idx_tmp]
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
