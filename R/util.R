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
        stop('nu not yet implemented')
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
        stop('nu not yet implemented')
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
        stop('nu not yet implemented')
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
        stop('nu not yet implemented')
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
        stop('nu not yet implemented')
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

        stop('nu not yet implemented')
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
        stop('nu not yet implemented')
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
#' @param rec_tau Should a starting value for the reciprocal of tau be given?
#' @param model_options List object containing the model options.
#' @param factor_start_range Factor to multiply the max/min/diagonal dimension of the bounding box to obtain a starting value for range. Default is 0.5.
#' @param type_start_range_bbox Which dimension from the bounding box should be used? The options are 'diag', the default, 'max' and 'min'.

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
                                  log_scale = FALSE,
                                  model_options = list(),
                                  rec_tau = TRUE,
                                  factor_start_range = 0.3,
                                  type_start_range_bbox = "diag"){

  check_graph(graph)

  type_start_range_bbox <- match.arg(type_start_range_bbox, 
                                   choices = c("diag", "max", "min"))

  


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

  if(is.null(model_options$start_range)){
    bounding_box <- graph$get_bounding_box(format = "sf")

    # Check if the bounding box object inherits from "bbox" (sf format)
  if (inherits(bounding_box, "bbox")) {
      # Extract coordinates from the bounding box
      min_x <- bounding_box["xmin"]
      max_x <- bounding_box["xmax"]
      min_y <- bounding_box["ymin"]
      max_y <- bounding_box["ymax"]

      # Create points in sf format with the appropriate CRS
      point_min_x <- sf::st_point(c(min_x, min_y)) |> sf::st_sfc(crs = sf::st_crs(bounding_box))
      point_max_x <- sf::st_point(c(max_x, min_y)) |> sf::st_sfc(crs = sf::st_crs(bounding_box))
      point_min_y <- sf::st_point(c(min_x, min_y)) |> sf::st_sfc(crs = sf::st_crs(bounding_box))
      point_max_y <- sf::st_point(c(min_x, max_y)) |> sf::st_sfc(crs = sf::st_crs(bounding_box))

      # Calculate the width and height
      width <- sf::st_distance(point_min_x, point_max_x)
      height <- sf::st_distance(point_min_y, point_max_y)

      if(type_start_range_bbox == "diag") {
        dimension_size <- sqrt(as.numeric(width)^2 + as.numeric(height)^2)/1000
      } else if(type_start_range_bbox == "max") {
        dimension_size <- max(as.numeric(width), as.numeric(height))/1000
      } else { # min case
        dimension_size <- min(as.numeric(width), as.numeric(height))/1000
      }

    } else {
      # If not sf format, assume it’s a standard list and compute Euclidean distances
      min_x <- bounding_box$min_x
      max_x <- bounding_box$max_x
      min_y <- bounding_box$min_y
      max_y <- bounding_box$max_y

      width <- max_x - min_x
      height <- max_y - min_y

      dimension_size <- switch(type_start_range_bbox,
                              "diag" = sqrt(width^2 + height^2),
                              "max" = max(width, height),
                              "min" = min(width, height))

    }

    prior.range.nominal <- dimension_size * factor_start_range
  } else{
    prior.range.nominal <- model_options$start_range
  }

  if(!is.null(model_options$fix_range)){
    prior.range.nominal <- model_options$fix_range
  }

  start_sigma <- NULL

  if(nu){
      start_sigma <- 1
      start_nu <- 1
      if(!is.null(model_options$start_nu)){
        start_nu <- model_options$start_nu
      }

  }

  if(!is.null(model_options$start_sigma)){
      start_sigma <- model_options$start_sigma
  }       

  if(!is.null(model_options[["fix_sigma"]])){
    start_sigma <- model_options[["fix_sigma"]]
  }

  if(!is.null(model_options[["fix_sigma_e"]])){
    start_sigma_e <- model_options[["fix_sigma_e"]]
  }

  if(!is.null(model_options$fix_nu)){
    start_nu <- model_options$fix_nu
  }



  if (model == "alpha1") {
    start_kappa <- sqrt(8 * 0.5) / prior.range.nominal
    #variance is sigma^2/2 kappa
    if(is.null(start_sigma)){
      if(data){
        start_sigma <- sqrt(2*start_kappa) * data_std
      } else{
        start_sigma <- 1
      }
    } 
    nu_tmp <- 0.5
    start_tau <- sqrt(gamma(nu_tmp) / (start_sigma^2 * start_kappa^(2 * nu_tmp) * (4 * pi)^(1 / 2) * gamma(nu_tmp + 1 / 2)))

  if(!is.null(model_options$start_tau)){
      start_tau <- model_options$start_tau
  }       

  if(!is.null(model_options$fix_tau)){
    start_tau <- model_options$fix_tau
  }

  if(!is.null(model_options$start_kappa)){
      start_kappa <- model_options$start_kappa
  }       

  if(!is.null(model_options$fix_kappa)){
    start_kappa <- model_options$fix_kappa
  }

  } else if (model == "alpha2") {
    start_kappa <- sqrt(8 * 1.5) / prior.range.nominal
    if(is.null(start_sigma)){
      if(data){
        #variance is sigma^2/(4 * kappa^3)
        # multiplying by 2 to help stabilize.
        start_sigma <- 2 * sqrt(4*start_kappa^3) * data_std
      } else{
        start_sigma <- 1
      }
    }
    nu_tmp <- 1.5
    start_tau <- sqrt(gamma(nu_tmp) / (start_sigma^2 * start_kappa^(2 * nu_tmp) * (4 * pi)^(1 / 2) * gamma(nu_tmp + 1 / 2)))   

  if(!is.null(model_options$start_tau)){
      start_tau <- model_options$start_tau
  }       

  if(!is.null(model_options$fix_tau)){
    start_tau <- model_options$fix_tau
  }

  if(!is.null(model_options$start_kappa)){
      start_kappa <- model_options$start_kappa
  }       

  if(!is.null(model_options$fix_kappa)){
    start_kappa <- model_options$fix_kappa
  }  

  } else if (model == "isoExp") {
    start_kappa <- sqrt(8 * 0.5) / prior.range.nominal
    if(is.null(start_sigma)){
      if(data){
        start_sigma <- data_std
      } else{
        start_sigma <- 1
      }
    }
    nu_tmp <- 0.5
    start_tau <- sqrt(gamma(nu_tmp) / (start_sigma^2 * start_kappa^(2 * nu_tmp) * (4 * pi)^(1 / 2) * gamma(nu_tmp + 1 / 2)))    

  if(!is.null(model_options$start_tau)){
      start_tau <- model_options$start_tau
  }       

  if(!is.null(model_options$fix_tau)){
    start_tau <- model_options$fix_tau
  }

  if(!is.null(model_options$start_kappa)){
      start_kappa <- model_options$start_kappa
  }       

  if(!is.null(model_options$fix_kappa)){
    start_kappa <- model_options$fix_kappa
  }  


  } else if (model == "GL1") {
    if(is.null(graph$Laplacian)) {
      graph$compute_laplacian()
    }
    h <- mean(graph$edge_lengths)
    k <- sqrt(8 * 0.5) / prior.range.nominal
    start_kappa <- exp(-k*h)/(1-exp(-2*k*h)) + 2*exp(-k*h) - 2

    if(is.null(start_sigma)){
      if(data){
        Q <- start_kappa^2*Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]]
        v <- rep(0,graph$nV)
        v[1] <- 1
        s2 <- solve(Q,v)[1]
        start_sigma <- data_std / sqrt(s2)
      } else{
        start_sigma <- 1
      }
    }
    nu_tmp <- 0.5
    start_tau <- sqrt(gamma(nu_tmp) / (start_sigma^2 * start_kappa^(2 * nu_tmp) * (4 * pi)^(1 / 2) * gamma(nu_tmp + 1 / 2)))    

  if(!is.null(model_options$start_tau)){
      start_tau <- model_options$start_tau
  }       

  if(!is.null(model_options$fix_tau)){
    start_tau <- model_options$fix_tau
  }

  if(!is.null(model_options$start_kappa)){
      start_kappa <- model_options$start_kappa
  }       

  if(!is.null(model_options$fix_kappa)){
    start_kappa <- model_options$fix_kappa
  }  


  } else if (model == "GL2") {
    if(is.null(graph$Laplacian)) {
      graph$compute_laplacian()
    }
    h <- mean(graph$edge_lengths)
    k <- sqrt(8 * 1.5) / prior.range.nominal
    start_kappa <- exp(-k*h)/(1-exp(-2*k*h)) + 2*exp(-k*h) - 2
    if(is.null(start_sigma)){
      if(data){
        Q <- start_kappa^2*Matrix::Diagonal(graph$nV, 1) + graph$Laplacian[[1]]
        v <- rep(0,graph$nV)
        v[1] <- 1
        s2 <- solve(Q %*% Q,v)[1]
        start_sigma <- data_std / sqrt(s2)
      } else{
        start_sigma <- 1
      }
    }
    nu_tmp <- 1.5
    start_tau <- sqrt(gamma(nu_tmp) / (start_sigma^2 * start_kappa^(2 * nu_tmp) * (4 * pi)^(1 / 2) * gamma(nu_tmp + 1 / 2)))    

  if(!is.null(model_options$start_tau)){
      start_tau <- model_options$start_tau
  }       

  if(!is.null(model_options$fix_tau)){
    start_tau <- model_options$fix_tau
  }

  if(!is.null(model_options$start_kappa)){
      start_kappa <- model_options$start_kappa
  }       

  if(!is.null(model_options$fix_kappa)){
    start_kappa <- model_options$fix_kappa
  }  


  } else {
    stop("wrong model choice")
  } 

  out_vec <- c()

  # reciprocal tau 

  if(like_format){
      if(is.null(model_options[["fix_sigma"]]) && is.null(model_options$fix_tau)){
        if(rec_tau){
          out_vec <- 1/start_tau
        } else{
          out_vec <- start_tau
        }
      }

      if(is.null(model_options$fix_range) && is.null(model_options$fix_kappa)){
        if(range_par){
          out_vec <- c(out_vec, prior.range.nominal)
        } else{
          out_vec <- c(out_vec, start_kappa)
        }
      }

      if(is.null(model_options$fix_nu)){
        if(nu){
          out_vec <- c(out_vec,start_nu)
        }
      }
      
      if(is.null(model_options[["fix_sigma_e"]])){
        if(is.null(model_options$start_sigma_e)){
          out_vec <- c(out_vec, 0.1 * data_std)
        } else{
          out_vec <- c(out_vec, model_options$start_sigma_e)
        }
      }
  } else{
      if(is.null(model_options[["fix_sigma_e"]])){
        if(is.null(model_options$start_sigma_e)){
          out_vec <- c(out_vec, 0.1 * data_std)
        } else{
          out_vec <- c(out_vec, model_options$start_sigma_e)
        }
      }

      if(is.null(model_options[["fix_sigma"]]) && is.null(model_options$fix_tau)){
        out_vec <- c(out_vec, 1/start_tau)
      }

      if(is.null(model_options$fix_range) && is.null(model_options$fix_kappa)){
        if(range_par){
          out_vec <- c(out_vec, prior.range.nominal)
        } else{
          out_vec <- c(out_vec, start_kappa)
        }
      }

      if(is.null(model_options$fix_nu)){
        if(nu){
          out_vec <- c(out_vec,start_nu)
        }
      }
  }

  out_fixed <- list()

  if(!is.null(model_options[["fix_sigma_e"]])){
    out_fixed <- c(out_fixed, fixed_sigma_e = start_sigma_e)
  }

  if(!is.null(model_options[["fix_sigma"]]) || !is.null(model_options$fix_tau)){
    if(rec_tau){
      out_fixed <- c(out_fixed, fixed_tau = 1/start_tau)
    } else{
      out_fixed <- c(out_fixed, fixed_tau = start_tau)
    }
  }

  if(!is.null(model_options$fix_range) || !is.null(model_options$fix_kappa)){
    out_fixed <- c(out_fixed, fixed_kappa = start_kappa)
  }

  if(!is.null(model_options$fix_nu)){
    out_fixed <- c(out_fixed, fixed_nu = start_nu)
  }

  if(log_scale){
    if(length(out_vec)>1){
      out_vec <- log(out_vec)
    }
    out_fixed <- lapply(out_fixed, log)  
  }

  out_list <- list(start_values = out_vec, fixed_values = out_fixed)

  return(out_list)
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
process_data_add_obs <- function(PtE, new_data, old_data, group_vector, suppress_warnings) {
  new_data[[".edge_number"]] <- PtE[, 1]
  new_data[[".distance_on_edge"]] <- PtE[, 2]

  # Ensure group_vector is initialized correctly
  if (is.null(group_vector)) {
    group_vector <- if (!is.null(old_data)) {
      rep(old_data[[".group"]][1], nrow(PtE))
    } else {
      rep(1, nrow(PtE))
    }
  }

  # Get unique group values
  group_val <- if (is.null(old_data)) {
    unique(group_vector)
  } else {
    unique(c(old_data[[".group"]], group_vector))
  }

  # Combine and order coordinates
  data_coords <- unique(rbind(
    data.frame(PtE1 = PtE[, 1], PtE2 = PtE[, 2]),
    if (!is.null(old_data)) data.frame(PtE1 = old_data[[".edge_number"]], PtE2 = old_data[[".distance_on_edge"]]) else NULL
  ))
  data_coords <- data_coords[order(data_coords$PtE1, data_coords$PtE2), ]

  # Expand coordinates for groups
  n_group <- length(group_val)
  data_coords <- data.frame(
    PtE1 = rep(data_coords$PtE1, n_group),
    PtE2 = rep(data_coords$PtE2, n_group),
    group = rep(group_val, each = nrow(data_coords))
  )
  data_coords[["idx"]] <- seq_len(nrow(data_coords))

  # Map new and old data to indices
  idx_new_entries <- merge(
    data.frame(PtE1 = PtE[, 1], PtE2 = PtE[, 2], group = group_vector),
    data_coords, by = c("PtE1", "PtE2", "group"), sort = FALSE
  )[["idx"]]

  if (!is.null(old_data)) {
    idx_old_entries <- merge(
      data.frame(PtE1 = old_data[[".edge_number"]], PtE2 = old_data[[".distance_on_edge"]], group = old_data[[".group"]]),
      data_coords, by = c("PtE1", "PtE2", "group"), sort = FALSE
    )[["idx"]]
  } else {
    idx_old_entries <- integer(0)
  }

  # Warn about conflicts if necessary
  if (!suppress_warnings && length(intersect(idx_old_entries, idx_new_entries)) > 0) {
    warning("Conflicting data detected. New data may overwrite existing data.")
  }

  # Combine data, prioritizing non-NA values from new_data
  list_result <- lapply(union(names(old_data), names(new_data)), function(col_name) {
    tmp <- rep(NA, nrow(data_coords))
    
    # Prioritize new_data values
    if (!is.null(new_data[[col_name]])) {
      tmp[idx_new_entries] <- new_data[[col_name]]
    }
    
    # Fill missing values with old_data where available
    if (!is.null(old_data[[col_name]])) {
      na_idx <- is.na(tmp[idx_old_entries])  # Check missing in new_data
      tmp[idx_old_entries[na_idx]] <- old_data[[col_name]][na_idx]
    }
    
    tmp
  })

  names(list_result) <- union(names(old_data), names(new_data))
  list_result[[".edge_number"]] <- data_coords$PtE1
  list_result[[".distance_on_edge"]] <- data_coords$PtE2
  list_result[[".group"]] <- data_coords$group

  return(list_result)
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

change_parameterization_graphlme <- function(#likelihood, 
nu, par, hessian, fix_vec
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

  if(all(!fix_vec)){
    new_observed_fisher <- t(grad_par[!fix_vec,!fix_vec]) %*% hessian[!fix_vec,!fix_vec] %*% (grad_par[!fix_vec,!fix_vec])
  } else{
    new_observed_fisher <- grad_par[!fix_vec,!fix_vec] * hessian[!fix_vec,!fix_vec] * (grad_par[!fix_vec,!fix_vec])
  }

  # No need to include the additional term as the gradient is approximately zero.
  # from some numerical experiments, the approximation without the additional term
  # seems to be better in general.

  inv_fisher <- tryCatch(solve(new_observed_fisher),
                         error = function(e) matrix(NA, nrow(new_observed_fisher),
                                                    ncol(new_observed_fisher)))

  std_err <- sqrt(diag(inv_fisher))

  std_err_tmp <- c(NA,NA)
  std_err_tmp[!fix_vec] <- std_err

  return(list(coeff = c(sigma, range), std_random = std_err_tmp))
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
  } else if(vertex_unit == "degree"){
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

compute_line_lengths <- function(edge, longlat, unit, crs, proj4string, which_longlat, vertex_unit, project_data, transform){
  if(!is.null(edge)){
      class(edge) <- setdiff(class(edge), "metric_graph_edge")
  }
    if(!longlat || project_data){
      fact <- process_factor_unit(vertex_unit, unit)
      return(compute_length(edge) * fact)
    } else if(which_longlat == "sf"){
      if(!is.null(edge)){
        if(!is.null(crs)){
          fact <- 1
        }
        linestring <- sf::st_sfc(sf::st_linestring(edge), crs = crs)
        # linestring <- sf::st_transform(linestring,  crs = 4326)        
        length <- sf::st_length(linestring)
        units(length) <- unit
        units(length) <- NULL
        return(length)
      } else{
        return(0)
      }
    } else{
      # Line <- sp::Line(edge)
      # Line <- sp::Lines(Line, ID = 1)
      # Line <- sp::SpatialLines(list(Line), proj4string =  proj4string)
      # Line <- sp::spTransform(Line, CRSobj = sp::CRS("+proj=longlat +datum=WGS84"))
      # Line <- sp::Line(sp::coordinates(Line))
      if(!transform){
        length <- sp::LineLength(edge, longlat = longlat)
      } else{
        Line <- sf::st_as_sf(as.data.frame(edge), coords = 1:2, crs = crs)
        Line <- sf::st_transform(Line, crs = 4326)
        Line <- sf::st_coordinates(Line) 
        length <- sp::LineLength(Line, longlat = longlat)
        units(length) <- "km"
        fact <- 1
        units(length) <- unit
        units(length) <- NULL
        return(length)
      }

      fact <- process_factor_unit(vertex_unit, unit)
      return(length * fact)
    }
}


#' @noRd 
#' 

compute_aux_distances <- function(lines, crs, longlat, proj4string, points = NULL, fact, which_longlat, length_unit, transform){
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
        if(transform){
          sf_points <- sf::st_transform(sf_points,  crs = 4326)
        }
        if(is.null(points)){
          dists <- sf::st_distance(sf_points, which = "Great Circle")
        } else{
          sf_p_points <- sf::st_as_sf(as.data.frame(points), coords = 1:2, crs = crs)
          if(transform){
            sf_p_points <- sf::st_transform(sf_p_points,  crs = 4326)
          }
          dists <- sf::st_distance(x = sf_points, y = sf_p_points, which = "Great Circle", by_element = TRUE)
        }

        units(dists) <- length_unit
        units(dists) <- NULL
    } else{
        sp_points <- sp::SpatialPoints(coords = lines, proj4string = proj4string) 
        if(transform){
          sp_points <- sp::spTransform(sp_points, CRSobj = sp::CRS("+proj=longlat +datum=WGS84"))
        }
        if(is.null(points)){
          dists <- sp::spDists(sp_points, longlat = TRUE) #* fact
        } else{
          sp_p_points <- sp::SpatialPoints(coords = points, proj4string = proj4string) 
          if(transform){
            sp_p_points <- sp::spTransform(sp_p_points, CRSobj = sp::CRS("+proj=longlat +datum=WGS84"))          
          }
          dists <- sp::spDists(x = sp_points, y=sp_p_points, longlat = TRUE, diagonal = TRUE) #* fact
        }
        units(dists) <- "km"
        units(dists) <- length_unit
        units(dists) <- NULL
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
#' @param compute_characteristics Should the characteristics of the graph be computed? If `NULL` it will be determined based on the size of the graph.
#' @param check_euclidean Check if the graph has Euclidean edges? If `NULL` it will be determined based on the size of the graph.
#' @param check_distance_consistency Check the distance consistency assumption?#' If `NULL` it will be determined based on the size of the graph.
#' @param ... not used.
#' @return An object of class \code{summary_graph_lme} containing information
#' about a *metric_graph* object.
#' @method summary metric_graph
#' @export
summary.metric_graph <- function(object, messages = FALSE, compute_characteristics = NULL, check_euclidean = NULL, check_distance_consistency = NULL, ...){
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
    if(!is.null(attr(x[[i]], "kirchhoff_weight"))){
      kw <- attr(x[[i]], "kirchhoff_weight")
      w_tmp <- attr(x[[i]], "weight")
      if(is.data.frame(w_tmp)){
        cat("Kirchhoff weight:", w_tmp[[kw]],"\n\n")
      } else{
        cat("Kirchhoff weight:", w_tmp,"\n\n")
      }
    }
    
    if(!is.null(attr(x[[i]], "directional_weight"))){
      dw <- attr(x[[i]], "directional_weight")
      w_tmp <- attr(x[[i]], "weight")
      if(is.data.frame(w_tmp)){
        cat("Directional weight:", w_tmp[[dw]],"\n\n")
      } else{
        cat("Directional weight:", w_tmp,"\n\n")
      }
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
  PtE <- cbind(attr(x, "id"), PtE)
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
    if(!is.null(attr(x, "kirchhoff_weight"))){
      kw <- attr(x, "kirchhoff_weight")
      w_tmp <- attr(x, "weight")
      if(is.data.frame(w_tmp)){
        cat("Kirchhoff weight:", w_tmp[[kw]],"\n\n")
      } else{
        cat("Kirchhoff weight:", w_tmp,"\n\n")
      }
    }    
    if(!is.null(attr(x, "directional_weight"))){
      dw <- attr(x, "directional_weight")
      w_tmp <- attr(x, "weight")
      if(is.data.frame(w_tmp)){
        cat("Directional weight:", w_tmp[[dw]],"\n\n")
      } else{
        cat("Directional weight:", w_tmp,"\n\n")
      }
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

# na.const <- function(x){
#   if(!any(is.na(x))){
#     return(x)
#   }
#   not_na <- which(!is.na(x))
#   min_nonna <- min(not_na)
#   max_nonna <- max(not_na)
#   if(min_nonna > 1){
#     x[1:(min_nonna-1)] <- x[min_nonna]
#   }
#   if(max_nonna < length(x)){
#     x[(max_nonna+1):length(x)] <- x[max_nonna]
#   }
#   return(x)
# }

na.const <- function(x) {
  # If no NA values, return as-is
  if(!any(is.na(x))) {
    return(x)
  }
  
  # Check if all values are NA
  if(all(is.na(x))) {
    return(x)  # or return a default value depending on your needs
  }
  
  not_na <- which(!is.na(x))
  min_nonna <- min(not_na)
  max_nonna <- max(not_na)
  
  # Safety check for vector creation
  if(min_nonna > 1 && (min_nonna - 1) < .Machine$integer.max) {
    x[1:(min_nonna-1)] <- x[min_nonna]
  }
  
  if(max_nonna < length(x)) {
    x[(max_nonna+1):length(x)] <- x[max_nonna]
  }
  
  return(x)
}


#' @noRd 

# Function factory to fix some variables, and return a new function to be passed to the likelihood
# func -> original function
# num_var -> dimension of the original parameter vector
# fix_vec -> boolean vector of size num_var with the variables to be fixed
# fix_val -> values to be fixed
# n_cov -> number of covariates

function_factory_fix_var <- function(func, fix_vec, num_var, fix_val, n_cov) {
  ret_fun <- function(theta){
    new_theta <- numeric(num_var)
    fix_vec_latent <- c(fix_vec, rep(FALSE,n_cov))
    new_theta[fix_vec_latent] <- fix_val
    new_theta[!fix_vec_latent] <- theta
    return(func(new_theta))
  }
  return(ret_fun)
}

#' @noRd 

# Get the starting values to be passed to optim when fixing variables

# start_values -> original starting values
# fix_vec -> boolean vec

start_values_fix <- function(start_values, fix_vec, n_cov){
  fix_vec_latent <- c(fix_vec, rep(FALSE,n_cov))
  return(start_values[!fix_vec_latent])
}

#' @noRd 

get_fixed_values <- function(start_values, fix_vec, n_cov){
  fix_vec_latent <- c(fix_vec, rep(FALSE,n_cov))
    return(start_values[!fix_vec_latent])
}


#' @noRd 

create_fix_vec_val <- function(fixed_values){
    if(is.null(fixed_values$fixed_sigma_e)){
      fix_vec <- FALSE
      fix_v_val <- NA
    } else{
      fix_vec <- TRUE
      fix_v_val <- fixed_values$fixed_sigma_e
    }

    if(is.null(fixed_values$fixed_tau)){
      fix_vec <- c(fix_vec, FALSE)
      fix_v_val <- c(fix_v_val, NA)
    } else{
      fix_vec <- c(fix_vec,TRUE)
      fix_v_val <- c(fix_v_val, fixed_values$fixed_tau)
    }

    if(is.null(fixed_values$fixed_kappa)){
      fix_vec <- c(fix_vec, FALSE)
      fix_v_val <- c(fix_v_val, NA)
    } else{
      fix_vec <- c(fix_vec,TRUE)
      fix_v_val <- c(fix_v_val, fixed_values$fixed_kappa)
    }

  return(list(fix_vec = fix_vec, fix_v_val = fix_v_val))
}


#' @noRd 

check_model_options <- function(model_options){
  if(length(model_options[["fix_tau"]]) > 1){
    stop("'fix_tau' must have length 1!")
  }
  if(length(model_options[["fix_sigma"]]) > 1){
    stop("'fix_sigma' must have length 1!")
  }
  if(length(model_options[["fix_sigma_e"]]) > 1){
    stop("'fix_sigma_e' must have length 1!")
  }
  if(length(model_options[["fix_kappa"]]) > 1){
    stop("'fix_kappa' must have length 1!")
  }  
  if(length(model_options[["fix_range"]]) > 1){
    stop("'fix_range' must have length 1!")
  }
  if(length(model_options[["fix_nu"]]) > 1){
    stop("'fix_nu' must have length 1!")
  }
  if(!is.null(model_options[["fix_sigma"]])){
    if(model_options[["fix_sigma"]] <= 0){
      stop("'fix_sigma' must be positive!")
    }
  }
  if(!is.null(model_options[["fix_tau"]])){
    if(model_options[["fix_tau"]] <= 0){
      stop("'fix_tau' must be positive!")
    }
  }  
  if(!is.null(model_options[["fix_kappa"]])){
    if(model_options[["fix_kappa"]] <= 0){
      stop("'fix_kappa' must be positive!")
    }
  }
  if(!is.null(model_options[["fix_range"]])){
    if(model_options[["fix_range"]] <= 0){
      stop("'fix_range' must be positive!")
    }
  }
  if(!is.null(model_options[["fix_nu"]])){
    if(model_options[["fix_nu"]] <= 0){
      stop("'fix_nu' must be positive!")
    }
  }  
  if(!is.null(model_options[["fix_sigma_e"]])){
    if(model_options[["fix_sigma_e"]] < 0){
      stop("'fix_sigma_e' must be non-negative!")
    }
  }
}

#' @noRd 

get_only_first <- function(vec){
  idx <- which(vec)
  vec <- rep(FALSE, length(vec))
  vec[idx[1]] <- TRUE
  return(vec)
}


#' @noRd 
# Create a map from vertices into reference edges

map_into_reference_edge <- function(graph, verbose = 0) {
  if (verbose == 2) {
    message("Creating a map from vertices into reference edges")
  }

  # Initialize the reference edge matrix
  ref_edge <- matrix(nrow = graph$nV, ncol = 2)

  idx_pos_0 <- match(seq_len(graph$nV), graph$E[, 1], nomatch = 0)
  idx_pos_1 <- match(seq_len(graph$nV), graph$E[, 2], nomatch = 0)

  # Determine which vertices map to the first column of graph$E
  ref_edge[, 1] <- ifelse(idx_pos_0 > 0, idx_pos_0, idx_pos_1)
  ref_edge[, 2] <- ifelse(idx_pos_0 > 0, 0, 1)

  return(ref_edge)
}

#' @noRd 
# Converts distance on edge equal 1 to distance on edge equal to 0

standardize_df_positions <- function(df, graph, edge_number = "edge_number", distance_on_edge = "distance_on_edge"){
  idx_pos1 <- which(df[[distance_on_edge]] == 1)
  idx_pos0 <- which(df[[distance_on_edge]] == 0)
  if(length(idx_pos1) + length(idx_pos0) == 0){
    return(df)
  }

  ref_edges <- graph$.__enclos_env__$private$ref_edges
  
  if(length(idx_pos1)>0){
    edge_num_pos1 <- df[[edge_number]][idx_pos1]
    vertices_pos1 <- graph$E[edge_num_pos1, 2]
    df[[edge_number]][idx_pos1] <- ref_edges[vertices_pos1,1]
    df[[distance_on_edge]][idx_pos1] <- ref_edges[vertices_pos1,2]
  }

  if(length(idx_pos0)>0){
    edge_num_pos0 <- df[[edge_number]][idx_pos0]
    vertices_pos0 <- graph$E[edge_num_pos0, 1]
    df[[edge_number]][idx_pos0] <- ref_edges[vertices_pos0,1]
    df[[distance_on_edge]][idx_pos0] <- ref_edges[vertices_pos0,2]
  }

  return(df)
}

#' Helper function to be used in the split_edge method
#' @noRd
fill_na_values_split_edge <- function(data) {
  # Sort by pos_edge to ensure interpolation occurs in order
  data <- data[order(data$pos_edge), ]

  # Perform linear interpolation for x and y
  data$x <- approx(data$pos_edge, data$x, data$pos_edge, method = "linear", rule = 2, ties = mean)$y
  data$y <- approx(data$pos_edge, data$y, data$pos_edge, method = "linear", rule = 2, ties = mean)$y

  # Filter to remove duplicates, keeping rows where is_t_values is TRUE or pos_edge is unique
  unique_rows <- !duplicated(data$pos_edge) | data$is_t_values
  data <- data[unique_rows & !is.na(data$pos_edge), ]

  return(data)
}


# Helper function for merging observations to be used with `add_observations()`
# strategies "remove", "merge", "average"
#' @noRd 

get_idx_within_merge_tolerance <- function(PtE, group_vector, aux_length, tolerance, dplyr = FALSE){
  if(is.null(group_vector)){
    group_vector <- rep(1, length(PtE[,1]))
  }

  if(!dplyr){
      # Loop over unique groups and edges
      selected_rows <- unlist(lapply(unique(group_vector), function(group) {
          # Get indices for the current primary group
          group_indices <- which(group_vector == group)
          # Further split by unique edges within this group
          unique_edges <- unique(PtE[group_indices, 1])
          unlist(lapply(unique_edges, function(edge) {
              # Get indices of rows for the current edge within the current group
              edge_indices <- group_indices[PtE[group_indices, 1] == edge]
              # Apply the filtering function
              filter_indices_by_tolerance(edge_indices, PtE[edge_indices, 2], tolerance)
          }))
      }))
  } else{
      group <- index <- edge <- NULL
      PtE_df <- as.data.frame(PtE)
      PtE_df$group <- group_vector
      colnames(PtE_df) <- c("edge", "position", "group")
      # Add an index column to keep track of original row numbers
      PtE_df$orig_index <- 1:nrow(PtE_df)
      # Group by both `group` and `edge` and apply the filtering function
      selected_indices <- PtE_df |>
          dplyr::group_by(group, edge) |>
          dplyr::group_modify(~ dplyr::tibble(index = filter_group_edge(.x, tolerance))) |>
          dplyr::pull(index)
  }
}


# Function to iteratively filter rows within a single group-edge subgroup
#' @noRd
filter_indices_by_tolerance <- function(edge_indices, positions, tolerance) {
    selected <- edge_indices  # Start with all indices in the subgroup

    while (TRUE) {
        diffs <- diff(positions[selected - min(selected) + 1])  # Calculate diffs on the current selection
        below_tolerance <- which(diffs < tolerance)
        
        if (length(below_tolerance) == 0) {
            # Stop if no diffs are below tolerance
            break
        }

        # Remove the first occurrence where diff < tolerance
        selected <- selected[-(below_tolerance[1] + 1)]
    }
    
    return(selected)  # Return the filtered indices for this subgroup
}


# Function to filter rows within each group-edge based on tolerance
#' @noRd
filter_group_edge <- function(df, tolerance) {
    # Start with all indices selected
    selected <- seq_len(nrow(df))
    positions <- df$position
    
    # Calculate initial differences
    diffs <- diff(positions)
    
    # While there are differences below tolerance
    while (any(diffs < tolerance)) {
        # Find the first position where diff is below tolerance
        first_below <- which(diffs < tolerance)[1]
        
        # Remove the second element in the violating pair
        selected <- selected[-(first_below + 1)]
        
        # Recalculate diffs only around the modified region
        if (first_below > 1) diffs[first_below - 1] <- positions[selected[first_below]] - positions[selected[first_below - 1]]
        diffs <- diffs[-first_below]
    }
    
    # Return the original indices of the selected rows
    return(df$orig_index[selected])
}

# Function to find merged indices for unselected rows only
#' @noRd
find_merged_indices_for_unselected <- function(selected_rows, total_rows) {
    # Identify unselected rows
    all_rows <- 1:total_rows
    unselected_rows <- setdiff(all_rows, selected_rows)
    
    # For each unselected row, find the nearest preceding selected row
    merged_indices <- sapply(unselected_rows, function(row) {
        max(selected_rows[selected_rows <= row])
    })
    
    return(merged_indices)
}

# Helper function to fill missing values using "merge" strategy
#' @noRd
fill_na_merge <- function(data, removed_merge, ref_idx, removed_indices) {
    for (col in names(data)) {
        # Check if the current entry has NA for the reference row
        if (is.na(data[[col]][ref_idx])) {
            # Find the first non-NA value in the removed observations for this column
            for (removed_idx in removed_indices) {
                if (!is.na(removed_merge[[col]][removed_idx])) {
                    # Fill the NA in data with the non-NA value from removed_merge
                    data[[col]][ref_idx] <- removed_merge[[col]][removed_idx]
                    break  # Stop after filling the first non-NA value
                }
            }
        }
    }
    return(data)
}

# Main function to apply the chosen merge strategy
#' @noRd
apply_merge_strategy <- function(data, removed_merge, merge_idx_map, ref_idx_merges, merge_strategy) {
    for (i in seq_along(ref_idx_merges)) {
        # Access the mapped index using the character version of ref_idx_merges[i]
        ref_idx <- as.integer(merge_idx_map[as.character(ref_idx_merges[i])])
        
        # Get the removed observations linked to this ref_idx
        removed_indices <- which(ref_idx_merges == ref_idx_merges[i])
        
        # Apply the merge or average strategy based on `merge_strategy`
        if (merge_strategy == "merge") {
            data <- fill_na_merge(data, removed_merge, ref_idx, removed_indices)
        } else if (merge_strategy == "average") {
            data <- fill_na_average(data, removed_merge, ref_idx, removed_indices)
        }
    }
    
    return(data)
}

# Helper function to fill values using "average" strategy
#' @noRd
fill_na_average <- function(data, removed_merge, ref_idx, removed_indices) {
    for (col in names(data)) {
        # Gather all values for averaging, including the reference index value
        values_to_average <- c(data[[col]][ref_idx], removed_merge[[col]][removed_indices])
        
        # Filter out NA values from values_to_average
        non_na_values <- values_to_average[!is.na(values_to_average)]
        
        if (length(non_na_values) > 0) {
            if (is.numeric(data[[col]][ref_idx])) {
                # Use the average of all non-NA values if the column is numeric
                data[[col]][ref_idx] <- mean(non_na_values)
            } else {
                # For non-numeric, just use the first available non-NA value
                data[[col]][ref_idx] <- non_na_values[1]
            }
        }
    }
    return(data)
}


# Function to build constraint matrix
#' @noRd
construct_directional_constraint_matrix <- function(E, nV, nE, alpha, V_indegree, V_outdegree, weight,
                                    DirectionalWeightFunction_out, DirectionalWeightFunction_in) {
  
  # Precompute out_edges and in_edges for each vertex
  out_edges_list <- split(seq_len(nrow(E)), E[, 1])
  in_edges_list <- split(seq_len(nrow(E)), E[, 2])

  # Calculate an upper bound on the number of elements in i_, j_, and x_
  nC <- sum((V_outdegree > 0 & V_indegree > 0) * V_outdegree * (1 + V_indegree) + 
            (V_indegree == 0) * (V_outdegree - 1)) * alpha

  # Initialize vectors to store the row indices (i_), column indices (j_), and values (x_) of the sparse matrix
  i_ <- integer(nC)
  j_ <- integer(nC)
  x_ <- numeric(nC)

  count_constraint <- 0
  count <- 0

  # Process vertices with both outdegree and indegree
  Vs <- which(V_outdegree > 0 & V_indegree > 0)
  for (v in Vs) {
    out_edges <- out_edges_list[[as.character(v)]]
    in_edges <- in_edges_list[[as.character(v)]]
    n_in <- length(in_edges)

    # Loop through each out edge and derivative level
    for (i in seq_along(out_edges)) {
      out_weight_values <- DirectionalWeightFunction_out(weight[out_edges[i]])
      in_weight_values <- DirectionalWeightFunction_in(weight[in_edges])
      
      for (der in seq_len(alpha)) {
        # Set row indices and column indices for the current out edge
        i_[count + 1] <- count_constraint + 1
        j_[count + 1] <- 2 * alpha * (out_edges[i] - 1) + der
        x_[count + 1] <- out_weight_values  # Apply out weight
        
        # Set row indices, column indices, and values for each in edge
        i_[count + seq(2, n_in + 1)] <- count_constraint + 1
        j_[count + seq(2, n_in + 1)] <- 2 * alpha * (in_edges - 1) + alpha + der
        x_[count + seq(2, n_in + 1)] <- in_weight_values  # Apply in weights
        
        count <- count + (n_in + 1)
        count_constraint <- count_constraint + 1
      }
    }
  }

  # Process vertices with indegree == 0
  Vs0 <- which(V_indegree == 0)
  for (v in Vs0) {
    out_edges <- out_edges_list[[as.character(v)]]

    if (length(out_edges) > 1) {
      for (i in 2:length(out_edges)) {
        for (der in seq_len(alpha)) {
          # Set indices and values for indegree == 0 vertices
          i_[count + 1:2] <- count_constraint + 1
          j_[count + 1] <- 2 * alpha * (out_edges[i] - 1) + der
          j_[count + 2] <- 2 * alpha * (out_edges[i - 1] - 1) + der
          x_[count + 1:2] <- c(1, -1)
          
          count <- count + 2
          count_constraint <- count_constraint + 1
        }
      }
    }
  }

  # Construct sparse matrix with the populated i_, j_, and x_
  C <- Matrix::sparseMatrix(
    i = i_[1:count],
    j = j_[1:count],
    x = x_[1:count],
    dims = c(count_constraint, 2 * alpha * nE)
  )

  return(C)
}
