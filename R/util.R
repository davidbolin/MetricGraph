#internal function for checking metric_graph inputs
gpgraph_check_graph <- function(graph)
{
  if (!inherits(graph, "GPGraph::graph")) {
    stop("The graph object is not a metric graph")
  }
  out <- list(has.mesh = FALSE,
              has.obs = FALSE)
  if(!is.null(graph$mesh)){
    out$has.mesh = TRUE
  }
  if(!is.null(graph$y) && !is.null(graph$PtE))
    out$has.data = TRUE
  return(out)
}


#' Create metric graphs for connected components of a SpatialLines object
#'
#' @param lines Object of class `SpatialLines`
#' @param by_length Sort the components by total edge length? If FALSE,
#' the components are sorted by the number of vertices.
#' @param only_largest if TRUE, only return the largest component.
#' Otherwise return an ordered list with the components (largest first)
#'
#' @return A `metric_graph` object created from the largest component, or a
#' list of `metric_graph` objects for all connected components
#' @export
#'
#' @examples
#' library(sp)
#' line1 <- Line(rbind(c(0, 0), c(1, 0)))
#' line2 <- Line(rbind(c(1, 0), c(2, 0)))
#' line3 <- Line(rbind(c(1, 1), c(2, 1)))

#' Lines <-  SpatialLines(list(Lines(list(line1), ID = "1"),
#'                            Lines(list(line2), ID = "2"),
#'                            Lines(list(line3), ID = "3")))
#' graphs <- graph_components(Lines, only_largest = FALSE)
#' p <- graphs[[1]]$plot(edge_color = "red")
#' graphs[[2]]$plot(p = p, edge_color = "blue")
graph_components <- function(lines, by_length = TRUE, only_largest = TRUE) {

  graph <- metric_graph$new(lines = lines)
  g <- graph(edges = c(t(graph$E)), directed = FALSE)
  igraph::E(g)$weight <- graph$edge_lengths
  components <- igraph::clusters(g, mode="weak")
  Graphs <- list()
  nc <- components$no
  for(k in 1:nc) {
    vert_ids <- igraph::V(g)[components$membership == k]
    edge_rem <- NULL
    for (i in 1:graph$nE) {
      if(!(graph$E[i, 1] %in% vert_ids) && !(graph$E[i, 2] %in% vert_ids))
        edge_rem <- c(edge_rem, i)
    }
    edge_keep <- setdiff(1:graph$nE, edge_rem)
    Graphs[[k]] = metric_graph$new(lines = Lines[edge_keep])
  }
  sizes <- components$csize
  lengths <- unlist(lapply(1:nc, function(x) sum(Graphs[[x]]$edge_lengths)))
  if(by_length) {
    reo <- sort(lengths, decreasing = TRUE, index.return = TRUE)$ix
  } else {
    reo <- sort(sizes, decreasing = TRUE, index.return = TRUE)$ix
  }
  g <- Graphs[reo]
  if(only_largest) {
    return(g[[1]])
  } else {
    return(g)
  }
}


#'
#' computes the covariance of free space (r'=0) neumann boundary
#' @param s (n x 1) location
#' @param t (n x 1) location
#' @export
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
#' @export
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
#' @export
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



