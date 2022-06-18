#'
#' simple dim corrector function
#' @param t (n x 1) where to evaluate phi (in [0,L])
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
        r[1,] <- matern.derivative(t, kappa=kappa, nu = nu, sigma = sigma , deriv = deriv)
        r[2,] <-  matern.derivative(t-L, kappa=kappa, nu = nu, sigma = sigma, deriv = deriv )

    }
    return(r)
}

#'
#' simple dim corrector function v2
#' @param t (n x 1) where to evaluate phi (in [0,L])
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
        r[1,] <- -matern.derivative(t, kappa=kappa, nu = 3/2, sigma = sigma , deriv=deriv+1)
        r[2,] <- -matern.derivative(t-L, kappa=kappa, nu = 3/2, sigma = sigma, deriv=deriv+1 )

        return(r)
    }else{
        stop('nu not yet implimented')
    }
}

#'
#' one dim corrector function
#' @param t (n x 1) where to evaluate phi (in [0,L])
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

#'
#' Transforms the boundary conditions to phi
#' phi'(t) = B phi(t)
#'
#'
phi.to.phid <- function(kappa, sigma, nu=3/2 , L=1)
{
    if(nu==3/2){
        B = matrix(0,4,4)
        B[1, 3] <- -1
        B[2, 4] <- -1
        B[3, 1] <-  kappa^2
        B[3, 3] <-  -2 * kappa
        B[4, 2] <-    kappa^2 #not sign change due t-L (multipcation of sign(x))
        B[4, 4] <-  2* kappa
    }else if(nu==5/2){
        alpha = nu + 0.5
        B = matrix(0,6,6)
        B[1, 3] <- -1
        B[2, 4] <- -1
        B[3, 5] <- -1
        B[4, 6] <- -1
        B[5, 1] <-  -choose(alpha,0) * kappa^alpha
        B[5, 3] <-  choose(alpha,1) * kappa^(alpha - 1)
        B[5, 5] <-  -choose(alpha, 2) *  kappa^(alpha - 2)
        B[6, 2] <-  kappa^alpha
        B[6, 4] <-  choose(alpha,1) * kappa^(alpha - 1)
        B[6, 6] <-  choose(alpha,2) *  kappa^(alpha - 2)
    }else{
        stop('nu not yet implimented')
    }
    return(B)

}


#' inverse corrector elements
#' builds the elements of the inverse of the corrector matrix
#' r(t,s)  =r_0(t-s) - phi(t) A phi(s)
#' where r_0 is stationary matern
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
corrector.invrese.e <- function(kappa, sigma, nu=3/2, L = 1){
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
        B12 <- matern.derivative(D,kappa=kappa,nu=3/2,sigma=sigma,1)
        B12[1,2] <-  -B12[1,2]
        B.element$B12 <- B12
        B22 <- -matern.derivative(D,kappa=kappa,nu=3/2,sigma=sigma,2)
        B.element$B22 <- B22
    }else{
        stop('nu not yet implimented')
    }
    return(B.element)
}

#' The corrector matrix, A, such that neumann is:
#' r(t,s)  =r_0(t-s) - phi(t) A phi(s)
#' where r_0 is stationary matern
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
corrector <- function(kappa, sigma, nu=3/2, L = 1){

    B <- corrector.invrese.e(kappa, sigma, nu, L)
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
#' The corrector matrix, A, such that free space neumann is:
#' r(t,s)  =r_0(t-s) - phi(t) A phi(s)
#' where r_0 is stationary matern
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
corrector.neumann.free <- function(kappa, sigma, nu=3/2, L = 1){
    if(nu==3/2){
        D <- matrix(c(0,L,L,0),2,2)
        B11 <- -matern.covariance(D,kappa=kappa,nu=3/2,sigma=sigma)
        B11[1,2] <- B11[2,1] <-  -B11[1,2]
    }else{
        stop('nu not yet implimented')
    }
    return(solve(B11))
}


#'
#' computes the covariance of free space (r'=0) neumann boundary
#' @param s (n x 1) location
#' @param t (n x 1) location
matern.neumann.free <- function(s, t, kappa, sigma, nu=3/2, L = 1, deriv = c(0,0)){

    if(nu==3/2){

    }else if(nu==1/2){

    }else{
        stop('nu not yet implimented')
    }
    D <- outer(s, t, "-")
    phi_t <- phi1(t, kappa, sigma, nu, L, deriv = deriv[2])
    phi_s <- phi1(s, kappa, sigma, nu, L, deriv = deriv[1])
    A <- corrector.neumann.free(kappa, sigma, nu, L)
    if(sum(deriv)==0){
        r0 <- matern.covariance(D,kappa=kappa,nu=nu,sigma=sigma)
    }else{
        r0 <- (-1)^(deriv[2]+2)*matern.derivative(D,kappa=kappa,nu=nu,sigma=sigma, deriv = sum(deriv))
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
matern.neumann.free2 <- function(s, t, C, kappa, sigma=1, nu=3/2, L = 1, deriv = c(0,0)){

    if(nu==3/2){

    }else{
        stop('nu not yet implimented')
    }
    D <- outer(s, t, "-")
    phi_t <- phi(t, kappa, sigma, nu, L, deriv = deriv[2])
    phi_s <- phi(s, kappa, sigma, nu, L, deriv = deriv[1])
    A <- corrector(kappa, sigma, nu, L)
    if(sum(deriv)==0){
        r0 <- matern.covariance(D,kappa=kappa,nu=nu,sigma=sigma)
    }else{
        r0 <- (-1)^(deriv[2]+2)*matern.derivative(D,kappa=kappa,nu=nu,sigma=sigma, deriv = sum(deriv))
    }
    r <- r0 - t(phi_s)%*%A%*%phi_t

    # extra term
    phi2_t <- phi_t[3:4,]
    phi2_s <- phi_s[3:4,]
    phi1_t <- phi_t[1:2,]
    phi1_s <- phi_s[1:2,]
    Ainv.e <- corrector.invrese.e(kappa, sigma, nu, L)
    B11.inv.B12 <- solve(Ainv.e$B11,t(Ainv.e$B12))
    Sigma_X1_tilde <- Ainv.e$B22 - t(Ainv.e$B12)%*%solve(Ainv.e$B11,Ainv.e$B12)
    Sigma_Xt_X1_tilde <- -t(phi2_t) + t(phi1_t)%*%B11.inv.B12
    Sigma_Xs_X1_tilde <- -t(phi2_s) + t(phi1_s)%*%B11.inv.B12
    Sigma_X1_tilde.inv <- solve(Sigma_X1_tilde)
    A2 <- Sigma_X1_tilde.inv%*%C%*%Sigma_X1_tilde.inv
    r <- r + Sigma_Xs_X1_tilde%*%A2%*%t(Sigma_Xt_X1_tilde)
    return(r)
}

#' computes the covariance of r(s,t) using the corrector matrix form
#' @param s (n x 1) location
#' @param t (n x 1) location
#' @param kappa  matern param
#' @param sigma  matern param
#' @param nu     shape param
#' @param L      interval length
#' @param deriv  derivative order s, t
matern.neumann.c <- function(s, t, kappa, sigma, nu=3/2, L = 1, deriv = c(0,0)){

    if(nu==1/2){
    }else if(nu==3/2){

    }else{
        stop('nu not yet implimented')
    }
    D <- outer(s, t, "-")
    phi_t <- phi(t, kappa, sigma, nu, L, deriv = deriv[2])
    phi_s <- phi(s, kappa, sigma, nu, L, deriv = deriv[1])
    A <- corrector(kappa, sigma, nu, L)
    if(sum(deriv)==0){
        r0 <- matern.covariance(D,kappa=kappa,nu=nu,sigma=sigma)
    }else{
        r0 <- (-1)^(deriv[2]+2)*matern.derivative(D,kappa=kappa,nu=nu,sigma=sigma, deriv = sum(deriv))
    }
    r <- r0 - t(phi_s)%*%A%*%phi_t
    return(r)
}



#' computes C matrix for beta=1
#' @param L edge length
#' @param kappa  matern parameter
#' @param sigma  matern parameter
#' @param nu  matern parameter
build.C.beta1 <- function(L, kappa, sigma, nu){
    C_0 <- matern.neumann.free(c(0,L),c(0,L),kappa,sigma=1, nu=3/2, L=L,deriv=c(1,1))
    return(sigma^2*solve(solve(C_0) -0.5*diag(2)/kappa^2))
}

#' computes matrices needed for graph for beta=1
#' @param P vertices
#' @param E edges
#' @param kappa  matern parameter
#' @param sigma  matern parameter
#' @param n  number of nodes to include per edge
#' @return a list with elements:
#'     Sigma : Covariance matrix of 'free' processes and derivatives on the edges
#'     A : Matrix that implements the Kirchoff constraints
#'     rep.ind : the indices of the repeated nodes
#'     p.ind : the indices in Sigma corresponding to the process
#'     d.ind : the indices in Sigma corresponding to the derivative
#'     P : matrix with the locations of the nodes where Sigma is evaluated
build.Sigma.beta1 <- function(P,E,kappa,sigma,n){

    #edge lengths
    L <- sqrt((P[E[,2],1]-P[E[,1],1])^2 + (P[E[,2],2]-P[E[,1],2])^2)

    #number of vertices
    n.v <- dim(P)[1]

    #number of edges
    n.e <- dim(E)[1]

    Sigma <- matrix(0,nrow=2*n*length(L),ncol=2*n*length(L))

    p.inds <- vector("list",n.v) #indices of process at vertices
    d.inds <- vector("list",n.v) #indices of derivative at vertices
    d.sign <- vector("list",n.v) #sign of derivative
    proc.index <- NULL #indices for process nodes
    deri.index <- NULL #indices for derivative nodes
    P.sigma <- NULL #coorinates of nodes

    #loop over edges and construct free covariances
    for(i in 1:n.e){
        C <- build.C.beta1(L[i], kappa=kappa, sigma=sigma, nu=3/2)
        x <- seq(from=0,to=L[i],length.out = n)
        P.sigma <- cbind(P.sigma, P[E[i,1],]+outer(P[E[i,2],]-P[E[i,1],],seq(from=0,to=1,length.out=n)))
        S1 <- matern.neumann.free2(x, x, C, kappa=kappa, sigma=sigma, nu=3/2, L = L[i])
        S2 <- matern.neumann.free2(x, x, C, kappa=kappa, sigma=sigma, nu=3/2, L = L[i], deriv = c(0,1))
        S3 <- matern.neumann.free2(x, x, C, kappa=kappa, sigma=sigma, nu=3/2, L = L[i], deriv = c(1,1))
        Sigma[(i-1)*2*n + (1:(2*n)),(i-1)*2*n + (1:(2*n))]  <- rbind(cbind(S1, S2), cbind(t(S2),S3))

        p.inds[[E[i,1]]] <- c(p.inds[[E[i,1]]],(i-1)*2*n + 1)
        p.inds[[E[i,2]]] <- c(p.inds[[E[i,2]]],(i-1)*2*n + n)

        d.inds[[E[i,1]]] <- c(d.inds[[E[i,1]]],(i-1)*2*n + n + 1)
        d.inds[[E[i,2]]] <- c(d.inds[[E[i,2]]],(i-1)*2*n + 2*n)

        d.sign[[E[i,1]]] <- c(d.sign[[E[i,1]]],1)
        d.sign[[E[i,2]]] <- c(d.sign[[E[i,2]]],-1)

        proc.index <- c(proc.index,(i-1)*2*n + (1:n))
        deri.index <- c(deri.index,(i-1)*2*n + n + (1:n))
    }

    #now build A matrix
    i <- j <- x <- NULL
    n.c = 0 #number of constraints
    for(k in 1:n.v){
        n.c <- n.c+1
        d <- length(p.inds[[k]]) #degree of vertex
        i <- c(i,rep(n.c,d))
        j <- c(j,d.inds[[k]])
        x <- c(x,d.sign[[k]])
        if(d>1){
            for(kk in 1:(d-1)){
                n.c <- n.c + 1
                i <- c(i,rep(n.c,2))
                j <- c(j,p.inds[[k]][kk:(kk+1)])
                x <- c(x,c(1,-1))
            }
        }
    }
    A <- Matrix::sparseMatrix(i=i,j=j,x=x,dims=c(n.c,2*n*n.e))

    #fix index vectors with the indices
    rep.ind <- NULL
    for(k in 1:dim(P)[1]){
        d <-length(p.inds[[k]])
        if(d>1){
            rep.ind <- c(rep.ind, p.inds[[k]][2:d],d.inds[[k]][2:d])
            proc.index <- setdiff(proc.index,p.inds[[k]][2:d])
            deri.index <- setdiff(deri.index,p.inds[[k]][2:d])
        }

    }

    return(list(Sigma=Sigma,
                A=A,
                rep.ind = rep.ind,
                p.ind = proc.index,
                d.ind = deri.index,
                P = t(P.sigma[,!duplicated(t(P.sigma))])))
}
