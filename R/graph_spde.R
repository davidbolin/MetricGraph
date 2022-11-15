

gpgraph_spde <- function(graph_object, alpha = 1, stationary_endpoints = "all",
 start_kappa = 1, start_sigma = 1, debug = FALSE){

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

  idx_ij <- idx_ij[idx_sub]
  idx_ij <- sort(idx_ij, index.return=TRUE)
  idx_ij <- idx_ij$ix

  idx_ij <- idx_ij - 1

  model <- do.call(
        'inla.cgeneric.define',
        list(model="inla_cgeneric_gpgraph_alpha1_model",
            shlib=paste0(gpgraph_lib, '/gpgraph_cgeneric_models.so'),
            n=as.integer(n.v), debug=debug,
            prec_graph_i = as.integer(i_),
            prec_graph_j = as.integer(j_),
            index_graph = as.integer(idx_ij),
            EtV2 = EtV2,
            EtV3 = EtV3,
            El = El,
            stationary_endpoints = as.integer(index),
            start_kappa = start_kappa,
            start_sigma = start_sigma))

}