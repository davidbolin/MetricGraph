
#' Function for plotting space-time covariances
#' 
#' Plots of marginal spatial and temporal covariances of spacetime model
#' 
#' @param graph metric_graph object
#' @param Q space-time precision matrix
#' @param Qs Spatial precision matrix as reference
#' @param t.ind time indices to plot covariances for
#' @param s.ind space indices to plot covariances for
#' @param t.shift time shifts to plot, i.e., covariances C(u(s(s.ind),t),u(*,t+t.shift))
#' are shown.
#' @param show.temporal Plot the marginal temporal covariances?
#' @param t vector with the time points for which Q is computed
#' @export
plot_spacetime_covariances <- function(graph,
                                       Q,
                                       Qs = NULL,
                                       t.ind,
                                       s.ind,
                                       t.shift=0,
                                       show.temporal = FALSE,
                                       t) {
  if(is.list(Q)) {
    L <- Q$L
    N <- dim(L)[1]
  } else {
    N <- dim(Q)[1]
  }

  n <- N/length(t)

  T <- N/n
  if(length(t.ind)>4)
    stop("max 4 curves allowed")
  if(s.ind > n)
    stop("too large space index")
  if(max(t.ind)>T-1)
    stop("too large time index")

  cols <- c("green", "cyan","blue","red")[(4-length(t.ind)+1):4]

  if(!is.null(Qs)){
    v <- rep(0,dim(Qs)[1]); v[s.ind] <- 1
    c.spatial <- solve(Qs,v)
    p <- graph$plot_function(as.vector(c.spatial), plotly = TRUE, support_width = 1, line_color = "black")
  }

  time.index <- n*(0:(T-1)) + s.ind
  ct <- matrix(0,nrow = length(t.ind),ncol = T)
  for(i in 1:length(t.ind)) {
    v <- rep(0,N)
    v[(t.ind[i]-1)*n+s.ind] <- 1
    if(is.list(Q)){
      if(Q$C.inv) {
        tmp <- solve(Q$L,solve(Q$C,solve(t(Q$L),v)))
      } else {
        tmp <- solve(Q$L,Q$C%*%solve(t(Q$L),v))
      }
    } else {
      tmp <- solve(Q,v)
    }
    ct[i,] <- tmp[time.index]
    for(j in 1:length(t.shift)) {
      ind <- ((t.ind[i]-t.shift[j]-1)*n+1):((t.ind[i]-t.shift[j])*n)
      c <- tmp[ind]
      if(length(t.shift)>1) {
        col <- cols[j]
      } else {
        col <- cols[i]
      }
      if(i == 1 && is.null(Qs)) {
        p <- graph$plot_function(as.vector(c), plotly = TRUE, support_width = 0, line_color = col)
      } else {
        p <- graph$plot_function(as.vector(c), plotly = TRUE, p = p, support_width = 0, line_color = col)
      }

    }
  }
  if(show.temporal){
    df <- data.frame(t=rep(t,length(t.ind)),y=c(t(ct)), i=rep(1:length(t.ind), each=length(t)))
    pt <- plotly::plot_ly(df, x = ~t, y = ~y, split = ~i, type = 'scatter', mode = 'lines')
    fig <- plotly::layout(plotly::subplot(p,pt), title = "Covariances",
                                    scene = list(domain=list(x=c(0,0.5),y=c(0,1))),
                                    scene2 = list(domain=list(x=c(0.5,1),y=c(0,1))))
    fig$x$layout <- fig$x$layout[grep('NA', names(fig$x$layout), invert = TRUE)]
  } else {
    fig <- p
  }

  print(fig)
  return(ct)
}
