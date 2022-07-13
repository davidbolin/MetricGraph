
#' Plots of graphs
#'
#' @param P Vertex matrix of graph
#' @param E Edge matrix of graph
#' @param f function to plot
#' @param nc number of colors to use
#' @param rgl Use RGL plot?
#' @param N.e Number of points at which the function is define on each edge
#' @param add Add plot to existing plot p (for RGL)?
#' @param p Plot to add results to
#' @param col color to use
#' @param size size of elements to plot
#' @param plot.vertices plot the vertices of the graph?
#'
#' @return The plot object
#' @export
plot_graph <- function(P,E,f=NULL,nc=100,rgl=FALSE,N.e=NULL,
                       add=FALSE,p=NULL,col=NULL,size=6,plot.vertices=TRUE){
  if(rgl==FALSE){
    plot(x=NULL,xlim=c(min(P[,1]),max(P[,1])),ylim=c(min(P[,2]),max(P[,2])))
    if(plot.vertices){
      points(P)
    }
    if(is.null(f)){
      for(i in 1:dim(E)[1]){
        xl <- c(P[E[i,1],1],P[E[i,2],1])
        yl <- c(P[E[i,1],2],P[E[i,2],2])
        lines(xl,yl)
      }
    } else{
      cols <-  turbo(nc)#viridis(nc)
      f.min <- min(f)
      f.max <- max(f)
      for(i in 1:dim(E)[1]){
        if(P[E[i,1],1]!= P[E[i,2],1] && P[E[i,1],2]!=P[E[i,2],2]){
          xl <- c(P[E[i,1],1],P[E[i,2],1])
          yl <- c(P[E[i,1],2],P[E[i,2],2])
          vl <- c(f[E[i,1]],f[E[i,2]])
          if(0){
            lines(xl,yl)
          } else {
            xy <- approx(xl,yl,n=10)
            vv <- approx(x=vl,y=NULL,n=10)$y
            lines(xy$x,xy$y,col=cols[100*round((vv - f.min)/(f.max-f.min),2)])
          }
        }
      }
    }
    if(!is.null(f) && plot.vertices){
      cols <-  turbo(nc)#viridis(nc)
      c <- cols[cut(f, breaks = nc)]
      points(P,col=c,pch=16)

    }
  } else {
    n.v <- dim(P)[1]
    n.e <- dim(E)[1]
    n <- (length(f)-n.v)/n.e
    data <- data.frame(x=NULL,y=NULL,z=NULL,i=NULL)
    for(i in 1:n.e){
      v1 <- E[i,1]
      v2 <- E[i,2]
      if(i==1){
        ind <- c(v1,n.v + (i-1)*n + 1:N.e[i],v2)
      } else {
        ind <- c(v1,n.v + sum(N.e[1:(i-1)]) + 1:N.e[i],v2)
      }
      x <- mesh$P[ind,1]
      y <- mesh$P[ind,2]
      z <- f[ind]
      li <- rep(i,length(x))
      df.i <- data.frame(x=x,y=y,z=z,i=li)
      data <- rbind(data,df.i)
    }
    colnames(data) <- data.frame("x","y","z","i")
    if(add==FALSE){
      p <- plot_ly(data=data, x = ~y, y=~x,z=~z)
      p <- p %>% add_trace(mode="lines",type="scatter3d",
                           line = list(width = size,color=col),split=~i,
                           showlegend=FALSE)
      start = dim(P)[1]+1
      for(i in 1:n.e){
        data <- data.frame(x = mesh$P[start:(start+mesh$n.e[i]-1),1],
                           y = mesh$P[start:(start+mesh$n.e[i]-1),2],
                           z = rep(0,mesh$n.e[i]))
        start <- start + mesh$n.e[i]
        p <- p %>% add_trace(data=data,x=~y,y=~x,z=~z,mode="lines",type="scatter3d",
                             line = list(width = size,color='rgb(0,0,0)'))
        data <- data.frame(x=P[,1],y=P[,2],z=rep(0,dim(P)[1]))
        p <- p %>% add_trace(data=data,x=~y,y=~x,z=~z, type="scatter3d",mode= "markers",
                             marker = list(width = size,color='rgb(0,0,0)'))

      }
      if(1){
        for(i in 1:dim(mesh$P)[1]){
          data <- data.frame(x = mesh$P[i,1]*c(1,1),
                             y = mesh$P[i,2]*c(1,1),
                             z = c(0,f[i]))
          p <- p %>% add_trace(data=data,x=~y,y=~x,z=~z,mode="lines",type="scatter3d",
                               line = list(width = size/2,color='rgba(0,0,0,0.1)'))
        }
      }
      p <- p %>% layout(showlegend=FALSE)
      p
      return(p)
    } else {
      p <- p %>% add_trace(data=data, x = ~y, y=~x,z=~z,mode="lines",type="scatter3d",
                           line = list(width = size,color=col),split=~i,showlegend=FALSE)
      p
      return(p)
    }
  }

}

plot.mesh <- function(mesh,...){
  plot.graph(mesh$P,mesh$E,...)
}
