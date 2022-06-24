#-----------------------------------------------------------------------------#
#                                                                             #
#                  QUALITY CONTROL STATISTICS IN R                            #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sanchez                                       #
#              Professor of Mathematic Department                             #
#              Escuela Politecnica Nacional, Ecuador                          #
#              miguel.flores@epn.edu.ec                                       #
#                                                                             #
#-----------------------------------------------------------------------------#
#-------------------------------------------------------------------------
# plot.mqcs
#-------------------------------------------------------------------------
##' Plot method for 'fdqcd' objects
##' 
##' Generic function for plotting Multivarite charts of object of class 'fdqcd' 
##' to perform statistical quality control.
##' 
##' @method plot fdqcd
##' @param x  Object fdqcd (pashe I)
##' @param y  Object fdqcd (monitoring)
##' @param title An overall title for the plot.
##' @param xlab A title for the x axis.
##' @param ylab A title for the y axis.
##' @param col  The color for curves.
##' @param ...  Arguments to be passed to or from methods.
##' @export
##' @examples 
##' library(qcr)
##' m <- 30
##' tt<-seq(0,1,len=m)
##' mu<-30 * tt * (1 - tt)^(3/2)
##' n0 <- 100
##' set.seed(12345)
##' mdata<-matrix(NA,ncol=m,nrow=n0)
##' sigma <- exp(-3*as.matrix(dist(tt))/0.9)
##' for (i in 1:n0) mdata[i,]<- mu+0.5*mvrnorm(mu = mu,Sigma = sigma )
##' fdchart <- fdqcd(mdata)
##' plot(fdchart,type="l",col="gray")

plot.fdqcd <- function(x, y = NULL,title=NULL, xlab=NULL, ylab=NULL, col=NULL,...)
  #.........................................................................                     
{
     if (inherits(x, "fdqcd")) {
      if(is.null(title)) title <- attributes(x)$data.name
      if(is.null(xlab)) xlab <- "t"
      if(is.null(ylab)) ylab <- "X(t)"
      if(is.null(col)) col <- "gray"
      
      data <- x[["data"]]
      tt <- x[["argvals"]]
      rtt <- x[["rangeval"]]
      names <- list(main=title,xlab=xlab,ylab=ylab) 
      x <- fdata(x$data,tt,rtt,names)
      plot.fdata(x,col=col,...)
      lines(func.mean(x),col=1,lwd=2) #original curve
      
      legend(x = min(tt), y = 0.99 * max(data), bty = "n",
             legend = c("Curves of Calibrating", 
                        "Median (Deepest)"), 
             lty = 1, 
             lwd = c(1,2),
             col = c(col,"black"), cex = 0.9, 
             box.col = 0) 
      if(!is.null(y)){
        if (any(class(y) == "fdqcd")) {
          data <- y[["data"]]
          tt <- y[["argvals"]]
          rtt <- y[["rangeval"]]
          names <- list(main=title,xlab=xlab,ylab=ylab) 
          y <- fdata(y$data,tt,rtt,names)
          lines(y,col="black",...)
          lines(func.mean(y),col=2,lwd=2) #original curve
          
          legend(x = max(tt)*0.8, y = 0.99 * max(data), bty = "n",
                 legend = c("Curves of Monitoring", 
                            "Median (Deepest)"), 
                 lty = 1, 
                 lwd = c(1,2),
                 col = c("black",2), cex = 0.9, 
                 box.col = 0) 
      }
        }
    }

  #.........................................................................
} # plot.fdqcd
#.........................................................................