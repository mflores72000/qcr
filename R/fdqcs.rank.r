#-----------------------------------------------------------------------------#
#                                                                             #
#                  QUALITY CONTROL STATISTICS IN R                            #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sanchez                                       #
#              Professor of the Mathematics Department                        #
#              Escuela Politecnica Nacional, Ecuador                          #
#              miguel.flores@epn.edu.ec                                       #
#                                                                             #
#-----------------------------------------------------------------------------#
#-------------------------------------------------------------------------
# RFD chart
#-------------------------------------------------------------------------
##' Function to plot rank functional data (RFD) - chart
##'
##' This function is used to compute statistics required by the RFD chart.
##' @param x   An R object (used to select the method). See details.
##' @param ... Arguments passed to or from methods.
##' @export
##' @references Flores, M.; Naya, S.; Fernández-Casal,R.; Zaragoza, S.; Raña, P.; Tarrío-Saavedra, J. 
##' Constructing a Control Chart Using Functional Data. Mathematics 2020, 8, 58.
##' @examples
##' \dontrun{
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
##' summary(fdchart)
##' plot(fdchart,type="l",col="gray")
##' out <- fddep$out
##' ## Outliers - State in Control
##' alpha <- 0.005
##' trim <- 0.1
##' while (length(out)>0) {
##'  mdata <- fddep$fdata$data[-out,]
##'  fddep <- fdqcs.depth(mdata,ns = alpha, trim=trim, plot=FALSE)
##'  out <- fddep$out
##'}
##' plot(fddep,title.fdata = "FD-State in Control",title.depth = "Depth")
##' # Ha
##' mu_a<- 30 * tt^(3/2) * (1 - tt)
##' n_a <- 50
##' set.seed(12345)
##' mdata_a<-matrix(NA,ncol=m,nrow=n_a)
##' for (i in 1:n_a) mdata_a[i,]<- mu_a+0.5*mvrnorm(mu = mu_a,Sigma = sigma )
##' fdchart_a <- fdqcd(mdata_a,"Curves Monitoring")
##' plot(fdchart_a)
##' plot(fdchart,fdchart_a,main="Phase II")
##' pashe2.chart <- fdqcs.rank(fdchart,fdchart_a)
##' plot(pashe2.chart,title.fdata = "FDA",title.rank = "Rank")
##' summary(pashe2.chart)
##' }

fdqcs.rank <- function(x, ...) {
  UseMethod("fdqcs.rank")
}

##' @rdname fdqcs.rank
##' @method fdqcs.rank fdqcd
## @inheritParams fdqcd
##' @param y   The set of new curves to evaluate the depth. fdqcd class object.
##' The set of reference curves respect to which the depth is computed. fdqcd class object.
##' @param func.depth Type of depth measure, by default depth.mode.
##' @param alpha      Quantile to determine the cutoff from the Bootstrap procedure.
##' @param plot       Logical value. If TRUE a RFD chart should be plotted.
##' @param trim       The percentage of the trimming.
##' @param draw.control It specifies the col, lty and lwd for objects: 
##' fdataobj, statistic, IN and OUT.
##' @export

fdqcs.rank.fdqcd <- function(x, y = x, 
                       func.depth = depth.FM,
                       alpha = 0.01, plot = TRUE, trim = 0.1, 
                       draw.control = NULL,...){
  
  tt = x[["argvals"]]
  rtt = x[["rangeval"]]
  namesx <- list(main=attributes(x)$data.name)
  namesy <- list(main=attributes(y)$data.name)
  fdataobjori <- fdata(x$data,tt,rtt,names = namesx)
  fdataobj <- fdata(y$data,tt,rtt,names = namesy)
  
  depthori <- func.depth(fdataobjori , fdataobjori, trim = trim,...)
  depori <- depthori$dep
  rankori <- apply(matrix(depori,nrow=length(depori)), 1, 
                   function(x) mean( depori <= x))
  
  depth <- func.depth(fdataobj, fdataobjori, trim = trim,...)
  dep <- depth$dep
  rank <- apply(matrix(dep,nrow=length(dep)), 1,
                function(x) mean( depori <= x))
  
  r <- c(rankori,rank)
  ind <- which(rank <= alpha)
  indori <- which(rankori <= alpha)
  
  if (plot) {
    ##Functional data
    par(mfrow=c(2,2))
    if (is.null(draw.control)) 
      draw.control = list(col = c("grey", "blue", "red", "red","green","black"), 
                          lty = c(1, 1, 1,3,1,3), lwd = c(1, 2, 2,2,2,2))
    if (is.null(draw.control$lwd)) 
      draw.control$lwd = c(1, 2, 2)
    if (is.null(draw.control$lty)) 
      draw.control$lty = c(1, 1, 1)
    if (is.null(draw.control$col)) 
      draw.control$col = c("grey", "blue", "red")
    # Calibrado
    plot(fdataobjori, lwd = draw.control$lwd[1], lty = draw.control$lty[1], 
         col = draw.control$col[1], main="Phase II: Fda Chart")
    
    fmin <- fdata(apply(fdataobjori[depthori$ltrim][["data"]],2,min), tt, rtt, 
                  names = list(main= "LCI (trim) - FDA"))
    lines(fmin, lwd = draw.control$lwd[4], 
          lty = draw.control$lty[4], col = draw.control$col[4])
    
    fmax <- fdata(apply(fdataobjori[depthori$ltrim][["data"]],2,max), tt, rtt, 
                  names = list(main= "LCS (trim) - FDA"))
    lines(fmax, lwd = draw.control$lwd[4], 
          lty = draw.control$lty[4], col = draw.control$col[4])
    lines(depthori$mtrim, lwd = draw.control$lwd[2], 
          lty = draw.control$lty[2], col = draw.control$col[2])
    lines(fdataobjori[depthori$lmed], lwd = draw.control$lwd[3], 
          lty = draw.control$lty[3], 
          col = draw.control$col[3])
    
    legend("topleft", bty = "n",
           legend = c("Curves of Calibrating", 
                      "Trim Mean","Median (Deepest)",paste("Envelope",(1-alpha)*100,"%")), 
           lty = draw.control$lty[1:4], 
           lwd = draw.control$lwd[1:4],
           col = draw.control$col[1:4], cex = 0.9, 
           box.col = 0)
    
    # Calibrado
    ylim <- range(c(x[["data"]]),
                  c(fdataobj[["data"]]))
    
    plot(fdataobjori, lwd = draw.control$lwd[1], lty = draw.control$lty[1], 
         col = draw.control$col[1],ylim = ylim, main="Phase II: Fda Chart")
    
    lines(fdataobj, lwd = draw.control$lwd[5], 
          lty = draw.control$lty[5], 
          col = draw.control$col[5])
    
    if(length(ind)>0){
      lines(fdataobj[ind,], lwd = draw.control$lwd[6], 
            lty = draw.control$lty[6], 
            col = draw.control$col[6])
    }
    
    
    lines(fmin, lwd = draw.control$lwd[4], 
          lty = draw.control$lty[3], col = draw.control$col[4])
    
    lines(fmax, lwd = draw.control$lwd[4], 
          lty = draw.control$lty[3], col = draw.control$col[4])
    
    legend("topleft", bty = "n",
           legend = c(paste("Envelope",(1-alpha)*100,"%"),
                      "Curves of Monitoring","Outliers"), 
           lty = c(draw.control$lty[3],draw.control$lty[5:6]), 
           lwd = draw.control$lwd[4:6],
           col = draw.control$col[4:6], cex = 0.9, 
           box.col = 0)
    
    
    ##Funciontal data
    if(length(indori)==0){
      plot(c(rankori,rank), type="b",pch=16, main = paste("Phase II: Rank Chart"),
           ylim = range(c(alpha,rankori)),ylab = "")
      abline(h = alpha, lty = 2, col = "red")
    }else{
      plot(c(rankori[-indori],rank), type="b",pch=16, main = paste("Phase II: Rank Chart"),
           ylim = range(c(alpha,rankori)), ylab = "")
      abline(h = alpha, lty = 2, col = "red")
    }
    
    plot(rank, type="b",pch=16, main = paste("Phase II: Rank Chart"), 
         ylim = range(c(alpha,rank)), ylab = "")
    abline(h = alpha, lty = 2, col = "red")
    par(mfrow=c(1,1))
  }
  
  result <- list(fdataobj =fdataobj,fdataobjori =fdataobjori,fmin=fmin,fmax=fmax,
                 rankori = rankori, depthori = depthori,
                 rank = rank, depth = depth, outliers = ind, outliersori = indori, alpha=alpha)
  oldClass(result) <- "fdqcs.rank"
  return(result)
} # fdqcs.rank
#.........................................................................

