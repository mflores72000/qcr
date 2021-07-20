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
# DFD chart
#-------------------------------------------------------------------------
##' Function to plot depth functional data (DFD) - chart
##'
##' This function is used to compute statistics required by the DFD chart.
##' @param x   an R object (used to select the method). See details.
##' @param ... arguments passed to or from methods.
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
##' plot.fdqcd(fdchart,type="l",col="gray")
##' set.seed(1234)
##' fddep <- fdqcs.depth(fdchart,plot = T)
##' plot(fddep,title.fdata = "Fdata",title.depth = "Depth")
##' summary(fddep)
##' }

fdqcs.depth <- function(x, ...) {
  UseMethod("fdqcs.depth")
}

##' @rdname fdqcs.depth
##' @method fdqcs.depth default
##' @inheritParams fdqcd
##' @param func.depth Type of depth measure, by default depth.mode.
##' @param nb         The number of bootstrap samples.
##' @param type       the method used to trim the data (trim or pond).
##' @param ns      Quantile to determine the cutoff from the Bootstrap procedure
##' @param plot       a logical value indicating that it should be plotted.
##' @param trim       The porcentage of the trimming.
##' @param smo        The smoothing parameter for the bootstrap samples.
##' @param draw.control ist that it specifies the col, lty and lwd for objects: fdataobj, statistic, IN and OUT.
##' @export

fdqcs.depth.default <- function(x, data.name=NULL,func.depth = depth.mode,nb=200,
                                type = c("trim","pond"),ns =  0.01, 
                                plot = TRUE, trim = 0.025, smo =0.05,
                                draw.control = NULL,...)
#.........................................................................
  {
  
  obj<-fdqcd(x= x, data.name = data.name)

  result<-fdqcs.depth.fdqcd(x = obj, func.depth = func.depth,nb=nb,
                           type = type,ns =  ns, 
                           plot = plot, trim = trim, smo =smo,
                           draw.control = draw.control, ...)

  return(result)
} # fdqcs.depth.default
#.........................................................................

##' @rdname  fdqcs.depth
##' @method fdqcs.depth fdqcd
##' @inheritParams fdqcs.depth.default
##' @export


fdqcs.depth.fdqcd <- function(x, func.depth = depth.mode,nb=200,
                              type = c("trim","pond"),ns =  0.01, 
                              plot = TRUE, trim = 0.025, smo =0.05,
                              draw.control = NULL,...){
  type <- match.arg(type)
  data <- x[["data"]]
  tt = x[["argvals"]]
  rtt = x[["rangeval"]]
  data.name = attributes(x)$data.name
  names <- list(main=data.name)
  x <- fdata(x$data,tt,rtt,names)
 
  
  if(type=="trim"){
    out <- outliers.depth.trim(x,nb = nb, ns =ns,
                               smo = smo, dfunc = func.depth,trim = trim)
    Dep <- out$Dep
    ind <- as.numeric(out$outliers[out$iteration==1])
  }else{
    out <- outliers.depth.pond(x,nb = nb, ns = ns,
                               smo = smo, dfunc = func.depth)
    Dep <- out$Dep
    ind <- as.numeric(out$outliers[out$iteration==1])
  }
  
  if(length(ind)>0) fdaenv <- x[-ind,] else fdaenv <- x
  
  LCL <- apply(fdaenv[["data"]], 2,min)
  fmin <- fdata(LCL, tt, rtt)
  
  UCL <- apply(fdaenv[["data"]], 2,max)
  fmax <- fdata(UCL, tt, rtt)
  
  med <- func.med.mode(x,trim=trim)
  
  if (plot) {
    par(mfrow=c(1,2))
    ##Funciontal data
    if (is.null(draw.control)) 
      draw.control = list(col = c("grey", "blue", "red"), 
                          lty = c(1, 1, 2), lwd = c(1, 2, 2))
    if (is.null(draw.control$lwd)) 
      draw.control$lwd = c(1, 2, 1)
    if (is.null(draw.control$lty)) 
      draw.control$lty = c(2, 1, 1)
    if (is.null(draw.control$col)) 
      draw.control$col = c("grey", "blue", "pink")
    
    plot(x, lwd = draw.control$lwd[1], lty = draw.control$lty[1], 
         col = draw.control$col[1])
    

    lines(fmin, lwd = draw.control$lwd[3], 
          lty = draw.control$lty[3], col = draw.control$col[3])

    lines(fmax, lwd = draw.control$lwd[3], 
          lty = draw.control$lty[3], col = draw.control$col[3])
    #    lines(depth$mtrim, lwd = draw.control$lwd[3], 
    #          lty = draw.control$lty[3], col = draw.control$col[3])
    
    lines(med, lwd = draw.control$lwd[2], 
          lty = draw.control$lty[2], 
          col = draw.control$col[2])
    
    
    if(length(ind)>0)
      lines(x[ind,], lwd = 2, 
            lty = 3, col = 1)
    
    legend(x = min(tt), y = 0.99 * max(data), bty = "n",
           legend = c("Curves of Calibrating", 
                      "Median (Deepest)", paste("Envelope",(1-ns)*100,"%"),
                      "Outliers"), 
           lty = c(1,1,2,2), 
           lwd = c(draw.control$lwd,draw.control$lwd[2],2),
           col = c(draw.control$col,"black"), cex = 0.9, 
           box.col = 0) 
    
    plot(Dep, type="b",pch=16, main = "Phase I: Depth Chart",
         ylim=c(min(Dep,out$quantile),max(Dep)),xlab="t",ylab="Depth")
    abline(h = out$quantile, lty = 2, col = "red")
    par(mfrow=c(1,1))
  }
  
  result <- list(fdata = x,Depth = Dep,LCL=out$quantile, out=ind,fmin=fmin,
                 fmax=fmax,fmed=med, ns=ns)
  
  oldClass(result) <- "fdqcs.depth"
  
  return(result)
} # fdqcs.depth.fdqcd
#.........................................................................

