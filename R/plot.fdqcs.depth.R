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
# plot.fdqcs.depth
#-------------------------------------------------------------------------
##' Plot method for 'fdqcs.depth' objects
##' 
##' Generic function for plotting charts of object of class 'fdqcs.depth' to perform statistical 
##' quality control.
##' @rdname  plot.fdqcs
##' @method plot fdqcs.depth
##' @param x  Object fdqcs.depth
##' @param title.fdata An overall title for the fdata plot.
##' @param title.depth An overall title for the depth plot. 
##' @param xlab A title for the x axis.
##' @param ylab A title for the y axis.
##' @param col The color for curves.
##' @param draw.control It specifies the col, lty and lwd for objects: 
##' fdataobj, statistic, IN and OUT.
##' @param ...  Arguments to be passed to or from methods.
##' @export
##' 
plot.fdqcs.depth <- function(x, title.fdata=NULL, title.depth=NULL,xlab=NULL, ylab=NULL, col=NULL,
                             draw.control=NULL, ...)
  #.........................................................................                     
{
  if (inherits(x, "fdqcs.depth")) {
    
    if(is.null(title.fdata)) title.fdata <- "Phase I: Fdata Chart"
    if(is.null(title.depth)) title.depth <- "Phase I: Depth Chart"
    if(is.null(xlab)) xlab <- "t"
    if(is.null(ylab)) ylab <- "X(t)"
    
    Dep <- x$Depth
    LCL <- x$LCL
    ind <- x$out
    fmin <- x$fmin
    fmax <- x$fmax
    fmed <- x$fmed
    tt = fmin$argvals
    ns <- x$ns
    x <- x$fdata
    
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
    
    plot(x,main = title.fdata, lwd = draw.control$lwd[1], lty = draw.control$lty[1], 
         col = draw.control$col[1], xlab = xlab,ylab=ylab)
    
    lines(fmin, lwd = draw.control$lwd[3], 
          lty = draw.control$lty[3], col = draw.control$col[3])
    lines(fmax, lwd = draw.control$lwd[3], 
          lty = draw.control$lty[3], col = draw.control$col[3])
    #    lines(depth$mtrim, lwd = draw.control$lwd[3], 
    #          lty = draw.control$lty[3], col = draw.control$col[3])
    lines(fmed, lwd = draw.control$lwd[2], 
          lty = draw.control$lty[2], 
          col = draw.control$col[2])
    
    
    if(length(ind)>0)
      lines(x[ind,], lwd = 2, 
            lty = 3, col = 1)
    
    legend(x = min(tt), y = 0.99 * max(x$data), bty = "n",
           legend = c("Curves of Calibrating", 
                      "Median (Deepest)", paste("Envelope",(1-ns)*100,"%"),
                      "Outliers"), 
           lty = c(1,1,2,2), 
           lwd = c(draw.control$lwd,draw.control$lwd[2],2),
           col = c(draw.control$col,"black"), cex = 0.9, 
           box.col = 0) 
    
    plot(Dep, type="b",pch=16, main = title.depth,
         ylim=c(min(Dep,LCL),max(Dep)),xlab=xlab,ylab="Depth")
    abline(h = LCL, lty = 2, col = "red")
    par(mfrow=c(1,1))

  #.........................................................................
}} # plot.fdqcs.depth
#.........................................................................

