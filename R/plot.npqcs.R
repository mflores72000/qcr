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
# plot.npqcs
#-------------------------------------------------------------------------
##' Plot method for 'npqcs' objects
##' 
##' Generic function for plotting Multivarite charts of object of class 'npqcs' to perform statistical 
##' quality control.
##' 
##' @method plot npqcs
##' @param x  Object npqcs (Multivarite Quality Control Statical)
##' @param title an overall title for the plot
##' @param subtitle a sub title for the plot
##' @param xlab a title for the x axis
##' @param ylab a title for the y axis
##' @param ylim the y limits of the plot
##' @param lim a logical value indicating that limits should be constant.
##' @param ...  arguments to be passed to or from methods.
##' @export


plot.npqcs <- function(x, title, subtitle, xlab, ylab, ylim, lim = TRUE, ...)
  #.........................................................................                     
{
  
  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)
  limits <- x$limits
  sample <- 1:length(x$statistics)

    plot(x$statistics ~ sample, type =  "n",pch  =  16, axes  =  FALSE, 
         main  =  title, sub  =  subtitle, xlab  = xlab, 
         ylab  =  ylab, ylim  =  ylim)    
  axis(1, at  =  sample, cex.axis  =  0.7)
  axis(2, cex.axis  =  0.7)
  
  axis(4, at = c(max(x$limits), min(x$limits)), 
       labels = c("CL",paste("LCL=",min(x$limits))), adj = 0, las = 1)
    
  
  rect(par("usr")[1],
       par("usr")[3],
       par("usr")[2],
       par("usr")[4],
       col  =  "#CCCCCC")
  box(col  =  "#CCCCCC")
  grid(col  =  "#EEEEEE")
  
  #x <-res.npqcs
      points(sample, x$statistics, pch = 3, cex = 0.8)
      lines(x$statistics ~ sample, type = "o", pch=20)      

      if (lim == TRUE){ 
        
        lcl <- x$limits[1]
        cl <- x$limits[2]
        
        abline(h  =  lcl, lwd  =  2, col  =  "red", lty = 2)
        abline(h  =  cl, lwd  =  2, col  =  "blue",lty = 2) 
      }else
        {
        lcl <- x$limits[-1]
        cl <- x$limits[1]
        
        lines(sample,  lcl, lwd  =  2, col  =  "red", lty = 2)
        abline(h  =  cl, lwd  =  2, col  =  "blue",lty = 2)
        }  
      
  
    beyond.limits <- x$violations

    
    points(x$statistics[beyond.limits]~sample[beyond.limits],
           col = "red",
           pch  =  19)
    
  par(oldpar)
  #.........................................................................
} # plot.npqcs
#.........................................................................

#-------------------------------------------------------------------------
# plot.npqcs.r
#-------------------------------------------------------------------------
##' @rdname  plot.npqcs
##' @method  plot npqcs.r
##' @export
##' 

plot.npqcs.r <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                         ylab  =  NULL, ylim  =  NULL, ...)
  #.........................................................................                     
{
#  names(x)
#    str(x)
  if(is.null(ylim)) 
    ylim <-  range(x$statistics, x$limits)
  
  if (is.null(title)) title <- expression(paste("Chart of control ", r," "))
  
  if (is.null(subtitle)) subtitle <- ""
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- expression(Rank)
  
  plot.npqcs(x, title , subtitle, xlab, ylab, ylim)
} #plot.npqcs.r
#.........................................................................                     



#-------------------------------------------------------------------------
# plot.npqcs.Q
#-------------------------------------------------------------------------
##' @rdname  plot.npqcs
##' @method  plot npqcs.Q
##' @export
##' 
#x<-data.npqcs.Q
plot.npqcs.Q <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                         ylab  =  NULL, ylim  =  NULL, ...)
  #.........................................................................                     
{
  #  names(x)
  #    str(x)
  if(is.null(ylim)) 
    ylim <-  range(x$statistics, x$limits)
  
  if (is.null(title)) title <- expression(paste("Chart of control ", Q," "))
  
  if (is.null(subtitle)) subtitle <- ""
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- expression(Q_n)
  
  plot.npqcs(x, title , subtitle, xlab, ylab, ylim)
} #plot.npqcs.Q
#.........................................................................                     



#-------------------------------------------------------------------------
# plot.npqcs.S
#-------------------------------------------------------------------------
##' @rdname  plot.npqcs
##' @method  plot npqcs.S
##' @export
##' 

#x <- res.npqcs
plot.npqcs.S <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                         ylab  =  NULL, ylim  =  NULL, ...)
  #.........................................................................                     
{
  #  names(x)
  #    str(x)
  if(is.null(ylim)) 
    ylim <-  range(x$statistics, x$limits)
  
  if (is.null(title)) title <- expression(paste("Chart of control ", S," "))
  
  if (is.null(subtitle)) subtitle <- ""
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- expression(S_n)
  
  plot.npqcs(x, title , subtitle, xlab, ylab, ylim,, lim = TRUE)
} #plot.npqcs.S
#.........................................................................                     
