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
#-----------------------------------------------------------------------------#
# Main function to create a 'qcs.hat.cpm' object
#-----------------------------------------------------------------------------#
##' Process capability index (estimate cpm)
##' 
##' Estimate  \code{"cpm"} using the method described by Kerstin Vannman(2001).
##' @aliases qcs.hat.cpm 
##' @param object qcs object of type \code{"qcs.xbar"} or \code{"qcs.one"}.
##' @param limits A vector specifying the lower and upper specification limits.
##' @param target A value specifying the target of the process. 
##' If it is \code{NULL}, the target is set at the middle value between specification limits.
##' @param mu A value specifying the mean of data.
##' @param std.dev A value specifying the within-group standard deviation.
##' @param nsigmas A numeric value specifying the number of sigmas to use.
##' @param k0 A numeric value. If the capacity index exceeds the \code{k} value, 
##' then the process is capable.
##' @param alpha The significance level (by default alpha=0.05).
##' @param n Size of the sample.
##' @param contour Logical value indicating whether contour graph should be plotted.
##' @param ylim The y limits of the plot.
##' @param ... Arguments to be passed to or from methods.
##' @export
##' @references 
##' Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd
##' ed, New York, John Wiley & Sons. \cr
##' Vannman, K. (2001). \emph{A Graphical Method to Control Process Capability}. Frontiers in Statistical Quality Control, 
##' No 6, Editors: H-J Lenz and P-TH Wilrich. Physica-Verlag, Heidelberg, 290-311.\cr
##' Hubele and Vannman (2004). \emph{The E???ect of Pooled and Un-pooled Variance Estimators on Cpm When Using Subsamples}.
##' Journal Quality Technology, 36, 207-222.\cr
##' @examples
##' library(qcr)
##' data(pistonrings) 
##' xbar <- qcs.xbar(pistonrings[1:125,],plot = TRUE)
##' mu <-xbar$center
##' std.dev <-xbar$std.dev
##' LSL=73.99; USL=74.01
##' qcs.hat.cpm(limits = c(LSL,USL),
##'            mu = mu,std.dev = std.dev,ylim=c(0,1))
##'qcs.hat.cpm(object = xbar, limits = c(LSL,USL),ylim=c(0,1))
qcs.hat.cpm <- function(object, limits = c(lsl = -3, usl = 3), 
                        target = NULL, mu = 0, std.dev = 1, nsigmas = 3, 
                        k0 = 1, alpha = 0.05, n = 50, contour =TRUE, ylim = NULL,...){
  
  if (!missing(object)){
    
    if (!inherits(object, "qcs"))
      stop("an object of class 'qcs' is required")
    if (!(object$type == "xbar" | object$type == "one")) 
      stop("Process Capability Analysis only available for charts type 
           \"qcs.xbar\" and \"qcs.one\" charts")
    x <- object[[3]][,1]
    mu <- object$center
    std.dev <- object$std.dev
  }
  if (nsigmas <= 0) 
    stop("nsigmas must be a value positive")
  
  if (alpha < 0 | alpha > 1) 
    stop("alpha must be a value between 0 and 1")
  
  if (length(limits)!=2)
    stop("specification limits must be two")
  lsl <- limits[1]
  usl <- limits[2]
  if (lsl>= usl)
    stop("lsl >= usl")
  if (!(is.numeric(usl) & is.finite(lsl)))
    lsl <- NA
  if (!(is.numeric(usl) & is.finite(lsl)))
    usl <- NA
  if (is.na(lsl) & is.na(usl))
    stop("invalid specification limits")
  
  if (is.null(target)) target <- mean(limits, na.rm = TRUE)
  
  if (is.na(lsl)) {
    if (target > usl)
      warning("target value larger than one-sided specification limit...")
  }
  if (is.na(usl)) {
    if (target < lsl)
      warning("target value smaller than one-sided specification limit...")
  }
  if (!is.na(lsl) & !is.na(usl)) {
    if (target < lsl || target > usl)
      warning("target value is not within specification limits...")
  }

  m <- (lsl+usl)/2
  d <- (usl-lsl)/2
  u <- 0
  v <- 1
  c_alpha <- k0 * sqrt(n/qchisq(alpha, n))
  if (any(1+v*((mu-target)/std.dev)^2<0)) 
    stop("sample size must be a value more than 50")
  
  ind <-  (d-u*abs(mu-m))/(nsigmas*std.dev*sqrt(1+v*((mu-target)/std.dev)^2))
  names(ind) <- c("Cpm")
  #theorical
  f.delta.t <- 1/(u + 3 * k0 * sqrt(v))
  delta.t <- seq(-f.delta.t, f.delta.t, length = 100)
  gamma.t <- sqrt((1 - u * abs(delta.t))^2/(9 * k0^2) - v * delta.t^2)  
  # empirical    
  f.delta.e <- 1/(u + 3 * c_alpha * sqrt(v))
  delta.e <- seq(-f.delta.e, f.delta.e, length = 100)
  gamma.e <- sqrt((1 - u * abs(delta.e))^2/(9 * c_alpha^2) - v * delta.e^2)  
  
  if (contour == TRUE){ 
    oldpar <- par(bg="#CCCCCC", mar = c(5, 4, 4, 3) + 0.1)
    point.delta <- round((mu-target)/d,4)
    point.gamma <- round(std.dev/d,4) 
    ymax <- point.gamma+0.2
    if (is.null(ylim)) ylim <- c(0,ymax)
    plot(delta.t, gamma.t, cex.lab = 0.7, cex.axis = 0.7, 
         type = "l", col = "blue", lwd = 2, ylim = ylim,
         axes = FALSE, xlab = "delta", ylab = "gamma",
         main = paste ("Contour plot:",names(ind)),...)
    
    rect(par("usr")[1],
         par("usr")[3],
         par("usr")[2],
         par("usr")[4],
         col  =  "white")
    box(col  =  "#CCCCCC")
    grid(col  =  "#CCCCCC")
    abline(h = 0)
    abline(v = 0)
    axis(1)
    axis(2)
    
    lines(delta.t, gamma.t, col = "black", lwd = 2,lty = 1)
    lines(delta.e, gamma.e, col = "blue", lwd = 2,lty = 2)
    points(x = point.delta,y = point.gamma,col="red",pch = 21, bg = "red", 
           lwd=2)
    legend("topleft", c("Theory", "Empirical"), 
           lwd = c(2,2), bty="n",lty = c(1,2), col =c("black","blue"), cex =0.7)
    par(oldpar)
  }
  
  result <- list(round(ind,4), 
              delta = round((mu-target)/d,4),
              gamma = round(std.dev/d,4),delta.t = delta.t,gamma.t = gamma.t,
              delta.e = delta.e,gamma.e = gamma.e)
  invisible(result)
}