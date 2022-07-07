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
# Main function to create a 'qcs.cp' object
#-----------------------------------------------------------------------------#
##' Process capability indices (parametric)
##' 
##' Calculates \eqn{Cp}{Cpk}, \eqn{Cpm}{Cpmk} using the formulation described by Kerstin Vannman(1995).
##' @aliases qcs.cp 
##' @param object qcs object of type \code{"qcs.xbar"} or \code{"qcs.one"}.
##' @param parameters A vector specifying the \code{u} and \code{v} parameters values. 
##' If \code{parameters} = c(u=0, v=0), the Cp indice is calculed; 
##' If \code{parameters} = c(u=1, v=0), the Cpk indice is calculed;
##' If \code{parameters} = c(u=0, v=1), the Cpm indice is calculed; 
##' If \code{parameters} = c(u=1, v=1), the Cpmk indice is calculed.  
##' @param limits A vector specifying the lower and upper specification limits.
##' @param target A value specifying the target of the process. 
##' If it is \code{NULL}, the target is set at the middle value between specification limits.
##' @param mu A value specifying the mean of data.
##' @param std.dev A value specifying the within-group standard deviation.
##' @param nsigmas A numeric value specifying the number of sigmas to use.
##' @param k A numeric value. If the capacity index exceeds the \code{k} value, 
##' then the process is capable.
##' @param contour Logical value indicating whether contour graph should be plotted.
##' @param ylim The 'y' limits of the plot.
##' @param ... Arguments to be passed to or from methods.
##' @export
##' @references 
##' Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd
##' ed, New York, John Wiley & Sons. \cr
##' Vannman, K (1995) \emph{A Unified Approach to Capability Indices}. Statitica Sinica,5,805-820.\cr
##' @examples
##' library(qcr)
##' data(pistonrings) 
##' xbar <- qcs.xbar(pistonrings[1:125,],plot = TRUE)
##' mu <-xbar$center
##' std.dev <-xbar$std.dev
##' LSL=73.99; USL=74.01
##'qcs.cp(parameters = c(0,0),limits = c(LSL,USL),
##'       mu = mu,std.dev = std.dev,ylim=c(0,1))
##' #calculating all the indices
##' qcs.cp(object = xbar,parameters = c(0,0), limits = c(LSL,USL),ylim=c(0,1))
##' qcs.cp(object = xbar,parameters = c(1,0), limits = c(LSL,USL),ylim=c(0,1))
##' qcs.cp(object = xbar,parameters = c(0,1), limits = c(LSL,USL),ylim=c(0,1))
##' qcs.cp(object = xbar,parameters = c(1,1), limits = c(LSL,USL),ylim=c(0,1))
qcs.cp <- function(object, parameters = c(u = 0,v = 0), limits = c(lsl = -3, usl = 3), 
                   target = NULL, mu = 0, std.dev = 1, nsigmas = 3, k = 1,contour =TRUE, ylim = NULL,...){
  
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
  u <- parameters[1]
  v <- parameters[2]
  
  ind <-  (d-u*abs(mu-m))/(nsigmas*std.dev*sqrt(1+v*((mu-target)/std.dev)^2))
  
  if (u == 0 & v == 0) {
    names(ind) <- c("Cp")
    delta <- seq(-1, 1, length = 100)
  }else{
    if (u == 1 & v == 0) names(ind) <- c("Cpk")
    if (u == 0 & v == 1) names(ind) <- c("Cpm")
    if (u == 1 & v == 1) names(ind) <- c("Cpmk")
    f.delta <- 1/(u + 3 * k * sqrt(v))
    delta <- seq(-f.delta, f.delta, length = 100)
  }
  gamma <- sqrt((1 - u * abs(delta))^2/(9 * k^2) - v * delta^2)  
  if (contour == TRUE){ 
    oldpar <- par(bg="#CCCCCC",mar = c(5, 4, 4, 3) + 0.1)
    point.delta <- round((mu-target)/d,4)
    point.gamma <- round(std.dev/d,4) 
    ymax <- point.gamma + 0.2
    if (is.null(ylim)) ylim <- c(0,ymax)
    plot(delta, gamma, cex.lab = 0.7, cex.axis = 0.7, 
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
    
    lines(delta, gamma, col = "blue", lwd = 2)
    points(x = point.delta,y = point.gamma,col="red",pch = 21, bg = "red", lwd=2)
    par(oldpar)
  }

  result <- c(round(ind,4), 
              delta = round((mu-target)/d,4),
              gamma = round(std.dev/d,4))
  return(result)
}
