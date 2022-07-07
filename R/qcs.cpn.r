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
#-----------------------------------------------------------------------------#
# Main function to create a 'qcs.cpn' object
#-----------------------------------------------------------------------------#
##' Process capability indices (Nonparametric)
##' 
##' Calculates \eqn{CNp}{CNpk}, \eqn{CNpm}{CNpmk} using the formulation 
##' described by Tong and Chen (1998).
##' @aliases qcs.cpn 
##' @param object qcs object of type \code{"qcs.xbar"} or \code{"qcs.one"}.
##' @param parameters A vector specifying the \code{u} and \code{v} parameters values. 
##' If \code{parameters} = c(u=0, v=0), the CNp indice is calculed; 
##' If \code{parameters} = c(u=1, v=0), the CNpk indice is calculed;
##' If \code{parameters} = c(u=0, v=1), the CNpm indice is calculed; 
##' If \code{parameters} = c(u=1, v=1), the CNpmk indice is calculed.  
##' @param limits A vector specifying the lower and upper specification limits.
##' @param q A vector specifying the lower and upper quantiles. These values are 
##' necessary, if \code{object} value is missing. 
##' @param target A value specifying the target of the process. 
##' If it is \code{NULL}, the target is set at the middle value between specification limits.
##' @param median A value specifying the median of data.
##' @param nsigmas A numeric value specifying the number of sigmas to use.
##' @param confidence A numeric value between 0 and 1 specifying the probabilities 
##' for computing the quantiles.
##' This values is used only when \code{object} values is provided. 
##' By default \code{confidence}=0.9973.
##' @export
##' @references 
##' Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd
##' ed, New York, John Wiley & Sons. \cr
##' Tong, L.I. and Chen, J.P. (1998), \emph{Lower confidence limits of process capability 
##' indices for nonnormal process distributions.} International Journal of Quality & Reliability Management, 
##' Vol. 15 No. 8/9, pp. 907-19.\cr
##' @examples
##' library(qcr)
##' ##' data(pistonrings) 
##' xbar <- qcs.xbar(pistonrings[1:125,],plot = TRUE)
##' x<-xbar$statistics[[1]]
##' LSL=73.99; USL=74.01
##' median <-median(x)
##' lq=as.numeric(quantile(x,probs=0.00135))
##' uq=as.numeric(quantile(x,probs=0.99865))
##' qcs.cpn(parameters = c(0,0),limits = c(LSL,USL),
##'        median = median, q=c(lq,uq))
### all capacibility indices
##' qcs.cpn(object = xbar,parameters = c(0,0), limits = c(LSL,USL))
##' qcs.cpn(object = xbar,parameters = c(1,0), limits = c(LSL,USL))
##' qcs.cpn(object = xbar,parameters = c(0,1), limits = c(LSL,USL))
##'qcs.cpn(object = xbar,parameters = c(1,1), limits = c(LSL,USL))
qcs.cpn <- function(object, parameters = c(u = 0,v = 0), limits = c(lsl = -3, usl = 3), 
                    q = c(lq = -3, uq = 3),
                    target = NULL, median = 0,  nsigmas = 3,confidence = 0.9973){
  
  if (!missing(object)){
    
    if (!inherits(object, "qcs"))
      stop("an object of class 'qcs' is required")
    if (!(object$type == "xbar" | object$type == "one")) 
      stop("Process Capability Analysis only available for charts type 
           \"qcs.xbar\" and \"qcs.one\" charts")
    q1<-(1-confidence)/2
    q2<-confidence+q1
    x <- object[[3]][,1]
    F2=as.numeric(quantile(x,probs=q2))
    F1=as.numeric(quantile(x,probs=q1))
    median <- median(x)
  }else{
    F1 <- q[1]
    F2 <- q[2]
  }
  
  if (nsigmas <= 0) 
    stop("nsigmas must be a value positive")
  
  confidence = 1- 2*pnorm(-nsigmas)
  
  std.dev <- (F2-F1)/6
  
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
  
  
  ind <-  (d-u*abs(median-m))/(nsigmas*std.dev*sqrt(1+v*((median-target)/std.dev)^2))
  
  if (u == 0 & v == 0) names(ind) <- c("CNp")
  if (u == 1 & v == 0) names(ind) <- c("CNpk")
  if (u == 0 & v == 1) names(ind) <- c("CNpm")
  if (u == 1 & v == 1) names(ind) <- c("CNpmk")
  
  result <-round(ind,4)
  
  return(result)
}