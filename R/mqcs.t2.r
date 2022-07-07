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
# t2 chart
#-------------------------------------------------------------------------
##' Function to plot t2 Hotelling chart
##'
##' This function is used to compute statistics required by the t2 of HOTELLING 
##' or Shewhart Multivariate chart.
##'
##' @param x   An object of class 'mqcd'
##' @param ... Arguments passed to or from methods.
##' @seealso
##' \code{\link{mqcd}}, \code{\link{mqcs}}
##' @export
## @references Montgomery, D.C. (2000)
##' @examples
##' 
##' ##
##' ##  Continuous data 
##' ##
##' library(qcr)
##' data(dowel1)
##' str(dowel1)
##' data.mqcd <- mqcd(dowel1)
##' res.mqcs <- mqcs.t2(data.mqcd)
##' summary(res.mqcs)
##' plot(res.mqcs, title =" Hotelling Control Chart for dowel1")
##'
##' data(archery1)
##' str(archery1)
##' data.mqcd <- mqcd(archery1)
##' res.mqcs <- mqcs.t2(data.mqcd)
##' summary(res.mqcs)
##' plot(res.mqcs, title =" Hotelling Control Chart for archery1")

mqcs.t2 <- function(x, ...) {
  UseMethod("mqcs.t2")
}

##' @rdname mqcs.t2
##' @method mqcs.t2 default
##' @inheritParams mqcd
##' @param limits A two-values vector specifying the control limits.
##' @param Xmv The mean vector. It is only specified for Phase II or when 
##' the parameters of the distribution are known.
##' @param S The sample covariance matrix. It is only used for Phase II or 
##' when the parameters of the distribution are known.
##' @param colm The number of samples (m) and it is only used in Hotelling 
##' control chart for Phase II.
##' @param alpha It is the the significance level (0.01 for default)
##' @param phase Allows to select the type of UCL to use. Only values of phase = 1 or 2 are allowed.
##' @param method The method employed to compute the covariance matrix
##' in the individual observation case. Two methods are used "sw" 
##' for compute according to (Sullivan,Woodall 1996a) and "hm" 
##' by (Holmes,Mergen 1993)
##' @param plot Logical value. If \code{TRUE} a t2 chart should be plotted. 
##' @author Edgar Santos-Fernandez
##' @export
##' 
mqcs.t2.default <- function(x, data.name = NULL, limits = NULL, Xmv = NULL, S = NULL,
                            colm = NULL, alpha = 0.01,
                            phase = 1, method = "sw", plot = FALSE, ...)
#.........................................................................
  {
  
  obj<-mqcd(data= x, data.name = data.name)

  result<-mqcs.t2.mqcd(x = obj, data.name = data.name, limits = NULL, Xmv = Xmv, 
                       S = S, colm = colm, alpha = alpha,
                       phase = phase, method = method, plot = plot, ...)

  return(result)
} # mqcs.t2.default
#.........................................................................

##' @rdname  mqcs.t2
##' @method mqcs.t2 mqcd
##' @inheritParams mqcs.t2.default
##' @export
##' 

#x <- mqcd(datos2)
#alpha <- 0.00135

mqcs.t2.mqcd <- function(x, limits = NULL, Xmv = NULL, S = NULL, colm = NULL, 
                         alpha = 0.01,
                         phase = 1, method = "sw", plot = FALSE, ...) 
#.........................................................................  
{
  
  if(is.null(x) || !inherits(x, "mqcd"))
    stop("data must be an objects of class (or extending) 'mqcd'")
  
#  if(!missing(Xmv))(phase <- 2)
  
  mqcs<-mqcs(x, method)
  if(is.null(Xmv)) Xmv <- mqcs$mean 
  if(is.null(S)) S <- mqcs$S
  if(is.null(colm)) colm <- nrow(x)
  x.jk <- mqcs$mean.jk
  
  p <- ncol(x) # quality characteristics
  m <- nrow(x) # number of samples or observations
  n <- dim(x)[3] # observations or sample size 
  
  
  statistics <- matrix(0,m,1)
  
  for (ii in 1 : m){
    statistics[ii,1] <- n * t(x.jk[ii,] - Xmv) %*% solve(S) %*% (x.jk[ii,] - Xmv)
  }
 
  
if (is.null(limits)){

  ifelse(n == 1, ifelse(phase == 1, 
                        ucl <- ((colm - 1) ^ 2) / colm * qbeta(1 - alpha,p / 2,(((2 * (colm - 1) ^ 2) / (3 * colm - 4) - p - 1) / 2)),
                        ucl<-((p * (colm + 1) * (colm - 1)) / ((colm ^ 2) - colm * p)) * qf(1 - alpha,p,colm - p)),
         ifelse(phase == 1, 
                ucl <- (p * (colm - 1) * (n - 1)) / (colm * n - colm - p + 1) * qf(1 - alpha,p,colm * n - colm - p + 1),
                ucl <- (p * (colm + 1) * (n - 1)) / (colm * n - colm - p + 1) * qf(1 - alpha,p,colm * n - colm - p + 1))
  )  
  
limits <- c(lcl = 0, ucl = ucl)

}
  
  violations <- which(statistics > limits[2])

  data.name <- attr(x, "data.name")
  result <- list(mqcd  =  x, type  =  "t2", statistics  =  statistics,
                 mean  =  Xmv, S  =  S, alpha = alpha,
                 limits  =  limits, data.name  =  data.name,
                 violations  =  violations)
  
  oldClass(result) <- c("mqcs.t2", "mqcs")
  
  if(plot) plot(result, ...)
  
  return(result)
#.........................................................................
} # mqcs.t2.mqcd
#.........................................................................

