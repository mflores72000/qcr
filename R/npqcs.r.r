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
# r chart
#-------------------------------------------------------------------------
##' Function to plot the r chart
##'
##' This function is used to compute statistics required by the r chart.
##' 
##' @param x  An object npqcd (Non parametric Quality Control Data)
##' @param ... Arguments passed to or from methods.
##' @export
##' @references Regina Liu (1995)
##' @examples
##' \dontrun{
##' library(qcr)
##' set.seed(356)
##' mu<-c(0,0)
##' Sigma<- matrix(c(1,0,0,1),nrow = 2,ncol = 2)
##' u <- c(2,2)
##' S <- matrix(c(4,0,0,4),nrow = 2,ncol = 2)
##' G <- rmvnorm(540, mean = mu, sigma = Sigma)
##' x<- rmvnorm(40,mean=u,sigma = S)
##' x <- rbind(G[501:540,],x)
##' M <- G[1:500,]
##' data.npqcd <- npqcd(x,M)
##' str(data.npqcd)
##' res.npqcs <- npqcs.r(data.npqcd,method = "Liu", alpha=0.025)
##' str(res.npqcs)
##' summary(res.npqcs)
##' plot(res.npqcs,title =" r Control Chart")  }

npqcs.r <- function(x, ...) {
  UseMethod("npqcs.r")
}

##' @rdname npqcs.r
##' @method npqcs.r default
##' @inheritParams npqcd
##' @param limits A two-value vector specifying the control limits lower and central.
##' @param method Character string which determines the depth function used. 
##' \code{method} can be "Tukey" (the default), "Liu", "Mahalanobis", "RP" Random Project or "LD" Likelihood depth.
##' @param alpha It is the significance level (0.01 for default)
##' @param plot Logical value. If TRUE a r chart should be plotted. 
##' @export

npqcs.r.default <- function(x, G, data.name = NULL, limits = NULL,
                            method = c("Tukey","Liu","Mahalanobis","RP","LD"), alpha = 0.01, plot = FALSE, ...)
#.........................................................................
  {
  
  method <- match.arg(method)
  
  obj <- npqcd(x = x , G = G, data.name = data.name)

  result <- npqcs.r.npqcd(x = obj, data.name = data.name, limits = NULL, method = method, alpha = alpha,
                       plot = plot, ...)

  return(result)
} # npqcs.r.default
#.........................................................................

##' @rdname  npqcs.r
##' @method npqcs.r npqcd
##' @export
##' 

# x <-data.npqcd
# method="Tukey"
# alpha=0.5
npqcs.r.npqcd <- function(x, data.name, limits = NULL, method = c("Tukey","Liu","Mahalanobis","RP","LD"), alpha = 0.01, plot = FALSE, ...) 
#.........................................................................  
{
  
  if(is.null(x) || !inherits(x, "npqcd"))
    stop("x must be an objects of class (or extending) 'npqcd'")
  
  if(dim(x[[1]])[3]!=1)
    stop("The dimension of object x must be valid for control chart r")
  
  method <- match.arg(method)
  
  if (is.null(limits)){
    central<-0.5 
    limits <- c(lcl = alpha, cl = central)  
  }else{
    limits <- limits
  }
  
  npqcs <- npqcs(x, method)
  statistics <- npqcs$rank.depth
  depth.data <- npqcs$depth.data
  violations <- which(statistics < alpha)  
  
  
  data.name <- attr(x, "data.name")
  result <- list(npqcd  =  x, type  =  "r", depth.data = depth.data, statistics  =  statistics, alpha = alpha,
                 limits  =  limits, data.name  =  data.name, method = method,
                 violations  =  violations)
  
  oldClass(result) <- c("npqcs.r", "npqcs")
  
  if(plot) plot(result, ...)
  
  return(result)
#.........................................................................
} # npqcs.r.npqcd
#.........................................................................