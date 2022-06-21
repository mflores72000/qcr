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
# Q chart
#-------------------------------------------------------------------------
##' Function to plot the Q chart
##'
##' This function is used to compute statistics required by the Q chart.
##' 
##' @param x  An object npqcd (Non parametric Quality Control Data)
##' @param ... Arguments passed to or from methods.
##' @export
##' @references Regina Liu (1995)
##' @examples
##' \dontrun{ 
##' ##
##' ##  Continuous data 
##' ##
##' library(qcr)
##' set.seed(12345)
##' mu<-c(0,0)
##' Sigma<- matrix(c(1,0,0,1),nrow = 2,ncol = 2)
##' u <- c(2,2)
##' S <- matrix(c(4,0,0,4),nrow = 2,ncol = 2)
##' G <- rmvnorm(540, mean = mu, sigma = Sigma)
##' x<- rmvnorm(40,mean=u,sigma = S)
##' x <- rbind(G[501:540,],x)
##' n <- 4 # samples
##' m <- 20  # measurements
##' k <- 2  # number of variables
##' x.a <- array(,dim=c(n,k,m))
##' for (i in 1:m){
##' x.a[,,i] <- x[(1+(i-1)*n):(i*n),] }
##' M <- G[1:500,]
##' data.npqcd <- npqcd(x.a,M)
##' str(data.npqcd)
##' res.npqcs <- npqcs.Q(data.npqcd,method = "Liu", alpha=0.025)
##' str(res.npqcs)
##' summary(res.npqcs)
##' plot(res.npqcs,title =" Q Control Chart")}

npqcs.Q <- function(x, ...) {
  UseMethod("npqcs.Q")
}

##' @rdname npqcs.Q
##' @method npqcs.Q default
##' @inheritParams npqcd
##' @param limits A two-value vector specifying the control limits lower and central.
##' @param method Character string which determines the depth function used. 
##' \code{method} can be "Tukey" (the default), "Liu", "Mahalanobis", "RP" Random Project or "LD" Likelihood depth.
##' @param alpha It is the significance level (0.01 for default)
##' @param plot Logical value. If TRUE a Q chart should be plotted. 
##' @export
 
npqcs.Q.default <- function(x, G, data.name = NULL, limits = NULL,
                            method = c("Tukey","Liu","Mahalanobis","RP","LD"), alpha = 0.01, plot = FALSE, ...)
#.........................................................................
  {
  
  method <- match.arg(method)
  
  obj <- npqcd(x = x , G = G, data.name = data.name)

  result <- npqcs.Q.npqcd(x = obj, data.name = data.name, limits = limits, method = method, alpha = alpha,
                       plot = plot, ...)

  return(result)
} # npqcs.Q.default
#.........................................................................

##' @rdname  npqcs.Q
##' @method npqcs.Q npqcd
##' @export

npqcs.Q.npqcd <- function(x, data.name, limits = NULL, method = c("Tukey","Liu","Mahalanobis","RP","LD"), alpha = 0.01, plot = FALSE, ...) 
#.........................................................................  
{
  
  if(is.null(x) || !inherits(x, "npqcd"))
    stop("x must be an objects of class (or extending) 'npqcd'")
  
  method <- match.arg(method)
  
  n <- dim(x[[1]])[1]
  m <- dim(x[[2]])[1]
  
  if(n == 1)
    stop("The dimension of object x must be valid for control chart Q")

  if (is.null(limits)){
    central<-0.5
    zalpha<-qnorm(1-alpha,0,1)
    if(n<5){
      lcl<-((factorial(n)*alpha)^(1/n))/n
    }
    if(n>=5){
      lcl<-0.5-zalpha*sqrt(1/12*((1/n)+(1/m)))
    }
    limits <- c(lcl = lcl, cl = central)  
  }else{
    limits <- limits
  }
  
  npqcs <- npqcs(x, method)
  statistics <- apply(npqcs$rank.depth,2,mean)
  violations <- which(statistics < lcl)  

  data.name <- attr(x, "data.name")
  result <- list(npqcd  =  x, type  =  "Q", statistics  =  statistics, alpha = alpha,
                 limits  =  limits, data.name  =  data.name, method = method,
                 violations  =  violations)
  
  oldClass(result) <- c("npqcs.Q", "npqcs")
  
  if(plot) plot(result, ...)
  
  return(result)
#.........................................................................
} # npqcs.Q.npqcd
#.........................................................................