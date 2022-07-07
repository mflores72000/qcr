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
# S chart
#-------------------------------------------------------------------------
##' Function to plot the S chart
##'
##' This function is used to compute statistics required by the Non Parametric S chart.
##' 
##' @param x  An object of class "npqcd".
##' @param ... Arguments passed to or from methods.
##' @export
##' @references Regina Liu (1995)
##' @examples
##' \dontrun{
##' ##
##' ##  Continuous data 
##' ##
##' set.seed(12345)
##' mu<-c(0,0)
##' Sigma<- matrix(c(1,0,0,1),nrow = 2,ncol = 2)
##' u <- c(2,2)
##' S <- matrix(c(4,0,0,4),nrow = 2,ncol = 2)
##' G <- rmvnorm(540, mean = mu, sigma = Sigma)
##' x<- rmvnorm(40,mean=u,sigma = S)
##' x.a <- rbind(G[501:540,],x)
##' M <- G[1:500,]
##' data.npqcd <- npqcd(x.a,M)
##' res.npqcs <- npqcs.S(data.npqcd,method = "Liu", alpha=0.05)
##' summary(res.npqcs)
##' plot(res.npqcs,title =" S Control Chart")}

npqcs.S <- function(x, ...) {
  UseMethod("npqcs.S")
}

##' @rdname npqcs.S
##' @method npqcs.S default
##' @inheritParams npqcd
##' @param limits A two-value vector specifying the control limits lower and central.
##' @param method Character string which determines the depth function used. 
##' \code{method} can be "Tukey" (the default), "Liu", "Mahalanobis", "RP" Random Project or "LD" Likelihood depth.
##' @param alpha It is the significance level (by default \code{alpha} = 0.01)
##' @param plot Logical value. If \code{TRUE} a S chart should be plotted.
##' @param standardize A logical value indicating data should be standardized. 
##' @export
 
npqcs.S.default <- function(x, G, data.name = NULL, limits =NULL,
                            method = c("Tukey","Liu","Mahalanobis","RP","LD"), alpha = 0.01, plot = FALSE, standardize = FALSE, ...)
#.........................................................................
  {
  
  method <- match.arg(method)
  
  obj <- npqcd(x = x , G = G, data.name = data.name)

  result <- npqcs.S.npqcd(x = obj, data.name = data.name, limits = limits, method = method, alpha = alpha,
                       plot = plot, standardize = standardize, ...)

  return(result)
} # npqcs.S.default
#.........................................................................

##' @rdname  npqcs.S
##' @method npqcs.S npqcd
##' @export 


npqcs.S.npqcd <- function(x, data.name, limits =NULL, method = c("Tukey","Liu","Mahalanobis","RP","LD"), alpha = 0.01, plot = FALSE, standardize = F, ...) 
#.........................................................................  
{
  
  if(is.null(x) || !inherits(x, "npqcd"))
    stop("x must be an objects of class (or extending) 'npqcd'")
  
#  if(dim(x[[1]])[3]==1)
#    stop("The dimension of object x must be valid for control chart S")
  method="Liu"
  #x <- data.npqcd
  #alpha=0.025
  method <- match.arg(method)
  m <- dim(x[[2]])[1]
  n <- 1:dim(x[[1]])[1]
  if (is.null(limits)){
    central<-0 
    zalpha<-qnorm(1-alpha,0,1)

    if(m<=30){
      lcl <- -(zalpha*sqrt((n^2)*((1/m)+(1/n))/12))
    }else{
      lcl <-  -zalpha
    }

    limits <- c(cl = central, lcl = lcl)  
  }else{
    limits <- limits
  }

  npqcs <- npqcs(x, method)
  statistics <- npqcs$rank.depth  
  
  b<-c(statistics-1/2)
  sumaacumulada<-cumsum(b)
  
  if(standardize){
    Sn <- sumaacumulada  
  } else {
    Sn<-sumaacumulada/sqrt((n^2)*((1/m)+(1/n))/12)
  }
  
  violations <- which(Sn < lcl)
  
  data.name <- attr(x, "data.name")
  result <- list(npqcd  =  x, type  =  "S", statistics  =  Sn, alpha = alpha,
                 limits  =  limits, data.name  =  data.name, method = method,
                 violations  =  violations)
  
  oldClass(result) <- c("npqcs.S", "npqcs")
  
  if(plot) plot(result, ...)
  
  return(result)
#.........................................................................
} # npqcs.S.npqcd
#.........................................................................