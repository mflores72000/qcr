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
# mcusum chart
#-----------------------------------------------------------------------------#
##' Function to plot mcusum chart
##'
##' This function is used to compute statistics required by the mcusum chart.
##'
##' @param x   An R object (used to select the method). See details.
##' @param ... Arguments passed to or from methods.
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
##' res.mqcs <- mqcs.mcusum(data.mqcd)
##' summary(res.mqcs)
##' plot(res.mqcs, title =" MCUSUM Control Chart for dowel1")

#-----------------------------------------------------------------------------#
# function mqcs.mcusum
#-----------------------------------------------------------------------------#
mqcs.mcusum <- function(x, ...) {
  UseMethod("mqcs.mcusum")
}

#-----------------------------------------------------------------------------#
# function mqcs.mcusum.default
#-----------------------------------------------------------------------------#
##' @rdname mqcs.mcusum
##' @method mqcs.mcusum default
##' @inheritParams mqcd
##' @param limits A two-values vector specifying the control limits.
##' @param Xmv Is the mean vector. It is only specified for Phase II or when 
##' the parameters of the distribution are known.
##' @param S Is the sample covariance matrix. It is only used for Phase II or 
##' when the parameters of the distribution are known.
##' @param k Is a constant used in MCUSUM chart. Frequently k = 0.5
##' @param h Is a constant used in MCUSUM chart. Usually h = 5.5
##' @param method Is the method employed to compute the covariance matrix
##' in the individual observation case. Two methods are used "sw" 
##' for compute according to (Sullivan,Woodall 1996a) and "hm" 
##' by (Holmes,Mergen 1993)
##' @param plot Logical value. If TRUE a mcusum chart should be plotted. 
##' @author Edgar Santos-Fernandez
##' @export
##' 
mqcs.mcusum.default <- function(x, data.name = NULL, limits = NULL, Xmv = NULL, S = NULL,
                            k = 0.5, h= 5.5,
                            method = "sw", plot = FALSE, ...)
#.........................................................................
  {
  
  obj <- mqcd(data= x, data.name = data.name)

  result <- mqcs.mcusum.mqcd(x = obj, Xmv = Xmv, limits = NULL,
                       S = S, k = k, h = h,
                       method = method, plot = plot, ...)

  return(result)
} # mqcs.mcusum.default
#.........................................................................

#-----------------------------------------------------------------------------#
# function mqcs.mcusum.mqcd
#-----------------------------------------------------------------------------#
##' @rdname  mqcs.mcusum
##' @method mqcs.mcusum mqcd
##' @inheritParams mqcs.mcusum.default
##' @export
##' 


mqcs.mcusum.mqcd <- function(x, limits = NULL, Xmv = NULL, S = NULL, 
                            k = 0.5, h = 5.5,
                            method = "sw", plot = FALSE, ...) 
#.........................................................................  
{
  
  if(is.null(x) || !inherits(x, "mqcd"))
    stop("data must be an objects of class (or extending) 'mqcd'")
    
  mqcs<-mqcs(x, method)
  if(is.null(Xmv)) Xmv <- mqcs$mean 
  if(is.null(S)) S <- mqcs$S
  x.jk <- mqcs$mean.jk
  
  p <- ncol(x) # quality characteristics
  m <- nrow(x) # number of samples or observations
  n <- dim(x)[3] # observations or sample size 
    
  statistics <- matrix(0,m,1)
  
  if (is.null(limits)) limits <- c(lcl = 0, ucl = h)
  
  dif <- sweep(x.jk,2,Xmv)
  s <- matrix(0,m,p)
  ci <- matrix(0,m,1)
  
  ci[1] <- sqrt(dif[1,] %*% solve((S / n)) %*% dif[1,])
  
  if (ci[1] > k) { s[1,] <- (s[1,] + dif[1,]) * (1 - k / ci[1])
  }else(s[1,] = matrix(0,ncol = p)) #compute s
  
  for (i in 2:m){ 
    ci[i,]=sqrt((s[i - 1,] + dif[i,]) %*% solve(S / n) %*% (s[i - 1,] + dif[i,])) 
    if (ci[i] > k){ s[i,] = (s[i - 1,] + dif[i,]) * (1 - k / ci[i])} 
    else {s[i,] = matrix(0,ncol = p)}
  }
  
  for (i in 1:m){
    statistics[i]=sqrt(s[i,]%*%solve((S/n))%*%(s[i,]))
  } 
  
  
  violations <- which(statistics > limits[2])
    

  data.name <- attr(x, "data.name")
  result <- list(mqcd  =  x, type  =  "mcusum", statistics  =  statistics,
                 mean  =  Xmv, S  =  S, k = k, h = h,
                 limits  =  limits, data.name  =  data.name,
                 violations  =  violations)
  
  oldClass(result) <- c("mqcs.mcusum", "mqcs")
  
  if(plot) plot(result, ...)
  
  return(result)
#.........................................................................
} # mqcs.mcusum.mqcd
#.........................................................................
