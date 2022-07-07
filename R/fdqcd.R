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
# Main function to create a 'fdqcd' object
#-----------------------------------------------------------------------------#
##' It creates a data object to be used in Functional Data Quality Control
##' 
##' Create an object of class 'fdqcd' to perform statistical quality control.
##' This object is used to plot Functional Data Control Charts.
##' 
##' @aliases fdqcd 
##' @param x Matrix of set cases with dimension (n x m), where 'n' is the number of curves 
##' and 'm' are the points observed in each curve.
##' @param data.name  A string that specifies the title displayed on the plots. 
##' If not provided it is taken from the name of the object \code{x}.
##' @param ... Arguments passed to or from methods.
##' @export
##' @examples
##' library(qcr)
##' m <- 30
##' tt<-seq(0,1,len=m)
##' mu<-30 * tt * (1 - tt)^(3/2)
##' n0 <- 100
##' set.seed(12345)
##' mdata<-matrix(NA,ncol=m,nrow=n0)
##' sigma <- exp(-3*as.matrix(dist(tt))/0.9)
##' for (i in 1:n0) mdata[i,]<- mu+0.5*mvrnorm(mu = mu,Sigma = sigma )
##' fdchart <- fdqcd(mdata)
##' plot(fdchart,type="l",col="gray")

fdqcd <- function(x, data.name = NULL,...)
  #.........................................................................
{

  if (!is.matrix(x) & !is.data.frame(x))
    stop("object must be a matrix or data.frame")
  
  p <- ncol(x) # quality characteristics
  m <- nrow(x) # number of samples or observations
  
  if (is.null(data.name))
    data.name <- " Phase I: Fda Chart"
  
  result <- fdata(x,...)
  
  attr(result, "data.name") <- data.name
  attr(result, "type.data") <- "Functional"
  
  oldClass(result) <- c("fdata", "fdqcd")
  
  return(result)
} # fdqcd
