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
# Main function to create a 'mqcd' object
#-----------------------------------------------------------------------------#
##' It creates a data object to be used in Multivariante Quality Control
##' 
##' Create an object of class 'mqcd' to perform statistical quality control.
##' This object is used to plot Multivariate Control Charts.
##' 
##' @aliases mqcd 
##' @param x A matrix, a data-frame or an array where it should contain data.
##' @param data.name  A string that specifies the title displayed on the plots. 
##' If not provided it is taken from the name of the object \code{x}.
##' @export
##' @examples
##' library(qcr)
##' data(dowel1)
##' str(dowel1)
##' data.mqcd <- mqcd(dowel1)
##' str(data.mqcd)

mqcd <- function(x, data.name = NULL)
  #.........................................................................
{

  if (!is.matrix(x) & !is.data.frame(x) & !is.array(x))
    stop("object must be a matrix or data.frame or array")
  
  p <- ncol(x) # quality characteristics
  m <- nrow(x) # number of samples or observations
  if (inherits(x, "matrix") || inherits(x, "data.frame")) {
      names <- colnames(x)    
      x <- array(data.matrix(x),c(m,p,1))
      colnames(x) <- names        
  }    
  n <- dim(x)[3] # observations or sample size 
  
  if (is.null(data.name))
    data.name <- " DATA"
  
  result <- x
  
  attr(result, "data.name") <- data.name
  attr(result, "type.data") <- "Multivariate"
  
  oldClass(result) <- c("mqcd", "array")
  
  return(result)
} # mqcd
#.........................................................................
