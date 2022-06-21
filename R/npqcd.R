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
# Main function to create a 'npqcd' object
#-----------------------------------------------------------------------------#
##' It creates a data object for Non Parametric Quality Control
##' 
##' It creates an object of class 'npqcd' to perform statistical quality control.
##' This object is used to plot Non Parametric Multivariate Control Charts.
##' 

##' @aliases npqcd 
##' @param x A matrix or data-frame or array which it should contain data. 
##' Dimension has to be the same as that of the observations.
##' @param G  The x as a matrix, data frame or list. If it is a matrix or data frame,
##' then each row is viewed as one multivariate observation.
##' @param data.name  A string that specifies the title displayed on the plots. 
##' If not provided it is taken from the name of the object x.
##' @export
##' @examples
##' 
##' library(qcr)
##' 
##' set.seed(356)
##' data <- matrix(rnorm(999), nc = 3)
##' x <-rexp(999,0.5) 
##' x <-matrix(x,ncol=3) 
##' data.npqcd <- npqcd(data,x)
##' str(data.npqcd)

npqcd <- function(x, G = NULL, data.name = NULL)
#.........................................................................
{
  if (!is.matrix(x) & !is.data.frame(x) & !is.array(x))
    stop("object must be a matrix or data.frame or array")
  
  if (is.null(G)) G <- x
  
  if (!is.matrix(G) & !is.data.frame(G) & !is.list(G))
    stop("object must be a matrix or data.frame or list")
  
  

  if (inherits(x, "matrix") || inherits(x, "data.frame")) {
    p <- ncol(x) # quality characteristics
    m <- nrow(x) # number of samples or observations
    names <- colnames(x)    
    x <- array(data.matrix(x),c(m,p,1))
    colnames(x) <- names        
  }   

  if(is.list(G))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  {
    m=length(G)
    n=length(G[[1]])
    for(i in 1:m){    
      if(length(G[,,i])!=n)
        stop("When using a list, each element must be a matrix of the same length.") 
    }
  }

    if (is.data.frame(G)) as.matrix(G)
    
    if (is.matrix(G)){
      n=dim(G)
      m=dim(x)[1:2]
      if(sum(n!=m)==2)
          stop("The object x must be of the same dimension than observations.")       
    }
      

  if (is.null(data.name))
     data.name <- "DATA"
   
  result <- list(x = x, G = G)
  
  attr(result, "data.name") <- data.name
  attr(result, "type.data") <- "Multivariate"
  
  oldClass(result) <- c("npqcd", "list")
  
   
  return(result)
  
} # npqcd
#.........................................................................
