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
# Main function to create a 'npqcs' object
#-----------------------------------------------------------------------------#
##' Statistical Quality Control Object
##' 
##' Create an object of class 'npqcs' to perform statistical quality control.
##' This function is used to compute statistics required to plot 
##' Non Parametric Multivariate Control Charts.
##' 
##' @aliases npqcs summary.npqcs print.npqcs
##' @param x  Object npqcd (Non Parametric Multivariante Quality Control Data)
##' @param method Character string which determines the depth function used. 
##' \code{method} can be "Tukey" (the default), "Liu", "Mahalanobis", "RP" Random Project or "LD" Likelihood depth.
##' @param ... Arguments passed to or from methods.
##' @export

npqcs <- function(x, method = c("Tukey","Liu","Mahalanobis", "RP", "LD"), ...)
  #.........................................................................  
  {
  if(is.null(x) || !inherits(x, "npqcd"))
    stop("x must be an objects of class (or extending) 'npqcd'")
  
  method <- match.arg(method)
  
  #npqcd <- data.npqcd
  npqcd <- x
  x <- npqcd[[1]]
  G <- npqcd[[2]]
  
  depth.data.x <- matrix(,nrow=dim(x)[1],ncol=dim(x)[3])
  
  for (i in 1:dim(x)[3]){
    depth.data.x[,i] <- f.depth(x[,,i], G, method)
  }
  
  depth.data.G <- f.depth(G, G, method)
  
  rank.depth.x <- vector (length = length(depth.data.x))
  depth.data.x <- as.vector(depth.data.x)

  for(i in 1:length(depth.data.x))  
  {  
    rank.depth.x[i]<-(sum(depth.data.G <= depth.data.x[i]))/length(depth.data.G) 
  }
  
  rank.depth.x <- matrix(rank.depth.x, nrow=dim(x)[1], ncol=dim(x)[3])
  depth.data.x <- matrix(depth.data.x, nrow=dim(x)[1], ncol=dim(x)[3])
  
  result <- list (npqcd = npqcd, depth.data = depth.data.x, rank.depth = rank.depth.x) 
  
  oldClass(result) <- c("npqcs")
  
  return(result)
} # npqcs
#.........................................................................

##' @export
##' @method print npqcs
print.npqcs <- function(x, ...) str(x,1)
#.........................................................................
##' @export
##' @method summary npqcs
summary.npqcs <- function(object, ...)
  #.........................................................................
{
  type <- object$type
  r <- object$statistics
  cat("\nSummary of group statistics:\n")
  print(summary(r))
  cat("\nNumber of quality characteristics: ", dim(object$npqcd[[1]])[2])
  cat("\nNumber of samples or observations: ", dim(object$npqcd[[1]])[1])
  cat("\nNumber of observations or sample size: ", dim(object$npqcd[[1]])[3])
  

  limits <- object$limits
  if (!is.null(limits)) 
  { cat("\nControl limits:", "\n") 
    print(limits)
  }
  

  if (length(object$violations)== 0){
    cat("\nNumber beyond limits: 0", "\n") 
  } 
  else {cat("\nBeyond limits of control:", "\n")
        print(object$statistics[object$violations])
  }
  
  invisible()
  #.........................................................................
} # summary.npqcs

f.depth <- function(x, G = NULL, method = c("Tukey","Liu","Mahalanobis", "RP", "LD")){
  
  if(is.null(G)) G <- x
  
  if (!is.matrix(G)) stop("object must be a matrix")  
  method <- match.arg(method)

  depth.data <- switch(method,"Tukey" =  {mdepth.TD(x,G)$dep},
                              "Liu" =   {mdepth.SD(x,G)$dep},
                              "Mahalanobis" = {mdepth.MhD(x,G)$dep},
                              "RP" = {mdepth.RP(x,G)$dep},
                              "LD" = {mdepth.LD(x,G)$dep})
  return (depth.data)
}



# npqcs.add function
#-------------------------------------------------------------------------
##' npqcs.add Add a matrix, data.frame or array object with a npqcs object
##' 
##' This function is used to join two objects of type matrix, data.frame or array and npqcs.
##' 
##' @param x   Object type npqcs
##' @export
##' 

npqcs.add <- function(x, ...){
  UseMethod("npqcs.add")
}

##' @rdname  npqcs.add 
##' @method npqcs.add default
##' @param value   Object type data.frame, matrix or array.
##' @param ...  Arguments to be passed to or from methods.
##' @export 


npqcs.add.default <- function(x, value, ...){
  
  
  if (!inherits(x, "npqcs"))
    stop("object must be npqcs")
  
  if (!is.matrix(value) & !is.data.frame(value) & !is.array(value))
    stop("object must be a matrix, data.frame or array")
  
  data <- x$npqcd$x
  limits <- x$limits
  data.name = x$data.name
  method = x$method
  type <- x$type
  G <- x$npqcd$G
  
  if (inherits(value, "matrix") || inherits(value, "data.frame")) {
    p <- ncol(value) # quality characteristics
    m <- nrow(value) # number of samples or observations
    names <- colnames(value)    
    value <- array(data.matrix(value),c(m,p,1))
    colnames(value) <- names        
  }   
  
  n1 <- dim(data)[3]
  n2 <- dim(value)[3]
  m1 <- dim(data)[1]
  m2 <- dim(value)[1]
  m <- m1 + m2
  k1 <- dim(data)[2]
  k2 <- dim(value)[2]
    
  if (n1 == n2) stop("The samples must be of the same dimension")
  if (k1 == k2) stop("The samples must be of the same quality characteristics")
  
  xx <- array(,dim = c(m,k1,n1))
  for (i in 1:n1 ){
    xx[,,i] <- cbind(data[,,i],value[,,i])     
  }
  
  z.npqcd <- npqcd(x = xx, G = G, data.name = data.name)
  
  z.npqcs <- switch(type, 
                  "r" = npqcs.r.npqcd(x = z.npqcd, limits = limits, method = method),
                  "Q" = npqcs.Q.npqcd(x = z.npqcd, limits = limits, method = method),
                  "S" = npqcs.S.npqcd(x = z.npqcd, limits = limits, method = method),
                  NULL)
  result <- z.npqcs
}


#-------------------------------------------------------------------------
# npstate.control
#-------------------------------------------------------------------------
##' Non parametric process state
##' 
##' This function removes observations from the sample which violates 
##' the rules of a process under control.
##' @aliases npstate.control
##' @param x  Object npqcd (Quality Control Statitical Non Parametric)
##' @param control A logical value indicating whether the initial sample comes 
##' from a process under control.
##' @export
##' @examples
##' \dontrun{
##' ##
##' ##  Continuous data 
##' ##
##' library(qcr)
##' set.seed(356)
##' mu<-c(0,0)
##' Sigma<- matrix(c(1,0,0,1),nrow = 2,ncol = 2)
##' mu <- c(2,2)
##' S <- matrix(c(4,0,0,4),nrow = 2,ncol = 2)
##' G <- rmvnorm(540, mean = mu, sigma = Sigma)
##' x<- rmvnorm(40,mean=mu,sigma = S)
##' x <- rbind(G[501:540,],x)
##' M <- G[1:500,]
##' data.npqcd <- npqcd(x,M)
##' str(data.npqcd)
##' res.npqcs <- npqcs.r(data.npqcd,method = "Liu", alpha=0.025)
##' str(res.npqcs)
##' summary(res.npqcs)
##' plot(res.npqcs)
##' new.npqcd <- npstate.control(x = res.npqcs)
##' res.npqcs <- npqcs.r(new.npqcd)
##' summary(res.npqcs)
##' plot(res.npqcs)  
##' }



npstate.control <- function(x, control = FALSE) 
  #.........................................................................  
{
  if (!inherits(x, "npqcs"))
    stop("an object of class 'npqcs' is required")
  if (length(x$violations)>0){
    ii<-x$violations  
    n <- dim(x$npqcd[[1]])[3]
    m <- dim(x$npqcd[[1]])[1]
    k <- dim(x$npqcd[[1]])[2]
    xx <- array(dim = c(m-length(ii),k,n))
    for (i in 1:n){ 
      xx[,,i] <- x$npqcd[[1]][,,i][-ii,]  
    }

    if (control == TRUE){ 
      G <- x$npqcd[[2]]
    }else{
      G <- x$npqcd[[2]][-ii,]
    }  
    
    result <- npqcd(x = xx, G, data.name = x$data.name)

  } else {
      cat("The process is under control")
  }
  
  
  oldClass(result) <- c("npqcd", "list")
  
  invisible(result)
  
} #npsate.control
#.........................................................................