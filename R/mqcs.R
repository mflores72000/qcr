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
# Main function to create a 'mqcs' object
#-----------------------------------------------------------------------------#
##'  It computes statistics to be used in Multivariante Quality Control
##' 
##' Create an object of class 'mqcs' to perform statistical quality control.
##' This function is used to compute statistics required to plot Multivariate Control Charts
##' 
##' @aliases mqcs summary.mqcs print.mqcs
##' 
##' @param x  Object mqcd (Multivariante Quality Control Data)
##' @param method Is the method employed to compute the covatiance matrix
##' in individual observation case. Two methods are used "sw" 
##' for compute according to (Sullivan,Woodall 1996a) and "hm" 
##' by (Holmes,Mergen 1993)
##' @param ... arguments passed to or from methods.
##' @export

mqcs <- function(x, method = "sw", ...)
  #.........................................................................  
  {

  if (!inherits(x, "mqcd"))
    stop("object must be mqcd")
  
  p <- ncol(x) # quality characteristics
  m <- nrow(x) # number of samples or observations
  n <- dim(x)[3] # observations or sample size 
  
  x.jk <- matrix(0,m,p)  
  x.jk <- apply(x,1:2,mean)
  
  Xmv <- colMeans(x.jk)
  S <- covariance(x,method = method)
  
  result <- list (mqcd = x, mean = Xmv, S = S, mean.jk = x.jk) 
  
  oldClass(result) <- c("mqcs")
  
  return(result)
} # mqcs
#.........................................................................

##' @export
##' @method print mqcs
print.mqcs <- function(x, ...) str(x,1)
#.........................................................................
##' @export
##' @method summary mqcs
summary.mqcs <- function(object, ...)
  #.........................................................................
{
  type <- object$type
  t2 <- object$statistics
  cat("\nSummary of group statistics:\n")
  print(summary(t2))
  cat("\nNumber of quality characteristics: ", ncol(object$mqcd))
  cat("\nNumber of samples or observations: ", nrow(object$mqcd))
  cat("\nNumber of observations or sample size: ", dim(object$mqcd)[3])
  
  center <- object$mean
  cat("\n\nMean Vector: \n", center)
  cat("\nCovariance Matrix:\n")
  S <-object$S
  print(S)
  
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
} # summary.mqcs



## ' Sample covariance
## ' 
## ' It allows to compute the sample covariance in presence of rational subgroups
## ' or for individuals according to (Sullivan,Woodall 1996) and (Holmes,Mergen
## ' 1993)
## ' 
## ' 
## ' @param x matrix or array of the quality characteristics.
## ' @param stat is the statistics
## ' @param method is the method used in individual observation case.
## ' @param \dots other parameters
## ' @note In individuals observation case (n = 1) use for default the
## ' (Sullivan,Woodall 1996) proposal
## ' @author Edgar Santos-Fernandez
## ' @references Holmes, D.S., Mergen, A.E.: Improving the performance of
## ' T-square control chart. Quality Engineering 5(4), 619-625 (1993)
## ' 
## ' Sullivan, J.H., Woodall, W.H.: A Comparison of Multivariate Quality Control
## ' Charts for Individual Observations. Journal of Quality Technology 28(4)
## ' (1996)
## ' @keywords ~kwd1 ~kwd2
## ' @examples
## ' 
## ' # individual case 
## ' data(dowel1)
## ' covariance(dowel1,method="sw")
## ' covariance(dowel1,method="hm")
 
covariance <- function(x, stat, method, ...){
    p <- ncol(x) # quality characteristics
    m <- nrow(x) # sample
    if (inherits(x, "matrix") || inherits(x, "data.frame"))
        x <- array(data.matrix(x),c(m,p,1))
    n <- dim(x)[3] # observations or sample size
    
    s.jk <- matrix(0,m,p ^ 2) # matrix of the covariances of observations
    SS <- matrix(0,m,1) # matrix of /S/ statistic 
    
    if(n > 1){
      arrays <- expand.grid(1:p,1:p)
      
      for (i in 1 : m){
        for(j in 1 : p ^ 2){
          s.jk[i,j] <- cov(x[i,arrays[j,1],],x[i,arrays[j,2],])
        }
      } 
      
      S <- matrix(colMeans(s.jk),p,p)
      
      for (ii in 1 : m){
        SS[ii] <- det(matrix(s.jk[ii,],p,p))
      }
      
      if(missing(stat)) (return(S))
      else (return(SS))
      
    }    
    
    if(n == 1){
      if(missing(method))(method="sw")
      
      if(method == "sw"){
        B <- matrix(0,p,p)
        w <- sweep(x,2,(apply(x,2,mean))) #compute de value minus the mean
        for(i in 1:m){
          B <- B + w[i,,] %*% t(w[i,,])
        }
        S <- s1 <- B/(m - 1)
      }
      
      if(method == "hm"){
        V <- matrix(0,m-1,p)
        for(i in 1:m-1){
          V[i,] <- x[i+1,,] - x[i,,]
        }
        S <- s2 <- .5 * t(V) %*% V / (m - 1)
      }
      
      
      return(S)
    }
    
    
  }


# mqcs.add function
#-------------------------------------------------------------------------
##' mqcs.add Add a matrix, data.frame or array object with a mqcs object
##' 
##' This function is used to join two objects of type matrix, data.frame or array and mqcs.
##' 
##' @param x   Object type mqcs
##' @export
##' 


mqcs.add <- function(x, ...){
  UseMethod("mqcs.add")
}

##' @rdname  mqcs.add 
##' @method mqcs.add default
##' @param value   Object type data.frame, matrix or array
##' @param ...  arguments to be passed to or from methods.
##' @export 


mqcs.add.default <- function(x, value, ...){
  
  if (!inherits(x, "mqcs"))
    stop("object must be mqcs")
  
  if (!is.matrix(value) & !is.data.frame(value) & !is.array(value))
    stop("object must be a matrix, data.frame or array")
  
  if (inherits(value, "matrix") || inherits(value, "data.frame")) {
    p <- ncol(value) # quality characteristics
    m <- nrow(value) # number of samples or observations
    names <- colnames(value)    
    value <- array(data.matrix(value),c(m,p,1))
    colnames(value) <- names        
  }   
  
  data <- x$mqcd
  limits <- x$limits
  data.name = x$data.name
  type <- x$type
  alpha <- x$alpha
  
  n1 <- dim(data)[3]
  n2 <- dim(value)[3]
  m1 <- dim(data)[1]
  m2 <- dim(value)[1]
  m <- m1 + m2
  k1 <- dim(data)[2]
  k2 <- dim(value)[2]
  
  if (n1 != n2) stop("The samples must be of the same dimension")
  if (k1 != k2) stop("The samples must be of the same quality characteristics")
  
  xx <- array(,dim = c(m,k1,n1))
  for (i in 1:n1 ){
    xx[,,i] <- rbind(data[,,i],value[,,i])     
  }
  
  z.mqcd <- mqcd(data =xx , data.name = data.name)
  
  mqcs.t2.mqcd(x = z.mqcd, limits = limits, alpha = alpha)
  
  z.mqcs <- switch(type, 
                    "t2" = mqcs.t2.mqcd(x = z.mqcd, limits = limits, alpha = alpha),
                    "mcusum" = mqcs.mcusum.mqcd(x = z.mqcd, limits = limits, alpha = alpha),
                    "mewma" = mqcs.mewma.mqcd(x = z.mqcd, limits = limits, alpha = alpha),
                    NULL)
  result <- z.mqcs
}


#-------------------------------------------------------------------------
# mstate.control
#-------------------------------------------------------------------------
##' Multivariate process state
##' 
##' This function removes observations from the sample which violates 
##' the rules of a process under control
##' @aliases mstate.control
##' @param x  Object mqcd (Multivariate Quality Control Statistical)
##' @param control a logical value indicating whether the initial sample comes from a process under control.
##' @export
##' @examples
##' 
##' ##
##' ##  Continuous data 
##' ##
##' library(qcr)
##' set.seed(356)
##' x <- matrix(rnorm(66),ncol=3)
##' x <- rbind(x,matrix(rexp(66,100),ncol=3))
##' dim(x)
##' x <-mqcd(x)
##' str(x)
##' x <-mqcs.mewma(x)
##' str(x)
##' plot(x)
##' data.mqcs <- mstate.control(x)
##' x <-mqcs.mewma(data.mqcs)
##' plot(x)

mstate.control <- function(x) 
  #.........................................................................  
{
  if (!inherits(x, "mqcs"))
    stop("an object of class 'mqcs' is required")
  
  if (length(x$violations)>0){
    ii<-x$violations  
    n <- dim(x$mqcd)[3]
    m <- dim(x$mqcd)[1]
    k <- dim(x$mqcd)[2]
    xx <- array(dim = c(m-length(ii),k,n))
    for (i in 1:n){ 
      xx[,,i] <- x$mqcd[,,i][-ii,]  
    }

    result <- mqcd(data = xx, data.name = x$data.name)
    
  } else {
    cat("The process is under control")
  }
  
  oldClass(result) <- c("mqcd","array")
  
  invisible(result)
  
} #msate.control
#.........................................................................