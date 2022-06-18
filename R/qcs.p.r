#-----------------------------------------------------------------------------#
#                                                                             #
#                  QUALITY CONTROL STATISTICS IN R                            #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sanchez                                       #
#              Professor of Mathematic Department                             #
#              Escuela Politecnica Nacional, Ecuador                          #
#              miguel.flores@epn.edu.ec                                       #
#                                                                             #
#-----------------------------------------------------------------------------#
#-------------------------------------------------------------------------
# p chart
#-------------------------------------------------------------------------
##' Function to plot Shewhart p chart
##'
##' This function is used to compute statistics required by the p chart.
##'
##' @param x   An R object (used to select the method). See details.
##' @param ... Arguments passed to or from methods.
##' @export
##' @examples
##' library(qcr)
##' data(orangejuice)
##' str(orangejuice)
##' attach(orangejuice)
##'
##' datos.qcd <- qcd(data = orangejuice, var.index = 1, sample.index = 2,
##'                 sizes = size, type.data = "atributte")
##'
##' res.qcs <- qcs.p(datos.qcd)
##' summary(res.qcs)
##' plot(res.qcs)
##'
##' datos.qcs <- qcs.p(orangejuice[trial,c(1,2)], sizes = orangejuice[trial,3])
##' plot(datos.qcs)

qcs.p <- function(x, ...) {
  UseMethod("qcs.p")
}

##' @rdname qcs.p
##' @method qcs.p default
##' @inheritParams qcd
##' @param center A value specifying the center of group statistics or the
##' ''target'' value of the process.
##' @param conf.nsigma  A numeric value used to compute control limits, specifying the
##' number of standard deviations (if \code{conf.nsigma} > 1) or the confidence level (if 0
##' < \code{conf.nsigma} < 1).
##' @param limits A two-values vector specifying control limits.
##' @param plot Logical value. If TRUE a p chart should be plotted.
##' @details
##' In the default method \code{qcs.p.default} parameter \code{x} is a matrix
##' or data-frame where it should contain data, index sample and, optionally, covariate(s).
##' @export
##' 
qcs.p.default <- function(x, var.index  =  1, sample.index  =  2,
                             covar.index  =  NULL, covar.names  =  NULL,
                             data.name = NULL,
                             sizes  =  NULL,
                             center = NULL,
                          conf.nsigma  =  3, limits = NULL, plot = FALSE, ...)
  {
  if (is.null(sizes)) 
    stop("sample sizes must be given for a attribute variable")


  obj<-qcd(data = x, var.index = var.index, sample.index = sample.index,
       covar.index = covar.index, covar.names = covar.names,
       data.name = data.name, sizes = sizes, type.data = "atributte")

  result<-qcs.p.qcd(x = obj, center = center, 
                    conf.nsigma = conf.nsigma, 
                    limits = limits, plot = plot)

  return(result)
}


##' @rdname  qcs.p
##' @method qcs.p qcd
##' @inheritParams qcs.p.default
##' @export
##' 
qcs.p.qcd <- function(x, center = NULL,
                         conf.nsigma  =  3, limits = NULL, plot = FALSE, ...) {
  #.........................................................................
  if(is.null(x) || !inherits(x, "qcd"))
    stop("data must be an objects of class (or extending) 'qcd'")
  sizes <- x$sizes
  type.data <- "atributte"
  
  qcs<-qcs(x = x$x, sample.index = x$sample, sizes = sizes, type  =  "p",
            center = center, 
           conf.nsigma = conf.nsigma, limits = limits, type.data = type.data)
  
  center <- qcs$center
  p <- qcs$statistics
  std.dev <- qcs$std.dev
  sizes <- qcs$sizes
  limits <- qcs$limits
  violations <- qcs$violations
  
  statistics <- data.frame(p)
  m <- length(x)
  sample <- x$sample

  if (m > 3) {
    new.x <- x[, -c(1, 2, length(x))]
    cov <- apply(new.x, 2, function(x) unlist(lapply(split(x, sample), unique)))
    statistics <- data.frame(p, cov)
  }
  
  row.names(statistics) <- unique(x$sample)
  data.name <- attr(x, "data.name")
  result <- list(qcd  =  x, type  =  "p", statistics  =  statistics,
                 center  =  center, std.dev  =  std.dev,
                 limits  =  limits, conf.nsigma  =  conf.nsigma,
                 sizes  =  sizes, data.name  =  data.name,
                 violations  =  violations)
  
  oldClass(result) <- c("qcs.p", "qcs")
  
  if(plot) plot(result, ...)
  
  return(result)
  #.........................................................................
} # qcs.p.qcd