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
#-----------------------------------------------------------------------------#
# Main function to create a 'qcs.pcr' object
#-----------------------------------------------------------------------------#
##' Process capability indices for a given dataset and distribution
##' 
##' Calculates the process capability indices cp, cpk, cpkL and cpkU for a given 
##' dataset and distribution. 
##' A histogram with a density curve is displayed along with the specification limits and a 
##' Quantile-Quantile Plot for the specified distribution.
##' @aliases qcs.pcr 
##' @param object qcs object of type \code{"qcs.xbar"} or \code{"qcs.one"}.
##' @param distribution Character string that represent the probability 
##' distribution of the data, such as: "normal", "beta", "chi-squared", 
##' "exponential", "f", "geometric", "lognormal", "log-normal", "logistic","t",
##' "negative binomial", "poisson", "weibull", "gamma".
##' @param limits A vector specifying the lower and upper specification limits.
##' @param target A value specifying the target of the process. 
##' If it is \code{NULL}, the target is set at the middle value between specification limits.
##' @param std.dev A value specifying the within-group standard deviation.
##' @param boxcox Logical value (by default \code{FALSE}). If \code{TRUE}, 
##' perform a Box-Cox transformation.
##' @param lambda A vector specifying or numeric value indicating lambda for the transformation.
##' @param confidence A numeric value between 0 and 1 specifying the nivel for 
##' computing the specification limits.
##' @param plot Logical value indicating whether graph should be plotted.
##' @param main Title of the plot.
##' @param ... Arguments to be passed to or from methods.
##' @export
##' @references 
##' Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd
##' ed, New York, John Wiley & Sons. \cr
##' @examples
##' library(qcr)
##' data(pistonrings) 
##' xbar <- qcs.xbar(pistonrings[1:125,],plot = TRUE)
##' limits = c(lsl = 73.99, usl = 74.01)
##' qcs.pcr(xbar, "normal", limits = limits) 
##' qcs.pcr(xbar, "weibull", limits = limits)
qcs.pcr <- function (object, distribution = c("normal","beta", "chi-squared", 
                                              "exponential", "f", "geometric", "lognormal", 
                                              "log-normal", "logistic","t",
                                              "negative binomial", "poisson", "weibull", "gamma"), 
                     limits = c(lsl = -3, usl = 3), 
                     target = NULL, std.dev = NULL, boxcox = FALSE, lambda = c(-5, 5), 
                     confidence = 0.9973, plot = TRUE, main = NULL, ...) 
{
  
  
  if (missing(object)){
    if (!inherits(object, "qcs"))
      stop("an object of class 'qcs' is required")
    if (!(object$type == "xbar" | object$type == "one")) 
      stop("Process Capability Analysis only available for charts type 
           \"qcs.xbar\" and \"qcs.one\" charts")
  }
  
  dis <- match.arg(distribution)
  x<-object$statistics
  numObs = sum(object$sizes)
  center = object$center
  std.dev = object$std.dev
  data.name = object$data.name

  parList = list(...)

  if (length(limits)!=2)
    stop("specification limits must be two")
  lsl <- limits[1]
  usl <- limits[2]
  if (lsl>= usl)
    stop("lsl >= usl")
  if (!(is.numeric(usl) & is.finite(lsl)))
    lsl <- NA
  if (!(is.numeric(usl) & is.finite(lsl)))
    usl <- NA
  if (is.na(lsl) & is.na(usl))
    stop("invalid specification limits")
  
  if (is.null(target)) target <- mean(limits, na.rm = TRUE)
  
  if (is.na(lsl)) {
    if (target > usl)
      warning("target value larger than one-sided specification limit...")
  }
  if (is.na(usl)) {
    if (target < lsl)
      warning("target value smaller than one-sided specification limit...")
  }
  if (!is.na(lsl) & !is.na(usl)) {
    if (target < lsl || target > usl)
      warning("target value is not within specification limits...")
  }
  m <- (lsl+usl)/2
  
  if (!is.numeric(lambda)) 
    stop("lambda needs to be numeric")
  
  if (any(x < 0) && any(dis == c("beta", "chi-squared", 
                                 "exponential", "f", "geometric", "lognormal", "log-normal", 
                                 "negative binomial", "poisson", "weibull", "gamma"))) 
    stop("choosen dis needs all values in x to be > 0")
  
  if (any(x > 1) && dis == "beta") 
    stop("choosen dis needs all values in x to be between 0 and 1")
  
  if (confidence < 0 | confidence > 1) 
    stop("conf.level must be a value between 0 and 1")
  
  confHigh = confidence + (1 - confidence)/2
  confLow = 1 - confidence - (1 - confidence)/2
  
  paramsList = vector(mode = "list", length = 0)
  estimates = vector(mode = "list", length = 0)
  varName = data.name
  
  if (boxcox) {
    dis = "normal"
    if (length(lambda) >= 2) {
      temp = boxcox(x[, 1] ~ 1, lambda = seq(min(lambda), 
                                             max(lambda), 1/10), plotit = FALSE)
      i = order(temp$y, decreasing = TRUE)[1]
      lambda = temp$x[i]
    }
    x = as.data.frame(x[, 1]^lambda)
  }  
  
  distWhichNeedParameters = c("weibull", "logistic", "gamma", 
                              "exponential", "f", "geometric", "chi-squared", "negative binomial", 
                              "poisson")
  
  
  qFun = .charToDistFunc(dis, type = "q")
  pFun = .charToDistFunc(dis, type = "p")
  dFun = .charToDistFunc(dis, type = "d")
  
  if (is.null(qFun) & is.null(pFun) & is.null(dFun)) 
    stop(paste(deparse(substitute(y)), "dis could not be found!"))
  
  fitList = vector(mode = "list", length = 0)
  fitList$x = x[, 1]
  fitList$densfun = dis
  
 
  fittedDistr = do.call(fitdistr, fitList)
  estimates = as.list(fittedDistr$estimate)
  paramsList = estimates
  
  paramsList = c(paramsList, .lfkp(parList, formals(qFun)))
  
  if (dis == "normal") {
    paramsList$mean = center
    paramsList$sd = std.dev
    estimates = paramsList
  }
  
  if (boxcox) {
    if (!is.null(lsl)) 
      lsl = lsl^lambda
    if (!is.null(usl)) 
      usl = usl^lambda
    if (!is.null(target)) 
      target = target^lambda
  }
  
  if (is.null(lsl) && is.null(usl)) {
    paramsList$p = confLow
    lsl = do.call(qFun, paramsList)
    paramsList$p = confHigh
    usl = do.call(qFun, paramsList)
  }
  
  
  if (!is.null(lsl) && !is.null(target) && target < lsl) 
    stop("target is less than lower specification limit")
  if (!is.null(usl) && !is.null(target) && target > usl) 
    stop("target is greater than upper specification limit")
  
  paramsList$p = c(confLow, 0.5, confHigh)
  paramsListTemp = .lfkp(paramsList, formals(qFun))
  qs = do.call(qFun, paramsListTemp)
  paramsListTemp = .lfkp(paramsList, formals(pFun))
  
  if (!is.null(lsl) && !is.null(usl)) 
    cp = (usl - lsl)/(qs[3] - qs[1])
  if (!is.null(usl)) {
    cpu = (usl - qs[2])/(qs[3] - qs[2])
    paramsListTemp$q = usl
    ppu = 1 - do.call(pFun, paramsListTemp)
  }
  
  if (!is.null(lsl)) {
    cpl = (qs[2] - lsl)/(qs[2] - qs[1])
    paramsListTemp$q = lsl
    ppl = do.call(pFun, paramsListTemp)
  }
  
  cpk = min(cpu, cpl)
  ppt = sum(ppl, ppu)
  
  if (plot == TRUE) {
    
    par.orig <- par(c("mar", "oma", "mfrow"))
    
    
    if (is.null(parList[["col"]])) 
      parList$col = "lightblue"
    if (is.null(parList[["border"]])) 
      parList$border = 1
    if (is.null(parList[["lwd"]])) 
      parList$lwd = 1
    if (is.null(parList[["cex.axis"]])) 
      parList$cex.axis = 1.5
    
    lineWidth = 1
    lineCol = "red" 
    lineType = "solid"
    specCol = "red3"
    specWidth = 1
    cex.text = 1.5 
    cex.val = 1
    cex.col = "darkgray"
    bounds.lty = 3 
    bounds.col = "red"
    xlim <- range(x[, 1], usl, lsl)
    xlim <- xlim + diff(xlim) * c(-0.2, 0.2)
    
    
    xVec <- seq(min(xlim), max(xlim), length = 200)
    dParamsList = .lfkp(paramsList, formals(dFun))
    dParamsList$x = xVec
    yVec = do.call(dFun, dParamsList)
    histObj <- hist(x[, 1], plot = FALSE)
    
    ylim <- range(histObj$density, yVec)
    ylim <- ylim + diff(ylim) * c(0, 0.05)
    
    par(bg="#CCCCCC",mar = c(0, 0, 0, 0) + 0.1)
    
    par(oma = c(2, 4, 7, 4) + 0.1)
    layout(matrix(c(1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 1, 3, 1, 
                    1, 1, 3), nrow = 4, byrow = TRUE))
    
    
    do.call(hist, c(list(x[, 1], freq = FALSE, xlim = xlim, 
                         ylim = ylim, main = "",col = "#EEEEEE")))
    
    tempList = parList
    tempList$col = "#EEEEEE"
    tempList$border = NULL
    
    abline(v = 0)
    abline(h = 0)
    lines(xVec, yVec, lwd = lineWidth, col = lineCol, lty = lineType)
    abline(v = usl, col = specCol, lwd = specWidth, lty = 5)
    abline(v = lsl, col = specCol, lwd = specWidth, lty = 5)
    abline(v = target, col = specCol, lwd = specWidth, lty = 5)
    
    do.call(grid,list(col = "#EEEEEE"))
    do.call(box,tempList)
    
    if (!is.null(lsl)) 
      axis(side = 3, at = lsl, labels = paste("LSL =", 
                                              format(lsl, digits = 3)), col = specCol)
    if (!is.null(usl)) 
      axis(side = 3, at = usl, labels = paste("USL =", 
                                              format(usl, digits = 3)), col = specCol)
    if (!is.null(lsl) && !is.null(usl)) 
      axis(side = 3, at = c(lsl, usl), labels = c(paste("LSL =", 
                                                        format(lsl, digits = 3)), paste("USL =", format(usl, 
                                                                                                        digits = 3))), col = specCol)
    if (!is.null(target)) 
      text(target, max(ylim), "TARGET", pos = 1, col = "black", 
           cex = cex.text)
    
    if (is.null(main)){
      if (boxcox) 
        main = paste("Process Capability using box cox transformation for", 
                     varName)
      else main = paste("Process Capability using", as.character(dis), 
                        "distribution for", varName)
    }
    
    par(mar = c(0, 0, 0, 0) + 0.1)
    index = 1:(length(estimates) + 3)
    title(main = main, outer = TRUE)
    plot(0:5, 0:5, type = "n", axes = FALSE, xlab = "", ylab = "", 
         main = "")
    
    names(x) = data.name
    adTestStats = .myADTest(x, dis)
    A = adTestStats$statistic
    p = adTestStats$p.value
    
    text(2.3, rev(index)[1], paste("n =",numObs), pos = 2, 
         cex = cex.val)
    
    text(2.3, rev(index)[2], paste("A =",format(A, digits = 3)), 
         pos = 2, cex = cex.val)
    
    
    if (!is.null(adTestStats$smaller) && adTestStats$smaller) 
      text(2.3, rev(index)[3], paste("p <", format(p, digits = 3)), 
           pos = 2, cex = cex.val)
    if (!is.null(adTestStats$smaller) && !adTestStats$smaller) 
      text(2.3, rev(index)[3], paste("p >=", format(p, 
                                                    digits = 3)), pos = 2, cex = cex.val)
    if (is.null(adTestStats$smaller)) 
      text(2.3, rev(index)[3], paste("p =", format(p, digits = 3)), 
           pos = 2, cex = cex.val)
    
    
    
    j = 1
    for (i in 3:(3 + length(estimates) - 1)) {
      try(text(2.3, rev(index)[i + 1], paste(names(estimates)[[j]],
                                             "=", 
                                             format(estimates[[j]], 
                                                    digits = 3)), 
               pos = 2, cex = cex.val), silent = TRUE)
      j = j + 1
    }
    
    do.call(box,tempList)
    #do.call(grid,list(col = "#EEEEEE"))
    
    
    do.call(qqPlot,list(x[, 1], y = dis, ylab = "", main = "", 
           axes = FALSE, bounds.lty = bounds.lty, 
           bounds.col = bounds.col))
    do.call(box,tempList)
    do.call(grid,list(col = "#EEEEEE"))
    on.exit(par(par.orig))
  }  
  
  if (!is.null(lsl)) {
    obsL = (sum(x < lsl)/length(x)) * 1e+06
  }
  else obsL = 0
  if (!is.null(usl)) {
    obsU = (sum(x > usl)/length(x)) * 1e+06
  }
  else obsU = 0
  
  obsT = obsL + obsU
  
  digits =4
  tab <- cbind(c(cp, cpl, cpu, cpk)) 
  
  rownames(tab) <- c("Cp", "Cp_l", "Cp_u", "Cp_k")
  
  colnames(tab) <- c("Value")
  
  cat("\nProcess Capability Analysis\n")
  cat("\nCall:\n", deparse(match.call()), "\n\n", sep = "")
  cat(paste(formatC("Number of obs = ", width = 16), formatC(numObs, 
                                                             width = 12, flag = "-"), formatC("Target = ", width = 10), 
            formatC(signif(target, digits = digits), 
                    flag = "-"), "\n", sep = ""))
  cat(paste(formatC("Center = ", width = 16), formatC(signif(center,digits = digits), width = 12, flag = "-"),
            formatC("LSL = ", width = 10),formatC(signif(lsl,digits = digits), flag = "-")), "\n", sep = "")
  
  cat(paste(formatC("StdDev = ", width = 16), formatC(signif(std.dev,digits = digits), width = 12, flag = "-"), 
            formatC("USL = ", width = 10), formatC(signif(usl, digits = digits), flag = "-")), "\n", sep = "")
  
  cat("\nCapability indices:\n\n")
  print(tab, digits = 4, na.print = "", print.gap = 2)
  cat("\n")
  cat("\nPPM:\n\n")
  cat(paste(formatC("Exp<LSL", width = 16), formatC(ppl* 1e+06, width = 12, flag = "-"), 
            formatC("Obs<LSL", width = 10), formatC(obsL* 1e+06, flag = "-")), "\n", sep = "")
  
  cat(paste(formatC("Exp>USL", width = 16), formatC(ppu* 1e+06, width = 12, flag = "-"), 
            formatC("Obs>USL", width = 10), formatC(obsU* 1e+06, flag = "-")), "\n", sep = "")
  
  cat(paste(formatC("Exp Total", width = 16), formatC(ppt* 1e+06, width = 12, flag = "-"), 
            formatC("Obs Total", width = 10), formatC(obsT* 1e+06, flag = "-")), "\n", sep = "")
  
  cat("\nTest:\n\n")
  print(adTestStats)
  
  
  
  result <- list(lambda = lambda, cp = cp, cpk = cpk, cpl = cpl, 
                 cpu = cpu, ppt = ppt, ppl = ppl, ppu = ppu, usl = usl, 
                 lsl = lsl, target = target,obsU = obsU, obsL = obsL, obsT = obsT)
  
  invisible(result)
}