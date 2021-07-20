# -------------------------------------
# This code is part of the qualityTools 1.55 package
# removed from the CRAN repository on 2021-02-23 as check problems were not corrected.
# Copyright Thomas Roth 2010-2016
# The complete package can be obtained from the archive: 
# https://cran.r-project.org/src/contrib/Archive/qualityTools/
# -------------------------------------
# Package: qualityTools
# Type: Package
# Title: Statistical Methods for Quality Science
# Version: 1.55
# Date: 2016-02-23
# Author: Thomas Roth
# Maintainer: Thomas Roth <thomas.roth@alumni.tu-berlin.de>
# Description: Contains methods associated with the Define, Measure, Analyze, 
# Improve and Control (i.e. DMAIC) cycle of the Six Sigma Quality Management 
# methodology.It covers distribution fitting, normal and non-normal process 
# capability indices, techniques for Measurement Systems Analysis especially 
# gage capability indices and Gage Repeatability (i.e Gage RR) and 
# Reproducibility studies, factorial and fractional factorial designs as well 
# as response surface methods including the use of desirability functions. 
# Improvement via Six Sigma is project based strategy that covers 5 phases: 
# Define - Pareto Chart; Measure - Probability and Quantile-Quantile Plots, 
# Process Capability Indices for various distributions and Gage RR Analyze i.e. 
# Pareto Chart, Multi-Vari Chart, Dot Plot; Improve - Full and fractional 
# factorial, response surface and mixture designs as well as the desirability 
# approach for simultaneous optimization of more than one response variable. 
# Normal, Pareto and Lenth Plot of effects as well as Interaction Plots; 
# Control - Quality Control Charts can be found in the 'qcc' package. 
# The focus is on teaching the statistical methodology used in the Quality 
# Sciences.
# License: GPL-2
# -------------------------------------

# PENDIENTE: documentacion -> internal

# @rdname qcr-internals
#' @keywords internal
qqPlot <- function (x, y, confbounds = TRUE, alpha, main, xlab, ylab, xlim, ylim, border = "red", bounds.col = "black", bounds.lty = 1,
                    start, ...) 
{
  DB = FALSE
  parList = list(...)
  if (is.null(parList[["col"]])) 
    parList$col = 1:2
  if (is.null(parList[["pch"]])) 
    parList$pch = 19
  if (is.null(parList[["lwd"]])) 
    parList$lwd = 1
  if (is.null(parList[["cex"]])) 
    parList$cex = 1
  
  #if (!require(MASS)) 
  #    stop("Package MASS needs to be installed!")
  
  if (class(x) == "distrCollection") {
    distList = x@distr
    for (i in 1:length(distList)) {
      d = distList[[i]]
      do.call(qqPlot, c(list(x = d@x, y = d@name), parList))
    }
    invisible()
  }
  if (missing(y)) 
    y = "normal"
  if(missing(alpha))
    alpha = 0.05
  if (alpha <=0 || alpha >=1) 
    stop(paste("alpha should be between 0 and 1!"))		
  if (missing(main)) 
    main = paste("Q-Q Plot for", deparse(substitute(y)), 
                 "distribution")
  if (missing(xlab)) 
    xlab = paste("Quantiles for", deparse(substitute(x)))
  if (missing(ylab)) 
    ylab = paste("Quantiles from", deparse(substitute(y)), 
                 "distribution")
  if (is.numeric(y)) {
    cat("\ncalling (original) qqplot from namespace stats!\n")
    return(stats::qqplot(x, y, ...))
  }
  qFun = NULL
  theoretical.quantiles = NULL
  xs = sort(x)
  distribution = tolower(y)
  distWhichNeedParameters = c("weibull", "logistic", "gamma", 
                              "exponential", "f", "geometric", "chi-squared", "negative binomial", 
                              "poisson")
  
  
  # new
  threeParameterDistr = c("weibull3", "lognormal3", "gamma3")                
  threeParameter = distribution %in% threeParameterDistr
  if(threeParameter) distribution = substr(distribution, 1, nchar(distribution)-1)
  # end new
  
  if (is.character(distribution)) {
    qFun = .charToDistFunc(distribution, type = "q")
    if (is.null(qFun)) 
      stop(paste(deparse(substitute(y)), "distribution could not be found!"))
  }
  theoretical.probs = ppoints(xs)
  
  xq = NULL
  yq = quantile(xs, prob = c(0.25, 0.75))
  dots <- list(...)
  if (TRUE) {
    if (DB) 
      print("TODO: Pass the estimated parameters correctly")
    fitList = .lfkp(parList, formals(qFun))
    fitList$x = xs
    fitList$densfun = distribution
    if (!missing(start)) 
      fitList$start = start
    if (DB) {
      print(fitList)
      print("Ende")
    }
    # new
    if(!threeParameter){
      fittedDistr = do.call(fitdistr, fitList)
      parameter = fittedDistr$estimate
      
      #save the distribution parameter#
      thethas = fittedDistr$estimate
      # save the cariance-covariance matrix
      varmatrix = fittedDistr$vcov
      # end of my code
      
      # new code for three parameter
    } else {
      parameter = do.call(paste(".",distribution, "3", sep = ""), list(xs) )    ####
      threshold = parameter$threshold
    }
    
    parameter = .lfkp(as.list(parameter), formals(qFun))
    params = .lfkp(parList, formals(qFun))
    parameter = .lfrm(as.list(parameter), params)
    parameter = c(parameter, params)
    theoretical.quantiles = do.call(qFun, c(list(c(theoretical.probs)), 
                                            parameter))
    
    # new
    if(!threeParameter){		
      # array containing names of the distributions, for which conf intervals can be computed
      confIntCapable = c("exponential", "log-normal", "logistic", "normal", "weibull", "gamma", "beta", "cauchy")
      getConfIntFun = .charToDistFunc(distribution, type = ".confint")
      # if possible, compute the conf intervals
      if(confbounds == TRUE){
        if(distribution %in% confIntCapable){
          confInt = getConfIntFun(xs, thethas, varmatrix, alpha)
        }
      }# end of my code
    }
    
    xq <- do.call(qFun, c(list(c(0.25, 0.75)), parameter))
    if (DB) {
      print(paste("parameter: ", parameter))
      print(xq)
    }
  }
  else {
    params = .lfkp(parList, formals(qFun))
    params$p = theoretical.probs
    theoretical.quantiles = do.call(qFun, params)
    params$p = c(0.25, 0.75)
    xq = do.call(qFun, params)
  }
  
  params = .lfkp(parList, c(formals(plot.default), par()))	
  
  if(!threeParameter){
    params$y = theoretical.quantiles
  }  else {
    params$y = theoretical.quantiles+threshold
  }
  params$x = xs
  params$xlab = xlab
  params$ylab = ylab
  params$main = main
  if (!(is.null(params$col[1]) || is.na(params$col[1]))) 
    params$col = params$col[1]
  if (!missing(xlim)) 
    params$xlim = xlim
  if (!missing(ylim)) 
    params$ylim = ylim
  params$lwd = 1
  do.call(plot, params)
  pParams = params
  pParams = .lfkp(pParams, list(x = 1, y = 1, col = 1, cex = 1))
  do.call(points, pParams)
  params = .lfkp(parList, c(formals(abline), par()))
  params$a = 0
  params$b = 1
  params$col = border
  do.call(abline, params)
  
  if(!threeParameter){
    # plot the confInt if available
    if(confbounds == TRUE){
      if(distribution %in% confIntCapable){
        params = .lfkp(parList, c(formals(lines), par()))	
        params$x = confInt[[3]]
        params$y = confInt[[1]]
        params$col = bounds.col
        params$lty = bounds.lty
        do.call(lines, params)
        
        params$x = confInt[[3]]
        params$y = confInt[[2]]
        params$col = bounds.col
        params$lty = bounds.lty
        do.call(lines, params)
      }
    } #end of my function
  }
  
  invisible(list(x = theoretical.quantiles, y = xs, int = params$a, 
                 slope = params$b))
}

# @rdname qcr-internals
#' @keywords internal
.confintbeta= function(xs, thethas, varmatrix, alpha) {

    th1= thethas[[1]]
    th2 =  thethas[[2]]
    
    prozent=ppoints(xs)

    perzentile=qbeta(prozent, th1, th2)
    
    h=1e-6
    dFdth1=(qbeta(prozent, th1, th2)-qbeta(prozent, th1+h, th2))/h
    dFdth2=(qbeta(prozent, th1, th2)-qbeta(prozent, th1, th2+h))/h
    
    Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
    zalpha=qnorm(1-alpha/2)
    halfwidth = zalpha*sqrt(Var)

    lci=perzentile-halfwidth
    uci=perzentile+halfwidth
	
	bounds = list(lci, uci, perzentile)


return (bounds) 
}

# @rdname qcr-internals
#' @keywords internal
.confintcauchy = function(xs, thethas, varmatrix, alpha) {

    th1= thethas[[1]]
    th2 =  thethas[[2]]
    
    prozent=ppoints(xs)

    perzentile=qcauchy(prozent, th1, th2)
    
    h=1e-6
    dFdth1=(qcauchy(prozent, th1, th2)-qcauchy(prozent, th1+h, th2))/h
    dFdth2=(qcauchy(prozent, th1, th2)-qcauchy(prozent, th1, th2+h))/h
    
    Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
    zalpha=qnorm(1-alpha/2)
    halfwidth = zalpha*sqrt(Var)

    lci=perzentile-halfwidth
    uci=perzentile+halfwidth
	
	bounds = list(lci, uci, perzentile)

return (bounds) 
}

# @rdname qcr-internals
#' @keywords internal
.confintexp=function(xs, thethas, varmatrix, alpha) {
  lambda=thethas[[1]]
  prozent=ppoints(xs)
  perzentile=qexp(prozent, lambda)
  logPerzentile = log(perzentile)
   zalpha=qnorm(1-alpha/2)
  halfwidth = zalpha*sqrt(varmatrix[[1, 1]]/lambda^2)
   lci = exp(logPerzentile - halfwidth);
   uci = exp(logPerzentile + halfwidth);

   bounds = list(lci, uci, perzentile)

return (bounds)

}

# @rdname qcr-internals
#' @keywords internal
.confintgamma= function(xs, thethas, varmatrix, alpha) {
  th1= thethas[[1]]
  th2 =  thethas[[2]]
  prozent=ppoints(xs)
  perzentile=qgamma(prozent, th1, th2)
  
    h=1e-6
    dFdth1=(qgamma(prozent, th1, th2)-qgamma(prozent, th1+h, th2))/h
    dFdth2=(qgamma(prozent, th1, th2)-qgamma(prozent, th1, th2+h))/h
    
    Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
    zalpha=qnorm(1-alpha/2)
    halfwidth = zalpha*sqrt(Var)

    lci=perzentile-halfwidth
    uci=perzentile+halfwidth
  
  bounds = list(lci, uci, perzentile)

return (bounds)
}

# @rdname qcr-internals
#' @keywords internal
.confintlnorm=function(xs, thethas, varmatrix, alpha){
    th1= thethas[[1]]
    th2 =  thethas[[2]]

prozent=ppoints(xs)

perzentile=qlnorm(prozent, th1, th2)

zp=qnorm(prozent)

varPerzentile = varmatrix[[1, 1]]+2*varmatrix[[1, 2]]*zp+varmatrix[[2, 2]]*zp*zp
 
   zalpha=qnorm(1-alpha/2)
lci=log(perzentile)-zalpha*sqrt(varPerzentile)
uci=log(perzentile)+zalpha*sqrt(varPerzentile)

bounds = list(exp(lci), exp(uci), perzentile)

return (bounds) 

}

# @rdname qcr-internals
#' @keywords internal
.confintlogis= function(xs, thethas, varmatrix, alpha) {

    th1= thethas[[1]]
    th2 =  thethas[[2]]
    
    prozent=ppoints(xs)


    perzentile=qlogis(prozent, th1, th2)
    
    h=1e-6
    dFdth1=(qlogis(prozent, th1, th2)-qlogis(prozent, th1+h, th2))/h
    dFdth2=(qlogis(prozent, th1, th2)-qlogis(prozent, th1, th2+h))/h
    
    Var = varmatrix[[1, 1]]*dFdth1^2 + 2*varmatrix[[1, 2]]*dFdth1*dFdth2 + varmatrix[[2, 2]]*dFdth2^2
    zalpha=qnorm(1-alpha/2)
    halfwidth = zalpha*sqrt(Var)

    lci=perzentile-halfwidth
    uci=perzentile+halfwidth
	
	bounds = list(lci, uci, perzentile)

return (bounds) 
}

# @rdname qcr-internals
#' @keywords internal
.confintnorm=function(xs, thethas, varmatrix, alpha){

   prozent=ppoints(xs)

zp=qnorm(prozent)
perzentile=qnorm(prozent, thethas[[1]], thethas[[2]])

varPerzentile = varmatrix[[1, 1]]+2*varmatrix[[1, 2]]*zp+varmatrix[[2, 2]]*zp*zp
   
zalpha=qnorm(1-alpha/2)
lci=perzentile-zalpha*sqrt(varPerzentile)
uci=perzentile+zalpha*sqrt(varPerzentile)

bounds = list(lci, uci, perzentile)
 
return (bounds) 
}

# @rdname qcr-internals
#' @keywords internal
.confintweibull= function(xs, thethas, varmatrix, alpha) {
    th1= thethas[[1]]
    th2 =  thethas[[2]]
     
    prozent=ppoints(xs)

perzentile=qweibull(prozent, th1, th2)
q=-log(1-prozent)
logPerzentile=log(perzentile)
logq=log(q)
dB=1/th2
dA=-1/(th1^2)

Var = varmatrix[[1, 1]]*(dA*logq)^2 + 2*varmatrix[[1, 2]]*dB*dA*logq + varmatrix[[2, 2]]*dB^2
zalpha=qnorm(1-alpha/2)
halfwidth = zalpha*sqrt(Var)


lci=exp(logPerzentile-halfwidth)
uci=exp(logPerzentile+halfwidth)

bounds = list(lci, uci, perzentile)

# print(data.frame(prozent, uci, perzentile, lci))

return (bounds)
}

# @rdname qcr-internals
#' @keywords internal
.gamma3 = function(data) {
     n=length(data)
     data=sort(data)

     pEmp= (seq(1:n)-0.5)/n

     weight = 1 / sqrt(pEmp*(1-pEmp))

     thld = .99*min(data)
     shape=1
     scale=1

   gammaEst = function(param) {
    return( sum(weight*(pgamma(data-param[3], shape = exp(param[1]), scale = exp(param[2]))-pEmp)^2) )
  }

     paramEst = optim(c(shape, scale, thld), gammaEst, method = "Nelder-Mead")
     paramEst = paramEst$par
     return(list(shape = exp(paramEst[1]), scale = exp(paramEst[2]), threshold = paramEst[3]))
}

# @rdname qcr-internals
#' @keywords internal
.lognormal3 = function(data) {

  n=length(data)
  data=sort(data)
  #compute the empirical cumulative distribution function of the data
  pEmp= (seq(1:n)-0.5)/n
   # will minimize the weighted sum of squared distances
   # so compute weights
  weight = 1 / sqrt(pEmp*(1-pEmp))

  # initial values for optimization
  thld = .99*min(data)
  mu0 = mean(log(data-thld))
  sigma0 = sd(log(data-thld))


  lnEst = function(param) {
    return( sum(weight*(plnorm(data-param[3], meanlog = param[1], sdlog = exp(param[2]))-pEmp)^2) )
  }

  logSigma0=log(sigma0)
  # optimize gammaEst using optim function
  paramEst = optim(c(mu0,logSigma0, thld), lnEst, method = "Nelder-Mead")
  param = paramEst$par

  return(list(meanlog = param[1], sdlog = exp(param[2]), threshold = param[3]))
}

# @rdname qcr-internals
#' @keywords internal
.weibull3 = function(x)
{
  if(any(x < 0))
    stop("x must be positive")

  n = length(x)
  x = sort(x)
  p = ((1:n)-0.5)/n
  interval = c(0.75*min(x), 0.9999*min(x))

  wb3RSquared = function(th)
  {
      return(summary(lm(log(x-th) ~ log(-log(1-p))))$r.squared)
  }

  th = (optimize(wb3RSquared, interval = interval, maximum = TRUE))$maximum

  lm.1 = lm(log(x-th) ~ log(-log(1-p)))
  estimates = list(shape = 1/coef(lm.1)[[2]], scale = exp(coef(lm.1)[[1]]), threshold = th)
  return(estimates)
}

# @rdname qcr-internals
#' @keywords internal
.charToDistFunc = function(distribution, type = "q") {                                                           ####   .CHARTODISTFUNC-FUNCTION
  fun = NULL
  if (identical("beta", distribution)) 
    fun = eval(parse(text = paste(type, "beta", sep = "")))
  if (identical("cauchy", distribution)) 
    fun = eval(parse(text = paste(type, "cauchy", sep = "")))
  if (identical("chi-squared", distribution)) 
    fun = eval(parse(text = paste(type, "chisq", sep = "")))
  if (identical("exponential", distribution)) 
    fun = eval(parse(text = paste(type, "exp", sep = "")))
  if (identical("f", distribution)) 
    fun = eval(parse(text = paste(type, "f", sep = "")))
  if (identical("geometric", distribution)) 
    fun = eval(parse(text = paste(type, "geom", sep = "")))
  if (identical("log-normal", distribution) || identical("lognormal", distribution))         ####
    fun = eval(parse(text = paste(type, "lnorm", sep = "")))
  if (identical("log-normal3", distribution) || identical("lognormal3", distribution))       ####
    fun = eval(parse(text = paste(type, "lnorm3", sep = "")))                              ####
  if (identical("logistic", distribution)) 
    fun = eval(parse(text = paste(type, "logis", sep = "")))
  if (identical("negative binomial", distribution)) 
    fun = eval(parse(text = paste(type, "nbinom", sep = "")))
  if (identical("normal", distribution)) 
    fun = eval(parse(text = paste(type, "norm", sep = "")))
  if (identical("poisson", distribution)) 
    fun = eval(parse(text = paste(type, "pois", sep = "")))
  if (identical("t", distribution)) 
    fun = eval(parse(text = paste(type, "t", sep = "")))
  if (identical("weibull", distribution)) 
    fun = eval(parse(text = paste(type, "weibull", sep = "")))
  if (identical("weibull3", distribution))                                                   ####
    fun = eval(parse(text = paste(type, "weibull3", sep = "")))                            ####
  if (identical("gamma", distribution)) 
    fun = eval(parse(text = paste(type, "gamma", sep = "")))
  if (identical("gamma3", distribution)) 
    fun = eval(parse(text = paste(type, "gamma3", sep = "")))
  return(fun)
}


# @rdname qcr-internals
#' @keywords internal
.lfrm = function(wholeList, filterList) {
  if (!is.list(wholeList)) 
    stop(paste(deparse(substitute(wholeList)), "is not a list!"))
  if (length(wholeList) == 0) 
    return(wholeList)
  if (!is.list(filterList)) 
    stop(paste(deparse(substitute(filterList)), "is not a list!"))
  if (length(filterList) == 0) 
    return(wholeList)
  logVec = lapply(names(wholeList), "%in%", names(filterList))
  filteredList = wholeList[!unlist(logVec)]
  return(filteredList)
}


# @rdname qcr-internals
#' @keywords internal
.lfkp = function(wholeList, filterList) {
  if (!is.list(wholeList)) 
    stop(paste(deparse(substitute(wholeList)), "is not a list!"))
  if (length(wholeList) == 0) 
    return(wholeList)
  if (!is.list(filterList)) 
    stop(paste(deparse(substitute(filterList)), "is not a list!"))
  if (length(filterList) == 0) 
    return(filterList)
  logVec = lapply(names(wholeList), "%in%", names(filterList))
  filteredList = wholeList[unlist(logVec)]
  return(filteredList)
} 

# @rdname qcr-internals
#' @keywords internal
.myADTest = function(x, distribution, ...) {                                                                     ####   .MYADTESTS-FUNCTION
  #require(MASS, quietly = TRUE)
  if (missing(distribution)) 
    distribution = "normal"
  data.name = names(x)
  if (is.data.frame(x)) 
    x = x[, 1]
  dots = list(...)
  parameter = NULL
  smaller = NULL
  pFun = NULL
  tableValue = FALSE
  A = 0
  x <- sort(x[complete.cases(x)])
  n = length(x)
  if (n < 8) 
    stop("sample size must be greater than 7")
  if (n > 40) 
    warning("sample size is greater than 40")
  if (is.character(distribution)) {
    pFun = .charToDistFunc(distribution, type = "p")
    distribution = tolower(distribution)
    if (is.null(pFun)) 
      stop(paste(deparse(substitute(distribution)), " is not supported!"))
  }
  else {
    pFun = match.fun(distribution)
  }
  #    if (identical(distribution, "log-normal")) {                               ####
  #        x = log(x)                                                             ####
  #        distribution = "normal"                                                ####
  #    }                                                                          ####
  if (length(dots) == 0) {
    fittedDistr = MASS::fitdistr(x, distribution)
    parameter = fittedDistr$estimate
    if (distribution == "normal") {
      parameter["mean"] = mean(x)
      parameter["sd"] = sd(x)
    }
    p = do.call(pFun, c(list(x), as.list(parameter)))
  }
  else {
    p = pFun(x, ...)
  }
  h = (2 * seq(1:n) - 1) * (log(p) + log(1 - rev(p)))
  A = -n - mean(h)
  AA = (1 + 0.75/n + 2.25/n^2) * A
  if (AA < 0.2) {
    pval <- 1 - exp(-13.436 + 101.14 * AA - 223.73 * AA^2)
  }
  else if (AA < 0.34) {
    pval <- 1 - exp(-8.318 + 42.796 * AA - 59.938 * AA^2)
  }
  else if (AA < 0.6) {
    pval <- exp(0.9177 - 4.279 * AA - 1.38 * AA^2)
  }
  else {
    pval <- exp(1.2937 - 5.709 * AA + 0.0186 * AA^2)
  }
  if (identical(distribution, "cauchy")) {
    pval = NA
  }
  if (identical(distribution, "beta")) {
    pval = NA
  }
  if (identical(distribution, "chi-squared")) {
    pval = NA
  }
  if (identical(distribution, "f")) {
    pval = NA
  }
  if (identical(distribution, "t")) {
    pval = NA
  }
  if (identical(distribution, "geometric")) {
    pval = NA
  }
  if (identical(distribution, "poisson")) {
    pval = NA
  }
  if (identical(distribution, "negative-binomial")) {
    pval = NA
  }
  if (identical(distribution, "weibull")) {
    AWei = A * (1 + 1/sqrt(n))
    tableValue = TRUE
    smaller = TRUE
    if (AWei < 0.474) {
      pval = 0.25
      smaller = FALSE
    }
    if (AWei >= 0.474) 
      pval = 0.25
    if (AWei >= 0.637) 
      pval = 0.1
    if (AWei >= 0.757) 
      pval = 0.05
    if (AWei >= 0.877) 
      pval = 0.025
    if (AWei >= 1.038) 
      pval = 0.01
  }
  if (identical(distribution, "exponential")) {
    AExp = A * (1 + 0.6/n)
    pval = NA
    if (0.95 < AExp) {
      pval = exp(0.731 - 3.009 * AExp + 0.15 * AExp^2)
    }
    if (0.51 < AExp & AExp < 0.95) {
      pval = exp(0.9209 - 3.353 * AExp + 0.3 * AExp^2)
    }
    if (0.26 < AExp & AExp < 0.51) {
      pval = 1 - exp(-6.1327 + 20.218 * AExp - 18.663 * AExp^2)
    }
    if (AExp < 0.26) {
      pval = 1 - exp(-12.2204 + 67.459 * AExp - 110.3 * AExp^2)
    }
  }
  if (identical(distribution, "logistic")) {
    ALogist = A * (1 + 0.25/n)
    tableValue = TRUE
    smaller = TRUE
    if (ALogist < 0.426) {
      pval = 0.25
      smaller = FALSE
    }
    if (ALogist >= 0.426) {
      pval = 0.25
    }
    if (ALogist >= 0.563) {
      pval = 0.1
    }
    if (ALogist >= 0.66) {
      pval = 0.05
    }
    if (ALogist >= 0.769) {
      pval = 0.025
    }
    if (ALogist >= 0.906) {
      pval = 0.01
    }
    if (ALogist >= 1.1) {
      pval = 0.005
    }
  }
  if (identical(distribution, "gamma")) {
    tableValue = TRUE
    gammaDF = data.frame(c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, Inf), c(0.486, 0.477, 0.475, 
                                                                        0.473, 0.472, 0.472, 0.471, 0.471, 0.471, 0.47, 0.47, 0.47), c(0.657, 0.643, 0.639, 0.637, 
                                                                                                                                       0.635, 0.635, 0.634, 0.633, 0.633, 0.632, 0.632, 0.631), c(0.786, 0.768, 0.762, 0.759, 
                                                                                                                                                                                                  0.758, 0.757, 0.755, 0.754, 0.754, 0.754, 0.753, 0.752), c(0.917, 0.894, 0.886, 0.883, 
                                                                                                                                                                                                                                                             0.881, 0.88, 0.878, 0.877, 0.876, 0.876, 0.875, 0.873), c(1.092, 1.062, 1.052, 1.048, 
                                                                                                                                                                                                                                                                                                                       1.045, 1.043, 1.041, 1.04, 1.039, 1.038, 1.037, 1.035), c(1.227, 1.19, 1.178, 1.173, 
                                                                                                                                                                                                                                                                                                                                                                                 1.17, 1.168, 1.165, 1.164, 1.163, 1.162, 1.161, 1.159))
    names(gammaDF) = c("m", 0.75, 0.9, 0.95, 0.975, 0.99, 0.995)
    critCheck <- gammaDF[min(which(gammaDF$m >= parameter["shape"])), 2:length(gammaDF)] > A
    if (any(critCheck)) {
      firPos <- min(which(critCheck))
    }
    else {
      firPos <- length(gammaDF)
    }
    if (firPos == 1) {
      pValue <- 1 - as.numeric(names(gammaDF)[2])
      pval = pValue
      pValue <- paste(">", pValue)
      smaller = FALSE
    }
    else {
      pValue <- 1 - as.numeric(names(gammaDF)[firPos])
      pval = pValue
      pValue <- paste("<=", pValue)
      smaller = TRUE
    }
  }
  out = list()
  out$data.name = data.name
  out$statistic = as.vector(data.frame(A = A))
  out$parameter = parameter
  out$p.value = as.vector(data.frame(p = pval))
  out$smaller = smaller
  out$tableValue = tableValue
  out$conf.int = NULL
  out$estimate = NULL
  temp = NULL
  if (is.character(distribution)) 
    temp = as.vector(distribution)
  else temp = deparse(substitute(distribution))
  names(temp) = "distribution"
  out$null.value = temp
  out$method = paste("Anderson Darling Test for", temp, "distribution")
  class(out) = "adtest"
  return(out)
}


# @rdname qcr-internals
#' @keywords internal
print.adtest = function(x, digits = 4, quote = TRUE, prefix = "", ...) {
  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data: ", x$data.name, "\n")
  out <- character()
  if (!is.null(x$statistic)) 
    out <- c(out, paste(names(x$statistic), "=", format(round(x$statistic, 4))))
  if (!is.null(x$parameter)) 
    out <- c(out, paste(names(x$parameter), "=", format(round(x$parameter, 3))))
  if (!is.null(x$p.value)) {
    fp <- format.pval(x$p.value, digits = digits)
    if (x$tableValue) {
      if (x$smaller) 
        out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("<=", fp)))
      else out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "=") fp else paste(">", fp)))
    }
    else {
      out <- c(out, paste("p-value", if (substr(fp, 1L, 1L) == "<") fp else paste("=", fp)))
    }
  }
  cat(strwrap(paste(out, collapse = ", ")), sep = "\n")
  cat("alternative hypothesis: ")
  if (!is.null(x$null.value)) {
    if (length(x$null.value) == 1) {
      cat("true", names(x$null.value), "is not equal to", x$null.value, "\n")
    }
    else {
      cat(x$alternative, "\nnull values:\n")
      print(x$null.value, ...)
    }
  }
  if (!is.null(x$conf.int)) {
    cat(format(100 * attr(x$conf.int, "conf.level")), "percent confidence interval:\n", format(c(x$conf.int[1L], x$conf.int[2L])), "\n")
  }
  if (!is.null(x$estimate)) {
    cat("sample estimates:\n")
    print(x$estimate, ...)
  }
  cat("\n")
  invisible(x)
} 