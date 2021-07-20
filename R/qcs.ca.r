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
# Main function to create a 'qcs.ca' object
#-----------------------------------------------------------------------------#
##' Capability Analysis
##' 
##' Calculates the process capability indices cp, cpk, cpL cpU, cpm, cpmk for a qcs object and normal distribution. 
##' Also, this function calculates confidence limits for \eqn{C_p}{C_p} using the method described by Chou et al. (1990). 
##' Approximate confidence limits for \eqn{C_{pl}}{C_pl}, \eqn{C_{pu}}{C_pu} and  \eqn{C_{pk}}{C_pk} are computed using the method in Bissell (1990). 
##' Confidence limits for \eqn{C_{pm}}{C_pm} are based on the method of Boyles (1991); this method is approximate and it assumes 
##' the target is midway between the specification limits.
##' Moreover, calculates the process capability indices cnp, cnpk, cnpm, cnpmk for a qcs object. 
##' A histogramm with a density curve is displayed along with the specification limits, a 
##' Quantile-Quantile Plot for the specified distribution and contour graph is plotted for estimate the indice cpm.
##' @aliases qcs.ca 
##' @param object qcs object of type \code{"qcs.xbar"} or \code{"qcs.one"}.
##' @param limits A vector specifying the lower and upper specification limits.
##' @param target A value specifying the target of the process. 
##' If is \code{NULL}, the target is set at the middle value bewteen specification limits.
##' @param std.dev A value specifying the within-group standard deviation.
##' @param nsigmas A numeric value specifying the number of sigmas to use.
##' @param confidence A numeric value between 0 and 1 specifying the probabilities for computing the quantiles.
##' This values is used only when object values is provided. The default value is 0.9973.
##' @param plot Logical value indicating whether graph should be plotted.
##' @param main Title of the plot.
##' @param ... Arguments to be passed to or from methods.
##' @export
##' @references 
##' Montgomery, D.C. (1991) \emph{Introduction to Statistical Quality Control}, 2nd
##' ed, New York, John Wiley & Sons. \cr
##' Tong, L.I. and Chen, J.P. (1998), \emph{Lower con???dence limits of process capability 
##' indices for nonnormal process distributions.} International Journal of Quality & Reliability Management, 
##' Vol. 15 No. 8/9, pp. 907-19.\cr
##' Vannman, K (1995) \emph{A Unified Approach to Capability Indices}. Statitica Sinica,5,805-820.\cr
##' Vannman, K. (2001). \emph{A Graphical Method to Control Process Capability}. Frontiers in Statistical Quality Control, 
##' No 6, Editors: H-J Lenz and P-TH Wilrich. Physica-Verlag, Heidelberg, 290-311.\cr
##' Hubele and Vannman (2004). \emph{The E???ect of Pooled and Un-pooled Variance Estimators on Cpm When Using Subsamples}.
##' Journal Quality Technology, 36, 207-222.\cr
##' @examples
##' library(qcr)
##' data(pistonrings) 
##' xbar <- qcs.xbar(pistonrings[1:125,],plot = TRUE)
##' LSL=73.99; USL=74.01
##' limits = c(lsl = 73.99, usl = 74.01)
##' qcs.ca(xbar, limits = limits)
qcs.ca <- function (object, 
                      limits = c(lsl = -3, usl = 3),
                      target = NULL, std.dev = NULL,nsigmas = 3,
                      confidence = 0.9973, plot = TRUE, main = NULL, ...)
 {

   if (missing(object)){
     if (!inherits(object, "qcs"))
       stop("an object of class 'qcs' is required")
     if (!(object$type == "xbar" | object$type == "one"))
       stop("Process Capability Analysis only available for charts type
            \"qcs.xbar\" and \"qcs.one\" charts")
   }

   x<-object$statistics
   n = sum(object$sizes)
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
   
   if (confidence < 0 | confidence > 1)
     stop("conf.level must be a value between 0 and 1")

   alpha <- 1 - confidence
   
   if (!is.null(lsl) && !is.null(usl)){
     cp <- (usl - lsl)/(2 * nsigmas * std.dev)
     cp.limits <- cp * sqrt(qchisq(c(alpha/2, 1 - alpha/2), n -
                                    1)/(n - 1))
     
     cpm <- cp/sqrt(1 + ((center - target)/std.dev)^2)
     
     df <- n * (1 + ((center - target)/std.dev)^2)/(1 + 2 * ((center -
                                                                target)/std.dev)^2)
     cpm.limits <- cpm * sqrt(qchisq(c(alpha/2, 1 - alpha/2),
                                     df)/df)
     
     cpmk <- qcs.cp(object,parameters = c(1,1),limits = c(lsl,usl),
                    contour = FALSE)[1]
     
     cnp <- qcs.cpn(object,parameters = c(0,0),limits = c(lsl,usl))[1]
     cnpk <- qcs.cpn(object,parameters = c(1,0),limits = c(lsl,usl))[1]
     cnpm <- qcs.cpn(object,parameters = c(0,1),limits = c(lsl,usl))[1]
     cnpmk <- qcs.cpn(object,parameters = c(1,1),limits = c(lsl,usl))[1]
     
     
     
   }else{
     cp <- NA
     cpm <- NA
     cpmk <- NA
     cnp <- NA
     cnpk <- NA
     cnpm <- NA
     cnpmk <- NA
    }
   
   if (!is.null(usl)) {
     cpu <- (usl - center)/(nsigmas * std.dev)
     cpu.limits <- cpu * (1 + c(-1, 1) * qnorm(confidence) *
                              sqrt(1/(9 * n * cpu^2) + 1/(2 * (n - 1))))
   }
   
   if (!is.null(lsl)) {
     cpl <- (center - lsl)/(nsigmas * std.dev)
     cpl.limits <- cpl * (1 + c(-1, 1) * qnorm(confidence) *
                              sqrt(1/(9 * n * cpl^2) + 1/(2 * (n - 1))))
   }
   
   cpk = min(cpu, cpl)
   cpk.limits <- cpk * (1 + c(-1, 1) * qnorm(1 - alpha/2) *
                            sqrt(1/(9 * n * cpk^2) + 1/(2 * (n - 1))))
   


   
   names(cp.limits) <- names(cpk.limits) <- names(cpm.limits) <- 
     c(paste(round(100 *alpha/2, 1), "%", sep = ""), 
       paste(round(100 * (1 - alpha/2), 1), "%", sep = ""))
   
   if (is.na(lsl))
     exp.LSL <- NA
   else {
     exp.LSL <- pnorm((lsl - center)/std.dev) * 100
     if (exp.LSL < 0.01)
       exp.LSL <- 0
   }
   
   if (is.na(usl))
     exp.USL <- NA
   else {
     exp.USL <- (1 - pnorm((usl - center)/std.dev)) * 100
     if (exp.USL < 0.01)
       exp.USL <- 0
   }
   
   obs.LSL <- sum(x < lsl)/n * 100
   obs.USL <- sum(x > usl)/n * 100

   tab <- cbind(c(cp, cpl, cpu, cpk, cpm), rbind(cp.limits,
                                                    cpl.limits, cpu.limits, cpk.limits, cpm.limits))
   rownames(tab) <- c("Cp", "Cp_l", "Cp_u", "Cp_k", "Cpm")
   colnames(tab) <- c("Value", names(cp.limits))
   
   tabn <- cbind(c(cnp, cnpk, cnpm,cnpmk))
   rownames(tabn) <- c("CNp", "CNpK", "CNpm", "CNpmk")
   colnames(tabn) <- c("Value")
   
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
     
     xlim <- range(x, usl,lsl, target, na.rm = TRUE)
     xlim <- xlim + diff(xlim) * c(-0.1, 0.1)
     xx <- seq(min(xlim), max(xlim), length = 250)
     dx <- dnorm(xx, center, std.dev)
     if (class(object)[1]=="qcs.xbar"){
        std.dev2 <- sd(object$qcd[,1])
        dx2<- dnorm(xx, center, std.dev2) 
     }
     h <- hist(x[, 1], plot = FALSE)
     ylim <- range(h$density, dx)
     ylim <- ylim + diff(ylim) * c(0, 0.05)
     xlim <- range(x[, 1], usl, lsl)
     xlim <- xlim + diff(xlim) * c(-0.2, 0.2)


     par(bg="#CCCCCC",mar = c(0, 0, 0, 0) + 0.1)

     par(oma = c(2, 4, 7, 4) + 0.1)
     layout(matrix(c(1, 1, 2, 3, 1, 1, 4, 5, 1, 1, 6, 7), 
                   nrow = 3, byrow = TRUE))

          do.call(hist, c(list(x[, 1], freq = FALSE, xlim = xlim,
                          ylim = ylim, main = "",col = "#EEEEEE")))

     tempList = parList
     tempList$col = "#EEEEEE"
     tempList$border = NULL

     abline(v = 0)
     abline(h = 0)
     lines(xx, dx, lwd =1 , col = "red", lty = 1)
     if (class(object)[1]=="qcs.xbar") 
       lines(xx, dx2, lwd = lineWidth, col = "blue", lty = lineType)
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

     if (is.null(main))
         main = paste("Process Capability for",data.name)
     
     if (class(object)[1]=="qcs.xbar") 
       legend("topright",legend = c("ST","LT"), 
              col=c(lineCol,"blue"),lty = c(1,2),bty="n")

#################
     par(mar = c(0, 0, 0, 0) + 0.1)
     title(main = main, outer = TRUE)
     plot(0:5, 0:5, type = "n", axes = FALSE, xlab = "", ylab = "",
          main = "")
     
     rect(par("usr")[1],
          par("usr")[3],
          par("usr")[2],
          par("usr")[4],
          col  =  "white")
     
     
     dis="normal"
     adTestStats = .myADTest(x, dis)
     A = adTestStats$statistic
     p = adTestStats$p.value

     text(3.2, 5, "Process Data", pos = 2,
          cex = cex.val, col = "blue")
     
     text(4, 4.5, paste("Sample n =",n), pos = 2,
          cex = cex.val)
     
     text(4, 4, paste("A =",format(A, digits = 3)),
          pos = 2, cex = cex.val)
     
     if (!is.null(adTestStats$smaller) && adTestStats$smaller)
       text(4, 3.5, paste("p <", format(p, digits = 3)),
            pos = 2, cex = cex.val)
     if (!is.null(adTestStats$smaller) && !adTestStats$smaller)
       text(4, 3.5, paste("p >=", format(p,
                                         digits = 3)), pos = 2, cex = cex.val)
     if (is.null(adTestStats$smaller))
       text(4, 3.5, paste("p =", format(p, digits = 3)),
            pos = 2, cex = cex.val)
     
     text(4, 3, paste("LSL =",format(lsl, digits = 3)), 
          pos = 2,cex = cex.val)

     text(4, 2.5, paste("Target =",format(target, digits = 3)), 
          pos = 2,cex = cex.val)
     
     text(4, 2, paste("USL =",format(usl, digits = 3)), 
          pos = 2,cex = cex.val)
     
     text(4, 1.5, paste("std.dev (ST) =",format(std.dev, digits = 3)), pos = 2,
          cex = cex.val)
     
     text(4, 1, paste("mean =",format(center, digits = 3)), 
          pos = 2,cex = cex.val)
     
    if (class(object)[1]=="qcs.xbar") 
     text(4, 0.5, paste("std.dev (LT) =",format(std.dev2, digits = 3)), pos = 2,
          cex = cex.val)

     
     

     qqPlot(x[, 1], y = dis, ylab = "", main = "",
            axes = FALSE, bounds.lty = bounds.lty,
            bounds.col = bounds.col)
     box()
     do.call(grid,list(col = "#EEEEEE"))
     
 #################
     par(mar = c(0, 0, 0, 0) + 0.1)
     title(main = main, outer = TRUE)
     plot(0:5, 0:5, type = "n", axes = FALSE, xlab = "", ylab = "",
          main = "")
     
     rect(par("usr")[1],
          par("usr")[3],
          par("usr")[2],
          par("usr")[4],
          col  =  "white")
     
     text(4.7, 5, " Parametric Capability Process (ST)", pos = 2,
          cex = cex.val, col = "blue")
 
     
     text(4, 4.2, paste("Cp =",format(cp, digits = 3)), pos = 2,
          cex = cex.val)
     
     text(4, 3.4, paste("Cpu =",format(cpu, digits = 3)),
          pos = 2, cex = cex.val)
     
     text(4, 2.6, paste("Cpl =",format(cpl, digits = 3)),
          pos = 2, cex = cex.val)
     
     text(4, 1.8, paste("Cpk =",format(cpk, digits = 3)), 
          pos = 2,cex = cex.val)
     
     text(4, 1, paste("Cpm =",format(cpm, digits = 3)), 
          pos = 2,cex = cex.val)
     
     text(4, 0.2, paste("Cpmk =",format(cpmk, digits = 3)), 
          pos = 2,cex = cex.val)
     
     plot(0:5, 0:5, type = "n", axes = FALSE, xlab = "", ylab = "",
          main = "")
     rect(par("usr")[1],
          par("usr")[3],
          par("usr")[2],
          par("usr")[4],
          col  =  "white")
     
     text(5, 5, " NonParametric Capability Process (ST)", pos = 2,
          cex = cex.val, col = "blue")

     
     text(4, 4, paste("CNp =",format(cnp, digits = 3)), pos = 2,
          cex = cex.val)
     
     text(4, 3, paste("CNpk =",format(cnpk, digits = 3)),
          pos = 2, cex = cex.val)
     
     text(4, 2, paste("CNpm =",format(cnpm, digits = 3)), 
          pos = 2,cex = cex.val)
     
     text(4, 1, paste("CNpmk =",format(cnpmk, digits = 3)), 
          pos = 2,cex = cex.val)
     
     
     
#################     
     
     par(mar = c(0, 0, 0, 0) + 0.1)
     title(main = main, outer = TRUE)
     plot(0:5, 0:5, type = "n", axes = FALSE, xlab = "", ylab = "",
          main = "")
     
     rect(par("usr")[1],
          par("usr")[3],
          par("usr")[2],
          par("usr")[4],
          col  =  "white")
     
     text(4, 5, " Performance - ST [%]", pos = 2,
          cex = cex.val, col = "blue")
     
     
     text(4, 4.2, paste("Exp<LSL =",format(exp.LSL, digits = 3)), pos = 2,
          cex = cex.val)
     
     text(4, 3.4, paste("Exp>USL =",format(exp.USL, digits = 3)),
          pos = 2, cex = cex.val)
     
     text(4, 2.6, paste("Exp Total =",format(exp.USL+exp.LSL, digits = 3)),
          pos = 2, cex = cex.val)
     
     text(4, 1.8, paste("Obs<LSL =",format(obs.LSL, digits = 3)), 
          pos = 2,cex = cex.val)
     
     text(4, 1, paste("Obs>USL =",format(obs.USL, digits = 3)), 
          pos = 2,cex = cex.val)
     
     text(4, 0.2, paste("Obs Total =",format(obs.USL+obs.LSL, digits = 3)), 
          pos = 2,cex = cex.val)
     
     par(mar = c(0, 0, 0, 0) + 0.2)
     hat.cpm <- qcs.hat.cpm(object, limits = limits,contour = FALSE)
     delta.t<- hat.cpm$delta.t
     gamma.t<-hat.cpm$gamma.t
     delta.e<-hat.cpm$delta.e
     gamma.e<-hat.cpm$gamma.e
     point.delta<-hat.cpm$delta
     point.gamma<-hat.cpm$gamma
     ylim<-c(0,point.gamma+0.2)
     xlim<-c(-max(delta.t)-0.2,max(delta.t)+0.2)
     plot(delta.t, gamma.t, type = "n", axes = FALSE, xlab = "", ylab = "",
          main = "",ylim=ylim,xlim=xlim)

     do.call(grid,list(col = "#EEEEEE"))
     abline(v = 0)
     lines(delta.t, gamma.t, col = "black", lwd = 2,lty = 1)
     lines(delta.e, gamma.e, col = "blue", lwd = 2,lty = 2)
     points(x = point.delta,y = point.gamma,col="red",pch = 21, bg = "red", 
            lwd=2)
     legend("topleft", c("Theory", "Empirical"), 
            lwd = c(2,2), bty="n",lty = c(1,2), col =c("black","blue"), cex =0.7)
     box()
     
     
     
     
     
#################     
     on.exit(par(par.orig))
   }

   digits = 4  
   cat("\nProcess Capability Analysis\n")
   cat("\nCall:\n", deparse(match.call()), "\n\n", sep = "")
   cat(paste(formatC("Number of obs = ", width = 16), formatC(n, 
                                                              width = 12, flag = "-"), formatC("Target = ", width = 10), 
             formatC(signif(target, digits = digits), 
                     flag = "-"), "\n", sep = ""))
   cat(paste(formatC("Center = ", width = 16), formatC(signif(center,digits = digits), width = 12, flag = "-"),
             formatC("LSL = ", width = 10),formatC(signif(lsl,digits = digits), flag = "-")), "\n", sep = "")
   
   cat(paste(formatC("StdDev = ", width = 16), formatC(signif(std.dev,digits = digits), width = 12, flag = "-"), 
             formatC("USL = ", width = 10), formatC(signif(usl, digits = digits), flag = "-")), "\n", sep = "")
   
   cat("\nParemetric Capability indices:\n\n")
   print(tab, digits = 4, na.print = "", print.gap = 2)
   cat("\n")
   
   
   cat("\nNon parametric Capability indices:\n\n")
   print(tabn, digits = 4, na.print = "", print.gap = 2)
   cat("\n")
   
   
   cat("\nPPM:\n\n")
   cat(paste(formatC("Exp<LSL", width = 16), formatC(exp.LSL* 1e+06, width = 12, flag = "-"), 
             formatC("Obs<LSL", width = 10), formatC(obs.LSL* 1e+06, flag = "-")), "\n", sep = "")
   
   cat(paste(formatC("Exp>USL", width = 16), formatC(exp.USL* 1e+06, width = 12, flag = "-"), 
             formatC("Obs>USL", width = 10), formatC(obs.USL* 1e+06, flag = "-")), "\n", sep = "")
   
   cat(paste(formatC("Exp Total", width = 16), formatC((exp.USL+exp.LSL)* 1e+06, width = 12, flag = "-"), 
             formatC("Obs Total", width = 10), formatC((obs.USL+obs.LSL)* 1e+06, flag = "-")), "\n", sep = "")
   
   cat("\nTest:\n\n")
   print(adTestStats)



   result <- list(cp = cp, cpk = cpk, cpl = cpl, cpu = cpu, cpm = cpm, cpmk = cpmk, 
                  cnp = cnp, cnpk = cnpk,cnpm = cnpm, cnpmk = cnpmk,
                  exp.LSL = exp.LSL* 1e+06, exp.USL = exp.USL* 1e+06, exp.T = (exp.LSL+exp.USL)* 1e+06, 
                  usl = usl, lsl = lsl, target = target,
                  obs.LSL = obs.LSL* 1e+06, obs.USL = obs.USL* 1e+06, obs.T = (obs.LSL+obs.USL)* 1e+06)

   invisible(result)
 }


