% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qcr.R
\docType{data}
\name{pistonrings}
\alias{pistonrings}
\title{Piston rings data}
\format{
A data frame with 200 observations on the following 3 variables.
\describe{ 
\item{diameter}{A numeric vector.}
\item{sample}{Sample ID}
\item{trial}{Trial sample indicator (\code{TRUE}/\code{FALSE})}
 }
}
\description{
This is a database that contains information on the inside 
diameter of piston rings. Piston rings are manufactured for an automobile 
engine that are produced through a forging process. This process is 
considered "in control" and is measured on 25 samples, each of size 5.
}
\examples{

data(pistonrings)
attach(pistonrings)
summary(pistonrings)
boxplot(diameter ~ sample)
plot(sample, diameter, cex=0.7)
lines(tapply(diameter,sample,mean))
detach(pistonrings)
}
\references{
Montgomery, D.C. (1991) \emph{Introduction to Statistical
Quality Control}, 2nd ed, New York, John Wiley & Sons, pp. 206--213
}
\keyword{datasets}
