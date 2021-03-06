% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.mqcs.R
\name{plot.mqcs}
\alias{plot.mqcs}
\alias{plot.mqcs.t2}
\alias{plot.mqcs.mcusum}
\alias{plot.mqcs.mewma}
\title{Plot method for 'mqcs' objects}
\usage{
\method{plot}{mqcs}(x, title, subtitle, xlab, ylab, ylim, ...)

\method{plot}{mqcs.t2}(
  x,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  ylim = NULL,
  ...
)

\method{plot}{mqcs.mcusum}(
  x,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  ylim = NULL,
  ...
)

\method{plot}{mqcs.mewma}(
  x,
  title = NULL,
  subtitle = NULL,
  xlab = NULL,
  ylab = NULL,
  ylim = NULL,
  ...
)
}
\arguments{
\item{x}{An Object of class 'mqcs' (Multivarite Quality Control Statical)}

\item{title}{An overall title for the plot.}

\item{subtitle}{A sub title for the plot.}

\item{xlab}{A title for the 'x' axis.}

\item{ylab}{A title for the 'y' axis.}

\item{ylim}{The 'y' limits of the plot.}

\item{...}{Arguments to be passed to or from methods.}
}
\description{
Generic function for plotting Multivarite charts of object of class 'mqcs' 
to perform statistical quality control.
}
\examples{
\dontrun{
## 
## Continuous data 
## 
data(dowel1) 
data.mqcd <- mqcd(dowel1)  
res.mqcs <- mqcs.mcusum(data.mqcd) 
plot(res.mqcs, title =" MCUSUM Control Chart ", subtitle="Database dowel1",
     xlab = "Observations", ylab = "MCUSUM", ylim = c(-1,6.5)) 

res1.mqcs <- mqcs.mewma(data.mqcd) 
plot(res1.mqcs, title =" MEWMA Control Chart", subtitle="Database dowel1",
     xlab = "Observations", ylab = "MEWMA", ylim = c(-1,10))
      
res2.mqcs <- mqcs.t2(data.mqcd)
plot(res2.mqcs, title =" Hotelling Control Chart",subtitle="Database dowel1",
     xlab = "Observations", ylab = "T2 Hotelling", ylim = c(-1,14))}
}
