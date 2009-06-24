\name{print.Daim}
\alias{print.Daim}
\alias{print.Daim.vector}
\title{Print Method for Daim Object}
\description{
Print object of class \code{Daim}.
}
\usage{
\method{print}{Daim}(x, digits=max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{x}{object of class \code{Daim}.}
  \item{digits}{a non-null value for \code{digits} specifies the minimum
    number of significant digits to be printed in values.}
  \item{\dots}{additional arguments.}
}
\details{
An object of class \code{Daim} is printed.
}
\seealso{                                              
  \code{\link{summary.Daim}}, \code{\link{Daim}}, \code{\link{plot.Daim}}   
}
\examples{
  library(Daim)
  data(Daim.data1)
  perform <- performDaim(Daim.data1$prob.oob, 
		Daim.data1$labels, Daim.data1$prob.app)
  perform
}                                                      
\keyword{print} \keyword{manip}
