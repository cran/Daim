\name{summary.Daim}
\alias{summary.Daim}
\alias{summary.Daim.vector}
\alias{summary.Daim.list}
\title{Summarizing a Daim Object}
\description{
  \code{summary} method for class \code{"Daim"}.
}
\usage{
\method{summary}{Daim}(object, \dots)
\method{summary}{Daim.vector}(object, \dots)
\method{summary}{Daim.list}(object, \dots)
}
\arguments{
  \item{object}{an object of class \code{Daim}.}
  \item{...}{further arguments passed to or from other methods.}
}
\seealso{
  \code{\link{Daim}}
}
\examples{
  library(Daim)
  data(Daim.data1)
  perform <- performDaim(Daim.data1$prob.oob, 
		Daim.data1$labels, Daim.data1$prob.app)
  perform
  summary(perform)

  data(Daim.data3)
  head(Daim.data3)

  M <- roc(Daim.data3[,2], Daim.data3[,1], labpos="pos")
  summary(M)

  M <- roc(Daim.data3[,-1], Daim.data3[,1], labpos="pos")
  summary(M)

}
\keyword{manip}
