\name{Daim.data3}
\alias{Daim.data3}
\docType{data}
\title{Data set: Artificial data for use with Daim}
  
\description{
  The Artificial data set.
}
\usage{data(Daim.data1)}
\format{A \code{data.frame} containing 4 variables and 196 observations.}

\examples{
  data(Daim.data3)
  head(Daim.data3)
  M <- roc(Daim.data3[,-1], Daim.data3[,1],labpos="pos")
  summary(M)
}
\keyword{datasets}
