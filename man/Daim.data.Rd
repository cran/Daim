\name{Daim.data 1 & 2}
\alias{Daim.data1}
\alias{Daim.data2}
\docType{data}
\title{Data sets: Artificial bootstrap data sets for use with Daim.}
  
\description{
  The data sets containing 100 sets of predictions, id's and
  corresponding labels were obtained from 100 bootstrap samples.
}
\usage{
data(Daim.data1)

data(Daim.data2)
}
\format{A three element list. The first element,
  \code{prob.oob}, is itself a matrix where rows are 
  observations and colums are the (bootstrap) samples. Each of these 100
  colums is a numerical prediction for each bootstrap sample.
  The second list entry, \code{prob.app} is a vector of 
  numerical apparent predictions.
  The third list entry, \code{labels} is a 100
  element list in which each element is a vector of true class
  labels corresponding to the predictions.}

\examples{
  data(Daim.data1)
  (perform <- performDaim(Daim.data1$prob.oob, Daim.data1$labels, Daim.data1$prob.app))
  plot(perform)
}
\keyword{datasets}
