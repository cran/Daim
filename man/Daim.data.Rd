\name{Daim.data1}
\alias{Daim.data1}
\alias{Daim.data2}
\docType{data}
\title{Data set: Artificial bootstrap data for use with Daim}
  
\description{
  The data set containing 100 sets of predictions, id's and
  corresponding labels were obtained from 100 bootstrap samples.
}
\usage{data(Daim.data1)}
\format{A three element list. The first element,
  \code{Daim.data1$prob.oob}, is itself a matrix where rows are 
  observations and colums are the (bootstrap) samples. Each of these 100
  colums is a numerical prediction for each bootstrap sample.
  The second list entry, \code{Daim.data1$prob.app} is a vector of 
  numerical apparent predictions.
  The third list entry, \code{Daim.data1$labels} is a 100
  element list in which each element is a vector of true class
  labels corresponding to the predictions.}

\examples{
  data(Daim.data1)
  perform <- performDaim(Daim.data1$prob.oob, Daim.data1$labels,
	Daim.data1$prob.app)
  perform
  plot(perform)
}
\keyword{datasets}
