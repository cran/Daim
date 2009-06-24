\name{roc}
\alias{roc}
\alias{roc.integer}
\alias{roc.numeric}
\alias{roc.matrix}
\alias{roc.data.frame}
\title{Compute a ROC curve}
\description{
  This function computes sensitivity and specificity for a variety of cut-points.
}
\usage{
\method{roc}{integer}(x, ...)
\method{roc}{numeric}(x, labels, labpos, thres=NULL, \dots)
\method{roc}{matrix}(x, labels, labpos, thres=NULL, \dots)
\method{roc}{data.frame}(x, ...)

}
\arguments{
  \item{x}{an object (\code{vector, matrix, data.frame}) used for prediction.}
  \item{labels}{a vector containing the true class labels. This can be a factor or character vector.}
  \item{labpos}{a character string of the variable \code{labels} that defines a "positive" event.}
  \item{thres}{a numeric vector with the cutoff values. By default, the \code{x} define the grid of cut-points.}
  \item{...}{additional parameters.}
}

\seealso{
  \code{\link{plot.Daim}}, \code{\link{auc.Daim}}
}
\examples{

  data(Daim.data3)

  M <- roc(Daim.data3[,2:5], Daim.data3$Gold, "pos")
  summary(M)
  plot(M,color=c("black","blue","green3","red"))

  roc.area(M)

}
\keyword{hplot} \keyword{classif}
