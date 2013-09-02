\name{performDaim}
\alias{performDaim}
\alias{performDaim.matrix}
\alias{performDaim.data.frame}
\title{Function to evaluate predictions}
\description{
  Perform the evaluation of predictions.
}
\usage{
\method{performDaim}{matrix}(x, labels, prediction=NULL,
            thres=seq(0,1,by=0.01),
            cutoff=0.5, labpos="1", \dots)
		
\method{performDaim}{data.frame}(x, labels, prediction=NULL,
            thres=seq(0,1,by=0.01),
            cutoff=0.5, labpos="1", \dots)
}
\arguments{
  \item{x}{a matrix or data frame containing the predictions.}
  \item{labels}{a vector containing the true class labels. It can be a factor or character vector.}
  \item{prediction}{a vector containing the prediction probability obtained by a model: see \code{\link{Daim}}.}
  \item{thres}{a numeric vector with the cutoff values.}
  \item{cutoff}{the cutoff value for error estimation. This can be a numeric value or a character string.
	If the \code{cutoff} set to:\cr
	\code{"loob"} - the optimal cut-point corresponding to \code{loob} estimation of the sensitivity and the specificity.\cr
	\code{"0.632"} - the optimal cut-point corresponding to \code{0.632} estimation of the sensitivity and the specificity.\cr 
	\code{"0.632+"} - the optimal cut-point corresponding to \code{0.632+} estimation of the sensitivity and the specificity.}
  \item{labpos}{a character string of the variable \code{labels} that defines a "positive" event.}
  \item{\dots}{additional parameters.}
}

\seealso{
  \code{\link{Daim}}, \code{\link{plot.Daim}}
}
\examples{
  data(Daim.data1)
  (perform <- performDaim(Daim.data1$prob.oob, Daim.data1$labels, Daim.data1$prob.app))
}
\keyword{classif} \keyword{models}
