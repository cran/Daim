\name{auc.Daim}
\alias{auc}
\alias{auc.Daim}
\title{The area under the ROC curve (AUC)}
\description{
  This function computs the area under an ROC curve.
}
\usage{
\method{auc}{Daim}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{Daim}.}
  \item{\dots}{additional parameters.}
}
\value{ a list with following components :
\item{auc.632p}{the .632+ estimation of the AUC.}
\item{auc.632}{the .632 estimation of the AUC.}
\item{auc.loob}{the LOOB estimation of the AUC.}
\item{auc.app}{the apparent estimation of the AUC.}
\item{auc.samples}{the AUC values for each bootstap sample.}
}
\seealso{
  \code{\link{auc.numeric}}, \code{\link{Daim}}, \code{\link{plot.Daim}}, \code{\link{performDaim}}
}
\examples{

  data(Daim.data1)
  perform <- performDaim(Daim.data1$prob.oob, Daim.data1$labels,
		Daim.data1$prob.app)

  auc(perform)
  roc.area(perform)

}
\keyword{classif} \keyword{manip}
