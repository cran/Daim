\name{auc}
\alias{auc.numeric}
\title{The area under the ROC curve (AUC)}
\description{
This function computs the area under an ROC curve.
}
\usage{
\method{auc}{numeric}(x, y, \dots)
}
\arguments{
  \item{x}{sensitivity.}
  \item{y}{specificity.}
  \item{\dots}{additional parameters.}
}
\value{
a scalar number
}
\seealso{
  \code{\link{auc.Daim}}, \code{\link{roc.area.Daim}}, \code{\link{Daim}}
}
\examples{

  data(Daim.data1)
  perform <- performDaim(Daim.data1$prob.oob, Daim.data1$labels,
		Daim.data1$prob.app)

  ####
  #### compute the .632+ estimation of the AUC.
  ####

  auc(perform$roc$sens632p,perform$roc$spec632p)

}
\keyword{classif} \keyword{manip}
