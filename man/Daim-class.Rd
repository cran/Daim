\name{Daim Class}
\docType{class}
\alias{Daim-class}
\title{Class "Daim"}

\description{Object includes the results of a performance evaluation 
from several cross-validation or bootstrapping runs.}

\section{Objects from the Class}{
Objects can be created by using the \code{Daim} or \code{performDaim} functions.
}


\section{Slots}{
  \describe{
  \item{\code{call}: }{the matched call.}
  \item{\code{formula}: }{the formula supplied.}
  \item{\code{method}: }{the list of control parameters.}
  \item{\code{err632p}: }{the \code{0.632+} estimation of the misclassification error.}
  \item{\code{err632}: }{the \code{0.632} estimation of the misclassification error.}
  \item{\code{errloob}: }{the \code{LOOB} estimation of the misclassification error.}
  \item{\code{errapp}: }{the apparent error.}
  \item{\code{sens632p}: }{the \code{0.632+} estimation of the sensitivity.}
  \item{\code{spec632p}: }{the \code{0.632+} estimation of the specificity.}
  \item{\code{sens632}: }{the \code{0.632} estimation of the sensitivity.}
  \item{\code{spec632}: }{the \code{0.632} estimation of the specificity.}
  \item{\code{sensloob}: }{the \code{LOOB} estimation of the sensitivity.}
  \item{\code{specloob}: }{the \code{LOOB} estimation of the specificity.}
  \item{\code{sensapp}: }{the apparent sensitivity.}
  \item{\code{specapp}: }{the apparent specificity.}
  \item{\code{roc}: }{a data frame with estimated values of sensitivity and specificity for a variety of cutoffs.}
  \item{\code{sample.roc}: }{a list in which each entry contains the values of the ROC curve of this special sample or cross-validation run.}
  \item{\code{sample.data}: }{a data frame with the results of this particular sample or cross-validation run.}
  }
}

\seealso{
  \code{\link{Daim}}, \code{\link{performDaim}}, \code{\link{auc.Daim}}, \code{\link{roc.area.Daim}}
}
\keyword{classes}
