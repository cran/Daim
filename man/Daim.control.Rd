\name{Daim.control}
\alias{Daim.control}
\title{Control parameters for the diagnostic accuracy of models.}
\description{
  Control of resampling methods.
}
\usage{
Daim.control(method="boot", number = 100, replace = TRUE, 
             boot.size = 1, k = 10, k.runs = 1,
             dependency = list(var = NULL, keep.id = FALSE))
}
\arguments{
  \item{method}{the resampling method:\cr 
  \code{boot} - bootstrapping, \code{cv} - k-fold cross-validation.}
  \item{number}{the number of bootstrap samples}
  \item{replace}{a logical indicating whether sampling of observations is 
                 done with or without replacement}
  \item{boot.size}{percentage of observations (\code{0 < boot.size < 1}) to draw without replacement (only relevant if \code{replace = FALSE}). In this case subsampling is performed.}
  \item{k}{the number of folds}
  \item{k.runs}{the number of runs of k-fold cross-validations}
  \item{dependency}{to handle paired data structures}
}
\references{ 

Breiman L. (1996), Out-Of-Bag Estimation. \emph{Technical Report}.\cr
\url{ftp://ftp.stat.berkeley.edu/pub/users/breiman/OOBestimation.ps.Z}.\cr

Adler W, Brenning A, Potapov S, Schmid M, Lausen B. (2011),\cr 
Ensemble classification of paired data. \emph{Computational Statistics and Data Analysis}.\cr
55(5):1933-1941.

}

\seealso{
  \code{\link{Daim}}
}
\examples{

  ###
  ### for bootstrap
  ###

  Daim.control(method="boot", number=100)

  ###
  ### for cross-validation
  ###

  Daim.control(method="cv", k=10, k.runs=10)

  ###
  ### for subagging or subsampling
  ###

  Daim.control(method="boot", number=100, replace=FALSE, boot.size=0.9)

}

\keyword{classif}
