\name{deLong.test}
\alias{deLong.test}
\title{DeLong test}
\description{
  A nonparametric test for comparing AUC of two or more correlated ROC curves.
}
\usage{
deLong.test(x, labels, labpos, ref=NULL, conf.level=0.95)
}
\arguments{
  \item{x}{at least 2 dimensional matrix containing the observations of numeric predictors.}
  \item{labels}{a vector containing the true class labels. It can be a factor or character vector.}
  \item{labpos}{a character string of the variable \code{labels} that defines a "positive" event.}
  \item{ref}{an index of reference predictors. Relevant only by pairwise comparisons.}
  \item{conf.level}{confidence level of the interval.}
}
\references{
	DeLong, E.R., DeLong, D.M., and Clarke-Pearson, D.L. (1988).\cr
	Comparing the Areas Under Two or More Correlated Receiver Operating Characteristics Curves: A Nonparametric Approach. 
	\emph{Biometrics}, \bold{44}, 837--845.\cr
}
\seealso{
  \code{\link{roc}}, \code{\link{Daim}}, \code{\link{plot.Daim}}, \code{\link{performDaim}}
}
\examples{


y <- rbinom(100, 1, .5)
x1 <- rnorm(100) + 1.5 * y
x2 <- rnorm(100) + .5 * y
x3 <- rnorm(100) + 2.5 * y
x <- data.frame(x1,x2,x3)
y <- as.factor(y)

(M1 <- deLong.test(x, labels=y, labpos="1"))



\dontrun{  
## just because of checking time on CRAN
  

data(iris)
(M2 <- deLong.test(iris[,1:4], labels=iris[,5], labpos="versicolor"))
}
}
\keyword{classif}
