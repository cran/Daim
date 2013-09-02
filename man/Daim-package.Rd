\name{Daim-package}
\alias{Daim-package}
\docType{package}
\title{Diagnostic accuracy of classification models}
\description{
Several functions for evaluating the accuracy of classification models. The package provides the following 
performance measures: repeated k-fold cross-validation, 0.632 and 0.632+ bootstrap estimation of the 
misclassification rate, sensitivity, specificity and AUC. If an application is computationally intensive, 
parallel execution can be used to reduce the computational effort.
}
\details{
\tabular{ll}{
Package: \tab Daim\cr
Type: \tab Package\cr
Version: \tab 1.1-0\cr
Date: \tab 2013-10-09\cr
License: \tab GPL (>= 2)\cr
}
See the help files for the following functions for more information:

\code{\link{Daim}}, \code{\link{performDaim}}, \code{\link{auc.Daim}}
}
\author{
Sergej Potapov, Werner Adler, Benjamin Hofner and Berthold Lausen

Maintainer: Sergej Potapov <sergej.potapov@gmail.com>
}
\references{
	Werner Adler and Berthold Lausen (2009).\cr
	Bootstrap Estimated True and False Positive Rates and ROC Curve.\cr
	\emph{Computational Statistics & Data Analysis}, \bold{53}, (3), 718--729.\cr
		
	Tom Fawcett (2006).\cr
	An introduction to ROC analysis.\cr
	\emph{Pattern Recognition Letters}, \bold{27}, (8).\cr

	Bradley Efron and Robert Tibshirani (1997).\cr
	Improvements on cross-validation: The.632+ bootstrap method.\cr
	\emph{Journal of the American Statistical Association}, \bold{92}, (438), 548--560.\cr

}

