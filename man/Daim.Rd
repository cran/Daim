\name{Daim}
\alias{Daim}
\title{
  Diagnostic accuracy of classification models.
}
\description{
Estimation of prediction error based on cross-validation (CV) or various bootstrap techniques.  
}
\usage{
Daim(formula, model=NULL, data=NULL, control = Daim.control(),			
		thres = seq(0,1,by=0.01), cutoff = 0.5, 
		labpos = "1", returnSample = FALSE,
		cluster = NULL, seed.cluster = NULL, \dots)
}
\arguments{
  \item{formula}{formula of the form \code{y ~ x1 + x2 + \ldots},
    where \code{y} must be a factor and \code{x1,x2,\ldots} are numeric or factor.}
  \item{model}{function. Modelling technique whose error rate is to be estimated. 
	The function \code{model} returns the predicted probability for each observation.}
  \item{data}{an optional data frame containing the variables
    	in the model (training data).}
  \item{control}{See \code{\link{Daim.control}}.}
  \item{thres}{a numeric vector with the cutoff values.}
  \item{cutoff}{the cutoff value for error estimation. This can be a numeric value or a character string.
	If the \code{cutoff} set to:\cr
	\code{".632"} - the estimated cut-point corresponding the \code{.632} estimation of the sensitivity and the specificity.\cr 
	\code{".632+"} - the estimated cut-point corresponding the \code{.632+} estimation of the sensitivity and the specificity.}
  \item{labpos}{a character string of the response variable that defines a "positive" event. 
	The labels of the "positive" events will be set to "pos" and other to "neg".}
  \item{returnSample}{a logical value for saving the data from each sample.}
  \item{cluster}{the name of the cluster, if parallel computing will be used.}
  \item{seed.cluster}{an integer value used as seed for the RNG.}
  \item{\dots}{additional parameters.}
}
\value{a list with the following components :
  \item{call}{the matched call.}
  \item{formula}{the formula supplied.}
  \item{method}{the list of control parameters.}
  \item{err632p}{the \code{.632+} estimation of the misclassification error.}
  \item{err632}{the \code{.632} estimation of the misclassification error.}
  \item{errloob}{the \code{LOOB} estimation of the misclassification error.}
  \item{errapp}{the apparent error.}
  \item{sens632p}{the \code{.632+} estimation of the sensitivity.}
  \item{spec632p}{the \code{.632+} estimation of the specificity.}
  \item{sens632}{the \code{.632} estimation of the sensitivity.}
  \item{spec632}{the \code{.632} estimation of the specificity.}
  \item{sensloob}{the \code{LOOB} estimation of the sensitivity.}
  \item{specloob}{the \code{LOOB} estimation of the specificity.}
  \item{sensapp}{the apparent sensitivity.}
  \item{specapp}{the apparent specificity.}
  \item{roc}{a data frame with estimated values of sensitivity and specificity for a variety of cutoffs.}
  \item{sample.roc}{a list in which each entry contains the values of the ROC curve of this special sample or cross-validation run.}
  \item{sample.data}{a data frame with the results of this particular sample or cross-validation run.}
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
\seealso{
  \code{\link{plot.Daim}}, \code{\link{performDaim}}, \code{\link{auc.Daim}}, \code{\link{roc.area.Daim}}
}
\examples{

  library(ipred)
  data(GlaucomaM)
  head(GlaucomaM)

  mylda <- function(formula,train,test){
  	model <- lda(formula,train)
  	predict(model,test)$posterior[,"pos"]
  }

  ACC <- Daim(Class~.,model=mylda,data=GlaucomaM,labpos="glaucoma")
  ACC
  summary(ACC)


  ####
  #### for parallel computing with snow cluster
  ####
 
  # library(snow)
  ### 
  ### create cluster with two slave nodes

  # cl <- makeCluster(2)

  ###
  ### Load used library on all slaves and execute the Daim in parallel
  ###

  # clusterEvalQ(cl, library(ipred))
  # ACC <- Daim(Class~.,model=mylda,data=GlaucomaM,labpos="glaucoma",cluster=cl)
  # ACC


  ####
  #### for parallel computing with multicore package 
  #### you need only to load this library
  ####

  # library(multicore)
  # ACC <- Daim(Class~.,model=mylda,data=GlaucomaM,labpos="glaucoma")
  # ACC


  library(randomForest)

  myRF <- function(formula,train,test){
  	model <- randomForest(formula,train)
  	predict(model,test,type="prob")[,"pos"]
  }

  ACC2 <- Daim(Class~.,model=myRF,data=GlaucomaM,labpos="glaucoma",
	control=Daim.control(number=25))
  ACC2
  summary(ACC2)


  ####
  #### for parallel computing with snow cluster
  ####

  # library(snow)
  ### 
  ### create cluster with two slave nodes

  # cl <- makeCluster(2)

  ###
  ### Load used library on all slaves and execute the Daim in parallel
  ###

  # clusterEvalQ(cl, library(randomForest))
  # ACC2 <- Daim(Class~.,model=myRF,data=GlaucomaM,labpos="glaucoma",cluster=cl)
  # ACC2

  ####
  #### for parallel computing with multicore package
  ####

  # library(multicore)
  # ACC2 <- Daim(Class~.,model=myRF,data=GlaucomaM,labpos="glaucoma")
  # ACC2

}
\keyword{classif} \keyword{models} \keyword{programming}
