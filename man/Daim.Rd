\name{Daim}
\alias{Daim}
\title{
  Diagnostic accuracy of classification models.
}
\description{
Estimation of misclassification rate, sensitivity, specificity and AUC based on 
cross-validation (CV) or various bootstrap techniques.  
}
\usage{
Daim(formula, model = NULL, data = NULL, control = Daim.control(),			
     thres = seq(0, 1, by = 0.01), cutoff = 0.5, labpos = "1", returnSample = FALSE,
     cluster = NULL, seed.cluster = NULL, multicore = FALSE, \dots)
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
  \item{cutoff}{the cutoff value for error estimation. This can be a numeric value or a character string:\cr
	\code{"cv"} - the optimal cut-point corresponding to \code{cv} estimation of the sensitivity and the specificity.\cr
	\code{"loob"} - the optimal cut-point corresponding to \code{loob} estimation of the sensitivity and the specificity.\cr
	\code{"0.632"} - the optimal cut-point corresponding to \code{0.632} estimation of the sensitivity and the specificity.\cr 
	\code{"0.632+"} - the optimal cut-point corresponding to \code{0.632+} estimation of the sensitivity and the specificity.}
  \item{labpos}{a character string of the response variable that defines a "positive" event. 
	The labels of the "positive" events will be set to "pos" and others to "neg".}
  \item{returnSample}{a logical value for saving the data from each sample.}
  \item{cluster}{the name of the cluster, if parallel computing is used.}
  \item{seed.cluster}{an integer value used as seed for the RNG.}
  \item{multicore}{a logical indicating whether multiple cores (if available) should be used for the computations.}
  \item{\dots}{additional parameters passed to \code{\link[parallel]{clusterApplyLB}} or \code{\link[parallel]{mclapply}}.}
}

\value{An object of class \code{\link{Daim-class}}.}

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

#############################
##      Evaluation of      ##
##           LDA           ##
#############################

library(TH.data)
library(MASS)
data(GlaucomaM)
head(GlaucomaM)

mylda <- function(formula, train, test){
  	model <- lda(formula, train)
  	predict(model, test)$posterior[,"pos"]
}
  
set.seed(1102013)
ACC <- Daim(Class~., model=mylda, data=GlaucomaM, labpos="glaucoma", 
            control=Daim.control(method="boot", number=50))
ACC
summary(ACC)

  
\dontrun{  
## just because of checking time on CRAN
  
  
  ####
  #### optimal cut point determination
  ####
  
  
  set.seed(1102013)
  ACC <- Daim(Class~., model=mylda, data=GlaucomaM, labpos="glaucoma", 
              control=Daim.control(method="boot", number=50), cutoff="0.632+")
  ACC
  summary(ACC)
  
  
  
  ####
  #### for parallel execution on multicore CPUs and computer clusters
  ####
  
  library(parallel)
  ### 
  ### create cluster with two slave nodes

  cl <- makeCluster(2)

  ###
  ### Load used package on all slaves and execute Daim in parallel
  ###

  clusterEvalQ(cl, library(ipred))
  ACC <- Daim(Class~., model=mylda, data=GlaucomaM, labpos="glaucoma", cluster=cl)
  ACC


  ####
  #### for parallel computing on multicore CPUs
  ####

  ACC <- Daim(Class~., model=mylda, data=GlaucomaM, labpos="glaucoma", multicore=TRUE)
  ACC
  
  
  
  
  
  #############################
  ##      Evaluation of      ##
  ##      randomForrest      ##
  #############################
  
  
  library(randomForest)

  myRF <- function(formula, train, test){
      model <- randomForest(formula, train)
  	  predict(model,test,type="prob")[,"pos"]
  }

  ACC2 <- Daim(Class~., model=myRF, data=GlaucomaM, labpos="glaucoma",
               control=Daim.control(number=50))
  ACC2
  summary(ACC2)
  
  
  ####
  #### optimal cut point determination
  ####
  
  
  set.seed(1102013)
  ACC2 <- Daim(Class~., model=myRF, data=GlaucomaM, labpos="glaucoma", 
              control=Daim.control(method="boot", number=50), cutoff="0.632+")
  summary(ACC2)
  
  
  
  ####
  #### for parallel execution on multicore CPUs and computer clusters
  ####
  
  
  library(parallel)
  ### 
  ### create cluster with two slave nodes

  cl <- makeCluster(2)

  ###
  ### Load used package on all slaves and execute Daim in parallel
  ###

  clusterEvalQ(cl, library(randomForest))
  ACC2 <- Daim(Class~., model=myRF, data=GlaucomaM, labpos="glaucoma", cluster=cl)
  ACC2

  ####
  #### for parallel computing on multicore CPUs
  ####

  ACC2 <- Daim(Class~., model=myRF, data=GlaucomaM, labpos="glaucoma", multicore=TRUE)
  ACC2
  }
}
\keyword{classif} \keyword{models} \keyword{programming}
