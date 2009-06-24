\name{roc.area}
\alias{roc.area}
\alias{roc.area.Daim}
\alias{roc.area.Daim.list}
\alias{roc.area.Daim.vector}
\title{Plot the area under the ROC curve}
\description{
  This function plots the ROC curve and fills the area under this curve.
}
\usage{
\method{roc.area}{Daim}(x, method=NULL, col="red", area.color=rgb(1,0,0,alpha=0.5), 
	xlab="False positive rate", ylab="True positive rate",
	density=NULL, angle=4, border=NULL, add=FALSE, \dots)
\method{roc.area}{Daim.list}(x, col="black", area.color=rgb(1,0,0,alpha=0.5),
	xlab="False positive rate", ylab="True positive rate", 
        main="ROC curves", density=NULL, angle=4, border=NULL, 
	add=FALSE, \dots)
\method{roc.area}{Daim.vector}(x, col="red", area.color=rgb(1,0,0,alpha=0.5),
	xlab="False positive rate", ylab="True positive rate", 
        main="ROC curve", density=NULL, angle=4, border=NULL, 
	add=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{Daim}, \code{Daim.list} or \code{Daim.vector}.}
  \item{method}{kind of the estimation of the ROC curve: '.632+', '.632', 'loob', 'cv', 'sample'}
  \item{col}{the color used to draw the ROC curve.}
  \item{area.color}{the color for filling the area.}
  \item{xlab}{a title for the x axis: see 'title'.}
  \item{ylab}{a title for the y axis: see 'title'.}
  \item{main}{a main title for the plot, see also \code{\link{title}}.}
  \item{density}{the density of shading lines, in lines per inch.  The default
          value of 'NULL' means that no shading lines are drawn. A zero
          value of 'density' means no shading nor filling whereas
          negative values (and 'NA') suppress shading (and so allow
          color filling).}
  \item{angle}{the slope of shading lines, given as an angle in degrees
          (counter-clockwise).}
  \item{border}{the color to draw the border.  The default, 'NULL', means to
          use 'par("fg")'.  Use 'border = NA' to omit borders.}
  \item{add}{logical specifying if roc-area should be added to an already
          existing plot.}
  \item{...}{graphical parameters can be given as arguments to 'plot'.}
}

\seealso{
  \code{\link{plot.Daim}}, \code{\link{auc.Daim}}
}
\examples{
  data(Daim.data1)
  data(Daim.data2)

  perform1 <- performDaim(Daim.data1$prob.oob, Daim.data1$labels, 
	Daim.data1$prob.app)
  perform2 <- performDaim(Daim.data2$prob.oob, Daim.data2$labels, 
	Daim.data2$prob.app)
  
  summary(perform1)
  summary(perform2)

  roc.area(perform2)
  roc.area(perform1, area.color=rgb(0,0,1,alpha=0.2), col="blue", add=TRUE)
  legend(0.7,0.2,c("Model-1","Model-2"), col=c("red","blue"),lty=1,bg="white")

  ####
  #### If your device don't support the semi-transparent colors use 
  #### the PDF Graphics Device.
  ####

  pdf("ROC-area.pdf",version=1.4)
  roc.area(perform2)
  roc.area(perform1, area.color=rgb(0,0,1,alpha=0.2), col="blue", add=TRUE)
  legend(0.7,0.2,c("Model 1","Model 2"), col=c("red","blue"),lty=1,bg="white")
  dev.off()

}
\keyword{hplot} \keyword{classif}
