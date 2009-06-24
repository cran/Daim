plot.Daim <- function(x, method=NULL, all.roc=FALSE, color="red", alpha=0.25,
				type="b", xlab="False positive rate",ylab="True positive rate", 
				main=NULL, add=FALSE, ...)
{		
	if(class(x)[2] != "cv"){
		if(is.null(method))
			method <- ".632+"
		meth <- charmatch(method,c(".632+",".632","loob","sample"),nomatch = 0)
		if(meth == 0)
			stop("\nThe value of 'method' must be one of '.632+', '.632', 'loob' or 'sample'\n")
		if(is.null(main) && !add)
			main <- paste("Method -",method)
		if(meth == 4)
			all.roc <- TRUE
		if(!add){
			plot(-1,-1, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab,
				main=main)
		}
		if(all.roc || meth == 4){
			for(i in 1:length(x$sample.roc)){
				lines(x$sample.roc[[i]],col=rgb(0,1,0,255*alpha,max=255), ...)
				}
		}
		if(meth == 1)
			lines(1-x$roc$spec632p,x$roc$sens632p,col=color,type=type, ...)
		else if(meth == 2)
			lines(1-x$roc$spec632,x$roc$sens632,col=color,type=type, ...)
		else if(meth == 3)
			lines(1-x$roc$specloob,x$roc$sensloob,col=color,type=type, ...)
		grid()
	}
	else{
		if(is.null(method))
			method <- "cv"
		meth <- charmatch(method,c("cv","sample"),nomatch = 0)
		if(meth == 0)
			stop("\nThe value of 'method' must be one of 'cv' or 'sample'\n")
		if(is.null(main) && !add)
			main <- paste("Method -",method)
		if(meth == 2)
			all.roc <- TRUE
		if(!add){
			plot(-1,-1, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab,
				main=main)
		}
		if(all.roc || meth == 2){
			for(i in 1:length(x$sample.roc)){
				lines(x$sample.roc[[i]],col=rgb(0,1,0,255*alpha,max=255), ...)
				}
		}
		if(meth == 1)
			lines(1-x$roc$specloob,x$roc$sensloob,col=color,type=type, ...)
		grid()
	}
}

