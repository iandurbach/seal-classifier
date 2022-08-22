##########
#
#		FUNCTIONS FOR TRANSITION MATRICES - MARKOV CHAIN - NETWORK GRAPH
#
#		Includes source code from package "markovchain"




padMatrix <- function(x, padrow=c(1,1), padcol=c(1,1), value=NA){
# Pad matrix by adding rows and columns
#	x: original matrix
#	arow: number of row to add above and below, respectively
#	acol: number of columns to add left and right, respectively
#	value: value used to fill in cells
	if (!is.matrix(x)) stop("'x' must be a matrix")
	if (padrow[1]>0) x <- rbind(matrix(value,nrow=padrow[1],ncol=ncol(x)),x)
	if (padrow[2]>0) x <- rbind(x,matrix(value,nrow=padrow[2],ncol=ncol(x)))
	if (padcol[1]>0) x <- cbind(matrix(value,nrow=nrow(x),ncol=padcol[1]),x)
	if (padcol[2]>0) x <- cbind(x,matrix(value,nrow=nrow(x),ncol=padcol[2]))
	return(x)
}


prettyMatrix <- function( x, empty=0, aspect='horiz', ... ){
# Make a matrix with automatic row & column numbers to make it pretty
#	empty: devault value
#	aspect: aspect ratio
	if (aspect=='horiz') {
		m <- floor(sqrt(length(x)))	# nrow
		n <- ceiling(length(x)/m)	# ncol
	} else {
		n <- floor(sqrt(length(x)))	# ncol		
		m <- ceiling(length(x)/n)	# nrow
	}
	if(length(x)<m*n) { x <- c(x, rep(empty,m*n-length(x))) }
	return(matrix(x,nrow=m,ncol=n,...))
}


expand.tmat <- function(x, X, fill=0){
#	Function to expand matrix 'x' to 'X'
	#	x: original matrix
	#	X: larger matrix
	#	fill: value to fill the large matrix. 0:ML (default), >0:Laplace, non-numeric: leave as is
	if (is.numeric(fill)) X[1:length(X)] <- fill
	for (i in 1:nrow(x)){
		X[rownames(x)[i],colnames(x)] <- X[rownames(x)[i],colnames(x)] + x[i,]
	}
	return(X)
}


plot.tmat <- function(X, grid=TRUE, col.grid=grey(0,alpha=0.1), ztransf, pt.cex=1.5,
					col=hsv(alpha=seq(0,1,0.1)), xlab, ylab, leg=NULL, inset=c(-0.1,0), ...){
#	Function to plot transition matrix
	require(scales)
	sleg <- leg
	if (!missing(ztransf)){
		X <- get(ztransf)(X)
		if(!is.null(sleg)) sleg <- get(ztransf)(leg)
	}
	par(xpd=FALSE)
	image( x=rescale(1:nrow(X)), y=rescale(1:nrow(X)), z=(t(X[nrow(X):1,])), 
				col=col, breaks=seq(0,1,length.out=length(col)+1), xaxt='n', yaxt='n', xlab="", ylab="", ...)
	if (grid){
		abline(h=rescale(1:nrow(X))+(diff(rescale(1:nrow(X)))/2)[1], col=col.grid)
		abline(v=rescale(1:nrow(X))+(diff(rescale(1:nrow(X)))/2)[1], col=col.grid)
	}
	axis(3, at=rescale(1:nrow(X)), labels=rownames(X), las=2, ...)
	axis(2, at=rescale(1:nrow(X)), labels=rev(rownames(X)), las=2, ...)
	if (!missing(xlab)) mtext(xlab, side=3, cex=2, padj=-3)
	if (!missing(ylab)) mtext(ylab, side=2, cex=2, padj=-3)
	if (!is.null(leg)){
		par(xpd=TRUE)
		legend( "topright", inset=inset, legend=leg, y.intersp=1.2, xjust=0,
					 pt.cex=pt.cex, pch=22, col=col.grid, pt.bg=hsv(alpha=sleg), title="Probability")
	} 	
}

#			Function to add a transition step to an existing transition matrix
#			(from haul out only to haul out + at sea)
addTransit <- function(Tm, Ds, Dt, prefx="at-sea ", check.prob=6){
	#	Tm: transition matrix
	#	Ds: vector of expected duration in each state (optional, default = 1)
	#	Dt: vector of expected duration in each transition phase (optional, default = 1)
	#	check.prob: checks that probabilities in Tm sum to one (with specified number of digits if > 0)
	#						no checks is performed if </= 0
	if ( missing(Tm) ) stop("'Tm' is missing")
	if ( check.prob>0 ){
			if ( !all((1-rowSums(Tm))<(1/10^check.prob)) ) stop("probabilities in 'Tm' do not sum to one")
	}
	n <- nrow(Tm)
	ho <- colnames(Tm)
	Tmod <- matrix(	0, nrow=n*2, ncol=n*2, dimnames=list(c(ho, paste(prefx,ho,sep="")),
					c(ho,paste(prefx,ho,sep=""))) )
	if ( missing(Ds) ) Ds <- rep(1, n)
	if ( missing(Dt) ) Dt <- rep(1, n)
	#		Probability of staying at haul out (continue haulout) <-- Ds
	diag(Tmod)[1:n] <- exp( log(0.5)/Ds )
	#		Probability of leaving haul out	(start trip)	<-- Ds
	diag(Tmod[1:n,(n+1):(n*2)]) <- 1 - exp( log(0.5)/Ds )
	#		Probability of staying at sea (continue trip) <-- Dt
	diag(Tmod[(n+1):(n*2),(n+1):(n*2)]) <- exp( log(0.5)/Dt )
	#		Probability of hauling out (end trip) <-- Dt
	Tmod[(n+1):(n*2),1:n] <- apply(Tm, 2, function(x) x*( 1 - exp( log(0.5)/Dt ) ))
	return(Tmod)
}







############		Following from package markovchain			##############


markovchainSequence <- function (n, markovchain, t0 = sample(markovchain@states, 1), 
														include.t0 = FALSE){
#	Sampler for univariate markov chains
	if (!(t0 %in% markovchain@states)) stop("Error! Initial state not defined")
	chain <- rep(NA,n)# CHANGED
	state <- t0
	for (i in 1:n) {
		rowProbs <- markovchain@transitionMatrix[which(markovchain@states == state), ]
		outstate <- sample(size = 1, x = markovchain@states, prob = rowProbs)
		chain[i] <- outstate #CHANGED
		state <- outstate
		}
	if (include.t0){
		out <- c(t0, chain)
	} else {
		out <- chain
	}
	return(out)
}



checkSequence<-function(object){
#	Check whether the subsequent states are included in the previous ones
	out<-TRUE
	if (dim(object)==1) return(out) #if the size of the list is one do
	for (i in 2:dim(object)){
		statesNm1<-states(object[[i-1]]) #evalutate mc n.1
		statesN<-states(object[[i]]) #evaluate mc n
		intersection<-intersect(statesNm1,statesN) #check the ibntersection
		if (setequal(intersection, statesNm1)==FALSE) { #the states at n-1 
			out<-FALSE
			break
		}
	}
	return(out)
}



rmarkovchain <- function(n,object,...){
#	Function to perform random sampling
	if (class(object)=="markovchain") out<-markovchainSequence(n=n, markovchain=object,...)
	if (class(object)=="markovchainList"){
		verify <-.checkSequence(object=object)
		if (!verify) warning("Warning: some states in the markovchain sequences are not contained in the following states!")
		iteration<-numeric()
		values<-character()
		for (i in 1:n){ #number of replicates
		#the first iteration may include initial state
			sampledValues<-markovchainSequence(n=1,markovchain=object[[1]],...)
			outIter<-rep(i,length(sampledValues))
			if (dim(object)>1){
				for (j in 2:dim(object)){
					pos2take<-length(sampledValues)
					newVals<-markovchainSequence(n=1,markovchain=object[[j]],t0=sampledValues[pos2take]) #the initial state is in the ending position of the mclist
					outIter<-c(outIter,i)
					sampledValues<-c(sampledValues,newVals)
				}
			}
			iteration<-c(iteration, outIter)
			values<-c(values, sampledValues)
		}
		out<-data.frame(iteration=iteration, values=values)
	}
	return(out)
}






