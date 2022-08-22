####################################################
#
#	Collection of functions for:
#	- Calculations
#	- Data manipulations
#	Use: source('~/Documents/R/funBasics.R', chdir = TRUE)
#	PC?: source.with.encoding('E:/R Code Archive/funBasics.R', encoding='UTF-8')
#



###############################
#
#	CALCULATIONS
#

ratio <- function( x1=NULL, x2=NULL, crossed=FALSE ){
  # Ratio of two vectors
  #	x1: numerator
  #	x2: denominator
  #	crossed: whether to calculate all possible combinations of x1 & x2 values: returns matrix x1 * x2
  #	- x1 and x2 must be multiples of each other; the shorter vector is recycled
  #	- if x2 is NULL, returns the ratio of values in x1
  if (is.null(x2)){
    out <- sapply(x1,function(x){x/x1})
  } else {
    if (crossed){
      out <- sapply(x1, function(x){x/x2})
    } else {
      out <- x1/x2
    }
  }
  return(out)
}


meanM <- function( x, w=NULL, na.rm=FALSE ){
  # Mean of matrices
  #	x: list of matrices with equal dimensions	
  #	w: vector of weights (default: all equal)
  if (is.null(w)){ w <- rep(1/length(x),x) }
  if (sum(w, na.rm=na.rm)!=1){ w <- w/sum(w, na.rm=na.rm) }
  x.mean <- x[[1]]*w[1]
  for (i in 2:length(x)){
    x.mean <- x.mean+ x[[i]]*w[i]
  }
  return(x.mean)
}


se <- function( x, na.rm=FALSE ){
  # Standard error of a sample
  if(na.rm) x <- na.exclude(x)
  return( sqrt( var(x, na.rm=na.rm) / length(x) ) )
}


ci <- function( x, p=0.95, tail=2, limit='both', na.rm=FALSE ){
  # Confidence interval of a sample mean based on the Normal distribution
  # (can add different distributions later)
  p <- 1 - (1-p)/tail; xmean <- mean(x, na.rm=na.rm); xsd <- sd(x, na.rm=na.rm)
  q1 <- qnorm( p=p, mean=xmean, sd=xsd, lower.tail=FALSE )
  q2 <- qnorm( p=p, mean=xmean, sd=xsd, lower.tail=TRUE )
  if (substr(tolower(limit),1,3)=='low') {
    return( q1 )
  } else if (substr(tolower(limit),1,2)=='up') {
    return( q2 )
  } else {
    return( c(q1,q2) )
  }
}


smode <- function( x, na.rm=FALSE ) {
  # Mode of a sample x
  if(na.rm) x <- na.exclude(x)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


hmean <- function( x, na.rm=FALSE ){
  # Harmonic mean
  #	often for rates and ratios
  #	means: Harmonic < Geometric < Arithmetic
  if(any(x==0)) stop('Zero values are not tolerated')
  return( 1/(mean(1/x,na.rm=na.rm)) )
}


gmean <- function(x, na.rm=FALSE){
  # Geometric mean
  #	often for exponential / lognormal data (positive numbers only)
  #	means: Harmonic < Geometric < Arithmetic
  if(any(x<=0)) stop('All values must be positive')
  exp( mean(log(x),na.rm=na.rm) )
}


gsd <- function( x, na.rm=FALSE ){
  #	Geometric SD
  if(any(x<=0)) stop('All values must be positive')
  return( exp( sqrt(var(x, na.rm=na.rm)) ) )
}


rmean <- function( x, intv=5, wts=NULL,... ){
  # Running mean / moving average
  #	calculates a running mean through entire series
  #	assumes data is in correct order
  #	see filter() and filtfilt() for something better
  #	x: vector of values
  #	intv: interval size for means
  #	wts: weights for each value (optional)
  y <- vector('numeric',length=length(x)-(intv-1))
  for (i in 1:length(y)){
    index <- i:(i+(intv-1))
    if (is.null(wts)) {
      y[i] <- mean( x[index], ... )
    } else { 
      y[i] <- weighted.mean( x[index], w=wts[index], ... )
    }
  }
  return(y)
}


cumave <- function( x ){
  #	Cumulative average
  sapply( 1:length(x), function(y)mean(x[1:y]) )
}


prod.logit <- function(p, vcv, alpha=0.05, ...){
  # Function to calculate mean, SE and confidence intervals for the product
  #	of random variables on a logit scale using the delta method
  #	p: vector of probabilities to multiply
  #	vcv: variance covariance matrix
  #	alpha: for confidence interval
  require(msm)
  n <- length(p)
  ylogit.se <- deltamethod(g=~log((paste("x",1:n,collapse="*", sep=""))/(1-(paste("x",1:n,collapse="*", sep="")))), 
                           mean=p, cov=vcv, ses=TRUE)
  ylogit.hat <- prod(qlogis(p))
  #   Calculate CI in logit scale and backtransform
  ylogit.uci <- ylogit.hat + qnorm(p=1-alpha/2)*ylogit.se
  ylogit.lci <- ylogit.hat - qnorm(p=1-alpha/2)*ylogit.se
  #   Estimates of mean, SE, lowerCI, upperCI
  p.hat <- plogis( c(estim=ylogit.hat, se=ylogit.se, lci=ylogit.lci, uci=ylogit.uci) )
  return(p.hat)
}


mean.inv <- function(x, vcv, alpha=0.05, ...){
  # Function to calculate mean, SE and confidence intervals for the mean
  #	of random variables on a inverse scale using the delta method
  #	Note: at the moment assumes a logit transformation
  #	p: vector of probabilities to multiply
  #	vcv: variance covariance matrix
  #	alpha: for confidence interval
  require(msm)
  n <- length(x)
  yinv.se <- deltamethod( 	g=formula(paste("~ ",n,"/(",paste("x",1:n,collapse="+", sep=""),")")), 
                           mean=x, cov=v, ses=TRUE )
  yinv.hat <- mean(1/x)
  #   Calculate CI in logit scale and backtransform
  yinv.uci <- yinv.hat + qnorm(p=1-alpha/2)*yinv.se
  yinv.lci <- yinv.hat - qnorm(p=1-alpha/2)*yinv.se
  #   Estimates of mean, SE, lowerCI, upperCI
  x.hat <- c( estim=yinv.hat, se=yinv.se, lci=yinv.lci, uci=yinv.uci )
  return(x.hat)
}


# #	See also package 'scales':::rescale() instead
# # Scale x to [0:1] using linear interpolation
# scale01 <- function ( x=NULL, xmin=NULL, xmax=NULL, inv=FALSE ){
# #	inv: logical to indicate reverse transformation (default=FALSE)
# #	Check parameters
# if ( is.null(x) ) { stop('x is NULL') }
# if ( is.null(xmin) ) { xmin <- min(x) }
# if ( is.null(xmax) ) { xmax <- max(x) }
# if (xmin>=xmax) { stop(paste('Invalid interval:',xmin,'-',xmax)) }
# #	Return scaled value
# if (inv) {
# if ( (x<0)||(x>1) ) {stop('\'x\' must be between zero & one')}
# return( x*(xmax-xmin) + xmin )
# } else {
# if ( (x<xmin)||(x>xmax) ) {
# print(paste('x =',x,'; xmin =',xmin,'; xmax =',xmax))
# stop(paste('\'x\' is not within [',xmin,':',xmax,']'))
# }
# return( (x-xmin)/(xmax-xmin) )
# }
# }


dsinexp <- function( x, a, b ){
  # Generate densities for a double bounded continuous distribution
  #	(Kumaraswamy,1976; see also Beta distribution)
  a*b*x^(a-1) * (1-x^a)^(b-1)
}
# Test
# x <- seq(0,1,0.01)
# plot(NULL, xlim=0:1, ylim=c(0,6), ylab='density')
# for(i in 2:11){points(x,dsinexp(x,i,2),col=rev(heat.colors(10))[i], type='l')}
# for(i in 2:11){points(x,dsinexp(x,2,pi*i),col=rev(heat.colors(10))[i], type='l')}

cumunique <- function(x,count=TRUE){
  #	Function to find (number of) unique values in a sequence
  #	x: sequence of values
  #	count: if TRUE, returns number of unique values (vector), else returns list of unique values
  index <- cbind(rep(1,length(x)),seq(1,length(x)))
  if (count){
    apply(index,1, function(y)length(unique(x[y[1]:y[2]])))
  } else {
    apply(index,1, function(y)unique(x[y[1]:y[2]]))
  }
}

cutQ <- function( x,n=2,type=8,dig=2, na.rm=FALSE, ... ){
  # Categorize variables based on quantiles
  #	x: numeric vector
  #	n: number of categories
  #	type:		1-3: discontinuous variables
  #				4-9: continuous variables
  #				'?quantile' for descriptions of types
  #	dig: number of digits for labels only
  #	Verify that data is numeric
  if(!all(is.numeric(x))) stop('all values of \'x\' must be numeric')
  #	Break points
  brk <- quantile(x,probs=seq(0,1,1/n),type=type)
  #	Output
  int <- cut(x, brk, include.lowest=TRUE, labels=FALSE)
  ltext <- gsub(' ', '', paste('[', round(brk[-length(brk)], digits=dig), '-',
                               round(brk[-1], digits=dig), '['))
  substring(ltext[length(ltext)],nchar(ltext[length(ltext)]))<-']'
  return( cbind(Categ=int, Labels=ltext[int] ) )
} # end function cutQ


which.maxa <- function( x, n=1, ties=FALSE ){
  # Return location of n maxima
  #	n: number of maxima to find
  #	ties: whether to return all locations equal to maxima
  #			- default: FALSE and returns only the first instance as which.max()
  if (!is.numeric(x)) {
    warning('coercing x to numeric')
    x <- as.numeric(x) }
  index <- order(x, 1:length(x), decreasing=TRUE)[1:n]
  if (ties){ 
    index <- sapply(index, FUN=function(y){which(x[y]==x)})
    index <- unique( as.vector(index) ) }
  return( index )
}
which.mina <- function( x, n=1, ties=FALSE ){
  # Return location of n minima
  #	n: number of minima to find
  #	ties: whether to return all locations equal to minima
  #			- default: FALSE and returns only the first instance as which.max() does
  if (!is.numeric(x)) {
    warning('coercing x to numeric')
    x <- as.numeric(x) }
  index <- order(x, 1:length(x), decreasing=FALSE)[1:n]
  if (ties){ 
    index <- sapply(index, FUN=function(y){which(x[y]==x)})
    index <- unique( as.vector(index) ) }
  return( index )
}


maxa <- function( x=NULL, n=1, ties=FALSE ){
  # Return n maxima
  return( x[which.maxa(x,n,ties)] )
}
mina <- function( x=NULL, n=1, ties=FALSE ){
  # Return n minima
  return( x[which.mina(x,n,ties)] )
}


which.ininterv <- function(pts, intervals,closed=c(TRUE,TRUE)){
  # Find indices of values inside an interval
  #	returns a list of indices
  #	pts: vector of values to assess
  #	intervals: table of intervals (columns: start & end)
  #	closed: should left and right intervals bounds be closed (default: TRUE for both)
  if (!is.matrix(intervals)) intervals <- as.matrix(intervals)
  if (ncol(intervals)!=2) stop("'intervals' must have 2 columns" )
  if (closed[1]){
    left <- sapply(intervals[,1], function(x){pts>=x})
  } else {
    left <- sapply(intervals[,1], function(x){pts>x})
  }
  if (closed[2]){
    right <- sapply(intervals[,2], function(x){pts<=x})
  } else {
    right <- sapply(intervals[,2], function(x){pts<x})
  }
  out <- apply(left&right, 2, which); names(out) <- rownames(intervals)
  return(out)
}


which.huginterv <- function(pts, intervals,closed=c(TRUE,TRUE)){
  # Find indices of values immediately preceding and immediately following an interval
  #	returns a matrix of indices (columns)
  #	pts: vector of values to assess
  #	intervals: table of intervals (columns: start & end)
  #	closed: should left and right intervals bounds be closed (default: TRUE for both)
  if (!is.matrix(intervals)) intervals <- as.matrix(intervals)
  if (ncol(intervals)!=2) stop("'intervals' must have 2 columns" )
  if (closed[1]){
    left <- sapply(intervals[,1], function(x){pts<x})
  } else {
    left <- sapply(intervals[,1], function(x){pts<=x})
  }
  if (closed[2]){
    right <- sapply(intervals[,2], function(x){pts>x})
  } else {
    right <- sapply(intervals[,2], function(x){pts>=x})
  }
  
  left <- apply(left, 2, function(x)if(any(x)){sum(x)}else{NA})
  right <- apply(right, 2, function(x)if(any(x)){min(which(x))}else{NA})
  out <- cbind(left,right); names(out) <- rownames(intervals)
  return(out)
}



which.hug <- function(pts, table){
  # Find indices of values immediately lower & higher than x
  #	returns a matrix of dimension (2, length(x))
  #	pts: numeric vector of values	 to match in table
  #	table: numeric vector of values to search * assumed ordered
  x.diff <- sapply(pts, function(x){x-table})
  index.hug <- apply(x.diff, 2, function(x){		# get the day before & day after
    if (all(x<0)){
      rep(which.max(x),2)
    } else if (all(x>0)) {
      rep(which.min(x),2)
    } else {
      c( max(which(x>=0)), min(which(x<=0)) )
    }
  })
  rownames(index.hug) <- c("below","above"); colnames(index.hug) <- pts
  return(index.hug)
}

hug <- function(pts, table){
  # Find values immediately lower and higher than x
  #	returns a matrix of dimension (2, length(x))
  #	pts: numeric vector of values	 to match in table
  #	table: numeric vector of values to search * assumed ordered
  matrix( table[which.hug(pts,table)], nrow=2, dimnames=list(c("below","above"),pts) )
}


peaks <- function( series, span=3, do.pad=TRUE ) {
  # Find peaks in time series data
  #	(Martin Maechler, Date: 25 Nov 2005)
  if((span <- as.integer(span)) %% 2 != 1) stop('\'span\' must be odd')
  s1 <- 1:1 + (s <- span %/% 2)
  if(span == 1) return(rep.int(TRUE, length(series)))
  z <- embed(series, span)
  v <- apply(z[,s1] > z[, -s1, drop=FALSE], 1, all)
  if(do.pad) {
    pad <- rep.int(FALSE, s)
    c(pad, v, pad)
  } else {
    return(v)
  }
}


peaksign <- function( series, span=3, do.pad=TRUE ){
  # Classify time series data (-1 / 0 / 1) if series[i] is ( trough / 'normal' / peak )
  #	(Martin Maechler, Date: 25 Nov 2005)
  if((span <- as.integer(span)) %% 2 != 1 || span == 1)
    stop('\'span\' must be odd and >= 3')
  s1 <- 1:1 + (s <- span %/% 2)
  z <- embed(series, span)
  d <- z[,s1] - z[, -s1, drop=FALSE]
  ans <- rep.int(0:0, nrow(d))
  ans[apply(d > 0, 1, all)] <- as.integer(1)
  ans[apply(d < 0, 1, all)] <- as.integer(-1)
  if(do.pad) {
    pad <- rep.int(0:0, s)
    return(c(pad, ans, pad))
  } else {
    return(ans)
  }
}




###############################
#
#	DATA MANIPULATIONS
#


odd <- function(x){ x[seq(1,length(x),2)] }
even <- function(x){ x[seq(2,length(x),2)] }


d2time <- function(x, sep=":"){
  # Convert days to H:M:S
  #	x: number of days in decimal format
  x <- x*24*60*60		# convert days to sec
  s <- round(x%%60)			# get remainder sec (/min)
  x <- (x-s)/60		# convert left over sec to min
  m <- round(x%%60)		# remainder min (/hr)
  x <- (x-m)/60		# convert min to hr
  h <- round(x)
  out <- paste(h,m,s,sep=sep)
  return(out)
}

d2intv <- function(x, units="days", digits=0){
  # Convert dates to intervals (days)
  #	x: vector of dates
  x <- difftime(tail(x,-1),head(x,-1), units=units)
  x <- c(0,cumsum(as.numeric(x)))
  x <- round(x, digits=digits)
  return( Intervals( cbind(head(x,-1),tail(x,-1)) ) )
}

cmerge <- function(x, index){
  # Sums the columns to be merged 
  #	(can expand to do other operations)
  #	x: matrix
  #	index: names of columns to be merged
  index <- unlist(index)
  if(length(index)<2){ stop('index: provide more than one column names to merge') }
  if(is.character(index)){
    index <- index[nchar(index)>0]
    if (any(index%in%colnames(x))){
      index <- index[index%in%colnames(x)]
      index <- match( index, colnames(x) )
    } else {
      index <- ''
    }
  }
  if (length(index)>1){
    x[,index[1]] <- rowSums(x[,index])
    return( x[,-index[-1]] )
  } else {
    return(x)
  }
}


substrN <- function( x, n, side='right' ){
  # Wrapper for substring() to extract last characters
  #	x: character			
  #	n: number of characters
  #	side: where to start (default = right because left side is easy with substr)
  if (tolower(substr(side,1,1))=='l') {
    return( substr(x, start=1, stop=n) )
  } else {
    return( substr(x, start=nchar(x)-n+1, stop=nchar(x)) )
  }
}

trimSpace <- function (x, at=3) {
  # Trim leading and/or trailing end of a character string
  #	x: character
  #	at: where to trim; 1:leading, 2:trailing, else both
  if (at==1){
    sub("^\\s+", "", x)
  } else	if (at==2){
    sub("\\s+$", "", x)
  } else {
    gsub("^\\s+|\\s+$", "", x)
  }
}

revChar <- function( x ){
  # Reverse a character string or number
  #	x: character
  options('scipen'=100)	# avoid exponential notation
  type <- is(x)[1]
  x <- unlist(strsplit(as.character(x), split=''))
  x <- paste(rev(x), collapse='')
  if(type=='character'){
    return(x)
  } else {
    return(as.numeric(x))
  }
}



###
# Function to check numbers
# - is.na(x): set to zero
# - x < 0: set to zero
# (left out checking for characters: not likely...)
# - int: check that value is an integer (whole organism)
#		- NULL: do not check
#		- 'down': round down
#		- 'up': round up
#		- 'round'
checkN <- function( x=NULL, int=NULL, val=0 ){
  if( any(is.na(x)) ) x[which(is.na(x))] <- val
  if( any(x<0) ) x[which(x<0)] <- val
  if( !is.null(int) ){
    int <- tolower(as.character(int))
    if( substr(int,1,1)=='d' ){ x<-floor(x)
    } else if( substr(int,1,1)=='u' ){ x<-ceiling(x)
    } else if( substr(int,1,1)=='r' ){ x<-round(x)
    }
  }
  return(x)
}


offdiag <- function( X=NULL, tri='upper', off=1 ){
  # Get indices of diagonals offcenter
  #	X: square matrix
  #	tri: which side of the diagonal to get (default='upper')
  #	off: number of row/col away from center diagonal (default=1)
  index <- which(diag(nrow(X))==1)
  if (substr(tri,1,1)=='u'){
    index <- index[-(1:off)] - off
  } else {
    index <- index[1:(length(index)-off)] + off
  }
  return(X[index])
}


diagz <- function( x ){
  # Get indices of diagonals in a non-square matrix with ncol = multiples of nrow
  #	for use in structured population models with
  # 	- more than one species
  # 	- species with equal nb of stages (can set dummy stages?)
  # See diag() for square matrices
  if (!is.matrix(x)) stop("'x' must be a matrix")
  n <- 1:nrow(x) + (0:(nrow(x)-1))*nrow(x)
  out <- sapply( n, FUN=function(n) {
    n + 0:(ncol(x)/nrow(x)-1)*nrow(x)^2 } )
  return(sort(out))
}


prob2freq<- function(X, X1, distr='uniform') {
  # Function to generate a matrix of random frequencies
  #	from a matrix population model
  #	X: a vector/matrix of a structured population
  #	X1: the projection matrix
  #	distr: family distribution (default=uniform)
  N <- matrix(0, nrow=length(X), ncol=length(X) )
  X2 <- rep(X, each=length(X) )	# will fill N by columns
  distr <- substr(tolower(distr),1,3)
  if (distr=='poi') {
    for (i in 1:length(X1)) {
      if (X2[i]<0) {
        stop('Negative probabilities are not valid')
      } else {
        N[i] <- sum( rpois(n=X2[i],lambda=X1[i]) )
      }
    }
  } else if (distr=='bin') {
    for (i in 1:length(X1)) {
      if (X2[i]<=0) {
        stop('Negative probabilities are not valid')
      } else {
        N[i] <- rbinom(n=1, size=X2[i], prob=X1[i] )
      }
    }
  } else if (distr=='uni') {
    for (i in 1:length(X1)) {
      if (X2[i]<=0) {
        stop('Negative probabilities are not valid')
      } else {
        N[i] <- sum( runif(n=X2[i]) < X1[i] )
      }
    }
  } else {
    stop('Invalid distribution')
  }
  return(N)
}


rate2freq <- function(x=NULL, rate=NULL, out='total', distr='Poisson'){
  # Convert rates to (randomly generated) frequencies
  #	- x : sampling unit (vector or matrix)
  #	- rate : rate of occurance (vector or matrix)
  #	- out : type of output, total freq or list of individual draws
  #	- distr : 'Poisson', the only one for now
  #	Notes : 
  #	- assumes x & rate are in correct order
  #	- if rate is shorter than x, values are recycled with warning
  out <- substr(tolower(out),1,3)
  if (length(x)>length(rate)) {
    rate <- rep( rate, ceiling(length(x)/length(rate)) )
    warning('The variable \'rate\' was shorter than \'x\' so its values were recycled.')
  }
  x.freq <- x
  for (i in 1:length(x)) {
    if (out=='tot') {
      x.freq[i] <- rpois(1, x[i]*rate[i])
    } else if (out=='ind') {
      x.freq[i] <- list(rpois(x[i], rate[i]))
    } else { stop('invalid \'out\' argument') }
  }
  return(x.freq)
}


# Convert probabilities <-> rates
# r: rate (per one time unit)
# p: probability of occurrence in time t
# Note: either p or t can be a vector of length > 1
#		if both have lengths > 1, the longer vector
#		must be a multiple of the shorter one. The
#		shorter will be recycled.
prob2rate <- function(p=NULL,t=1){
  # convert a probability to a rate
  return( -log(1-p)/t )
}
rate2prob <- function(r=NULL,t=1){
  # convert a rate to a probability
  return( 1-exp(-r*t) )	
}


calcBeta.ab <- function(x, q=exp(1)){
  # Compute alpha & beta for the Beta distribution using:
  #	x: mode of distribution
  #	q: shape parameter
  if (x >= 0.5){
    a <- q
    b <- (a-1)/x -a +2	# shape parameter beta = f(xi, a)
  } else {
    b <- q				# shape parameter alpha (~ autocorrelation)
    a <- (2*x-b*x-1)/(x-1)	# shape parameter beta = f(xi, a)
  }
  return(c(a,b))
}
# TEST Plot density for illustration
# layout(matrix(1:4,2, byrow=TRUE))
# maxT <- 45			# maximum temperature
# minT <- 5			# minimum temperature
# T <- seq(minT,maxT,0.1)		# range of temperatures
# x <- scale01(seq(5,45,0.1),0,45)	# scaled range of temperatures [0:1]
# Ti <- 35					# temperature at time i
# xi <- (Ti-minT)/(maxT-minT)	# scaled temperature at time i (mode of distribution)
# xi <- 0.8
# q <- exp(4)
# B <- calcBeta.ab(xi,q)
# plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab='Scaled Temperature', 
# 	ylab='Density', yaxt='n')#, main='Density of Beta distribution')
#	mtext(paste('(mode = ', signif(xi,2), '; q = ', signif(q,2), ')', sep=''))
# abline(v=xi, lty=2, col='blue')
# points(x, dbeta(x,B[1],B[2])/max(dbeta(x,B[1],B[2])), type='l', col='blue')


sample.col <- function( X=NULL, size=NULL, w=NULL ){
  # Sample columns of a (frequency) table
  #	checks for sampling exceeding frequencies
  #	for single column (vector), use sample()
  #	X: frequency table to sample with multiple columns
  #	size: sample size, length(size) = ncol(table)
  #	w: not implemented, extra weighting
  if (is.null(dim(X))) stop('invalid dimensions for \'table\'')
  if (ncol(X)!=length(size)) stop('ncol(X) not equal to length(size)')
  X1 <- X
  X1[1:length(X1)] <- 0
  if(sum(X)==0||sum(size)==0){
    warning('zero probability / sample size so no sample was taken')
  } else {
    for (i in 1:length(size)){
      if (size[i]>0){
        ok <- FALSE
        while(!ok){
          if (size[i]==sum(X[,i])){
            x <- X[,i]
          } else if (size[i]<0.5*sum(X[,i])){
            x <- sample( x=1:nrow(X), size=size[i], prob=X[,i], replace=TRUE )
            x <- table(x)
          } else {
            x <- sample( x=1:nrow(X), size=sum(X[,i])-size[i], prob=X[,i], replace=TRUE )
            x <- table(x)
            x <- X[,i] - x
          }
          ok <- all( x <= X[,i] )
          ok <- all( x >= 0 )
        }
        X1[,i] <- X1[,i] + x
      }
    }
  }
  return(X1)
}


# varHyp <- function(mean, sd){
# #	Propagate the variance of x1 & x2 to the resulting vector/hypothenuse
# #		(Faber 2012 Statistics and Probability Theory: p.135)
# #		mean: vector for the means of x1 & x2
# #		sd: vector for the standard deviations of x1 & x2
# #		NOTE: Assumes no covariance between x1 & x2

# }



