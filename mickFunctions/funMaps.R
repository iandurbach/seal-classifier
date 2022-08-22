############################
#
#	FUNCTIONS FOR WORKING WITH MAPS AND SPATIAL DATA
#


#	Load required packages
#	Windows: require(Rtools)
require('sp')
require('Rcpp')

is.spatial <- function(x){ length(grep("Spatial",class(x)))>0 }

map.axes <- function(x=TRUE,y=TRUE,box=TRUE, ...){
#	Function to add map axes with formatting
#	(modified from map:::map.axes; "map" is not available for R v3.1.3)
	if (x) axis(1, ...)
	if (y) axis(2, ...)
	if (box) box()
}


map.circle <- function(sp, r=1, npt=100){
#	Function to generate a circle (SpatialPolygon) with a given radius (m)
#	sp: SpatialPoint object or matrix of coordinates (in UtM)
#	r: radius of circle in metres
#	npt: number of points to draw circle
	if ("Spatial"%in%is(sp)){
		proj <- proj4string(sp)
		coord <- coordinates(sp)
	} else {
		coord <- sp
	}
	x <- coord[,1] + ( r * sin(seq(0,2*pi, length=npt)) )
	y <- coord[,2] + ( r * cos(seq(0,2*pi, length=npt)) )
	poly <- Polygons(srl=list(Polygon(cbind(x,y))), ID=paste("radius_",r,sep=""))
	if ("Spatial"%in%is(sp)){
		return(SpatialPolygons(Srl=list(poly), proj4string=proj))
	} else {
		return(poly)
	}
}


utmzone <- function(coord){
#	Get the UTM zone for a coordinate (longitude, latitude)
#	coord: lon,lat coordinates; latitudes between -80 and 84
	if((coord[2]< -80) | (coord[2]> 84)) stop("UTM is limited to latitudes between -80 and +84")
	breaks <- seq(-80,84,8)
	breaks[length(breaks)] <- 84
	zoney <- LETTERS[3:24][-c(7,13)][sum(coord[2]>breaks)]
	zonex <- (floor((coord[1] + 180)/6) %% 60) + 1
	#	Deal with exceptions
	if (zoney=="V" & zonex==31) if(coord[1]>3) zonex <- 32
	if (zoney=="X" & zonex==32) zonex <- ifelse(coord[1]<9,31,33)
	if (zoney=="X" & zonex==34) zonex <- ifelse(coord[1]<21,33,35)
	if (zoney=="V" & zonex==36) if(coord[1]>33) zonex <- 37
	paste(zonex, zoney, sep="")
}



rm.holes <- function(Poly){
#	Function to remove holes inside a polygon (from package 'wild1' v1.09)
    is.hole <- lapply(Poly@Polygons,function(P)P@hole)
    is.hole <- unlist(is.hole)
    polys <- Poly@Polygons[!is.hole]
    Poly <- Polygons(polys,ID=Poly@ID)
    return(Poly)
}


gClip <- function(shp, bb){
#	Clip polygons to area within a bounding box
	#	shp: shapefile containing spatialPolygons
	#	bb: bounding box, a matrix of coordinates for vertices or a SpatialPolygon
	if (class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
	else b_poly <- as(extent(bb), "SpatialPolygons")
	gIntersection(shp, b_poly, byid = T)
}


spDistsPr <- function(pts, longlat=TRUE){
#	Function to get distances (km) between pairs of coordinates
	#	pts: matrix of coordinates, where each row is a pair of coordinates
	#		- columns 1:2 are LON,LAT of 1st point
	#		- columns 3:4 are LON,LAT of 2nd point
	#	longlat: whether coordinates are in longitude & latitude (default=TRUE)
	if (!is.matrix(pts)) pts <- as.matrix(pts)
	if (nrow(pts)==1){
		return(spDists(pts, longlat=longlat))
	} else {
		dists <- rep(NA, nrow(pts))
		index <- !apply(pts,1,function(x)any(is.na(x)))
		dists[index] <- apply( pts[index,], 1, function(x){spDistsN1(pts=matrix(x[1:2], ncol=2), 
						pt=matrix(x[3:4], ncol=2), longlat=longlat)} )
		return(dists)
	}
}


spDistsSeq <- function(pts, lag=1, longlat=TRUE ){
#	Calculate distances between sequential coordinates
	#	pts: matrix of coordinates ('LON','LAT')
	#	lag: interval size (default=1)
	#	longlat: whether coordinates are in longitude & latitude (default=TRUE)
	if (!is.matrix(pts)) pts <- as.matrix(pts)
	pts <- cbind(head(pts,-1),pts[-1,])
	dists <- c(NA, spDistsPr(pts, longlat=longlat))
	return(dists)
}


which.ininterv <- function(x, intervals,closed=c(TRUE,TRUE)){
#	Function to find indices of values inside an interval; returns a list of indices
#	x: vector of values to assess
#	intervals: table of intervals (columns: start & end)
#	closed: should left and right intervals bounds be closed (default: TRUE for both)
	if (!is.matrix(intervals)) intervals <- matrix(intervals, ncol=2)
	if (ncol(intervals)!=2) stop("'intervals' must have 2 columns" )
	if (any(intervals[,1]>intervals[,2])) stop("some intervals are invalid")
	if (closed[1]){
		left <- sapply(intervals[,1], function(int){x>=int})
	} else {
		left <- sapply(intervals[,1], function(int){x>int})
	}
	if (closed[2]){
		right <- sapply(intervals[,2], function(int){x<=int})
	} else {
		right <- sapply(intervals[,2], function(int){x<int})
	}
	out <- apply(left&right, 2, which); names(out) <- rownames(intervals)
	return(out)
}


which.huginterv <- function(x, intervals,closed=c(TRUE,TRUE)){
#	Function to find indices of values immediately preceding and immediately following an interval
#	returns a matrix of indices (columns)
#	x: vector of values to assess
#	intervals: table of intervals (columns: start & end)
#	closed: should left and right intervals bounds be closed (default: TRUE for both)
	if (!is.matrix(intervals)) intervals <- matrix(intervals, ncol=2)
	if (ncol(intervals)!=2) stop("'intervals' must have 2 columns" )
	if (closed[1]){
		left <- sapply(intervals[,1], function(int){x<int})
	} else {
		left <- sapply(intervals[,1], function(int){x<=int})
	}
	if (closed[2]){
		right <- sapply(intervals[,2], function(int){x>int})
	} else {
		right <- sapply(intervals[,2], function(int){x>=int})
	}
	left <- apply(left, 2, function(y)if(any(y)){sum(y)}else{NA})
	right <- apply(right, 2, function(y)if(any(y)){min(which(y))}else{NA})
	out <- cbind(left,right); names(out) <- rownames(intervals)
	return(out)
}



fitdistrXY <- function(vcov, n=10000, family="lognormal", longlat=TRUE){
#	Find variance of location error (hypotenuse) assuming covariance is negligeable
	#	vcov: variance-covariance matrix between lon & lat coordinates
	#	n: size of sample for numerical solution (default=10000)
	#	family: distribution family for the error (at the moment, only "normal" or "lognormal")
	#	longlat: whether coordinates are in lon-lat
	require(mvtnorm, quietly=TRUE); require(MASS, quietly=TRUE)
	mode(vcov) <- "numeric"
	if (sum(vcov)==0){
		out <- c(meanlog=NA,sdlog=0)
	} else {
		xy.boot <- rmvnorm(n=n, mean=c(0,0), sigma=vcov)
		d <- spDistsN1(pts=xy.boot, pt=c(0,0), longlat=longlat)	
		out <- fitdistr(x=d, densfun=family)$estimate
	}
	return(out)
}



pErrLoc <- function(xy, vcov, q, lower.tail=FALSE, n=5000, longlat=FALSE){
#	Calculate probability of error distance
	#	xy: matrix of coordinates (lon, lat)
	#	vcov: variance-covariance matrix between lon & lat coordinates
	#	q: threshold value
	#	lower.tail: whether to give probability of lower tail (default=FALSE)
	#	n: size of sample for numerical solution (default=10000)
	#	longlat: whether coordinates are in lon-lat
#	* Consider fitting Rayleigh distribution
	require(mvtnorm, quietly=TRUE)
	mode(vcov) <- "numeric"
	if (sum(vcov)==0){
		return(0)
	} else {
		xy.boot <- rmvnorm(n=n, mean=as.vector(xy), sigma=vcov)
		d <- spDistsN1(pts=xy.boot, pt=xy, longlat=longlat)	
		if (lower.tail){
			return(sum(d<q)/length(d))
		} else {
			return(sum(d>q)/length(d))
		}
	}
}


qErrLoc <- function(xy, vcov, p=0.95, n=5000, longlat=FALSE){
#	Calculate quantile of error distance
	#	xy: matrix of coordinates (lon, lat)
	#	vcov: variance-covariance matrix between lon & lat coordinates
	#	p: quantile
	#	n: size of sample for numerical solution
	#	longlat: whether coordinates are in lon-lat
#	* Consider fitting Rayleigh distribution
	require(mvtnorm, quietly=TRUE)
	mode(vcov) <- "numeric"
	if (sum(vcov)==0){
		return(0)
	} else {
		mode(xy) <- "numeric"
		xy.boot <- rmvnorm(n=n, mean=as.vector(xy), sigma=vcov)
		d <- spDistsN1(pts=xy.boot, pt=xy, longlat=longlat)	
		return(quantile(d, probs=p))
	}
}



#	Function to find intersection between two lines
#	modified from http://stackoverflow.com/questions/20519431/finding-point-of-intersection-in-r
#	(~7X faster than R code)
#	(think you can code? vectorized algorithm at https://software.intel.com/en-us/courseware/249568)
	# line1: vertices of line1 [x,y] (2x2 matrix) in UTM
	# line2: vertices of line2 [x,y] (2x2 matrix) in UTM
	# Returns the coordinates of the intersection [x,y] or if not intersecting: [NA,NA]
cppFunction('NumericVector intersectLinesRcpp(NumericMatrix line1, NumericMatrix line2) {
// Declare variables
	NumericVector xpoint(2);
	double x1 = line1(0,0);
	double x2 = line1(1,0);
	double x3 = line2(0,0);
	double x4 = line2(1,0);
	double y1 = line1(0,1);
	double y2 = line1(1,1);
	double y3 = line2(0,1);
	double y4 = line2(1,1);
	double denom;
	double ua;
	double ub;

	xpoint[0] = NA_REAL;
	xpoint[1] = NA_REAL;

	denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1);
// Assume parallel lines do not cross
	// TO DO: Check whether parallel lines overlap
	if ( abs(denom) > 1.0e-10 ) {
		ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom;
		ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom;
		if (ua>=0 && ua<=1 && ub>=0 && ub<=1) {
			xpoint[0] = x1 + ua*(x2-x1);
			xpoint[1] = y1 + ua*(y2-y1);
		}		
	}
// Return intersection point or NAs
	return(xpoint);
}')

intersectLinePoly <- function(ln, sl, longlat=FALSE){
#	Find where a line crosses other lines
	#	ln: vertices for one line (2x2 matrix) in UTM 
	#	sl: list of vertices for polygons (nx2 matrix) in UTM
	#	longlat: whether coordinates are longitude & latitude. Only FALSE for now
	if (longlat) stop("Long-lat not currently implemented; use UTM coordinates")
	#	Check whether there is overlap in at least one dimension
	sl.lim <- t(sapply(sl,function(x)apply(x,2,range)))		# Exclude:
	index <- which( !apply(min(ln[,1])>sl.lim[,1:2],1,all)	# polygons to the left
			 & !apply(max(ln[,1])<sl.lim[,1:2],1,all)		# polygons to the right
			 & !apply(min(ln[,2])>sl.lim[,3:4],1,all)		# polygons to the below
			 & !apply(max(ln[,2])<sl.lim[,3:4],1,all) )	# polygons to the above
	#	Find cross points (if any) for each polygon
	cpts <- lapply(index, function(y){
		sl.j <- cbind(head(sl[[y]],-1), sl[[y]][-1,])
		tmp <- t(apply(sl.j, 1, function(x){
				intersectLinesRcpp(ln, matrix(x,nrow=2,byrow=TRUE)) }) )
		tmp <- tmp[!is.na(tmp[,1]),]
		return(tmp)
		})
	# if none cross, return the first point of the segment 'ln'
	# (this can be due to imprecision when defining the segment)
	if (any(sapply(cpts,length)>0)){
		cpts <- cpts[sapply(cpts,length)!=0]
	} else {
		cpts <- list(ln[1,])
	}
	return(cpts)
}



#	Function to find point in line (sl) nearest to the point (pt)
#	pt: coordinates of reference point [x,y] in UTM
#	sl: vertices of the line [x,y] (2x2 matrix) in UTM
cppFunction('NumericVector nearestPtLineRcpp(NumericVector point, NumericMatrix line) {
// Declare variables
	NumericVector npoint(2);
	double x;
	double y;
	double m1;
	double m2;
// In case of horizontal line
	if ( (line(1,0)-line(0,0)) == 0 ) {
		x = line(0,0);
		y = point[1];
// In case of vertical line
	} else if ( (line(1,1)-line(0,1)) == 0 ) {
		x = point[0];
		y = line(0,0);
// In case of diagonal line
	} else {
		m1 = (line(1,1)-line(0,1)) / (line(1,0)-line(0,0));
		m2 = -1/m1;
		x = (m1*line(0,0) - m2*point[0] + point[1] - line(0,1)) / (m1 - m2);
		y = m2*(x-point[0]) + point[1];
	}
// Add avoid point beyond vertices
	if ( x < min(line(_,0)) ) {
		x = min(line(_,0));
		if ( y < min(line(_,1)) ) {
			y = min(line(_,1));
		} else {
			y = max(line(_,1));
		}
	} else if ( x > max(line(_,0)) ) {
		x = max(line(_,0));
		if ( y < min(line(_,1)) ) {
			y = min(line(_,1));
		} else {
			y = max(line(_,1));
		}
	}
// Return nearest point
	npoint[0] = x;
	npoint[1] = y;		
	return(npoint);
}')



#	See also maptools:::snapPointsToLines, which uses gDistance to find nearest line and then snaps to it
nearestPtPoly <- function(pt, sl, longlat=FALSE, ...){
#	Function to find point of Polygon nearest a given point
	#	pt: coordinates of point [x,y] in UTM
	#	sl: list of vertices for polygons (nx2 matrix) in UTM
	#	longlat: whether coordinates are longitude & latitude. Only FALSE for now
	if (longlat) stop("longlat not currently implemented")
	#	Find nearest point on each line of the Polygon
	npts <- vector(mode="list", length=length(sl))
	#	Loop through each polygon
	for (j in 1:length(sl)){
		# Loop through each segment of the Polygon
		npts[[j]] <- matrix(NA, nrow=nrow(sl[[j]])-1, ncol=2)
		for (i in 1:(nrow(sl[[j]])-1)){
			npts[[j]][i,] <- unlist(nearestPtLineRcpp(pt, sl[[j]][i:(i+1),]))
		}
	}
	npts <- do.call(rbind,npts)
	npts <- npts[which.min(spDistsN1(pts=npts, pt=pt, longlat=FALSE)),]
	return(npts)
}



aveCoord <- function(coords, family, fpar, longlat=FALSE){
#	Find average coordinates with error by maximum likelihood
#	assume distribution of error is the same at any angle
#	requires the function nlogLik.coord()
	#	coords: matrix of coordinates of points to average (x,y in UTM)
	#	family: character value for the family of error distribution
	#	fpar: matrix of parameters for error distribution; values are recycled
	if (!is.matrix(coords)) coords <- as.matrix(coords)
	if (!is.matrix(fpar)) fpar <- as.matrix(fpar)
	if (length(family)>1)warnings("only the first value of 'family' is used")
	if (length(unique(coords[,1]))==1 & length(unique(coords[,2]))==1){
		return(coords[1,])
	} else {
		out <- optim( par=colMeans(coords), fn=nlogLik.coord, coords=coords, 
					family=family, fpar=fpar, longlat=longlat)
		# out <- optim( par=colMeans(coords), fn=nlogLik.coord, coords=coords,family=family,fpar=fpar,
					  # method="L-BFGS-B", lower=apply(coords,2,min), upper=apply(coords,2,max))
		if (out$convergence!=0) warnings("possible convergence problems")
		return(out$par)
	}
}

nlogLik.coord <- function(x, coords, family, fpar, longlat=FALSE){
#	Negative log likelihood function to average coordinates: aveCoord()
	#	x: x and y coordinates of location to optimise
	#	coords: matrix of coordinates of points to average (x,y in UTM)
	#	family: character value for the family of error distribution
	#	fpar: matrix of parameters for error distribution; values are recycled
	fn <- switch(substr(tolower(family),1,1), n = pnorm, l = plnorm, g = pgamma, b = pbeta)
	d <- spDistsN1(pts=coords, pt=x, longlat=longlat)
	d <- ifelse(d>0,d,min(d[d>0])/1000)	#	avoid problems with zero distances 
	nlogLik <- -sum(fn(d,fpar[,1],fpar[,2], lower.tail=FALSE, log.p=TRUE))
	return(nlogLik)
}



velocity <- function(pts, tsec, longlat=FALSE){
#	Calculate velocity in m/s
	#	pts: matrix of x & y UTM coordinates (lon, lat)
	#	tsec: time variable in seconds
	#	Returns change in x and y coordinates (m/s)
	if(longlat) stop("Coordinates must be in UTM")
	tsec <- c(NA, diff(tsec))
	dx <- c(NA,diff(pts[,1]))*1000/tsec
	dy <- c(NA,diff(pts[,2]))*1000/tsec
	return(cbind(dx,dy))
}


extentPtsPolys <- function(pts, sp, tol.p=0, tol.km=0, step="max", longlat=FALSE){
#	Function to define the extent of Polygons with bounding boxes within the extent points
#	Returns a list of: extent of pts, extent of polygons included, condition for stopping function
#	TO DO:	maybe use convex polygons instead of rectangular bounding boxes
#			allow tolerance to include Polygons nearby
	#	pts: SpatialPoints (UTM coordinates) or matrix of x & y coordinates
	#	sp: SpatialPolygons (UTM coordinates) or list of vertices for each polygon
	#	tol.p, tol.km: increase extent of pts (proportion/km). Not yet implemented
	#	step: number of steps to find adjacent polygons (positive integer):
	#		- step=1: finds polygons near pts
	#		- step=2: adds polygons near previous polygons
	#		- (default="max"): recurse the process until there are no more nearby polygons
	#	longlat: whether coordinates are longitude & latitude. Only FALSE for now
	if (longlat) stop("'longlat' not currently implemented")
	proj <- NA
#	Extract extent of Polygons
	if (is.spatial(sp)){
		proj <- proj4string(sp)
		sp <- lapply(sp@polygons[[1]]@Polygons,coordinates)
	}
	sp.lim <- t(sapply(sp,function(x)apply(x,2,range)))
#	Extract extent of pts
	if (is.spatial(pts)){
		if(proj4string(pts)!=proj) stop("'pts' and 'sp' have different projection")
		pts <- coordinates(pts)
	}
	ext <- c(range(pts[,1]),range(pts[,2]))
#	Check tolerance
	if (tol.p!=0 | tol.km!=0) stop("tolerance 'tol.p' and 'tol.km' not yest implemented")
#	Recursive search for nearby polygons (overlapping bounding boxes)
	ext.prev <- ext
	cond <- TRUE
	i <- 0
	while(cond==TRUE){	
		i <- i+1
		index.lim <- which(	!( (ext.prev[1]>sp.lim[,2]) | (ext.prev[2]<sp.lim[,1])
					| (ext.prev[3]>sp.lim[,4]) | (ext.prev[4]<sp.lim[,3]) ) )
		ext.poly <- c(range(sp.lim[index.lim,1:2]),range(sp.lim[index.lim,3:4]))
#	Check whether to keep looping
		if (length(index.lim)==nrow(sp.lim)){cond <- "all polygons available"}
		if (is.numeric(step)){if(i>=step) cond <- paste("reached max step",step)	}
		if (all(ext.poly==ext.prev)){cond <- paste("found extent at step",i)}
		ext.prev <- ext.poly
	}
#	Return output
	return(list(pts=ext,poly=ext.poly,cond=cond))
}




interpolSL <- function(x, sl, xout, longlat=FALSE, ...){
#	interpolate along a SpatialLines object
	# x: numeric vector for start & end of track
	# sl: coordinates from a SpatialLine along which to interpolate
	# xout: point at which to interpolate
	# longlat: whether SpatialLines are in Lon Lat (only FALSE at the moment)
	# ...: other parameters for the approx function
	if (is.spatial(sl)) sl <- coordinates(sl)
	dcum <- c(0,cumsum(offdiag(spDists(sl))))
	i <- approx(x=x, y=c(0,max(dcum)), xout=xout, ...)$y
	out <- cbind(approx(hug(i,dcum),sl[which.hug(i,dcum),1], xout=i)$y,
				approx(hug(i,dcum),sl[which.hug(i,dcum),2], xout=i)$y)
	return(out)
}





