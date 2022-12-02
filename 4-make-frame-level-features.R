# Calculates frame-level features and adds these to the target.nfo data frame
# 
# frame-level = movement, size/shape

# Note: In 2021 Tritech provided their own calculations of derived variables in target.nfo, which is 
#       read in by script '1-extract-....R'. These variables can also be recalculated from a combination
#       of a few input variables from target.nfo (vx, vy, step) and the pixel intensity matrices.
#       These should give the same answers but don't, need to ask Tritech or someone why.

library(parallel)
library(raster)
library(dplyr)
source('mickFunctions/funBasics.R', chdir = TRUE)
source("mickFunctions/funMatrix.R", chdir = TRUE)
mc <- 2 # ifelse(detectCores()>4, detectCores()-2, detectCores())
options(scipen=999)	# avoid scientific notation (long sequences of numbers were rounded)

load("output/targets-w2022.Rdata")

####
#	Calculate movement variables
####

#	distance, angle
frame.nfo$dist <- sqrt(frame.nfo$vx^2 + frame.nfo$vy^2)
frame.nfo$radi <- atan2(frame.nfo$vx,frame.nfo$vy)
#		Calculate movement variables by ID
frame.nfo.j <- split(frame.nfo, frame.nfo$targ_id)
#	delta distance
frame.nfo$ddist <- unlist(mclapply(frame.nfo.j,function(x){
  c(NA,diff(x$dist)/diff(x$step)) }, mc.cores=mc))
#	delta angle
frame.nfo$dradi <- unlist(mclapply(frame.nfo.j,function(x){
  c(NA,diff(x$radi)/diff(x$step)) }, mc.cores=mc))
#	delta velocity
frame.nfo$dvec <- unlist(mclapply(frame.nfo.j,function(x){
  c(NA,offdiag(as.matrix(dist(x[,c("vx","vy")])))) },mc.cores=mc))

# matrix size (assumes padded to same size)
size <- max(unlist(sapply(target.pix,function(x)sapply(x,dim))))

#		Calculate target shape and pixel variables
#	Prepare target raster
# Temporary fix because of duplicate path numbers
tmp <- lapply(names(target.pix), function(j){
  print(j)
  # j="1825580036"	# test with known duplicates
  # j="1748370037"	# test with file producing NA's
  # j="1748370038"	# test with file producing NA's
  x <- target.pix[[j]]
  #names(x) <- frame.nfo.j[[j]]$step		# correct any duplicate step numbers
  # tmp <- mclapply(target.pix, function(x){
  out.mat <- matrix(NA,nrow=length(x),ncol=18, dimnames=list(names(x),
                                                             c(	"area","perim","length","pa.ratio","la.ratio","lp.ratio",
                                                                "pixmin","pixmax","pixmed","pixave","pixsd","pixrng",
                                                                "pixpmin","pixpmax","pixpmed","pixpave","pixpsd","pixprng")))
  for (i in 1:length(x)){
    r <- raster(x[[i]]) #; image(r)
    values(r)[values(r)<max(values(r))/3] <- NA
    rmask <- values(boundaries(r, directions=4))
    xyperim <- coordinates(boundaries(r, directions=4))[!is.na(rmask),]
    out.mat[i,] <- c(
      #	Area, Perimeter, AP ratio
      area=sum(!is.na(rmask)),
      perim=sum(rmask, na.rm=TRUE),
      length=max(dist(xyperim))*size,
      pa=sum(rmask, na.rm=TRUE)/sum(!is.na(rmask)),
      la=max(dist(xyperim))*size/sum(!is.na(rmask)),
      lp=max(dist(xyperim))*size/sum(rmask, na.rm=TRUE),
      #	Whole target pixel: min, max, median, mean, sd, range
      pixmin=min(values(r)[!is.na(rmask)], na.rm=TRUE),
      pixmax=max(values(r)[!is.na(rmask)], na.rm=TRUE),
      pixmed=median(values(r)[!is.na(rmask)], na.rm=TRUE),
      pixave=mean(values(r)[!is.na(rmask)], na.rm=TRUE),
      pixsd=sd(values(r)[!is.na(rmask)], na.rm=TRUE),
      pixrng=diff(range(values(r)[!is.na(rmask)], na.rm=TRUE)),
      #	Perimeter target pixel: median, mean, sd, range
      ppmin=min(values(r)[rmask==1], na.rm=TRUE),
      ppmax=max(values(r)[rmask==1], na.rm=TRUE),
      ppmed=median(values(r)[rmask==1], na.rm=TRUE),
      ppave=mean(values(r)[rmask==1], na.rm=TRUE),
      ppsd=sd(values(r)[rmask==1], na.rm=TRUE),
      pprng=diff(range(values(r)[rmask==1], na.rm=TRUE)) )
  }
  return(out.mat)	# return(cbind(id=j,step=names(x),out.mat))	# for checking that rows match
})
tmp <- do.call(rbind,tmp)
# paste(tmp[,1],tmp[,2]) == paste(frame.nfo$id,frame.nfo$step)		# Check that rows match
# index <- which(apply(tmp,1,function(x)any(is.na(x))))	# Check for NA's
frame.nfo <- cbind(frame.nfo,tmp)
rm(tmp); gc()

# # going to use this version in subsequent analyses but save in another name so can switch between
# # versions if desired
# frame.nfo.FromPIMs <- frame.nfo

#		Save data
save(frame.nfo, target.nfo, target.pix, file="output/targets-w2022.RData")
save(frame.nfo, file="output/targets-w2022-framenfo.RData")

