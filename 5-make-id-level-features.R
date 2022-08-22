# Calculates ID-level features and makes a new df target.var containing these
#
# Some ID-level features are aggregated frame-level features over frames, some are
# new and can't be calculated per frame (spectral features)

library(dplyr)
require(parallel); require(unbalanced); require(kernlab)
require(McSpatial); require(raster)
require(animation); require(scales)
source('mickFunctions/funBasics.R', chdir = TRUE)
source('mickFunctions/funMatrix.R', chdir = TRUE)
mc <- ifelse(detectCores()>4, detectCores()-2, detectCores())
options(scipen=999)	# avoid scientific notation (long sequences of numbers were rounded)

load("output/targets-w2022-hq.RData")

# # can create target.var from either target.nfo.FromPIMs or target.nfo.Tritech
# frame.nfo <- frame.nfo
# 
# # make 2021 look like 2015
# frame.nfo <- frame.nfo %>% mutate(id = targetID, valid = targetType) %>%
#   mutate(class = 0) %>% # for now make all class = non-seal; will fix up later
#   arrange(id) %>% group_by(id) %>% mutate(step = row_number()) %>% ungroup()
#   
#	SUMMARISE FEATURES BY TARGET

#		Create data frame
target.var <- data.frame(targ_id=unique(frame.nfo$targ_id),
                         class=tapply(frame.nfo$class, frame.nfo$targ_id, unique),
                         nsteps=tapply(frame.nfo$step, frame.nfo$targ_id, max), stringsAsFactors=FALSE)

#		Calculate basic summary statistics
frame.nfo.j <- split(frame.nfo, frame.nfo$targ_id)
tmp <- mclapply(frame.nfo.j, function(x){ 
  out <- as.vector(apply( x[,c("dist", "radi", "ddist", "dradi", "dvec",
                               "area", "perim", "length", "pa.ratio", "la.ratio", "lp.ratio",
                               "pixmin", "pixmax", "pixmed", "pixave", "pixsd", "pixrng", 
                               "pixpmin", "pixpmax", "pixpmed", "pixpave", "pixpsd", "pixprng") ], 2, function(y){
                                 c(	min(y,na.rm=TRUE), max(y,na.rm=TRUE),
                                    median(y,na.rm=TRUE), mean(y,na.rm=TRUE),
                                    sd(y,na.rm=TRUE), diff(range(y,na.rm=TRUE))) }))
  out <- c(	out, apply(x[,c("dist", "radi", "ddist", "dradi", "dvec")],2, 
                       function(y)sum(y==0, na.rm=TRUE)/length(y)) )
  return(out)
}, mc.cores=mc)
tmp <- do.call(rbind,tmp)
colnames(tmp) <- c(  as.vector(sapply(c("dist", "radi", "ddist", "dradi", "dvec", 
                                        "area", "perim", "length", "pa.ratio", "la.ratio", "lp.ratio",
                                        "pixmin", "pixmax", "pixmed", "pixave", "pixsd", "pixrng", 
                                        "pixpmin", "pixpmax", "pixpmed", "pixpave", "pixpsd", "pixprng"), 
                                      function(x){ paste(c("min","max","med","ave","sd","rng"),x , sep="." )} ) ),
                     paste("pzero",c("dist", "radi", "ddist", "dradi", "dvec"),sep=".")  )
target.var <- cbind(target.var, tmp)

####
#		Calculate spectral features
####

#	Methods
#		Autoregressive time series model: spec.ar()
#		- order of AR fitted by AIC
#		- Warning about validity; see help
#		Fast Fourier transform (+ modified Daniell smooth): spec.pgram()
#		- remove trend ~ stationary process
#		- specify spans of moving average for smoothing: affects minimum length of time series
#	Output
#		npks: number of peaks to output
#		Frequency of npks peaks
#		Amplitude of npks peaks
#		Distribution of spectral density ~ Chi square
frame.nfo.j <- split(frame.nfo, frame.nfo$targ_id)
spans <- c(7,7)	# smoothing but larger values exclude short time series. smooth peaks at (7,7)
npks <- 2
tmp <- mclapply(frame.nfo.j, function(x){ 
  if (nrow(x)<(sum(spans)+1)){
    out <- rep(0,2*npks+2)
  } else {
    out <- as.vector( apply( x[,c("dist", "radi", "ddist", "dradi", "dvec",
                                  "area", "perim", "length", "pa.ratio", "la.ratio", "lp.ratio",
                                  "pixmin", "pixmax", "pixmed", "pixave", "pixsd", "pixrng") ], 2, function(y){
                                    out.spec <- rep(0,2*npks+2)
                                    names(out.spec) <- c(paste(rep(c("pkfreq","pkampl"),each=npks), rep(1:npks,2), sep=""), "npk", "spdens")
                                    per <- spec.pgram(na.exclude(y), spans=spans, detrend=TRUE, na.action=na.omit, plot=FALSE)
                                    index <- na.exclude(which(peaks(per$spec)))
                                    if (length(index)>0){
                                      out.spec["npk"] <- length(index)
                                      out.spec["spdens"] <- per$df
                                      index <- index[na.exclude(which.maxa(per$spec[index], npks))]
                                      out.spec[1:length(index)] <- per$freq[index]
                                      out.spec[npks+1:length(index)] <- per$spec[index]
                                    }
                                    return(out.spec)
                                  }) )
  }
  return(out)
}, mc.cores=mc)
tmp <- do.call(rbind,tmp)
colnames(tmp) <- as.vector(sapply( c("dist", "radi", "ddist", "dradi", "dvec",
                                     "area", "perim", "length", "pa.ratio", "la.ratio", "lp.ratio",
                                     "pixmin", "pixmax", "pixmed", "pixave", "pixsd", "pixrng"), function(x){
                                       paste(c("pkfreq1","pkfreq2","pkampl1","pkampl2","pkn","spdens"),x, sep="." )} ) )
target.var <- cbind(target.var,tmp)

# # outcomment as needed
# target.var.FromPIMs <- target.var
# target.var.Tritech <- target.var

# variable names

#	Movement (minima almost all zero so skip)
vmov <- c("nsteps", "pzero.dist", "pzero.radi", "pzero.dvec", # propZero ddist & dradi same as dvec
          "med.dist", "med.radi", "med.ddist", "ave.dradi", "med.dvec", #  "med.dradi" is constant so use mean
          "sd.dist", "sd.radi", "sd.ddist", "sd.dradi", "sd.dvec")

#	Shape (skip size, already in Tritech algorithm)
vshp <- c("med.length", "med.perim", "med.pa.ratio", "med.la.ratio",
          "sd.length", "sd.perim", "sd.pa.ratio", "sd.la.ratio",
          "min.length", "min.perim", "min.pa.ratio", "min.la.ratio",
          "max.length", "max.perim", "max.la.ratio") 	# "max.pa.ratio" constant at 1, so exclude

#	Pixels
vpix <- c("med.pixmin", "med.pixmax", "med.pixmed", "med.pixsd", 
          "min.pixmax", "min.pixmed", "min.pixsd",
          "max.pixmin", "max.pixmax", "max.pixmed", "max.pixsd",
          "sd.pixmin", "sd.pixmax", "sd.pixmed", "sd.pixsd")

#	Movement spectrum
vmspec <- c("pkfreq1.dist", "pkfreq2.dist", "pkampl1.dist", "pkampl2.dist", "pkn.dist", "spdens.dist",
            "pkfreq1.dradi", "pkfreq2.dradi", "pkampl1.dradi", "pkampl2.dradi", "pkn.dradi", "spdens.dradi", 
            "pkfreq1.dvec", "pkfreq2.dvec", "pkampl1.dvec", "pkampl2.dvec", "pkn.dvec", "spdens.dvec")

#	Shape spectrum
vsspec <- c("pkfreq1.area", "pkfreq2.area", "pkampl1.area", "pkampl2.area", "pkn.area", "spdens.area", 
            "pkfreq1.length", "pkfreq2.length", "pkampl1.length", "pkampl2.length", "pkn.length", "spdens.length", 
            "pkfreq1.pa.ratio", "pkfreq2.pa.ratio", "pkampl1.pa.ratio", "pkampl2.pa.ratio", "pkn.pa.ratio", "spdens.pa.ratio",
            "pkfreq1.la.ratio", "pkfreq2.la.ratio", "pkampl1.la.ratio", "pkampl2.la.ratio", "pkn.la.ratio", "spdens.la.ratio")

#	Pixels spectrum
vpspec <- c("pkfreq1.pixmin", "pkfreq2.pixmin", "pkampl1.pixmin", "pkampl2.pixmin", "pkn.pixmin", "spdens.pixmin", 
            "pkfreq1.pixmax", "pkfreq2.pixmax", "pkampl1.pixmax", "pkampl2.pixmax", "pkn.pixmax", "spdens.pixmax",
            "pkfreq1.pixmed", "pkfreq2.pixmed", "pkampl1.pixmed", "pkampl2.pixmed", "pkn.pixmed", "spdens.pixmed", 
            "pkfreq1.pixsd", "pkfreq2.pixsd", "pkampl1.pixsd", "pkampl2.pixsd", "pkn.pixsd", "spdens.pixsd")

vars <- c(vmov,vshp,vpix,vmspec,vsspec,vpspec)

# # choose which one to use in next steps (can overwrite later if needed, got both)
# target.var <- target.var.FromPIMs

target.var <- target.var %>% left_join(target.nfo %>% dplyr::select(targ_id, rawclass, class_conf), by = "targ_id")
target.nfo <- target.var

#		Save data
save(frame.nfo, target.nfo, target.pix, 
     file="output/targets-w2022-hq.RData")

