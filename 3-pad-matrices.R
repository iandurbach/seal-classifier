# pads matrices

source('mickFunctions/funBasics.R', chdir = TRUE)
source('mickFunctions/funMatrix.R', chdir = TRUE)

load("output/targets-w2022.Rdata")

# target.pix.nopad <- target.pix

# pad as much as needed
size <- max(unlist(sapply(target.pix,function(x)sapply(x,dim))))
target.pix <- lapply(target.pix, function(x){
  lapply(x,function(y){
    pad <- cbind(floor((size-dim(y))/2), ceiling((size-dim(y))/2))
    return(padMatrix(y,padrow=pad[1,],padcol=pad[2,],value=0))
  })
})
names(target.pix) <- sapply(strsplit(names(target.pix),"/"),tail,1)
target.pix <- target.pix[order(names(target.pix))]

save(target.pix, frame.nfo, target.nfo, file = "output/targets-w2022.Rdata")
