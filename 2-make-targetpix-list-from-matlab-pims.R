# Reads in pixel intensity matrices stored in a Matlab struct. These will have been 
# created by running the script 'matchTargetMatrices.m' beforehand.
# Must have run '1-extract-IDs...' to create target_nfo file.

library(R.matlab)
library(dplyr)
library(readxl)
library(stringr)
library(lubridate)
library(tidyr)
library(moveHMM)

matfiles <- list.files(path = "data/mat", pattern = ".mat")

dat1 <- read_excel("data/Sonar_notes_123.xlsx", sheet = "Part1", col_types = c("date", rep("guess", 10)))
dat2 <- read_excel("data/Sonar_notes_123.xlsx", sheet = "Part2", col_types = c("date", rep("guess", 10)))
dat3 <- read_excel("data/Sonar_notes_123.xlsx", sheet = "Part3", col_types = c("date", rep("guess", 10)))

dat <- rbind(dat1,dat2,dat3)

# remove rows with no filenames
dat <- dat %>% filter(!is.na(filename))

# add recording date
dat <- dat %>% mutate(date = ymd(str_extract(dat$filename, "[0-9].{9}")))

# files with multiple targets have multiple rows, need to autofill down for these
dat <- dat %>% fill(site, start_time, end_time, class)

# select filenames by quality (or some other filter)
mydat <- dat %>% filter(quality >= 3)
myfiles <- paste0(mydat$filename, ".mat")
myfiles <- intersect(matfiles, myfiles)
length(myfiles)
mydat <- mydat %>% filter(filename %in% str_remove_all(myfiles, ".mat"))

target.pix <- vector(mode = "list", length = length(myfiles))
frame.nfo <- data.frame(targ_id = as.character(), frame_id = as.integer(), imgindex = as.integer(), 
                        x = as.numeric(), y = as.numeric())
targetIDs <- c()
minpix <- 300
maxpix <- 0
for(i in 1:length(myfiles)){
  print(i)
  x <- readMat(paste0("data/mat/",myfiles[i])) 
  targetID <- str_extract(myfiles[i], "[0-9].*[0-9]")
  nframes <- dim(x$sealtrackdat[[2]])[3]
  # number of non-NaN entries (NaN x means a 0x0 pixel intensity matrix)
  usable_frames <- !is.nan(unlist(x$sealtrackdat[[2]][,,1:nframes]["x",]))
  usable_frames_id <- (1:nframes)[usable_frames]
  nusframes <- length(usable_frames_id)
  target.pix[[i]] <- vector(mode = "list", length = nusframes)
  cnt <- 0
  for(j in 1:nusframes){
    cnt <- cnt + 1
    frame.nfo <- rbind(frame.nfo, 
                       data.frame(targ_id = targetID, frame_id = j, imgindex = x$sealtrackdat[[2]][,,j]$imgindex, 
                                  x = x$sealtrackdat[[2]][,,j]$x, y = x$sealtrackdat[[2]][,,j]$y))
    target.pix[[i]][[cnt]] <- x$sealtrackdat[[2]][,,j]$targetSS
    minpix <- min(minpix, as.numeric(target.pix[[i]][[cnt]]), na.rm = TRUE)
    maxpix <- max(maxpix, as.numeric(target.pix[[i]][[cnt]]), na.rm = TRUE)
    names(target.pix[[i]][[cnt]]) <- as.character(usable_frames_id[j])
  }
}
names(target.pix) <- as.character(myfiles)

length(target.pix)
length(target.pix[[1]])
image(target.pix[[1]][[10]])

# MIGHT CHAGNE THIS LATER
# change all NAs to frame minimum
for(i in 1:length(target.pix)){
  for(j in 1:length(target.pix[[i]])){
    target.pix[[i]][[j]][is.na(target.pix[[i]][[j]])] <- min(target.pix[[i]][[j]], na.rm = TRUE)
  }
}


# recode to 0-256 not -128,128
hist(target.pix[[1]][[22]])
for(i in 1:length(target.pix)){
  for(j in 1:length(target.pix[[i]])){
    target.pix[[i]][[j]] <- ifelse(target.pix[[i]][[j]] < 0, target.pix[[i]][[j]] + 256, target.pix[[i]][[j]])
  }
}
hist(target.pix[[1]][[2]])

# frame.nfo vx, vy calculation
frame.nfo <- frame.nfo %>% group_by(targ_id) %>%
  mutate(vx = x - lag(x), vy = y - lag(y)) %>%
  mutate(step = frame_id) %>% ungroup()

# target.nfo
target.nfo <- data.frame(id = as.character(myfiles), rawclass = mydat$class)
target.nfo <- target.nfo %>% 
  mutate(class = ifelse(str_detect(rawclass, "seal"), "seal", "nonseal")) %>%
  separate(col = rawclass, into = c("class_conf", "class_cat"), sep = " ", remove = FALSE, fill = "left") %>%
  mutate(class_cat = str_replace_all(class_cat, "\\?", "unknown")) %>%
  mutate(class_conf = str_replace_all(class_conf, "poss.", "possible")) %>%
  mutate(class_conf = replace_na(class_conf, "likely"))

target.nfo <- target.nfo %>% mutate(targ_id = str_extract(id, "[0-9].*[0-9]"))

save(target.pix, frame.nfo, target.nfo, file = "output/targets-w2022-hq.Rdata")

