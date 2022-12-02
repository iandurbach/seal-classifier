# generate training/validation/test splits for all models

library(dplyr)

load("output/targets-w2022.RData")

#######
#######  single train/val/test sample
#######

# a much faster cleaner way to split
set.seed(222)
table(target.nfo$quality)
tr <- target.nfo %>% group_by(class,quality) %>% slice_sample(prop = 0.6) %>% ungroup()
tt <- target.nfo %>% filter(!(targ_id %in% tr$targ_id))
vl <- tt %>% group_by(class,quality) %>% slice_sample(prop = 0.25) %>% ungroup()
te <- tt %>% filter(!(targ_id %in% vl$targ_id))

training_ids_q1234 <- tr$targ_id
validation_ids_q1234 <- vl$targ_id
test_ids_q1234 <- te$targ_id

length(training_ids_q1234)
length(validation_ids_q1234)
length(test_ids_q1234)

# 
# # s_id for seals 
# s_id_seals <- unique(unlist(target.nfo[target.nfo$class == 1, "targ_id"]))
# # s_id for any non-seals 
# s_id_nonseals <- unique(unlist(target.nfo[target.nfo$class == 0, "targ_id"]))
# 
# set.seed(222)
# 
# # randomly shuffle IDs
# s_id_seals_t <- sample(s_id_seals)
# s_id_nonseals_t <- sample(s_id_nonseals)
# 
# # allocate to tr/v/te
# n_s_train <- round(length(s_id_seals_t) * .7, 0)
# n_s_val <- round(length(s_id_seals_t) * .1, 0)
# n_ns_train <- round(length(s_id_nonseals_t) * .7, 0)
# n_ns_val <- round(length(s_id_nonseals_t) * .1, 0)
# 
# training_ids_q1234 <- c(s_id_seals_t[1:n_s_train], s_id_nonseals_t[1:n_ns_train])
# validation_ids_q1234 <- c(s_id_seals_t[(1+n_s_train):(n_s_train+n_s_val)], s_id_nonseals_t[(1+n_ns_train):(n_ns_train+n_ns_val)])
# test_ids_q1234 <- c(s_id_seals_t[(1+n_s_train+n_s_val):length(s_id_seals_t)], 
#               s_id_nonseals_t[(1+n_ns_train+n_ns_val):length(s_id_nonseals_t)])

# shuffle again
training_ids_q1234 <- sample(training_ids_q1234)
validation_ids_q1234 <- sample(validation_ids_q1234)
test_ids_q1234 <- sample(test_ids_q1234)

### DO FOR QUALITY > 1 ONLY

target.nfo.q234 <- target.nfo %>% filter(quality > 1)

# s_id for seals 
s_id_seals <- unique(unlist(target.nfo.q234[target.nfo.q234$class == 1, "targ_id"]))
# s_id for any non-seals 
s_id_nonseals <- unique(unlist(target.nfo.q234[target.nfo.q234$class == 0, "targ_id"]))

# allocate IDs to t/v/te

set.seed(222)

# randomly shuffle IDs
s_id_seals_t <- sample(s_id_seals)
s_id_nonseals_t <- sample(s_id_nonseals)

# allocate to tr/v/te
n_s_train <- round(length(s_id_seals_t) * .7, 0)
n_s_val <- round(length(s_id_seals_t) * .1, 0)
n_ns_train <- round(length(s_id_nonseals_t) * .7, 0)
n_ns_val <- round(length(s_id_nonseals_t) * .1, 0)

training_ids <- c(s_id_seals_t[1:n_s_train], s_id_nonseals_t[1:n_ns_train])
validation_ids <- c(s_id_seals_t[(1+n_s_train):(n_s_train+n_s_val)], s_id_nonseals_t[(1+n_ns_train):(n_ns_train+n_ns_val)])
test_ids <- c(s_id_seals_t[(1+n_s_train+n_s_val):length(s_id_seals_t)], 
                    s_id_nonseals_t[(1+n_ns_train+n_ns_val):length(s_id_nonseals_t)])

# shuffle again
training_ids <- sample(training_ids)
validation_ids <- sample(validation_ids)
test_ids <- sample(test_ids)

# split into downriver and upriver
downriver_ids <- target.nfo[target.nfo$traveldir == "downriver", "targ_id"]
upriver_ids <- target.nfo[target.nfo$traveldir == "upriver", "targ_id"]
other_ids <- target.nfo[target.nfo$traveldir == "other", "targ_id"]

training_ids_up <- intersect(training_ids, upriver_ids)
training_ids_down <- intersect(training_ids, downriver_ids)
training_ids_other <- intersect(training_ids, other_ids)

validation_ids_up <- intersect(validation_ids, upriver_ids)
validation_ids_down <- intersect(validation_ids, downriver_ids)
validation_ids_other <- intersect(validation_ids, other_ids)

test_ids_up <- intersect(test_ids, upriver_ids)
test_ids_down <- intersect(test_ids, downriver_ids)
test_ids_other <- intersect(test_ids, other_ids)

training_ids_q1234_up <- intersect(training_ids_q1234, upriver_ids)
training_ids_q1234_down <- intersect(training_ids_q1234, downriver_ids)
training_ids_q1234_other <- intersect(training_ids_q1234, other_ids)

validation_ids_q1234_up <- intersect(validation_ids_q1234, upriver_ids)
validation_ids_q1234_down <- intersect(validation_ids_q1234, downriver_ids)
validation_ids_q1234_other <- intersect(validation_ids_q1234, other_ids)

test_ids_q1234_up <- intersect(test_ids_q1234, upriver_ids)
test_ids_q1234_down <- intersect(test_ids_q1234, downriver_ids)
test_ids_q1234_other <- intersect(test_ids_q1234, other_ids)

save(list = ls(pattern = "(training_ids)|(validation_ids)|(test_ids)"), 
     file = "output/ids_traintestval_1ttv.Rdata")


#######
#######  K-fold
#######

#### DO FOR ALL QUALITY CLASSES

# make train/test split

# s_id for seals 
s_id_seals <- unique(unlist(target.nfo[target.nfo$class == 1, "targ_id"]))
# s_id for any non-seals 
s_id_nonseals <- unique(unlist(target.nfo[target.nfo$class == 0, "targ_id"]))

# allocate IDs for k-fold CV

set.seed(222)

# randomly shuffle IDs
s_id_seals_t <- sample(s_id_seals)
s_id_nonseals_t <- sample(s_id_nonseals)

# allocate IDs to folds
k <- 10 
n_seal_per_fold <- floor(length(s_id_seals)/k)
n_nonseal_per_fold <- floor(length(s_id_nonseals)/k)
seals_fold_id <- rep(1:k, n_seal_per_fold + 1)[1:length(s_id_seals)]
nonseals_fold_id <- rep(1:k, n_nonseal_per_fold + 1)[1:length(s_id_nonseals)]

test_ids <- list()
training_ids <- list()
validation_ids <- list()
for(i in 1:k){
  test_ids[[i]] <- c(s_id_seals_t[seals_fold_id == i], s_id_nonseals_t[nonseals_fold_id == i])
  nontest_ids_seal <- s_id_seals_t[seals_fold_id != i]
  nontest_ids_nonseal <- s_id_nonseals_t[nonseals_fold_id != i]  
  training_ids[[i]] <- c(sample(nontest_ids_seal, size = 0.9 * length(nontest_ids_seal)),
                         sample(nontest_ids_nonseal, size = 0.9 * length(nontest_ids_nonseal)))
  validation_ids[[i]] <- setdiff(c(nontest_ids_seal, nontest_ids_nonseal), training_ids[[i]])
  # shuffle again
  test_ids[[i]] <- sample(test_ids[[i]])
  training_ids[[i]] <- sample(training_ids[[i]])
  validation_ids[[i]] <- sample(validation_ids[[i]])
}

training_ids_q1234 <- training_ids
validation_ids_q1234 <- validation_ids
test_ids_q1234 <- test_ids

### DO FOR QUALITY > 1 ONLY

target.nfo.q234 <- target.nfo %>% filter(quality > 1)

# make train/test split

# s_id for seals 
s_id_seals <- unique(unlist(target.nfo[target.nfo$class == 1, "targ_id"]))
# s_id for any non-seals 
s_id_nonseals <- unique(unlist(target.nfo[target.nfo$class == 0, "targ_id"]))

# allocate IDs for k-fold CV

set.seed(222)

# randomly shuffle IDs
s_id_seals_t <- sample(s_id_seals)
s_id_nonseals_t <- sample(s_id_nonseals)

# allocate IDs to folds
k <- 10 
n_seal_per_fold <- floor(length(s_id_seals)/k)
n_nonseal_per_fold <- floor(length(s_id_nonseals)/k)
seals_fold_id <- rep(1:k, n_seal_per_fold + 1)[1:length(s_id_seals)]
nonseals_fold_id <- rep(1:k, n_nonseal_per_fold + 1)[1:length(s_id_nonseals)]

test_ids <- list()
training_ids <- list()
validation_ids <- list()
for(i in 1:k){
  test_ids[[i]] <- c(s_id_seals_t[seals_fold_id == i], s_id_nonseals_t[nonseals_fold_id == i])
  nontest_ids_seal <- s_id_seals_t[seals_fold_id != i]
  nontest_ids_nonseal <- s_id_nonseals_t[nonseals_fold_id != i]  
  training_ids[[i]] <- c(sample(nontest_ids_seal, size = 0.9 * length(nontest_ids_seal)),
                         sample(nontest_ids_nonseal, size = 0.9 * length(nontest_ids_nonseal)))
  validation_ids[[i]] <- setdiff(c(nontest_ids_seal, nontest_ids_nonseal), training_ids[[i]])
  # shuffle again
  test_ids[[i]] <- sample(test_ids[[i]])
  training_ids[[i]] <- sample(training_ids[[i]])
  validation_ids[[i]] <- sample(validation_ids[[i]])
}

save(training_ids_q1234, validation_ids_q1234, test_ids_q1234,
     training_ids, validation_ids, test_ids, file = "output/ids_traintestval_kfold.Rdata")
