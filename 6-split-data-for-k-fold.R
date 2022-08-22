# generate training/validation/test splits for all models

library(dplyr)

load("output/targets-w2022-hq.RData")

# drop one very long track
target.nfo <- target.nfo %>% group_by(id) %>% mutate(nf = n()) %>% ungroup()
max(target.nfo$nf)
target.nfo <- target.nfo %>% filter(nf < 300)

target.nfo$class <- as.factor(target.nfo$class)

# # add HMM stuff
# target.nfo.HMM <- target.nfo %>% dplyr::select(ID = id, x, y) %>% as.data.frame()
# target.nfo.HMM <- target.nfo.HMM %>% mutate(x = x + runif(nrow(target.nfo.HMM), -0.0001, 0.0001),
#                                             y = y + runif(nrow(target.nfo.HMM), -0.0001, 0.0001))
# target.nfo.HMM <- prepData(target.nfo.HMM, type = "UTM")
# target.nfo$hmmstep <- target.nfo.HMM$step
# target.nfo$hmmangle <- target.nfo.HMM$angle

target.nfo$class <- ifelse(target.nfo$class == "seal", 1, 0)

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
}

save(training_ids, validation_ids, test_ids, file = "output/ids_traintestval.Rdata")
