library(dplyr)
#library(keras)
#use_virtualenv("myenv", required = TRUE) 

# load("output/targets-w2022.RData")
load("output/lrcn-raw-inputs-frames-from-mid.Rdata")
load("output/ids_traintestval_1ttv.Rdata")

# # if you want to make a smaller dataset for testing, else comment out
# set.seed(123)
# training_ids <- sample(orig_training_ids, 20)
# validation_ids <- sample(orig_validation_ids, 10)
# test_ids <- sample(orig_test_ids, 10)
 
# train_data <- target.nfo[target.nfo$targ_id %in% training_ids, ] %>% dplyr::select(targ_id, class) 
# validation_data <- target.nfo[target.nfo$targ_id %in% validation_ids, ] %>% dplyr::select(targ_id, class) 
# test_data <- target.nfo[target.nfo$targ_id %in% test_ids, ] %>% dplyr::select(targ_id, class) 

xt_row_train <- unlist(seq_df_filt[seq_df_filt$targ_id %in% training_ids_up, "xtL_row"])
xt_row_validation <- unlist(seq_df_filt[seq_df_filt$targ_id %in% validation_ids_up, "xtL_row"])
xt_row_test <- unlist(seq_df_filt[seq_df_filt$targ_id %in% test_ids_up, "xtL_row"])

x_train <- xtL[xt_row_train,,,] 
y_train <- as.matrix(seq_df_filt$class[xt_row_train], ncol = 1)

x_validation <- xtL[xt_row_validation,,,] 
y_validation <- as.matrix(seq_df_filt$class[xt_row_validation], ncol = 1)

x_test <- xtL[xt_row_test,,,] 
y_test <- as.matrix(seq_df_filt$class[xt_row_test], ncol = 1)

save(training_ids, validation_ids, test_ids, x_train, x_validation, x_test, y_train, y_validation, y_test, file = "output/lrcn-inputs.Rdata")

# downsample here if you want
for(i in 1:4){
  selvec <- rep(FALSE, 4)
  selvec[i] <- TRUE
  xt_row_train_t <- xt_row_train[selvec]
  xt_row_validation_t <- xt_row_validation[selvec]
  xt_row_test_t <- xt_row_test[selvec]
  
  x_train <- xtL[xt_row_train_t,,,] 
  y_train <- as.matrix(seq_df_filt$class[xt_row_train_t], ncol = 1)
  
  x_validation <- xtL[xt_row_validation_t,,,] 
  y_validation <- as.matrix(seq_df_filt$class[xt_row_validation_t], ncol = 1)
  
  x_test <- xtL[xt_row_test_t,,,] 
  y_test <- as.matrix(seq_df_filt$class[xt_row_test_t], ncol = 1)
  
  save(training_ids, validation_ids, test_ids, x_train, x_validation, x_test, y_train, y_validation, y_test, file = paste0("output/lrcn-inputs-",i,".Rdata"))
}

# # Redefine  dimension of train/test inputs
# x_train <- array_reshape(x_train, c(dim(x_train), 1))
# x_validation <- array_reshape(x_validation, c(dim(x_validation), 1))
# x_test <- array_reshape(x_test, c(dim(x_test), 1))
# 
# # convert to binary class matrices
# y_train <- to_categorical(y_train, 2)
# y_validation <- to_categorical(y_validation, 2)
# y_test <- to_categorical(y_test, 2)
# 
# # Transform RGB values into [0,1] range
# x_train <- x_train / max(x_train)
# x_validation <- x_validation / max(x_train)
# x_test <- x_test / max(x_train)
# 
# # depending on what you choose above
# #save(training_ids, validation_ids, test_ids, x_train, x_validation, x_test, y_train, y_validation, y_test, file = "output/lrcn-inputs.Rdata")
# save(training_ids, validation_ids, test_ids, x_train, x_validation, x_test, y_train, y_validation, y_test, file = "output/lrcn-inputs-proc-all.Rdata")
