# for extensions see https://tensorflow.rstudio.com/guide/keras/examples/conv_lstm/
# https://blogs.rstudio.com/ai/posts/2020-12-17-torch-convlstm/

library(keras)
library(dplyr)
library(purrr)

load('output/targets-w2022-hq.RData')
load("output/ids_traintestval.Rdata")

n_mat_per_target <- target.pix %>% map_depth(1, length) %>% unlist() 
n_mat <- sum(n_mat_per_target) # number of matrices
dim_mat <- dim(target.pix[[1]][[1]]) # assume all matrices have same dim
xt <- array(rep(0, n_mat * dim_mat[1] * dim_mat[2]), dim = c(n_mat, dim_mat[1], dim_mat[2]))
cnt <- 0
for(i in 1:length(target.pix)){
  for(j in 1:length(target.pix[[i]])){
    cnt <- cnt + 1
    xt[cnt,,] <- target.pix[[i]][[j]]
  }
}

seq_df <- data.frame(xt_row = 1:nrow(xt),
                     s_id = rep(1:length(target.pix), n_mat_per_target)) %>%
  group_by(s_id) %>% mutate(f_id = row_number()) %>% ungroup() %>%
  left_join(target.var %>% dplyr::select(targ_id, class) %>% mutate(s_id = row_number()), by = "s_id")

# CNN

num_classes <- 2

# Input image dimensions
img_rows <- 105
img_cols <- 105
input_shape <- c(img_rows, img_cols, 1)

# Define model
model1 <- keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu',
                input_shape = input_shape) %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') %>% 
  layer_max_pooling_2d(pool_size = c(2, 2)) %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_flatten() %>% 
  layer_dense(units = 128, activation = 'relu') %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_dense(units = num_classes, activation = 'softmax')

model2 <- keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu',
                input_shape = input_shape) %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') %>% 
  layer_max_pooling_2d(pool_size = c(2, 2)) %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_flatten() %>% 
  layer_dense(units = 128, activation = 'relu') %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_dense(units = num_classes, activation = 'softmax')

model3 <- keras_model_sequential() %>%
  layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = 'relu',
                input_shape = input_shape) %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = 'relu') %>% 
  layer_max_pooling_2d(pool_size = c(2, 2)) %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_flatten() %>% 
  layer_dense(units = 128, activation = 'relu') %>% 
  layer_dropout(rate = 0.25) %>% 
  layer_dense(units = num_classes, activation = 'softmax')

summary(model1)
summary(model2)
summary(model3)

run_cnn <- function(seq_df, target.nfo, training_ids, validation_ids, test_ids, cnn_model, learn_rate = 0.001, batch_size = 64, epochs = 50, run_id = 0, hyppar_id = 0, mod_only = FALSE){
  
  train_data <- target.nfo[target.nfo$id %in% training_ids, ] %>% dplyr::select(id, class) %>% mutate(run_id = run_id, hyppar_id = hyppar_id)
  validation_data <- target.nfo[target.nfo$id %in% validation_ids, ] %>% dplyr::select(id, class) %>% mutate(run_id = run_id, hyppar_id = hyppar_id)
  test_data <- target.nfo[target.nfo$id %in% test_ids, ] %>% dplyr::select(id, class) %>% mutate(run_id = run_id, hyppar_id = hyppar_id)
  
  xt_row_train <- unlist(seq_df[seq_df$id %in% training_ids, "xt_row"])
  xt_row_validation <- unlist(seq_df[seq_df$id %in% validation_ids, "xt_row"])
  xt_row_test <- unlist(seq_df[seq_df$id %in% test_ids, "xt_row"])
  
  x_train <- xt[xt_row_train,,] 
  y_train <- as.matrix(target.nfo$class[xt_row_train], ncol = 1)
  
  x_validation <- xt[xt_row_validation,,] 
  y_validation <- as.matrix(target.nfo$class[xt_row_validation], ncol = 1)
  
  x_test <- xt[xt_row_test,,] 
  y_test <- as.matrix(target.nfo$class[xt_row_test], ncol = 1)
  
  # Redefine  dimension of train/test inputs
  x_train <- array_reshape(x_train, c(nrow(x_train), img_rows, img_cols, 1))
  x_validation <- array_reshape(x_validation, c(nrow(x_validation), img_rows, img_cols, 1))
  x_test <- array_reshape(x_test, c(nrow(x_test), img_rows, img_cols, 1))
  
  # Transform RGB values into [0,1] range
  x_train <- x_train / max(x_train)
  x_validation <- x_validation / max(x_train)
  x_test <- x_test / max(x_train)
  
  # Convert class vectors to binary class matrices
  y_train <- to_categorical(y_train, num_classes)
  y_validation <- to_categorical(y_validation, num_classes)
  y_test <- to_categorical(y_test, num_classes)
  
  # Callbacks
  es <- list(callback_early_stopping(monitor = "val_acc", 
                                     min_delta = 0.1, 
                                     patience = 20, 
                                     mode = "auto"),
             callback_model_checkpoint(filepath = paste0("best_model_", run_id, ".h5"), 
                                       monitor = "val_loss", 
                                       save_best_only = TRUE, 
                                       mode = "auto"))
  
  # cnn
  # Compile model
  cnn_model %>% compile(
    loss = "binary_crossentropy",
    optimizer = optimizer_adam(lr = learn_rate),
    metrics = c('accuracy')
  )
  
  cnn_model %>% fit(
    x_test, y_test,
    batch_size = batch_size,
    epochs = epochs,
    validation_data = list(x_validation, y_validation),
    callbacks = es
  )
  
  print(paste0("saved best_model_", run_id, ".h5"))
  
  if(!mod_only) {
    
    best_model <- load_model_hdf5(paste0("best_model_", run_id, ".h5"))
    
    # add frame-level predicted probabilities and classes
    train_data$predclass <- best_model %>% predict_classes(x_train, verbose = 0) 
    
    validation_data$predprob <- best_model %>% predict_proba(x_validation, verbose = 0) %>% .[, 2]
    validation_data$predclass <- best_model %>% predict_classes(x_validation, verbose = 0) 
    
    test_data$predprob <- best_model %>% predict_proba(x_test, verbose = 0) %>% .[, 2]
    test_data$predclass <- best_model %>% predict_classes(x_test, verbose = 0) 
    
    res <- list(train_data = train_data, validation_data = validation_data, test_data = test_data)
    
  } else { res <- list(mod = cnn_model) }
  
  return(res)
}


run_id <- 1:10
cnn1 <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_cnn, seq_df = seq_df, target.nfo = target.nfo, cnn_model = model1, learn_rate = 0.001, epochs = 50, batch_size = 64, hyppar_id = 1)
cnn2 <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_cnn, seq_df = seq_df, target.nfo = target.nfo, cnn_model = model1, learn_rate = 0.0005, epochs = 50, batch_size = 64, hyppar_id = 2)
cnn3 <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_cnn, seq_df = seq_df, target.nfo = target.nfo, cnn_model = model2, learn_rate = 0.001, epochs = 50, batch_size = 64, hyppar_id = 3)
cnn4 <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_cnn, seq_df = seq_df, target.nfo = target.nfo, cnn_model = model2, learn_rate = 0.0005, epochs = 50, batch_size = 64, hyppar_id = 4)
cnn5 <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_cnn, seq_df = seq_df, target.nfo = target.nfo, cnn_model = model3, learn_rate = 0.001, epochs = 50, batch_size = 64, hyppar_id = 5)
cnn6 <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_cnn, seq_df = seq_df, target.nfo = target.nfo, cnn_model = model3, learn_rate = 0.0005, epochs = 50, batch_size = 64, hyppar_id = 6)

# separate out models (take space) and rest of results, and save
all_res <- c(cnn1, cnn2, cnn3, cnn4, cnn5, cnn6)
temp_train <- all_res %>% map_depth(1, "training_data")
temp_validation <- all_res %>% map_depth(1, "validation_data")
temp_test <- all_res %>% map_depth(1, "test_data")
cnn_preds <- list(train_data = temp_train, validation_data = temp_validation, test_data = temp_test)
rm(all_res, temp_train, temp_validation, temp_test)

save(training_ids, validation_ids, test_ids, cnn_preds, file = "cnn_predictions.Rdata")
