# for extensions see https://tensorflow.rstudio.com/guide/keras/examples/conv_lstm/
# https://blogs.rstudio.com/ai/posts/2020-12-17-torch-convlstm/

library(dplyr)
library(keras)
use_virtualenv("myenv", required = TRUE) 

#load('output/targets-w2022.RData')
#load("output/ids_traintestval.Rdata")
# load("output/lrcn-inputs-tinyex.Rdata")  # testing
#load("output/lrcn-inputs.Rdata")        # full run
load("output/lrcn-raw-inputs-frames-from-mid.Rdata")
load("output/ids_traintestval_1ttv.Rdata")

num_classes <- 2

# Input image dimensions
img_rows <- 24
img_cols <- 24
input_shape <- c(img_rows, img_cols, 1)

# Define model
model1 <- keras_model_sequential() %>% 
  #time_distributed(layer_dense(units = 128, kernel_initializer = initializer_zeros(), input_shape = c(40, 1))) %>%
  time_distributed(layer_conv_2d(input_shape = c(5,img_rows,img_cols,1), filters = 32, kernel_size = c(3,3), activation = "relu", 
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu",
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu", 
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_flatten()) %>%
  layer_lstm(64, return_sequences = FALSE, dropout = 0.25) %>%
  layer_dense(2, activation = "softmax")

model2 <- keras_model_sequential() %>% 
  #time_distributed(layer_dense(units = 128, kernel_initializer = initializer_zeros(), input_shape = c(40, 1))) %>%
  time_distributed(layer_conv_2d(input_shape = c(5,img_rows,img_cols,1), filters = 32, kernel_size = c(3,3), activation = "relu", 
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu",
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu", 
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu",
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu", 
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_flatten()) %>%
  layer_lstm(64, return_sequences = FALSE, dropout = 0.25) %>%
  layer_dense(2, activation = "softmax")

model3 <- keras_model_sequential() %>% 
  #time_distributed(layer_dense(units = 128, kernel_initializer = initializer_zeros(), input_shape = c(40, 1))) %>%
  time_distributed(layer_conv_2d(input_shape = c(5,img_rows,img_cols,1), filters = 32, kernel_size = c(3,3), activation = "relu", 
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu",
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu", 
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu",
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_conv_2d(filters = 64, kernel_size = c(3,3), activation = "relu", 
                                 padding = "same", kernel_initializer = 'glorot_uniform', kernel_regularizer=regularizer_l2(0.001))) %>%
  time_distributed(layer_flatten()) %>%
  layer_lstm(64, return_sequences = TRUE, dropout = 0.25) %>%
  layer_lstm(64, return_sequences = FALSE, dropout = 0.25) %>%
  layer_dense(2, activation = "softmax")

# Callbacks
es <- list(callback_early_stopping(monitor = "val_accuracy", 
                                   min_delta = 0.1, 
                                   patience = 10, 
                                   mode = "auto"),
           callback_model_checkpoint(filepath = paste0("best_lrcn.h5"), 
                                     monitor = "val_loss", 
                                     save_best_only = TRUE, 
                                     mode = "auto"))


lrcn_model <- model1

# Compile model
lrcn_model %>% compile(
  loss = "binary_crossentropy", 
  optimizer = optimizer_adam(learning_rate = 0.001),
  metrics = c('accuracy')
)

lrcn_model %>% fit(
  x_train, y_train,
  batch_size = 64,
  epochs = 5,
  validation_data = list(x_validation, y_validation),
  callbacks = es
)

best_model <- load_model_hdf5(paste0("best_lrcn.h5"))

# add frame-level predicted probabilities and classes
train_predprob <- best_model %>% predict(x_train, verbose = 0) %>% .[, 2]
train_predclass <- best_model %>% predict(x_train) %>% k_argmax() %>% as.numeric()

validation_predprob <- best_model %>% predict(x_validation, verbose = 0) %>% .[, 2]
validation_predclass <- best_model %>% predict(x_validation) %>% k_argmax() %>% as.numeric()

test_predprob <- best_model %>% predict(x_test, verbose = 0) %>% .[, 2]
test_predclass <- best_model %>% predict(x_test) %>% k_argmax() %>% as.numeric()

save(train_predprob, train_predclass, validation_predprob, validation_predclass,
     test_predprob, test_predclass, file = "output/lrcn_predictions.Rdata")
