
library(dplyr)
library(randomForest)
library(purrr)
library(stringr)
#library(unbalanced)

load("output/targets-w2022.RData")
load("output/ids_traintestval_1ttv.Rdata")
load("output/targets-w2022-framenfo.Rdata")
load("output/targetnfo-only.Rdata")

frame.nfo$class <- as.factor(frame.nfo$class) # frame-level
target.nfo$class <- as.factor(target.nfo$class) # target-level

# variables for frame-level
rf_vars <- c("dist", "radi", "ddist", "dradi", "dvec",
             "area", "perim", "length", "pa.ratio", "la.ratio", "lp.ratio",
             "pixmin", "pixmax", "pixmed", "pixave", "pixsd", "pixrng")
# variables for target-level
varnames <- names(target.nfo)
rf_vars_tgt <- varnames[str_detect(varnames, paste(rf_vars, collapse = "|"))]

# remove variables with zero variance or high correlations with existing vars
vars_already_in <- c()
for(i in rf_vars_tgt){
  any_na <- sum(is.na(target.nfo[,i]))
  zero_var <- ifelse(sd(target.nfo[,i], na.rm = TRUE) < 0.00001, TRUE, FALSE)
  corr_with_already_in <- ifelse(max(cor(target.nfo[,i], target.nfo[,vars_already_in], use = "pairwise.complete.obs")) > 0.9, TRUE, FALSE)
  if(!zero_var & !corr_with_already_in){  vars_already_in <- c(vars_already_in, i) }
}
rf_vars_tgt <- vars_already_in

# add step variable so can use run_rf function as is (gets dropped anyway) for target-level model
target.nfo <- target.nfo %>% mutate(step = 999)

#
#frame.nfo$id <- frame.nfo$targ_id
#
rf_vars <- rf_vars[rf_vars != "lp.ratio"]
frame.nfo.com <- frame.nfo[, c("targ_id", "class", "step", rf_vars, "traveldir")]
frame.nfo.com <- frame.nfo.com[complete.cases(frame.nfo.com), ]

frame.nfo <- frame.nfo.com %>% mutate_at(4:19, scale) %>% mutate_at(4:19, as.numeric)
summary(frame.nfo)

# frame.nfo = frame.nfo
# training_ids = training_ids[[1]]
# validation_ids = validation_ids[[1]]
# test_ids = test_ids[[1]]
# rf_vars = rf_vars
# mtry = 4
# run_id = 0
# hyppar_id = 0
# SMOTE = 0
# mod_only = FALSE

run_rf <- function(frame.nfo, training_ids, validation_ids, test_ids, rf_vars, mtry, run_id = 0, hyppar_id = 0, SMOTE = 0, mod_only = FALSE){
  
  train_data <- frame.nfo[frame.nfo$targ_id %in% training_ids, ] %>% dplyr::select(targ_id, step, class, rf_vars) %>% mutate(run_id = run_id, hyppar_id = hyppar_id)
  validation_data <- frame.nfo[frame.nfo$targ_id %in% validation_ids, ] %>% dplyr::select(targ_id, step, class, rf_vars) %>% mutate(run_id = run_id, hyppar_id = hyppar_id)
  test_data <- frame.nfo[frame.nfo$targ_id %in% test_ids, ] %>% dplyr::select(targ_id, step, class, rf_vars) %>% mutate(run_id = run_id, hyppar_id = hyppar_id)
  
  # missing values for some variables at step 1, discard
  train_data <- train_data %>% filter(step != 1) 
  validation_data <- validation_data %>% filter(step != 1) 
  test_data <- test_data %>% filter(step != 1) 
  
  # SMOTE the training fold if desired
  if(SMOTE == 1){
    temp_train_data <- train_data %>% dplyr::select(-targ_id, -step, -run_id)
    train_data <- ubSMOTE(X = temp_train_data %>% dplyr::select(-class), Y = temp_train_data$class, k = 5, perc.over = 100, perc.under = 200)
    train_data <- cbind(class = train_data$Y, train_data$X)
  }
  
  f <- formula(paste("class ~",paste(rf_vars, collapse="+")))
  
  rf <- randomForest(f, 
                     data = train_data, 
                     mtry = mtry, # default for classification is sqrt of no of parameters
                     importance = FALSE, localImp = FALSE, 
                     ntree = 500)
  
  if(!mod_only){
    
    # add frame-level predicted probabilities and classes
    train_data$predprob <- predict(rf, type = "prob")[,2]
    train_data$predclass <- predict(rf, type = "class")
    
    validation_data$predprob <- predict(rf, type = "prob", newdata = validation_data)[,2]
    validation_data$predclass <- predict(rf, type = "class", newdata = validation_data)
    
    test_data$predprob <- predict(rf, type = "prob", newdata = test_data)[,2]
    test_data$predclass <- predict(rf, type = "class", newdata = test_data)
    
    res <- list(mod = rf, train_data = train_data, validation_data = validation_data, test_data = test_data)
    
  } else { res <- list(mod = rf) }
  
  return(res)
  
}

table(frame.nfo.com$traveldir, frame.nfo.com$class)

# frame-level
r1 <- run_rf(frame.nfo.com, training_ids_q1234, validation_ids_q1234, test_ids_q1234, rf_vars = rf_vars, mtry = 4)
xx <- table(r1$test_data$class, r1$test_data$predclass)
xx
sum(diag(xx))/sum(xx)

r1d <- run_rf(frame.nfo.com %>% filter(traveldir == "downriver"), training_ids_q1234_down, validation_ids_q1234_down, test_ids_q1234_down, rf_vars = rf_vars, mtry = 4)
xx <- table(r1d$test_data$class, r1d$test_data$predclass)
xx
sum(diag(xx))/sum(xx)

r1u <- run_rf(frame.nfo.com %>% filter(traveldir == "upriver"), training_ids_q1234_up, validation_ids_q1234_up, test_ids_q1234_up, rf_vars = rf_vars, mtry = 4)
xx <- table(r1u$test_data$class, r1u$test_data$predclass)
xx
sum(diag(xx))/sum(xx)

# target-level
yy<-r1$test_data %>% group_by(targ_id) %>% 
  summarize(predclass = mean(predclass == "1"), class = first(class)) %>%
  ungroup() %>% mutate(predclass = ifelse(predclass > 0.5, "1", "0"))
yy <- yy %>% left_join(target.nfo %>% dplyr::select(targ_id, quality, site, traveldir))

yy %>% group_by(quality) %>% summarize(acc = sum(predclass == class)/n(), n = n())
yy %>% group_by(traveldir) %>% summarize(acc = sum(predclass == class)/n(), n = n())
yy2 <- yy %>% mutate(site2 = case_when(
  site == "Ardoe" ~ "Ardoe",
  site == "watering hole" ~ "watering hole",
  site == "Waterside Fm." ~ "Waterside Fm.",
  TRUE ~ "Other"
)) %>%
  group_by(site2) %>% summarize(acc = sum(predclass == class)/n(), n = n())
yy2
table(yy$site)
xx<-table(yy$predclass, yy$class)
xx
sum(diag(xx))/sum(xx)


# target-level downriver
yy<-r1d$test_data %>% group_by(targ_id) %>% 
  summarize(predclass = mean(predclass == "1"), class = first(class)) %>%
  ungroup() %>% mutate(predclass = ifelse(predclass > 0.45, "1", "0"))
yy <- yy %>% left_join(target.nfo %>% dplyr::select(targ_id, quality, site, traveldir))

yy %>% group_by(quality) %>% summarize(acc = sum(predclass == class)/n(), n = n())
yy %>% group_by(traveldir) %>% summarize(acc = sum(predclass == class)/n(), n = n())

# target-level upnriver
yy<-r1u$test_data %>% group_by(targ_id) %>% 
  summarize(predclass = mean(predclass == "1"), class = first(class)) %>%
  ungroup() %>% mutate(predclass = ifelse(predclass > 0.5, "1", "0"))
yy <- yy %>% left_join(target.nfo %>% dplyr::select(targ_id, quality, site, traveldir))

yy %>% group_by(quality) %>% summarize(acc = sum(predclass == class)/n(), n = n())
yy %>% group_by(traveldir) %>% summarize(acc = sum(predclass == class)/n(), n = n())


yy2 <- yy %>% mutate(site2 = case_when(
  site == "Ardoe" ~ "Ardoe",
  site == "watering hole" ~ "watering hole",
  site == "Waterside Fm." ~ "Waterside Fm.",
  TRUE ~ "Other"
)) %>%
  group_by(site2) %>% summarize(acc = sum(predclass == class)/n(), n = n())
yy2
table(yy$site)
xx<-table(yy$predclass, yy$class)
xx
sum(diag(xx))/sum(xx)



# target-level
yy<-r1$test_data %>% group_by(targ_id) %>% 
  slice_head(n = 5) %>%
  summarize(predclass = mean(predclass == "1"), class = first(class)) %>%
  ungroup() %>% mutate(predclass = ifelse(predclass > 0.5, "1", "0"))
yy %>% summarize(acc = sum(predclass == class)/n(), n = n(),
                 fpr = sum(predclass == 1 & class == 0)/sum(predclass == 1),
                 fnr = sum(predclass == 0 & class == 1)/sum(class == 1))

yy <- yy %>% left_join(target.nfo %>% dplyr::select(targ_id, quality, site))


# frame-level
r2 <- run_rf(frame.nfo.com %>% filter(traveldir == "downriver"), training_ids, validation_ids, test_ids, rf_vars = rf_vars, mtry = 4)
xx <- table(r1$test_data$class, r1$test_data$predclass)
sum(diag(xx))/sum(xx)

yy<-r1$test_data %>% group_by(targ_id) %>% 
  summarize(predclass = mean(predclass == "1"), class = first(class)) %>%
  ungroup() %>% mutate(predclass = ifelse(predclass > 0.7, "1", "0"))
xx<-table(yy$predclass, yy$class)
xx
sum(diag(xx))/sum(xx)

r1$mod$importance
varImpPlot(r1$mod)



# most NB vars only
r2 <- run_rf(frame.nfo.com, training_ids, validation_ids, test_ids, rf_vars = c("dist", "radi", "pa.ratio"), mtry = 4)
xx <- table(r2$test_data$class, r2$test_data$predclass)
xx
sum(diag(xx))/sum(xx)

yy<-r2$test_data %>% group_by(targ_id) %>% 
  summarize(predclass = mean(predclass == "1"), class = first(class)) %>%
  ungroup() %>% mutate(predclass = ifelse(predclass > 0.25, "1", "0"))
xx<-table(yy$predclass, yy$class)
xx
sum(diag(xx))/sum(xx)

# run once to test
# frame-level
r1 <- run_rf(frame.nfo.com, training_ids[[1]], validation_ids[[1]], test_ids[[1]], rf_vars = rf_vars, mtry = 4)
# target-level
# r1 <- run_rf(target.nfo, training_ids[[1]], validation_ids[[1]], test_ids[[1]], rf_vars = rf_vars_tgt, mtry = 4)
# target-level with SMOTE
r1 <- run_rf(target.nfo, training_ids[[1]], validation_ids[[1]], test_ids[[1]], rf_vars = rf_vars_tgt, mtry = 4, SMOTE = 1)

## frame-level model

frame.nfo.red <- frame.nfo %>% filter(site %in% c("Waterside Fm."))

# run all
run_id <- 1:length(training_ids)
rf_mtry4_res <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_rf, frame.nfo = frame.nfo.red, rf_vars = rf_vars, mtry = 4, hyppar_id = 1)
rf_mtry8_res <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_rf, frame.nfo = frame.nfo.red, rf_vars = rf_vars, mtry = 8, hyppar_id = 2)
# <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_rf, frame.nfo = frame.nfo, rf_vars = rf_vars, mtry = 16, hyppar_id = 3)

# separate out models (take space) and rest of results, and save
all_res <- c(rf_mtry4_res, rf_mtry8_res)
rf_mod <- all_res %>% map_depth(1, "mod")
temp_train <- all_res %>% map_depth(1, "train_data")
temp_validation <- all_res %>% map_depth(1, "validation_data")
temp_test <- all_res %>% map_depth(1, "test_data")
rf_preds <- list(train_data = temp_train, validation_data = temp_validation, test_data = temp_test)
rm(all_res, temp_train, temp_validation, temp_test)

save(training_ids, validation_ids, test_ids, rf_preds, file = "output/rf_predictions_waterside.Rdata")
#save(training_ids, validation_ids, test_ids, rf_mod, file = "output/rf_models.Rdata")

# rf_mtry4_res <- list()
# for(i in 1:k){rf_mtry4_res[[i]] <- run_rf(frame.nfo, training_ids[[i]], validation_ids[[i]], test_ids[[i]], rf_vars = rf_vars, mtry = 4)}

## target-level model

# run all
run_id <- 1:length(training_ids)
mtry_base <- round(sqrt(length(rf_vars_tgt)))
rf_mtry1_res <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_rf, frame.nfo = target.nfo, rf_vars = rf_vars_tgt, mtry = mtry_base, hyppar_id = 1)
rf_mtry2_res <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_rf, frame.nfo = target.nfo, rf_vars = rf_vars_tgt, mtry = mtry_base * 2, hyppar_id = 2)
rf_mtry3_res <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_rf, frame.nfo = target.nfo, rf_vars = rf_vars_tgt, mtry = mtry_base * 3, hyppar_id = 3)

# separate out models (take space) and rest of results, and save
all_res <- c(rf_mtry1_res, rf_mtry2_res, rf_mtry3_res)
rf_mod <- all_res %>% map_depth(1, "mod")
temp_train <- all_res %>% map_depth(1, "train_data")
temp_validation <- all_res %>% map_depth(1, "validation_data")
temp_test <- all_res %>% map_depth(1, "test_data")
rf_preds <- list(train_data = temp_train, validation_data = temp_validation, test_data = temp_test)
rm(all_res, temp_train, temp_validation, temp_test)

save(training_ids, validation_ids, test_ids, rf_preds, file = "output/runs/rf_tgtlevel_predictions.Rdata")

## target-level model with SMOTE

# run all
run_id <- 1:length(training_ids)
mtry_base <- round(sqrt(length(rf_vars_tgt)))
rf_mtry1_res <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_rf, frame.nfo = target.nfo, rf_vars = rf_vars_tgt, mtry = mtry_base, hyppar_id = 1, SMOTE = 1)
rf_mtry2_res <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_rf, frame.nfo = target.nfo, rf_vars = rf_vars_tgt, mtry = mtry_base * 2, hyppar_id = 2, SMOTE = 1)
rf_mtry3_res <- pmap(list(training_ids, validation_ids, test_ids, run_id), .f = run_rf, frame.nfo = target.nfo, rf_vars = rf_vars_tgt, mtry = mtry_base * 3, hyppar_id = 3, SMOTE = 1)

# separate out models (take space) and rest of results, and save
all_res <- c(rf_mtry1_res, rf_mtry2_res, rf_mtry3_res)
rf_mod <- all_res %>% map_depth(1, "mod")
temp_train <- all_res %>% map_depth(1, "train_data")
temp_validation <- all_res %>% map_depth(1, "validation_data")
temp_test <- all_res %>% map_depth(1, "test_data")
rf_preds <- list(train_data = temp_train, validation_data = temp_validation, test_data = temp_test)
rm(all_res, temp_train, temp_validation, temp_test)

save(training_ids, validation_ids, test_ids, rf_preds, file = "output/runs/rf_tgtlevelSMOTE_predictions.Rdata")


####

yy<-test_data %>% group_by(targ_id) %>% summarize(predclass2 = mean(predclass == "1"),
                                             class = first(class)) %>%
  ungroup() %>% mutate(predclass = ifelse(predclass2 > 0.65, "1", "0"))
xx<-table(yy$predclass, yy$class)
xx
sum(diag(xx))/sum(xx)


