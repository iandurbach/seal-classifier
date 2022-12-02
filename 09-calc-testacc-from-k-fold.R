
library(dplyr)
library(purrr)

options(dplyr.summarise.inform = FALSE)

find_best_hyperpar = function(res, frame_threshold_df, vote_threshold = 0.5, prob_threshold = 0.5, agg = "vote", return_per_fold = FALSE, maxvote = 1){
  
  # combine datasets across runs (have run_id to identify)
  #x <- res %>% map_depth(1, dataname) %>% map_df(bind_rows)
  x <- res %>% map_df(bind_rows)
  
  # add in best frame_threshold for each hyppar_id and use that to compute predicted class
  x <- x %>% left_join(frame_threshold_df, by = "hyppar_id") %>%
    mutate(predclass = ifelse(predprob > frame_threshold, 1, 0))
  
  # add percentage of frames predicted as seal
  # maxvote == 1 unless predclass is the sum of votes across multiple models; then maxvote = n models
  x <- x %>% group_by(hyppar_id, run_id, id) %>% 
    mutate(percseal = sum(predclass)/(n() * maxvote)) %>% ungroup() %>%
    mutate(class = as.character(class))
  
  # assess individual level accuracy by vote > threshold
  x <- x %>% group_by(hyppar_id, run_id, id) %>%
    summarize(class = first(class),
              percseal = sum(predclass)/(n()*maxvote),
              meanprob = (prod(predprob))^(1/n())) %>% ungroup() %>%
    mutate(tgtpred_vote = ifelse(percseal > vote_threshold, "1", "0"),
           tgtpred_prob = ifelse(meanprob > prob_threshold, "1", "0")) 
  
  acc_per_fold <- x %>% group_by(hyppar_id, run_id) %>% 
    summarize(acc_vote = sum(class == tgtpred_vote)/n(),
              acc_prob = sum(class == tgtpred_prob)/n()) %>% ungroup()
  
  acc <- acc_per_fold %>% group_by(hyppar_id) %>% 
    summarize(acc_vote = mean(acc_vote), acc_prob = mean(acc_prob)) %>% ungroup()
  
  if(agg == "vote"){
    res <- data.frame(hyppar_id = acc$hyppar_id[which.max(acc$acc_vote)], 
                      threshold = vote_threshold, acc = max(acc$acc_vote))
  } else if (agg == "prob") {
    res <- data.frame(hyppar_id = acc$hyppar_id[which.max(acc$acc_prob)], 
                      threshold = prob_threshold, acc = max(acc$acc_prob)) 
  } else { stop("Check agg argument") }
  
  if(return_per_fold == TRUE){ res <- acc }
  
  res <- res %>% left_join(frame_threshold_df, by = "hyppar_id")
  
  return(res)
  
}

get_frame_accuracy = function(res, frame_threshold = 0.5){
  
  # combine datasets across runs (have run_id to identify)
  #x <- res %>% map_depth(1, dataname) %>% map_df(bind_rows)
  x <- res %>% map_df(bind_rows)
  
  x <- x %>% mutate(predclass = ifelse(predprob > frame_threshold, 1, 0))
  
  x <- x %>% group_by(hyppar_id, run_id) %>% summarize(acc_frame = sum(class == predclass)/n()) %>% ungroup()
  
  x <- x %>% group_by(hyppar_id) %>% summarize(acc_frame = mean(acc_frame)) %>% ungroup() %>% mutate(frame_threshold = frame_threshold)
  
  return(x)
  
}

get_accuracy = function(res, hp, frame_threshold, threshold, agg = "vote", maxvote = 1, use_existing_predclass = FALSE){
  
  # combine datasets across runs (have run_id to identify)
  #x <- res %>% map_depth(1, dataname) %>% map_df(bind_rows)
  x <- res %>% map_df(bind_rows)
  
  x <- x %>% dplyr::filter(hyppar_id == hp, threshold == threshold)
  
  if(!use_existing_predclass){ x <- x %>% mutate(predclass = ifelse(predprob > frame_threshold, 1, 0)) }
  
  # add percentage of frames predicted as seal
  x <- x %>% group_by(hyppar_id, run_id, id) %>% 
    mutate(percseal = sum(predclass)/(n()*maxvote)) %>% ungroup() %>%
    mutate(class = as.character(class))
  
  # assess individual level accuracy by vote > threshold
  x <- x %>% group_by(hyppar_id, run_id, id) %>%
    summarize(class = first(class),
              percseal = sum(predclass)/(n()*maxvote),
              meanprob = (prod(predprob))^(1/n())) %>% ungroup() %>%
    mutate(tgtpred_vote = ifelse(percseal > threshold, "1", "0"),
           tgtpred_prob = ifelse(meanprob > threshold, "1", "0")) 
  
  acc_per_fold <- x %>% group_by(hyppar_id, run_id) %>% 
    summarize(acc_vote = sum(class == tgtpred_vote)/n(),
              precision_vote = sum(class == "1" & tgtpred_vote == "1") / sum(class == "1"),
              recall_vote = sum(class == "1" & tgtpred_vote == "1") / sum(tgtpred_vote == "1"),
              f1_vote = (2 * precision_vote * recall_vote) / (precision_vote + recall_vote),
              fpr_vote = sum(class == "0" & tgtpred_vote == "1") / sum(class == "0"),
              fnr_vote = sum(class == "1" & tgtpred_vote == "0") / sum(class == "1"),
              fomr_vote = sum(class == "1" & tgtpred_vote == "0") / sum(tgtpred_vote == "0"),
              acc_prob = sum(class == tgtpred_prob)/n(),
              precision_prob = sum(class == "1" & tgtpred_prob == "1") / sum(class == "1"),
              recall_prob = sum(class == "1" & tgtpred_prob == "1") / sum(tgtpred_prob == "1"),
              f1_prob = (2 * precision_prob * recall_prob) / (precision_prob + recall_prob),
              fpr_prob = sum(class == "0" & tgtpred_prob == "1") / sum(class == "0"),
              fnr_prob = sum(class == "1" & tgtpred_prob == "0") / sum(class == "1"),
              fomr_prob = sum(class == "1" & tgtpred_prob == "0") / sum(tgtpred_prob == "0")) %>% ungroup()
  
  acc <- acc_per_fold %>% group_by(hyppar_id) %>% 
    summarize_all(mean, na.rm = TRUE) %>% ungroup() %>% dplyr::select(-run_id)
  
  if(agg == "vote"){
    res <- data.frame(hyppar_id = hp, frame_threshold = frame_threshold, threshold = threshold, 
                      acc = acc$acc_vote, precision = acc$precision_vote, recall = acc$recall_vote, f1 = acc$f1_vote,
                      fpr = acc$fpr_vote, fnr = acc$fnr_vote, fomr = acc$fomr_vote)
  } else if (agg == "prob") {
    res <- data.frame(hyppar_id = hp, frame_threshold = frame_threshold, threshold = threshold,
                      acc = acc$acc_prob, precision = acc$precision_prob, recall = acc$recall_prob, f1 = acc$f1_prob,
                      fpr = acc$fpr_prob, fnr = acc$fnr_prob, fomr = acc$fomr_prob)
  } else { stop("Check agg argument") }
  
  return(res)
  
}

make_var_numeric <- function(x){
  xt <- x
  if (is.character(x)) { xt <- as.numeric(x) } else if (is.factor(x)) { xt <- as.numeric(as.character(x)) }
  return(xt)
}

get_all_accuracies = function(predictions, method, opt_predclass = FALSE, nmod = 1, noopt_hyppar = 1, use_existing_predclass = FALSE){
  
  # convert any factor predclass variables to numeric
  predictions$validation_data <- map(predictions$validation_data, mutate, across(contains("predclass"), make_var_numeric))
  predictions$test_data <- map(predictions$test_data, mutate, across(contains("predclass"), make_var_numeric))
  
  # optimize frame threshold
  ft_df_opt <- map_dfr(seq(0.3,0.7,by=0.05), get_frame_accuracy, res = predictions$validation_data)
  ft_df_opt <- ft_df_opt %>% group_by(hyppar_id) %>% arrange(desc(acc_frame)) %>% slice(1) %>% ungroup()
  ft_df_noopt <- ft_df_opt %>% mutate(frame_threshold = 0.5)
  if(opt_predclass){ ft_df <- ft_df_opt } else { ft_df <- ft_df_noopt }
  
  # optimize vote threshold
  vote_hp <- map_dfr(seq(0.3,0.7,by=0.05), find_best_hyperpar, res = predictions$validation_data, frame_threshold_df = ft_df, prob_threshold = 0.5, agg = "vote", maxvote = nmod)
  vote_hp <- vote_hp %>% arrange(desc(acc)) %>% slice(1)
  # get test accuracy at optimized vote threshold
  ft <- as.numeric(ft_df[ft_df$hyppar_id == vote_hp$hyppar_id, "frame_threshold"])
  vote_testacc <- get_accuracy(predictions$test_data, hp = vote_hp$hyppar_id, frame_threshold = ft, threshold = vote_hp$threshold, agg = "vote", maxvote = nmod, use_existing_predclass = use_existing_predclass)
  # no optimization
  vote_testacc_noopt <- get_accuracy(predictions$test_data, hp = noopt_hyppar, frame_threshold = 0.5, threshold = 0.5, agg = "vote", use_existing_predclass = use_existing_predclass)
  
  # optimize prob threshold (frame_threshold doesn't matter here since working with frame predprobs)
  prob_hp <- map_dfr(seq(0.3,0.7,by=0.05), find_best_hyperpar, res = predictions$validation_data, frame_threshold_df = ft_df, vote_threshold = 0.5, agg = "prob", maxvote = nmod)
  prob_hp <- prob_hp %>% arrange(desc(acc)) %>% slice(1)
  # get test accuracy at optimized vote threshold
  prob_testacc <- get_accuracy(predictions$test_data, hp = prob_hp$hyppar_id, frame_threshold = 0.5, threshold = prob_hp$threshold, agg = "prob", maxvote = nmod, use_existing_predclass = use_existing_predclass)
  # no optimizati
  prob_testacc_noopt <- get_accuracy(predictions$test_data, hp = noopt_hyppar, frame_threshold = 0.5, threshold = 0.5, agg = "prob", maxvote = nmod, use_existing_predclass = use_existing_predclass)
  
  df1 <- data.frame(method = method, 
                    aggregation = c("vote", "prob", "vote", "prob"),
                    optimized = c("optimized", "optimized", "no opt.", "no opt."))
  df <- cbind(df1, rbind(vote_testacc, prob_testacc, vote_testacc_noopt, prob_testacc_noopt))
  
  return(df)
  
}


# Model 2a
load("output/rf_predictions_waterside.Rdata")
acc_rf_tgtlevel <- get_all_accuracies(predictions = rf_preds, "rf_tgt", noopt_hyppar = 2)
acc_rf_tgtlevel

# ad hoc analysis, fix up later (2022-09-20)

rf_preds

xt <-  rf_preds$test_data[[1]] %>% dplyr::select(id, class, run_id, hyppar_id, predprob, predclass)
xt <- xt %>% left_join(target.nfo %>% dplyr::select(targ_id, quality.x, site), by = c("id" = "targ_id"))
for(i in 2:20){
 xtt <- rf_preds$test_data[[i]] %>% dplyr::select(id, class, run_id, hyppar_id, predprob, predclass)
 xtt <- xtt %>% left_join(target.nfo %>% dplyr::select(targ_id, quality.x, site), by = c("id" = "targ_id"))
 xt <- rbind(xt, xtt)
}

xt3 <- xt %>% group_by(run_id, hyppar_id, id, site, quality.x, class) %>% 
  summarize(predseal = sum(predclass == "1"), prednoseal = sum(predclass == "0")) %>%
  mutate(predclass = ifelse(predseal > prednoseal, "1", "0")) %>% ungroup()

site_acc <- xt3 %>% filter(hyppar_id == 2) %>% group_by(hyppar_id, site) %>% summarize(accuracy = sum(class == predclass)/n())
qual_acc <- xt3 %>% filter(hyppar_id == 2) %>% group_by(hyppar_id, quality.x) %>% summarize(accuracy = sum(class == predclass)/n())
waterside_qual_acc <- xt3 %>% filter(site == "Waterside Fm.") %>% group_by(hyppar_id, quality.x) %>% summarize(accuracy = sum(class == predclass)/n())

res1 <- xt3 %>% filter(hyppar_id == 2) %>% group_by(hyppar_id, site, quality.x) %>% summarize(accuracy = sum(class == predclass)/n()) %>% ungroup()

library(ggplot2)
res1 %>% filter(site %in% c("Ardoe", "Friars", "watering hole", "Waterside Fm."),
                quality.x %in% c("1", "2", "3", "4")) %>%
  mutate(quality = as.numeric(quality.x)) %>%
  ggplot(aes(x = site, y = accuracy, fill = as.factor(quality))) + geom_bar(stat = "identity",position = 'dodge') +
  scale_fill_viridis_d(name = "quality") + theme_classic()

# Model 2b
load("output/runs/rf_tgtlevelSMOTE_predictions.Rdata")
acc_rf_tgtlevelSMOTE <- get_all_accuracies(predictions = rf_preds, "rf_tgt_SMOTE", noopt_hyppar = 2)
acc_rf_tgtlevelSMOTE

# Model 3
load("output/runs/svm_predictions.Rdata")
acc_svm <- get_all_accuracies(predictions = svm_preds, "svm", noopt_hyppar = 2, opt_predclass = FALSE)
acc_svm

# Model 4
load("output/runs/rf_predictions.Rdata")
acc_rf <- get_all_accuracies(predictions = rf_preds, "rf", noopt_hyppar = 1, opt_predclass = FALSE)
acc_rf

# Model 5
load("output/runs/mlp_predictions.Rdata")
acc_mlp <- get_all_accuracies(predictions = mlp_preds, "mlp", noopt_hyppar = 1, opt_predclass = FALSE)
acc_mlp

# Model 6
load("output/runs/cnn_predictions.Rdata")
acc_cnn <- get_all_accuracies(predictions = cnn_preds, "cnn", noopt_hyppar = 5, opt_predclass = FALSE)
acc_cnn

# Model 7
# for ensemble models need to specify how to calculate a joint predprob and predclass
load("output/runs/ensemble_cnn_rf_predictions.Rdata")
ens_preds$validation_data[[1]] <- ens_preds$validation_data[[1]] %>% mutate(predclass = rowSums(across(contains("predclass")))) %>% dplyr::select(-m1_predclass, -m2_predclass)
ens_preds$test_data[[1]] <- ens_preds$test_data[[1]] %>% mutate(predclass = rowSums(across(contains("predclass")))) %>% dplyr::select(-m1_predclass, -m2_predclass)
ens_preds$validation_data[[1]] <- ens_preds$validation_data[[1]] %>% mutate(predprob = sqrt(m1_predprob * m2_predprob)) %>% dplyr::select(-m1_predprob, -m2_predprob)
ens_preds$test_data[[1]] <- ens_preds$test_data[[1]] %>% mutate(predprob = sqrt(m1_predprob * m2_predprob)) #%>% dplyr::select(-m1_predprob, -m2_predprob)
#ens_preds$validation_data[[1]] <- ens_preds$validation_data[[1]] %>% mutate(predclass = ifelse(predprob > 0.5, 1, 0))
#ens_preds$test_data[[1]] <- ens_preds$test_data[[1]] %>% mutate(predclass = ifelse(predprob > 0.5, 1, 0))
acc_ens <- get_all_accuracies(predictions = ens_preds, "ens_cnn_rf1", nmod = 2, noopt_hyppar = 5)
acc_ens

# Model 8
load("output/runs/cnn_mlp_predictions.Rdata")
acc_cnn_mlp <- get_all_accuracies(predictions = cnn_mlp_preds, "cnn_mlp", noopt_hyppar = 1)
acc_cnn_mlp

# Model 9
load("output/runs/lrcn_predictions.Rdata")
acc_lrcn <- get_all_accuracies(predictions = lrcn_preds, "lrcn", noopt_hyppar = 1, opt_predclass = FALSE)
acc_lrcn


all_res <- rbind(acc_svm_tgtlevel, acc_rf_tgtlevel, acc_svm, acc_rf, acc_mlp, acc_cnn, acc_ens, acc_cnn_mlp, acc_lrcn)

all_res_opt <- all_res %>% filter(optimized == "optimized") %>% group_by(method) %>% arrange(desc(acc)) %>% slice(1) %>% ungroup()
all_res_noopt <- all_res %>% filter(optimized == "no opt.") %>% group_by(method) %>% arrange(desc(acc)) %>% slice(1) %>% ungroup()

write.csv(all_res_opt, file = "output/runs/all_accuracy_results_opt.csv")
write.csv(all_res_noopt, file = "output/runs/all_accuracy_results_noopt.csv")

# # # for testing, not run
# xx <- find_best_hyperpar(rf_preds$validation_data, vote_threshold = 0.5, prob_threshold = 0.5, agg = "vote")
# xx$hyppar_id <- 2
# xx2 <- get_accuracy(rf_preds$test_data, hp = xx$hyppar_id, threshold = xx$threshold, agg = "vote")
# xx2
# 
# xx
# xx2


#####
library(ggplot2)
library(viridis)
library(ggthemes)
library(RColorBrewer)

display.brewer.all(n = NULL, type = "qual", select = NULL, colorblindFriendly = TRUE)

# load("output/runs/svm_tgtlevel_predictions.Rdata") #m1
# load("output/runs/svm_tgtlevelSMOTE_predictions.Rdata") #m1
# load("output/runs/rf_tgtlevel_predictions.Rdata") #m2
# load("output/runs/rf_tgtlevelSMOTE_predictions.Rdata") #m2
load("output/runs/svm_predictions.Rdata") #m3
# load("output/runs/rf_predictions.Rdata") #m4
# load("output/runs/mlp_predictions.Rdata") #m5
load("output/runs/cnn_predictions.Rdata") #m6
load("output/runs/ensemble_cnn_rf_predictions.Rdata") #m7
ens_preds$validation_data[[1]] <- ens_preds$validation_data[[1]] %>% mutate(predclass = rowSums(across(contains("predclass")))) %>% dplyr::select(-m1_predclass, -m2_predclass)
ens_preds$test_data[[1]] <- ens_preds$test_data[[1]] %>% mutate(predclass = rowSums(across(contains("predclass")))) %>% dplyr::select(-m1_predclass, -m2_predclass)
ens_preds$validation_data[[1]] <- ens_preds$validation_data[[1]] %>% mutate(predprob = sqrt(m1_predprob * m2_predprob)) %>% dplyr::select(-m1_predprob, -m2_predprob)
ens_preds$test_data[[1]] <- ens_preds$test_data[[1]] %>% mutate(predprob = sqrt(m1_predprob * m2_predprob)) #%>% dplyr::select(-m1_predprob, -m2_predprob)
load("output/runs/cnn_mlp_predictions.Rdata") #m8
load("output/runs/lrcn_predictions.Rdata") #m9

fpr_fnr_m3 <- map_dfr(seq(0, 1, by = 0.05), get_accuracy, res = svm_preds$test_data, hp = 2, frame_threshold = 0.5, agg = "prob", maxvote = 1) %>% mutate(model = "SVM")
fpr_fnr_m6 <- map_dfr(seq(0, 1, by = 0.05), get_accuracy, res = cnn_preds$test_data, hp = 5, frame_threshold = 0.5, agg = "vote", maxvote = 1) %>% mutate(model = "CNN")
fpr_fnr_m7 <- map_dfr(seq(0, 1, by = 0.05), get_accuracy, res = ens_preds$test_data, hp = 5, frame_threshold = 0.5, agg = "vote", maxvote = 1) %>% mutate(model = "CNN + RF")
fpr_fnr_m8 <- map_dfr(seq(0, 1, by = 0.05), get_accuracy, res = cnn_mlp_preds$test_data, hp = 1, frame_threshold = 0.5, agg = "vote", maxvote = 1) %>% mutate(model = "Mixed-input CNN")
fpr_fnr_m9 <- map_dfr(seq(0, 1, by = 0.05), get_accuracy, res = lrcn_preds$test_data, hp = 1, frame_threshold = 0.5, agg = "vote", maxvote = 1) %>% mutate(model = "CNN + RNN")

fpr_fnr <- rbind(fpr_fnr_m3, fpr_fnr_m6, fpr_fnr_m7, fpr_fnr_m8, fpr_fnr_m9)

fpr_fnr <- fpr_fnr %>% mutate(acc_alpha = pmax((acc - 0.75) / (max(acc) - 0.75),0))

fpr_fnr$model <- factor(fpr_fnr$model, levels = c("CNN + RNN", "CNN + RF", "Mixed-input CNN", "CNN", "SVM"))

p1 <- fpr_fnr %>% ggplot(aes(x = fnr*100, y = fpr*100, alpha = acc_alpha, color = model)) + geom_path() + geom_point() +
  scale_colour_colorblind() + scale_alpha(guide = "none") + 
  xlab("False negative rate (%)") + ylab("False positive rate (%)") +
  coord_equal() +
  theme_bw(base_size = 12) + theme(panel.grid = element_blank()) +
  theme(legend.position="bottom", legend.title = element_blank(),
        legend.text=element_text(size=10), legend.spacing.y = unit(0, 'cm')) + 
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

p1 <- p1 + annotate("text", x = 3.5, y = 3, label = "B") + 
  annotate("text", x = 3.5, y = 8, label = "C") +
  annotate("text", x = 11.5, y = 6, label = "A")

ggsave("docs/fpr-fnr.png", p1, width=5.5, height=5, dpi = 300)

