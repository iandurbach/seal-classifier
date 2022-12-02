library(dplyr)
library(purrr)
library(stringr)

load("output/targets-w2022.RData")

# # there are some duplicated rows in here, check this
# which(duplicated(target.nfo$targ_id))
# target.nfo$targ_id[which(duplicated(target.nfo$targ_id))]
# target.nfo <- target.nfo %>% filter(!duplicated(targ_id))

# set number of previous frames to include here (including current one)
L <- 5

n_mat_per_target <- target.pix %>% map_depth(1, length) %>% unlist() 
include_target <- n_mat_per_target >= L

# remove targets with fewer than L frames
target.pix <- target.pix[include_target]
include_ids <- names(target.pix)
include_ids <- str_extract_all(include_ids,"\\d.*\\d")
target.nfo <- target.nfo %>% filter(targ_id %in% include_ids)

# make array with TIMs
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
  group_by(s_id) %>% mutate(f_id = row_number(), frames = n()) %>% 
  ungroup() %>%
  left_join(target.nfo %>% dplyr::select(targ_id, class) %>% mutate(s_id = row_number()), by = "s_id")

# how many sequences of length L can you make from each sequence
n_Lseq_per_target <- n_mat_per_target - L + 1
n_Lseq <- data.frame(s_id = 1:length(n_Lseq_per_target), n_Lseq = n_Lseq_per_target)
seq_df <- seq_df %>% left_join(n_Lseq, by = "s_id") 

# if more than maxL possible sequences, limit to maxL, taking the middlemost sequences (ignoring start and end bits of the sequence)
maxL <- 200
seq_df <- seq_df %>% group_by(s_id) %>% 
  mutate(f_id = ifelse(n_Lseq > maxL, f_id - round((n_Lseq - maxL)/2), f_id)) %>% 
  ungroup() %>%
  mutate(n_Lseq = pmin(n_Lseq, maxL))

n_seqL <- sum(pmin(n_Lseq_per_target, maxL))
xtL <- array(rep(0, n_seqL  * L * dim_mat[1] * dim_mat[2]), dim = c(n_seqL, L, dim_mat[1], dim_mat[2]))
cnt <- 0
for(i in 1:length(target.pix)){
  print(i)
  n_Lseq_i <- min(n_Lseq_per_target[i], maxL)
  for(j in 1:n_Lseq_i){
    cnt <- cnt + 1
    xt_ind_to_get <- seq_df %>% filter(s_id == i, f_id == j) %>% dplyr::select(xt_row) %>% unlist() %>% as.numeric()
    for(k in 1:L){
      xtL[cnt, k,,] <- xt[xt_ind_to_get + k - 1,,]
    }
  }
}

# training/test split so that individual targets do not appear in both
seq_df_filt <- seq_df %>% filter(f_id >= 1, f_id <= n_Lseq) %>% mutate(xtL_row = row_number())
y <- as.matrix(seq_df_filt$class, ncol = 1)

# # rename objects here with frame number appended
# xtL_f5 <- xtL
# seq_df_filt_f5 <- seq_df_filt
# y_f5 <- y
# 
# # check that things have been correctly done
# x1 <- xtL_f5[8,1,,]
# x2 <- target.pix[[2]][[1]]
# sum(x1!=x2)

#load("output/runs/lrcn-inputs.Rdata")
# save
save(xtL, seq_df_filt, y, file = "output/lrcn-raw-inputs-frames-from-mid.Rdata")
#save(xtL_f5, seq_df_filt_f5, y_f5, xtL_f9, seq_df_filt_f9, y_f9, file = "output/runs/lrcn-inputs.Rdata")