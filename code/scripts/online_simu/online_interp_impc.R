setwd("CONSISTENT_FDP/code")
library(lubridate)
library(rjson)
library(OnlineSuperUnif)
library(DiscreteFDR)
source("scripts/online_bounds.R")
source("scripts/utils.R")

#------------------ Load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config_files/online_interp_impc.json") ##

#------------------- Params ----------------------------------------------------

m = param.list$m
gamma = gamma_sequence(param.list$gamma_type, m, param.list$q_1)
delta <- param.list$delta

alpha_range <- seq(param.list$alpha$begin, param.list$alpha$end, param.list$alpha$by)

#----------------- Pvalues -----------------------------------------------------

Male_df <- read.csv("../xp_data/real_data/my_impc_data/impc_male_df.csv")
Male_test <- fisher.pvalues.support(Male_df, alternative = "greater", input = "noassoc")

for (alpha in alpha_range) {
  w0 = alpha / 2
  
  # ------------------ Procedure --------------------------------------------------
  lord_proc_male <- lord_OnlineSuperUnif(alpha, w0, Male_test$raw[1:m], gamma[1:m])
  
  rej_male <- numeric(m)
  rej_male[lord_proc_male$rej] = 1
  rej_male <- cumsum(rej_male)
  
  #-------------------------------------------------------------------------------
  Freed_bounds <- sapply(1:m, online_Freedman, lord_proc_male$cv, delta, Male_test$raw[1:m], return_Vk = TRUE)
  KR_bounds <- sapply(1:m, online_KR, lord_proc_male$cv, delta, Male_test$raw[1:m], return_Vk = TRUE)
  KR_U_bounds <- sapply(1:m, online_KR_U, lord_proc_male$cv, delta, Male_test$raw[1:m], return_Vk = TRUE)
  
  KR_U_bounds <- unlist_KR_U(KR_U_bounds)
  
  Freed_bounds_interp <- c()
  KR_bounds_interp <- c()
  KR_U_bounds_interp <- c()
  
  for (t in 1:m) {
    
    Freed_bounds_interp <- c(Freed_bounds_interp, get_interpolated_FDP_bound(t, Freed_bounds[1:t], rej_male[1:t])$fdp_interp)
    KR_bounds_interp <- c(KR_bounds_interp, get_interpolated_FDP_bound(t, KR_bounds[1:t], rej_male[1:t])$fdp_interp)
    KR_U_bounds_interp <- c(KR_U_bounds_interp, get_interpolated_FDP_bound(t, KR_U_bounds$KR_U_bound[1:t], rej_male[1:t])$fdp_interp)
  }
  
  
  #-------------------- create data frame and save it ----------------------------
  
  df_bound <- data.frame(KR_bounds_interp, KR_U_bounds_interp, Freed_bounds_interp)
  names(df_bound)[1] = "KR"
  names(df_bound)[2] = "KR_U"
  names(df_bound)[3] = "Freedman"
  file_name_df_bound = gsub(" " , "", paste("../xp_data/online/interp/", gsub(" ", "_", paste("data_bound_interp_impc_alpha", as.character(alpha), now(), sep="_")), ".csv"))
  write.csv(df_bound, file_name_df_bound, row.names = FALSE)
}
