setwd("CONSISTENT_FDP/code")
library(lubridate)
# library(ggplot2)
library(rjson)
source("scripts/topk_bounds.R")
source("scripts/utils.R")

#------------------ load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config_files/xp_topk_interp_vary_m_and_alpha.json") ##
#-------------------------------------------------------------------------------

nb_run <- param.list$nb_run
nb_bound <- 5

# mu_1_range <- seq(param.list$mu_1$begin, param.list$mu_1$end, param.list$mu_1$by)
alpha_range = seq(param.list$alpha$begin, param.list$alpha$end, param.list$alpha$by)
m_range <- c(10^(param.list$m_range$begin:param.list$m_range$end))

mu_1 <- param.list$mu_1
pi0 = param.list$pi0
delta = param.list$delta
alpha = param.list$alpha
adaptive = param.list$adaptive
beta = param.list$beta

m_column = rep(rep(m_range, each = nb_run), nb_bound)

for (alpha in alpha_range) {
  print(alpha)
  k_hat_vec <- c()
  
  Simes_bound_vec <- c()
  DKW_bound_vec <- c()
  KR_bound_vec  <- c()
  Wellner_bound_vec <- c()
  hybrid_KR_Well_bound_vec <- c()
  
  for (m in m_range) {
    print(m)
    m1 = (1 - pi0) * m 

    if(beta != 0) {
      pi0 = 1 - 0.5*(m**(-beta))
      # mu_1 <- sqrt(2 * log(m))
      mu_1 = 10
    }
    for (i in 1:nb_run) {
      
      if(param.list$DU) pvalues <- generate_DU_pvalues(m, pi0)
      else {
        data__ <- generate_pvalues(m, mu_1, pi0)
        pvalues <- data__$raw
      }
      
      if (adaptive) {
        hatpi0_Simes = m0_hat_Simes(pvalues, delta) 
        hatpi0_DKW = m0_hat_DKW(pvalues, delta) 
        hatpi0_KR = m0_hat_KR(pvalues, delta) 
        hatpi0_Well = m0_hat_Well(pvalues, delta) 
      }
      else hatpi0_Simes = hatpi0_DKW = hatpi0_KR = hatpi0_Well = NULL
    
      k_hat <- BH_k_hat(alpha, pvalues)
      k_hat_vec <- c(k_hat_vec, k_hat)
      
      if(k_hat != 0) {
        
        Simes_bound <- sapply(1:k_hat, topk_Simes, delta, pvalues, adaptive, hatpi0_Simes, return_Vk = TRUE)
        DKW_bound <- sapply(1:k_hat, topk_DKW, delta, pvalues, adaptive, hatpi0_DKW, return_Vk = TRUE)
        KR_bound <- sapply(1:k_hat, topk_KR, delta, pvalues, adaptive, hatpi0_KR, return_Vk = TRUE)
        Wellner_bound <- sapply(1:k_hat, topk_Wellner, delta, pvalues, adaptive, hatpi0_Well, return_Vk = TRUE)
        hybrid_KR_Well_bound <- sapply(1:k_hat, topk_hybrid_KR_Well, delta, pvalues, adaptive, hatpi0_KR, hatpi0_Well, return_Vk = TRUE)
        
        Simes_bound_vec <- c(Simes_bound_vec, get_interpolated_FDP_bound(k_hat, Simes_bound)$fdp_interp)
        DKW_bound_vec <- c(DKW_bound_vec, get_interpolated_FDP_bound(k_hat, DKW_bound)$fdp_interp)
        KR_bound_vec <- c(KR_bound_vec, get_interpolated_FDP_bound(k_hat, KR_bound)$fdp_interp)
        Wellner_bound_vec <- c(Wellner_bound_vec, get_interpolated_FDP_bound(k_hat, Wellner_bound)$fdp_interp)
        hybrid_KR_Well_bound_vec <- c(hybrid_KR_Well_bound_vec, get_interpolated_FDP_bound(k_hat, hybrid_KR_Well_bound)$fdp_interp)
      }

      else {
        Simes_bound_vec <- c(Simes_bound_vec, 0)
        DKW_bound_vec <- c(DKW_bound_vec, 0)
        KR_bound_vec <- c(KR_bound_vec, 0)
        Wellner_bound_vec <- c(Wellner_bound_vec, 0)
        hybrid_KR_Well_bound_vec <- c(hybrid_KR_Well_bound_vec, 0)
      }

    }
  }

  #-------------------- create data frame and save it ----------------------------
  bound = rep(c("DKW", "KR", "SIMES", "WELLNER", "HYBRID"), each = nb_run * length(m_range))
  df_bound <- data.frame(bound, c(DKW_bound_vec, KR_bound_vec, Simes_bound_vec, Wellner_bound_vec, hybrid_KR_Well_bound_vec), m_column)
  names(df_bound)[1] = "bound"
  names(df_bound)[2] = "FDP_bar"
  names(df_bound)[3] = "m"
  if(beta != 0) file_name_df_bound = gsub(" " , "", paste("../xp_data/topk/vary_m_alpha/interp/", gsub(" ", "_", paste("data_interp_bound_sparse_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  if(adaptive) file_name_df_bound = gsub(" " , "", paste("../xp_data/topk/vary_m_alpha/interp/adapt/", gsub(" ", "_", paste("data_interp_bound_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  else file_name_df_bound = gsub(" " , "", paste("../xp_data/topk/vary_m_alpha/interp/", gsub(" ", "_", paste("data_interp_bound_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  write.csv(df_bound, file_name_df_bound, row.names = FALSE)

  df_rej <- data.frame(k_hat_vec, m_column)
  names(df_rej)[1] = "k_hat"
  names(df_rej)[2] = "m"
  if(beta != 0) file_name_df_rej = gsub(" " , "", paste("../xp_data/topk/vary_m_alpha/interp/", gsub(" ", "_", paste("data_nb_rej_sparse_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  if(adaptive) file_name_df_rej = gsub(" " , "", paste("../xp_data/topk/vary_m_alpha/interp/adapt/", gsub(" ", "_", paste("data_nb_rej_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  else file_name_df_rej = gsub(" " , "", paste("../xp_data/topk/vary_m_alpha/interp/", gsub(" ", "_", paste("data_nb_rej_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  write.csv(df_rej, file_name_df_rej, row.names = FALSE)

  
} 