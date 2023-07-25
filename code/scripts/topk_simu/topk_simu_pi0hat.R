setwd("CONSISTENT_FDP/code")
library(lubridate)
# library(ggplot2)
library(rjson)
source("scripts/topk_bounds.R")

#------------------ load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config_files/xp_pi0hat_vary_m_and_alpha.json") ##
#-------------------------------------------------------------------------------

nb_run <- param.list$nb_run
nb_bound <- 4

alpha_range = seq(param.list$alpha$begin, param.list$alpha$end, param.list$alpha$by)
m_range <- c(10^(param.list$m_range$begin:param.list$m_range$end))

mu_1 <- param.list$mu_1
pi0 = param.list$pi0
alpha = param.list$alpha
delta = param.list$delta

m_column = rep(rep(m_range, each = nb_run), nb_bound)

for (alpha in alpha_range) {
  
  hatpi0_Simes <- c()
  hatpi0_DKW <- c()
  hatpi0_KR  <- c() 
  hatpi0_Well <- c()
  
  for (m in m_range) {
    print(m)
    for (i in 1:nb_run) {
      
      data__ <- generate_pvalues(m, mu_1, pi0)
      pvalues <- data__$raw
      
      
      hatpi0_Simes = c(hatpi0_Simes, m0_hat_Simes(pvalues, delta) / m)
      hatpi0_DKW = c(hatpi0_DKW, m0_hat_DKW(pvalues, delta) / m)
      hatpi0_KR = c(hatpi0_KR, m0_hat_KR(pvalues, delta) / m)
      hatpi0_Well = c(hatpi0_Well, m0_hat_Well(pvalues, delta) / m)
      
      # else hatpi0_Simes = hatpi0_DKW = hatpi0_KR = hatpi0_Well = NULL
      # 
      # 
      # k_hat <- BH_k_hat(alpha, pvalues)
      
      # Simes_bound <- c(Simes_bound, topk_Simes_BH(alpha, delta, adaptive, hatpi0_Simes))
      # DKW_bound <- c(DKW_bound, topk_DKW_BH(alpha, delta, m, k_hat, adaptive, hatpi0_DKW))
      # KR_bound <- c(KR_bound, topk_KR_BH(alpha, delta, k_hat, adaptive, hatpi0_KR))
      # Wellner_bound <- c(Wellner_bound, topk_Wellner_BH(alpha, delta, m, k_hat, adaptive, hatpi0_Well))
      # hybrid_KR_Well_bound <- c(hybrid_KR_Well_bound, topk_hybrid_KR_Well_BH(alpha, delta / 2, m, k_hat, adaptive, hatpi0_KR, hatpi0_Well))
    }
  }
  
  
  #-------------------- create data frame and save it ----------------------------
  estim = rep(c("DKW", "KR", "SIMES", "WELLNER"), each = nb_run * length(m_range))
  df_bound <- data.frame(estim, c(hatpi0_DKW, hatpi0_KR, hatpi0_Simes, hatpi0_Well), m_column)
  names(df_bound)[1] = "bound"
  names(df_bound)[2] = "pi_hat"
  names(df_bound)[3] = "m"
  file_name_df_bound = gsub(" " , "", paste("../xp_data/topk/pi0hat/", gsub(" ", "_", paste("data_pi0hat_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  write.csv(df_bound, file_name_df_bound, row.names = FALSE)
  
}