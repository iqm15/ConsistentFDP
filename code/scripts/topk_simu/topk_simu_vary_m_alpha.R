setwd("CONSISTENT_FDP/code")
library(lubridate)
# library(ggplot2)
library(rjson)
source("scripts/topk_bounds.R")

#------------------ load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config_files/xp_topk_vary_m_and_alpha.json") ##
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

m_column = rep(rep(m_range, each = nb_run), nb_bound)

for (alpha in alpha_range) {
  
  k_hat_vec <- c()
  
  Simes_bound <- c()
  DKW_bound <- c()
  KR_bound  <- c()
  Wellner_bound <- c()
  hybrid_KR_Well_bound <- c()
  
  for (m in m_range) {
    for (i in 1:nb_run) {
      
      if(param.list$DU) pvalues <- generate_DU_pvalues(m, pi0)
      else {
        data__ <- generate_pvalues(m, mu_1, pi0)
        pvalues <- data__$raw
      }
      
      if (adaptive) {
        hatpi0_Simes = m0_hat_Simes(pvalues, delta) / m
        hatpi0_DKW = m0_hat_DKW(pvalues, delta) / m
        hatpi0_KR = m0_hat_KR(pvalues, delta) / m
        hatpi0_Well = m0_hat_Well(pvalues, delta) / m
      }
      else hatpi0_Simes = hatpi0_DKW = hatpi0_KR = hatpi0_Well = NULL
    
      
      k_hat <- BH_k_hat(alpha, pvalues)
      
      k_hat_vec <- c(k_hat_vec, k_hat)
      
      Simes_bound <- c(Simes_bound, topk_Simes_BH(alpha, delta, adaptive, hatpi0_Simes))
      DKW_bound <- c(DKW_bound, topk_DKW_BH(alpha, delta, m, k_hat, adaptive, hatpi0_DKW))
      KR_bound <- c(KR_bound, topk_KR_BH(alpha, delta, k_hat, adaptive, hatpi0_KR))
      Wellner_bound <- c(Wellner_bound, topk_Wellner_BH(alpha, delta, m, k_hat, adaptive, hatpi0_Well))
      hybrid_KR_Well_bound <- c(hybrid_KR_Well_bound, topk_hybrid_KR_Well_BH(alpha, delta / 2, m, k_hat, adaptive, hatpi0_KR, hatpi0_Well))
    }
  }
  
  
  #-------------------- create data frame and save it ----------------------------
  bound = rep(c("DKW", "KR", "SIMES", "WELLNER", "HYBRID"), each = nb_run * length(m_range))
  df_bound <- data.frame(bound, c(DKW_bound, KR_bound, Simes_bound, Wellner_bound, hybrid_KR_Well_bound), m_column)
  names(df_bound)[1] = "bound"
  names(df_bound)[2] = "FDP_bar"
  names(df_bound)[3] = "m"
  file_name_df_bound = gsub(" " , "", paste("../xp_data/topk/vary_m_alpha/", gsub(" ", "_", paste("DU_data_bound_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  write.csv(df_bound, file_name_df_bound, row.names = FALSE)
  
  df_rej <- data.frame(k_hat_vec, m_column)
  names(df_rej)[1] = "k_hat"
  names(df_rej)[2] = "m"
  file_name_df_rej = gsub(" " , "", paste("../xp_data/topk/vary_m_alpha/", gsub(" ", "_", paste("DU_data_nb_rej_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  write.csv(df_rej, file_name_df_rej, row.names = FALSE)
  
  #-------------------- make boxplot and save it  ----------------------------
  
  # plot <- ggplot(df_bound, aes(x=m, y=FDP_bar, color=bound)) +
  #         geom_boxplot() +
  #         facet_wrap(~m, scale="free") +
  #         theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  #         geom_hline(yintercept = alpha) + 
  #         labs(x = "m", y = "FDP_bar")
  # plot_name = gsub(" " , "", paste("../xp_plot/topk/", gsub(" ", "_", paste("data_bound_vary_m_and_mu", as.character(mu_1),  now(), sep="_")), ".eps"))
  # ggsave(plot_name, width = 16, height = 10, units = "cm", device=cairo_ps)
  
} 