setwd("CONSISTENT_FDP/code")
library(lubridate)
library(rjson)
source("scripts/preorder_bounds.R")

#------------------ load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config_files/xp_preorder_vary_m_and_alpha.json") ##
#-------------------------------------------------------------------------------

nb_run <- param.list$nb_run
nb_bound <- 3

m_range <- c(10^(param.list$m_range$begin:param.list$m_range$end))
alpha_range = seq(param.list$alpha$begin, param.list$alpha$end, param.list$alpha$by)

b <- param.list$b
pi0 = param.list$pi0
beta = param.list$beta
pi1 = 1 - pi0
mu_1 = param.list$mu_1
f1 = param.list$f1
delta = param.list$delta
lambda = param.list$lambda
knockoff = param.list$knockoff

z = param.list$z

m_column = rep(rep(m_range, each = nb_run), nb_bound + 1)

for (alpha in alpha_range) {
  s = param.list$s
  if(!knockoff) {
    
    s = s * alpha
  }
  
  
  k_hat_vec <- c()
  r_hat_vec <- c()
  
  Freed_bound <- c()
  KR_bound  <- c()
  KR_U_bound <- c()
  KR_U_a <- c()
  
  for (m in m_range) {
    
    if (beta != 0) mu_1 <- sqrt(2 * log(m))
    
    print(m)

    
    for (i in 1:nb_run) {
      
      if(!knockoff) {
       
        data_ <- generate_preorder_data(m, mu_1, b, pi1, beta)
        # data_ <- generate_preorder_data(round(m**(1-beta)), mu_1, b, pi1, 0)
        pvalues <- data_$pvalues
      }
      else {
        pvalues <- generate_dummy_knockoff(m, z, f1, beta)
        
        # alt_probas_new <- 1/2 + pmax(0, 1/2 * ((z - ((1:m) * (m**(beta - 1)))) / (1 - z)) )
        # alt_pvalues <- sample(c(0.5, 1), m, TRUE, c(f1, 1 - f1))
        # pvalues <- numeric(m) + (rbinom(m, 1, alt_probas_new) * alt_pvalues)
        # nb_null <- sum(pvalues == 0)
        # pvalues[which(pvalues == 0)] = 1 - 0.5 * rbinom(nb_null, 1, 1/2)
      }

      proc <- LF_k_hat(alpha, s, lambda, pvalues)
      k_hat <- proc$k_hat
      nb_rej <- proc$nb_rej
      
      k_hat_vec <- c(k_hat_vec, k_hat)
      r_hat_vec <- c(r_hat_vec, nb_rej)

      if(k_hat != 0) {
        Freed_bound <- c(Freed_bound, preorder_freedman_LF(alpha, delta, s, lambda, k_hat, pvalues))
        KR_bound <- c(KR_bound, preorder_KR_LF(alpha, delta, s, lambda, k_hat, pvalues))
        KR_U_bound_list <- preorder_KR_U_LF(alpha, delta, s, lambda, k_hat, pvalues)
        
        KR_U_bound <- c(KR_U_bound, KR_U_bound_list$bound)
        KR_U_a <- c(KR_U_a, KR_U_bound_list$a)
      }
      
      else {
        Freed_bound <- c(Freed_bound, 0)
        KR_bound <- c(KR_bound, 0)
        
        KR_U_bound <- c(KR_U_bound, 0)
        KR_U_a <- c(KR_U_a, 0)
      }
      
      
    }
  }
  
  
  #-------------------- create data frame and save it ----------------------------
  
  bound = rep(c("FREEDMAN", "KR", "KR-U", "KR_U_a"), each = nb_run * length(m_range))
  df_bound <- data.frame(bound, c(Freed_bound, KR_bound, KR_U_bound, KR_U_a), m_column)
  names(df_bound)[1] = "bound"
  names(df_bound)[2] = "FDP_bar"
  names(df_bound)[3] = "m"
  file_name_df_bound = gsub(" " , "", paste("../xp_data/preorder/vary_m_alpha/", gsub(" ", "_", paste("data_bound_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  write.csv(df_bound, file_name_df_bound, row.names = FALSE)


  df_rej <- data.frame(k_hat_vec, r_hat_vec, m_column)
  names(df_rej)[1] = "k_hat"
  names(df_rej)[2] = "r_hat"
  names(df_rej)[3] = "m"
  file_name_df_rej = gsub(" " , "", paste("../xp_data/preorder/vary_m_alpha/", gsub(" ", "_", paste("data_nb_rej_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  write.csv(df_rej, file_name_df_rej, row.names = FALSE)
  
  #-------------------- make boxplot and save it  ----------------------------
  
  # plot <- ggplot(df_bound, aes(x=m, y=FDP_bar, color=bound)) +
  #   geom_boxplot() +
  #   facet_wrap(~m, scale="free") +
  #   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  #   geom_hline(yintercept = alpha) + 
  #   labs(x = "m", y = "FDP_bar")
  # plot_name = gsub(" " , "", paste("../xp_plot/preorder", gsub(" ", "_", paste("data_bound_vary_m_and_mu", as.character(mu_1),  now(), sep="_")), ".eps"))
  # ggsave(plot_name, width = 16, height = 10, units = "cm", device=cairo_ps)
  
} 





