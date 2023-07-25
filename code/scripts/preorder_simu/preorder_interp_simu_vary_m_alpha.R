setwd("CONSISTENT_FDP/code")
library(lubridate)
library(rjson)
source("scripts/preorder_bounds.R")
source("scripts/utils.R")

#------------------ load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config_files/xp_preorder_interp_vary_m_and_alpha.json") ##
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
to_interp = param.list$interpolate

m_column = rep(rep(m_range, each = nb_run), nb_bound)

for (alpha in alpha_range) {
  s = param.list$s
  if(!knockoff) {
    
    s = s * alpha
  }
  
  
  k_hat_vec <- c()
  r_hat_vec <- c()
  
  Freed_bound_vec <- c()
  KR_bound_vec  <- c()
  KR_U_bound_vec <- c()
  
  for (m in m_range) {
    if (beta != 0) mu_1 <- sqrt(2 * log(m))

    print(m)
    
    for (i in 1:nb_run) {
      
      if(!knockoff) {
        data_ <- generate_preorder_data(m, mu_1, b, pi1, beta)
        pvalues <- data_$pvalues
      }
      else {
        pvalues <- generate_dummy_knockoff(m, z, f1, beta)
      }
      
      proc <- LF_k_hat(alpha, s, lambda, pvalues)
      k_hat <- proc$k_hat
      nb_rej <- proc$nb_rej

      k_hat_vec <- c(k_hat_vec, k_hat)
      r_hat_vec <- c(r_hat_vec, nb_rej)
  
      if(k_hat != 0) {

        nb_rej_vec <- proc$nb_rej_vec

        Freed_bound <- sapply(1:k_hat, preorder_freedman, s, lambda, delta, pvalues, return_Vk = TRUE)
        KR_bound <- sapply(1:k_hat, preorder_KR, s, lambda, delta, pvalues, return_Vk = TRUE)
        KR_U_bound <- sapply(1:k_hat, preorder_KR_U, s, lambda, delta, pvalues, return_Vk = TRUE)
        
        unlist_KR_U_ <- unlist_KR_U(KR_U_bound)
        KR_U_a_vec <- unlist_KR_U_$KR_U_best_a

        Freed_bound_vec <- c(Freed_bound_vec, get_interpolated_FDP_bound(k_hat, Freed_bound, nb_rej_vec)$fdp_interp)
        KR_bound_vec <- c(KR_bound_vec, get_interpolated_FDP_bound(k_hat, KR_bound, nb_rej_vec)$fdp_interp)
        KR_U_bound_vec <- c(KR_U_bound_vec, get_interpolated_FDP_bound(k_hat, unlist_KR_U_$KR_U_bound, nb_rej_vec)$fdp_interp)
      }

      else {
        Freed_bound_vec <- c(Freed_bound_vec, 0)
        KR_bound_vec <- c(KR_bound_vec, 0)
        KR_U_bound_vec <- c(KR_U_bound_vec, 0)
      }

    }
  }
  
  #-------------------- create data frame and save it ----------------------------
  
  bound = rep(c("FREEDMAN", "KR", "KR-U"), each = nb_run * length(m_range))
  df_bound <- data.frame(bound, c(Freed_bound_vec, KR_bound_vec, KR_U_bound_vec), m_column)
  names(df_bound)[1] = "bound"
  names(df_bound)[2] = "FDP_bar"
  names(df_bound)[3] = "m"
  file_name_df_bound = gsub(" " , "", paste("../xp_data/preorder/vary_m_alpha/interp/", gsub(" ", "_", paste("data_bound_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  write.csv(df_bound, file_name_df_bound, row.names = FALSE)
  
  
  df_rej <- data.frame(k_hat_vec, r_hat_vec, m_column)
  names(df_rej)[1] = "k_hat"
  names(df_rej)[2] = "r_hat"
  names(df_rej)[3] = "m"
  file_name_df_rej = gsub(" " , "", paste("../xp_data/preorder/vary_m_alpha/interp/", gsub(" ", "_", paste("data_nb_rej_vary_m_and_alpha", as.character(alpha),  now(), sep="_")), ".csv"))
  write.csv(df_rej, file_name_df_rej, row.names = FALSE)
  
  
} 





