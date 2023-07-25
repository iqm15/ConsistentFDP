setwd("CONSISTENT_FDP/code")
library(lubridate)
# library(ggplot2)
library(rjson)
source("scripts/preorder_bounds.R")

#------------------ load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config_files/xp_preorder_vary_alpha.json") ##
#-------------------------------------------------------------------------------

nb_run <- param.list$nb_run
nb_bound <- 3

alpha_range = seq(param.list$alpha$begin, param.list$alpha$end, param.list$alpha$by)

b <- param.list$b
m <- param.list$m
pi0 = param.list$pi0
beta = pram.list$beta
delta = param.list$delta
pi1 = 1 - pi0
mu_1 = param.list$mu_1
lambda = param.list$lambda
knockoff = param.list$knockoff
z = param.list$z * m 
s = param.list$s


alpha_column = rep(rep(alpha_range, each = nb_run), nb_bound)
  
Freed_bound <- c()
KR_bound  <- c()
KR_U_bound <- c()
  
for (alpha in alpha_range) {
  
  if(!knockoff) {
    s = s * alpha
  }

  for (i in 1:nb_run) {
    
    if(!knockoff) {
      data_ <- generate_preorder_data(m, mu_1, b, pi1, beta)
      pvalues <- data_$pvalues
    }
    else {
      pvalues <- generate_dummy_knockoff(m, z)
    }
      
    proc <- LF_k_hat(alpha, s, lambda, pvalues)
    k_hat <- proc$k_hat
    nb_rej <- proc$nb_rej
      
    Freed_bound <- c(Freed_bound, preorder_freedman_LF(alpha, delta, s, lambda, k_hat, pvalues))
    KR_bound <- c(KR_bound, preorder_KR_LF(alpha, delta, s, lambda, k_hat, pvalues))
    KR_U_bound <- c(KR_U_bound, preorder_KR_U_LF(alpha, delta, s, lambda, k_hat, pvalues, nb_rej))
      
  }
}
  
  
#-------------------- create data frame and save it ----------------------------
bound = rep(c("FREEDMAN", "KR", "KR-U"), each = nb_run * length(alpha_range))
df_bound <- data.frame(bound, c(Freed_bound, KR_bound, KR_U_bound), alpha_column)
names(df_bound)[1] = "bound"
names(df_bound)[2] = "FDP_bar"
names(df_bound)[3] = "alpha"
file_name_df_bound = gsub(" " , "", paste("../xp_data/preorder/vary_alpha/", gsub(" ", "_", paste("data_bound_vary_alpha",  now(), sep="_")), ".csv"))
write.csv(df_bound, file_name_df_bound, row.names = FALSE)
  
#-------------------- make boxplot and save it  ----------------------------
  
  # plot <- ggplot(df_bound, aes(x=m, y=FDP_bar, color=bound)) +
  #   geom_boxplot() +
  #   facet_wrap(~m, scale="free") +
  #   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  #   geom_hline(yintercept = alpha) + 
  #   labs(x = "m", y = "FDP_bar")
  # plot_name = gsub(" " , "", paste("../xp_plot/preorder", gsub(" ", "_", paste("data_bound_vary_alpha_and_mu", as.character(mu_1),  now(), sep="_")), ".eps"))
  # ggsave(plot_name, width = 16, height = 10, units = "cm", device=cairo_ps)
  
