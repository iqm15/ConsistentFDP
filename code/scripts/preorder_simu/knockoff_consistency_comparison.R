setwd("CONSISTENT_FDP/code")
library(knockoff)
library(lubridate)
library(knockoff)
library(rjson)
source("scripts/preorder_bounds.R")
source("scripts/utils.R")

#------------------ load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config_files/knockoff_consistency_comparison.json") ##
#-------------------------------------------------------------------------------

nb_methods = 6
delta = param.list$delta
proc_s = param.list$proc_s
lambda = param.list$lambda
# alpha = param.list$alpha

alpha_column <- rep(rep(param.list$alpha_grid, each = param.list$nb_run), nb_methods)

load(paste0(param.list$tuning_param_path))

if(delta != 0.05) { 
  Bernoulli_mat <- get_Bernoulli_mat(num_trial=200, B=100000)
  
  k_list <- list( get_k_two_step(v_list[[1]], Bernoulli_mat, delta, step1_size_vec=c(1, 0.5, 0.1, 0.05, 0.01), step2_size_vec=c(50, 25, 10, 5, 1))$k_list[[3]],
                  get_k_two_step(v_list[[2]], Bernoulli_mat, delta, step1_size_vec=c(1, 0.5, 0.1, 0.05, 0.01), step2_size_vec=c(50, 25, 10, 5, 1))$k_list[[3]],
                  get_k_two_step(v_list[[3]], Bernoulli_mat, delta, step1_size_vec=c(1, 0.5, 0.1, 0.05, 0.01), step2_size_vec=c(50, 25, 10, 5, 1))$k_list[[3]],
                  get_k_two_step(v_list[[4]], Bernoulli_mat, delta, step1_size_vec=c(1, 0.5, 0.1, 0.05, 0.01), step2_size_vec=c(50, 25, 10, 5, 1))$k_list[[3]]
  )
}

p = param.list$p        # number of variables
n = param.list$n         # number of observations

rho = param.list$rho
Sigma = toeplitz(rho^(0:(p-1)))


for (amplitude in param.list$amplitude_grid){
  

  
  
  for (s in param.list$sparsity_grid) {
    
    KR_bound_interp <- c()
    KR_U_bound_interp <- c()
    FDP_JKI_A <- c()
    FDP_JKI_B <- c()
    FDP_JKI_C <- c()
    FDP_JKI_D <- c()
    
    # matrix <- matrix(data = NA, nrow = nb_methods * param.list$nb_run, ncol = p, byrow = FALSE)
    
    k = ceiling(s*p)
    
    for (alpha in param.list$alpha_grid) {
    
      for (r in 1:param.list$nb_run){
        # set.seed(2022)
        nonzero = sample(p, k) # number of variables with nonzero coefficients
        beta = amplitude * (1:p %in% nonzero) / sqrt(n)
        
        ######### Generate data
        X = matrix(rnorm(n*p),n) %*% chol(Sigma)
        y = X %*% beta + rnorm(n)
        
        ######### get original W using knockoff method
        W_ori <- knockoff.filter(X, y, knockoffs = create.fixed, statistic = stat.glmnet_lambdasmax, fdr=0.1)$statistic
        W_list <- preprocess_W_func(W_ori)  # Get the preprocessed W: |W_1|>...>|W_m|>0
        W <- W_list[[1]]
        W_order <- W_list[[2]]
        
        #### get the true FDP
        true_index <- mapply(function(i, W_order ){return(which(W_order==i))}, which(beta!=0), MoreArgs=list(W_order))
        FDP_true <- get_FDP_true(W, true_index)
        
        # get binary pvalues to compute the bounds
        sorted_signed_knockoff <-  preprocess_W_func(W_ori)$W_sort # get_sorted_signed_stat(W_ori)
        binary_pvalues <- 0.5 * (sorted_signed_knockoff > 0) + 1 * (sorted_signed_knockoff <= 0)
        
        
        # Run LF procedue 
        proc <- LF_k_hat(alpha, proc_s, lambda, binary_pvalues)
        k_hat <- proc$k_hat
        nb_rej <- proc$nb_rej
        
        if(k_hat != 0) {
        
        # compute base bounds and fdp estimate
        # FDP_hat <- preorder_FDP_hat(proc_s, lambda, binary_pvalues)
        KR_bound <- sapply(1:k_hat, preorder_KR, proc_s, lambda, delta, binary_pvalues, return_Vk = TRUE)
        KR_U_bound <- sapply(1:k_hat, preorder_KR_U, proc_s, lambda, delta, binary_pvalues, return_Vk = TRUE)
        
        # interpolate base bounds
        unlist_KR_U_bound <- unlist_KR_U(KR_U_bound)
        # nb_rej <- cumsum(binary_pvalues <= proc_s)
        
        # KR_bound_interp <- c()
        # KR_U_bound_interp <- c()
        
        # for (j in 1:length(binary_pvalues)) {
        KR_bound_interp <- c(KR_bound_interp, get_interpolated_FDP_bound(k_hat, KR_bound, proc$nb_rej_vec)$fdp_interp)
        KR_U_bound_interp <- c(KR_U_bound_interp, get_interpolated_FDP_bound(k_hat, unlist_KR_U_bound$KR_U_bound, proc$nb_rej_vec)$fdp_interp)
        # }
        
        # compute Li and Goeman KJI bounds
        
        FDP_JKI_A <- c(FDP_JKI_A, get_FDP_KJI(W, k_vec=k_list[[1]], v_vec=v_list[[1]])[k_hat])
        FDP_JKI_B <- c(FDP_JKI_B, get_FDP_KJI(W, k_vec=k_list[[2]], v_vec=v_list[[2]])[k_hat])
        FDP_JKI_C <- c(FDP_JKI_C, get_FDP_KJI(W, k_vec=k_list[[3]], v_vec=v_list[[3]])[k_hat])
        FDP_JKI_D <- c(FDP_JKI_D, get_FDP_KJI(W, k_vec=k_list[[4]], v_vec=v_list[[4]])[k_hat])
        
        }
        
        else {
          KR_bound_interp <- c(KR_bound_interp, 0)
          KR_U_bound_interp <- c(KR_U_bound_interp, 0)
          FDP_JKI_A <- c(FDP_JKI_A, 0)
          FDP_JKI_B <- c(FDP_JKI_B, 0)
          FDP_JKI_C <- c(FDP_JKI_C, 0)
          FDP_JKI_D <- c(FDP_JKI_D, 0)
          
        }
      }
    }
    
    bound = rep(c("KR", "KR-U", "KJI-A", "KJI-B", "KJI-C", "KJI-D"), each = param.list$nb_run * length(param.list$alpha_grid))
    print(length(bound))
    print(length(c(KR_bound_interp, KR_U_bound_interp, FDP_JKI_A, FDP_JKI_B, FDP_JKI_C, FDP_JKI_D)))
    df_bound <- data.frame(bound, c(KR_bound_interp, KR_U_bound_interp, FDP_JKI_A, FDP_JKI_B, FDP_JKI_C, FDP_JKI_D), alpha_column)
    names(df_bound)[1] = "bound"
    names(df_bound)[2] = "FDP_bar"
    names(df_bound)[3] = "alpha"
    file_name_df_bound = gsub(" " , "", paste("../xp_data/preorder/LGcompar/", gsub(" ", "_", paste("consistency_comparison_amplitude", as.character(amplitude), "sparsity", as.character(s), now(), sep="_")), ".csv"))
    write.csv(df_bound, file_name_df_bound, row.names = FALSE)
  }
}     
      
      
      
      
      
      
    #   # matrix containing bound for each method and each round of MC simulation
    #   if (length(binary_pvalues) != p ) {
    #     KR_bound_interp <- c(KR_bound_interp, rep(0, p - length(binary_pvalues)))
    #     KR_U_bound_interp <- c(KR_U_bound_interp, rep(0, p - length(binary_pvalues)))
    #     FDP_JKI_A <- c(FDP_JKI_A, rep(0, p - length(binary_pvalues)))
    #     FDP_JKI_B <- c(FDP_JKI_B, rep(0, p - length(binary_pvalues)))
    #     FDP_JKI_C <- c(FDP_JKI_C, rep(0, p - length(binary_pvalues)))
    #     FDP_JKI_D <- c(FDP_JKI_D, rep(0, p - length(binary_pvalues)))
    #     FDP_true <- c(FDP_true, rep(0, p - length(binary_pvalues)))
    #   }
    #   
    #   matrix[7 * (r-1) + 1, ] <- KR_bound_interp
    #   matrix[7 * (r-1) + 2, ] <- KR_U_bound_interp
    #   matrix[7 * (r-1) + 3, ] <- FDP_JKI_A
    #   matrix[7 * (r-1) + 4, ] <- FDP_JKI_B
    #   matrix[7 * (r-1) + 5, ] <- FDP_JKI_C
    #   matrix[7 * (r-1) + 6, ] <- FDP_JKI_D
    #   matrix[7 * (r-1) + 7, ] <- FDP_true
    #   
    # }
    # 
    # 
    # 
    # method = rep(c("KR", "KR-U", "KJI_A", "KJI_B", "KJI_C", "KJI_D", "FDP_true"), param.list$nb_run)
    # 
    # # convert it into dataframe
    # # df <- data.frame(cbind(method, matrix))
    # 
    # # get mean value over all MC simu rounds for each method
    # bounds_means <- matrix(data = NA, nrow = nb_methods, ncol = p, byrow = FALSE)
    # 
    # for (l in 2:(p+1)) {
    #   
    #   df <- data.frame(cbind(method, matrix))
    #   df <- data.frame(df$method, as.numeric(df[, l]))
    #   names(df)[2] = "bound"
    #   
    #   bounds_means[, l-1] <- aggregate(.~df.method, data=df, mean)$bound
    #   
    # }
    # 
    # data_mean <- data.frame(cbind(c("FDP_true", "KJI_A", "KJI_B", "KJI_C", "KJI_D", "KR", "KR-U"), bounds_means))
    # 
    # file_name = gsub(" " , "", paste("../xp_data/preorder/LGcompar/", gsub(" ", "_", paste("comparison_amplitude", as.character(amplitude), "sparsity", as.character(s), now(), sep="_")), ".csv"))
    # 
    # write.csv(data_mean, file_name, row.names = FALSE)
    
