get_interpolated_FDP_bound <- function(index, init_Vk_bound_vec, nb_rej_vec = NULL) {
  
  # card <- cumsum(rep(1, length(init_Vk_bound_vec)))
  # d_interp <- rep(0, length(init_Vk_bound_vec))
  # interp_index <- rep(0, length(init_Vk_bound_vec))
  
  # for (i in 1:length(init_Vk_bound_vec)) {
  #   d_interp[i] <- max(card[1:i] - init_Vk_bound_vec[1:i])
  #   interp_index[i] <- which.max(card[1:i] - init_Vk_bound_vec[1:i])
  # }
  
  # if in the topk setting, card(R_k) = k
  if(is.null(nb_rej_vec)) {
    d_interp <- max(1:index - init_Vk_bound_vec)
    interp_index <- which.max(1:index - init_Vk_bound_vec)
    card = index
    }
  
  # if in the preorder or online setting
  else {
    d_interp <- max(nb_rej_vec - init_Vk_bound_vec)
    interp_index <- which.max(nb_rej_vec - init_Vk_bound_vec)
    card = nb_rej_vec[index]
  } 
  
  fdp_interp = min((card - d_interp) / card, 1)
  
  output <- list(fdp_interp = fdp_interp, d_interp = d_interp, interp_index = interp_index)
  return(output)
}


get_interpolated_FDP_bound_bis = function(init_Vk_bound_vec, num_discoveries_vec){
  
  max_false_positives = pmin(num_discoveries_vec, init_Vk_bound_vec)

  min_true_positives = num_discoveries_vec - max_false_positives
  
  # interpolation
  max_false_positives = rev(cummin(rev(max_false_positives)))
  min_true_positives = cummax(min_true_positives)
  max_false_positives = pmin(max_false_positives, num_discoveries_vec - min_true_positives)
  
  max_false_positives <- pmin(max_false_positives, num_discoveries_vec)
  interp_FDP_bound <- max_false_positives / pmax(1, num_discoveries_vec)
  
  return(interp_FDP_bound)
}
 
unlist_KR_U <- function(KR_U_bound) {
  
  KR_U_bound_new <- c()
  KR_U_best_a <- c()
  for (i in 1:(length(KR_U_bound)/2.)) {
    
    # if(length(KR_U_bound[[i]]) == 1) {
    #   # if (is.na(KR_U_bound[[i]])) {}
    #   # else 
    #   KR_U_bound_new <- c(KR_U_bound_new, KR_U_bound[[i]])
    #   KR_U_best_a <- c(KR_U_best_a, 1)
    # }
    # 
    # else {
    KR_U_bound_new <- c(KR_U_bound_new, KR_U_bound[, i]$bound)
    KR_U_best_a <- c(KR_U_best_a, KR_U_bound[, i]$a)
    # } 
  }
  output <- list(KR_U_bound = KR_U_bound_new, KR_U_best_a = KR_U_best_a)
  return(output)
}

exp_scale <- function(x) {
  
  # if(round(exp(x), digits = 2) == 0) return(signif(exp(x), 2))
  # else return(round(exp(x), digits = 2))
  return(signif(exp(x), 2))
}

# form Li and Goeman
preprocess_W_func <- function(W){    
  W_order <- order(abs(W), decreasing=TRUE)
  W_sort <- W[W_order]
  # delete zeros
  zero_index <- which(W_sort==0)
  if(length(zero_index)>0){
    W_sort <- W_sort[-zero_index]
    W_order <- W_order[-zero_index]
  }
  
  # break the tie if there is any without changing the order
  W_sort_abs <- abs(W_sort)
  for (i in 1:length(W_sort_abs)) {
    temp <- W_sort_abs[i]
    if(sum(W_sort_abs==temp) >= 2){     # if there is tie
      tie_index <- which(W_sort_abs==temp)
      first_index <- tie_index[1]
      last_index <- tail(tie_index,1)
      
      if(last_index != length(W_sort_abs)){
        max_value <- W_sort_abs[last_index] - W_sort_abs[last_index+1]
        W_sort_abs[tie_index] <- W_sort_abs[tie_index]- seq(0, max_value/2,length.out=length(tie_index))
      }else{ # if the tie contains the last element
        max_value <- W_sort_abs[first_index-1] - W_sort_abs[first_index]
        W_sort_abs[tie_index] <- W_sort_abs[tie_index] + seq(max_value/2,0,length.out=length(tie_index))
      }
    } # end if
  } # end for
  
  W_sort <- sign(W_sort) * W_sort_abs
  return(list(W_sort=W_sort,W_order=W_order))
}


get_Bernoulli_mat <- function(num_trial, B){
  ### generate the Bernoulli matrix: B x num_trial
  return(matrix(sample(c(-1,1), B*num_trial, replace=TRUE),nrow=B, ncol=num_trial))
}

get_k_two_step <- function(v_vec, Bernoulli_mat, alpha, step1_size_vec=c(1,0.5,0.1,0.05,0.01), 
                           step2_size_vec=c(50,25,10,5,1)){
  if(length(v_vec)<2){stop("when length(v_vec)=1, can calculate k based on JS, no need to use this function!")}
  c_alpha <- -log(alpha)/log(2-alpha) 
  
  # based on the simulated Bernoulli trials, get 'stat_matrix' using 'v_vec' (B x m)
  Test_stat_mat <- mapply(function(i, Bernoulli_mat_input, v_vec_input){
    return(apply(Bernoulli_mat_input, 1, get_NBstat_simu, v_vec_input[i]))},
    1:length(v_vec), MoreArgs=list(Bernoulli_mat, v_vec))
  
  # find the smallest 'j' s.t. 'j-c_j+1=max(v_vec)': to obtain 'k' for 'max(v_vec)'
  j_max <- max(v_vec)
  c_temp <- floor(c_alpha*(1+j_max)/(1+c_alpha)) + 1
  v_max_temp <- j_max - c_temp + 1
  while (v_max_temp < max(v_vec)) {
    j_max <- j_max+1
    c_temp <- floor(c_alpha*(1+j_max)/(1+c_alpha)) + 1
    v_max_temp <- j_max - c_temp + 1
  }
  
  # Step0: get 'k^{raw}' and 'j^*' based on 'v_vec' (see definition of 'k^{raw}' in the paper)
  c_raw <- floor(c_alpha*(1+(1:j_max))/(1+c_alpha)) + 1
  b_raw <- (1:j_max) - c_raw + 1
  j_star <- mapply(function(i, b_input){
    return(which(b_input==i)[1])}, v_vec, MoreArgs=list(b_raw))
  k_raw <- pmin(c_raw[j_star], ncol(Bernoulli_mat)+1) # k-1 must be less than number of trials
  prob_raw <- 1 - mean(apply(t(t(Test_stat_mat) >= k_raw), 1, sum) > 0)  # original probability: didn't exhaust the alpha level.
  
  ### Two-step approach
  ### Step 1
  Step1_list <- get_k_updateC(Test_stat_mat, j_star, alpha, step1_size_vec)
  k_step1 <- pmin(Step1_list[[1]], ncol(Bernoulli_mat)+1)
  prob_step1 <- Step1_list[[2]]   # take a look: more or less exhaust the alpha level.
  c_step1 <- Step1_list[[3]]      # compare to the previous c(alpha)
  
  ### Step 2
  Step2_list <- get_k_greedy(Test_stat_mat, k_step1, v_vec, alpha, step2_size_vec)
  k_step2 <- pmin(Step2_list[[1]], ncol(Bernoulli_mat)+1)
  prob_step2 <- Step2_list[[2]]   # take a look: more or less exhaust the alpha level.
  
  return(list(k_list=list(k_raw, k_step1, k_step2), 
              #b_list=list(b_raw, b_step1, b_step2), 
              c_step1,
              prob_list=list(prob_raw, prob_step1, prob_step2)))
}

############################### Other functions
##### get NB test statistic based on one sequence of Bernoulli trials
get_NBstat_simu <- function(vec, para){
  if(sum(vec<0) >= para){return( sum(vec[1: which(vec<0)[para] ] > 0) )}
  if(sum(vec<0) < para){return( sum(vec>0) )}
}

##### step 1 update
get_k_updateC <- function(Test_stat_mat, j_star, alpha, step_size_vec){
  c_alpha <- -log(alpha)/log(2-alpha) 
  
  for (i in 1:length(step_size_vec)) {
    step_size <- step_size_vec[i]
    
    k_temp <- floor(c_alpha*(1+j_star)/(1+c_alpha)) + 1
    prob <- 1 - mean(apply(t(t(Test_stat_mat) >= k_temp), 1, sum) > 0)  
    
    c_alpha_temp <- c_alpha
    while (prob>=1-alpha) {
      c_alpha <- c_alpha_temp
      c_alpha_temp <- c_alpha_temp - step_size
      
      k_temp <- floor(c_alpha_temp*(1+j_star)/(1+c_alpha_temp)) + 1
      prob <- 1 - mean(apply(t(t(Test_stat_mat) >= k_temp), 1, sum) > 0)   
    }
  }
  
  k_step1 <- floor(c_alpha*(1+j_star)/(1+c_alpha)) + 1
  prob <- 1 - mean(apply(t(t(Test_stat_mat) >= k_step1), 1, sum) > 0)
  return(list(k_step1, prob, c_alpha))
}

######## step 2 update
get_k_greedy <- function(Test_stat_mat, k_vec, v_vec, alpha, step_size_vec){
  for (i in 1:length(k_vec)) {
    for (j in 1:length(step_size_vec)) {
      step_size <- step_size_vec[j]
      
      if(k_vec[i]-step_size < v_vec[i] ){next} # if step size is too large
      
      Prob_temp <- 1 - mean(apply(t(t(Test_stat_mat) >= k_vec), 1, sum) > 0)
      
      k_greedy_temp <- k_vec
      while (Prob_temp >= 1-alpha) {
        k_vec <- k_greedy_temp
        k_greedy_temp[i] <- k_greedy_temp[i] - step_size
        Prob_temp <- 1 - mean(apply(t(t(Test_stat_mat) >= k_greedy_temp), 1, sum) > 0)
      }
    }
  }
  
  prob <- 1 - mean(apply(t(t(Test_stat_mat) >= k_vec), 1, sum) > 0)
  return(list(k_vec, prob))
}



