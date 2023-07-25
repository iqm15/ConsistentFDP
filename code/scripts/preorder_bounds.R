library(Matrix)
# library("Rlab")

generate_preorder_data <- function(m, mu1, b, pi1, beta) {
  alt_proba <- pi1 * exp(-b * (1:m)*m**(beta - 1)) * (b / (1 - exp(-b)))
  alt_loc <- rbinom(m, 1, alt_proba)
  
  data <- rnorm(n = m) + mu1 * alt_loc

  pvalues <- pnorm(data, lower.tail = FALSE)
  output <- list(data = data, pvalues = pvalues)


  return(output)
  
}

generate_dummy_knockoff <- function(m, z, f1, beta=0) {
  # alt_probas <- 1/2 + pmax(0, 1/2 * ((z - ((1:m) * m**beta)) / (z-1)) )
  alt_probas <- 1/2 + pmax(0, 1/2 * ((z - ((1:m) * (m**(beta - 1)))) / (1 - z)))
  alt_pvalues <- sample(c(0.5, 1), m, TRUE, c(f1, 1 - f1))
  binary_pvalues <- numeric(m) + (rbinom(m, 1, alt_probas) * alt_pvalues)
  nb_null <- sum(binary_pvalues == 0)
  binary_pvalues[which(binary_pvalues == 0)] = 1 - 0.5 * rbinom(nb_null, 1, 1/2)
  
  return(binary_pvalues)
}


preorder_FDP_hat <- function(s, lambda, ordered_pvalues) {
  
  num <- cumsum(ordered_pvalues > lambda)
  denom <- pmax(1, cumsum(ordered_pvalues <= s))
  output <- list(fdp_hat = pmin(1, (s / (1 - lambda)) * ((1 + num) / denom)), Rk = denom)
  
  return(output)
}


LF_k_hat <- function(alpha, s, lambda, ordered_pvalues) {
  
  m <- length(ordered_pvalues)
  FDP_hats <- preorder_FDP_hat(s, lambda, ordered_pvalues)$fdp_hat
  
  if (sum(FDP_hats <= alpha) >= 1) {
    k_hat <- max(which(FDP_hats <= alpha))
    nb_rej <- sum(ordered_pvalues[1:k_hat] <= s)
    nb_rej_vec <- cumsum(ordered_pvalues[1:k_hat] <= s)
    output <- list(k_hat = k_hat, nb_rej = nb_rej, nb_rej_vec = nb_rej_vec)
    return(output)
    }
  else {
    k_hat <- 0
    nb_rej <- 0 
    output <- list(k_hat = k_hat, nb_rej = nb_rej)
    return(output)
  }
  
}

eps_function <- function(x, delta) {
  res <- log((1 + (pi**2 / 6)) / delta) + 2 * log(1 + max(0, log(x)/log(2)))
  return(res)
}

# from Li and Goeman
get_FDP_true <- function(W, true_index){
  if(sum(W==0)>0 | length(unique(W))!=length(W) | !identical(W[order(abs(W), decreasing=TRUE)], W)){
    stop("The input W might have ties or zeros or not sorted!")}
  
  p <- length(W)
  true_num_disc <- rep(0,p)
  for (i in 1:p) {
    count_temp <- 0
    for (j in 1:i) {
      if(j %in% true_index & W[j]>0) {count_temp <- count_temp+1}
    }
    true_num_disc[i] <- count_temp
  }
  
  num_disc <- cumsum(W>0)
  FDP_true <- (num_disc-true_num_disc)/pmax(num_disc,1)
  return(FDP_true)
}

#-------------------------------------------------------------------------------------------  
preorder_freedman <- function(k, s, lambda, delta, ordered_pvalues, return_Vk = FALSE) {
  
  beta = s * (1 + (min(s, lambda) / (1 - lambda)))
  eps_ <- eps_function(beta * k, delta)
  cte = 2 * sqrt(eps_) * sqrt(beta * k) + eps_
  
  nb_above_lambda <- sum(ordered_pvalues[1:k] > lambda)
  nb_below_s <- sum(ordered_pvalues[1:k] <= s)
  
  B <- s / (1 - lambda)
  
  
  if (return_Vk) return(B * nb_above_lambda + cte)
  
  else return(min(1, (B * nb_above_lambda + cte) / max(1, nb_below_s)))
  
}

preorder_KR <- function(k, s, lambda, delta, ordered_pvalues, a=1, return_Vk = FALSE) {
  
  B <- s / (1 - lambda)
  cte <- (log(1 / delta)) / (a * log(1 + ((1 - delta**(B/a)) / B )))
  
  nb_above_lambda <- sum(ordered_pvalues[1:k] > lambda)
  nb_below_s <- sum(ordered_pvalues[1:k] <= s)
  
  if (return_Vk) return (Vk = cte * (a + B * nb_above_lambda))
  
  else return(min(1, cte * ((a + B * nb_above_lambda) / max(1, nb_below_s))))
  
}


preorder_KR_optim <- function(a_max, k, s, lambda, delta, ordered_pvalues, return_Vk = FALSE) {

  
  # nb_below_s <- sum(ordered_pvalues[1:k] <= s)
  # a <- round(a)
  # delta <- delta / ((pi**2 / 6) * a**2)
  # B <- s / (1 - lambda)
  # cte <- (log(1 / delta)) / (a * log(1 + ((1 - delta**(B/a)) / B )))
  
  kapa = (pi**2 / 6)
  B = s / (1 - lambda)
  a_seq <- 1:a_max
  delta_a_seq <- delta / (kapa * a_seq**2)
  cte_seq <- log(1 / delta_a_seq) / (a_seq * log(1 + ((1 - delta_a_seq**(B / a_seq)) / B)))
  nb_above_lambda <- sum(ordered_pvalues[1:k] > lambda)
  bound_seq <- cte_seq * (a_seq + B * nb_above_lambda)
  
  if (return_Vk) return(bound_seq)
  else {
    nb_below_s <- sum(ordered_pvalues[1:k] <= s)
    return(pmin(1, bound_seq / max(1, nb_below_s)))
  }
}


preorder_KR_U <- function(k, s, lambda, delta, ordered_pvalues, return_Vk = FALSE) {
  
  nb_below_s <- sum(ordered_pvalues[1:k] <= s)
  
  if(nb_below_s <= 1) {
    delta <-  (6 * delta) / (pi**2)
    bound <- preorder_KR(k, s, lambda, delta, ordered_pvalues, return_Vk = return_Vk)
    best_a <- 1
  }  
  else {
    bound_seq <- preorder_KR_optim(nb_below_s, k, s, lambda, delta, ordered_pvalues, return_Vk)
  # optim <- optimize(preorder_KR_optim, c(1, nb_below_s), k, s, lambda, delta,
  #                   ordered_pvalues, return_Vk, lower = 1, upper = nb_below_s, maximum = FALSE)
  # 
  # bound <- optim$objective
  # best_a <- optim$minimum
    
  bound <- min(bound_seq)
  best_a <- which.min(bound_seq)
  }
  
  output <- list(bound = bound, a = round(best_a))
  return(output)
  
}

# from Li and Goeman
get_FDP_KJI <- function(W, k_vec, v_vec){
  if(sum(W==0)>0 | length(unique(W))!=length(W) | !identical(W[order(abs(W), decreasing=TRUE)],W)){
    stop("The input W might have ties or zeros or not sorted!")}
  m <- length(k_vec)
  S_list <- lapply(1:m, function(i_S, W, v_vec){
    candidate_set <- which(cumsum(W<0)==v_vec[i_S])
    ifelse(length(candidate_set)==0, threshold <- min(abs(W)), threshold <- abs(W[candidate_set[1]]))
    S <- which(W>=threshold)
    return(S)
  }, W, v_vec)
  
  ###
  p <- length(W)
  FDP_bound_vec <- rep(1,p)
  for (i in 1:p) {
    R <- which(W >= abs(W[i]))
    
    if(length(R)==0){FDP_bound_vec[i] <- 0; next}  # this means all W are negative, so discovery set is empty, so FDP is 0.
    
    FDP_k_temp <- mapply(function(j, p, R, S_list, k_vec){
      S_temp <- S_list[[j]]
      return( min(length(R), k_vec[j]-1+length(setdiff(R,S_temp))) / max(1,length(R)) )
    }, 1:m, MoreArgs=list(p, R, S_list, k_vec))
    
    FDP_bound_vec[i] <- min(FDP_k_temp)
  }
  return(FDP_bound_vec)
}

#-------------------------------------------------------------------------------------------

preorder_freedman_LF <- function(alpha, delta, s, lambda, k_hat, ordered_pvalues, return_Vk = FALSE) {
  beta = s * (1 + (min(s, lambda) / (1 - lambda)))
  eps_ <- eps_function(beta * k_hat, delta)
  cte = 2 * sqrt(eps_) * sqrt(beta * k_hat) + eps_
  
  nb_below_s <- sum(ordered_pvalues[1:k_hat] <= s)
  
  if(return_Vk) return(alpha * max(1, nb_below_s) + cte)
  else return(min(1, alpha + (cte / max(1, nb_below_s))))
}

preorder_KR_LF <- function(alpha, delta, s, lambda, k_hat, ordered_pvalues, a=1, return_Vk = FALSE) {
  B <- s / (1 - lambda)
  cte <- (log(1 / delta)) / (a * log(1 + ((1 - delta**(B/a)) / B )))
  
  nb_below_s <- sum(ordered_pvalues[1:k_hat] <= s)
  
  if(return_Vk) return(cte * (alpha * max(1, nb_below_s) + 1))
  else return(min(1, cte * (alpha + (1 / max(1, nb_below_s)))))
  
}


preorder_KR_LF_optim <- function(a_max, alpha, delta, s, lambda, k_hat, ordered_pvalues, return_Vk = FALSE) {
  
  # a <- round(a)
  # delta <- delta / ((pi**2 / 6) * a**2)
  # B <- s / (1 - lambda)
  # cte <- (log(1 / delta)) / (a * log(1 + ((1 - delta**(B/a)) / B )))
  # nb_below_s <- sum(ordered_pvalues[1:k_hat] <= s)
  # 
  # if(return_Vk) return(cte * (alpha * max(1, nb_below_s) + a ))
  # return(min(1, cte * (alpha + (a / max(1, nb_below_s)))))
  
  kapa = (pi**2 / 6)
  B = s / (1 - lambda)
  a_seq <- 1:a_max
  delta_a_seq <- delta / (kapa * a_seq**2)
  cte_seq <- log(1 / delta_a_seq) / (a_seq * log(1 + ((1 - delta_a_seq**(B / a_seq)) / B)))
  # nb_above_lambda <- sum(ordered_pvalues[1:k_hat] > lambda)
  nb_below_s <- sum(ordered_pvalues[1:k_hat] <= s)
  
  if (return_Vk) {
    bound_seq <- cte_seq * (alpha * max(1, nb_below_s) + a_seq)
    return(bound_seq)
  }
  else {
    
    bound_seq <- cte_seq * (alpha + (a_seq / max(1, nb_below_s)))
    return(pmin(1, bound_seq))
  }
  
}

preorder_KR_U_LF <- function(alpha, delta, s, lambda, k_hat, ordered_pvalues, return_Vk = FALSE) {
  
  nb_rej <- sum(ordered_pvalues[1:k_hat] <= s)
  
  if(nb_rej <= 1) {
    delta <-  (6 * delta) / (pi**2)
    bound <- preorder_KR_LF(alpha, delta, s, lambda, k_hat, ordered_pvalues, return_Vk = return_Vk)
    best_a <- 1 
  
  }  
  else {
    bound_seq <- preorder_KR_LF_optim(nb_rej, alpha, delta, s, lambda, k_hat, pvalues, return_Vk)
    
    # optim <- optimize(preorder_KR_LF_optim, c(1, nb_rej), alpha, delta, s, lambda, k_hat,
    #                 ordered_pvalues, return_Vk, lower = 1, upper = nb_rej, maximum = FALSE)
    # bound <- optim$objective
    # best_a <- optim$minimum
    
    bound <- min(bound_seq)
    best_a <- which.min(bound_seq)

  }
    output <- list(bound = bound, a = round(best_a))
    return(output)
  
}


#-------------------------------------------------------------------------------------------
get_sorted_signed_stat <- function(knockoff_stat){
  abs_knockoffstat <- abs(knockoff_stat)
  pi_knockoffstat <- order(abs_knockoffstat, decreasing = TRUE)
  sort_knockoffstat <- sort(abs_knockoffstat, decreasing = TRUE)
  
  init_sgn <- (knockoff_stat <= 0) * (-1) + (knockoff_stat > 0) * 1
  permutation_mat <- as.matrix(sparseMatrix(seq_along(pi_knockoffstat), pi_knockoffstat, x=1))
  sort_sgn <- as.vector(permutation_mat %*% init_sgn)
  
  return(sort_knockoffstat * sort_sgn)
}

knockoff <- function(alpha, knockoff_stat) {
  abs_knockoffstat <- abs(knockoff_stat)
  pi_knockoffstat <- order(abs_knockoffstat, decreasing = TRUE)
  sort_knockoffstat <- sort(abs_knockoffstat, decreasing = TRUE)
  
  init_sgn <- (knockoff_stat <= 0) * (-1) + (knockoff_stat > 0) * 1
  
  permutation_mat <- as.matrix(sparseMatrix(seq_along(pi_knockoffstat), pi_knockoffstat, x=1))
  sort_sgn <- as.vector(permutation_mat %*% init_sgn)
  
  pos_neg_ratio <- (1 + cumsum((sort_knockoffstat * sort_sgn) < 0)) / pmax(1, cumsum((sort_knockoffstat * sort_sgn) > 0))
  k_hat = length(which(pos_neg_ratio <= alpha))
  
  return(k_hat)
  
}

