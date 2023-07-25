eps <- function(x, delta) {
  return(log((1 + (pi**2 / 6)) / delta) + 2 * log(1 + max(0, log(x)/log(2))))
}

log_bar <- function(delta, a) {
  return((log(1 / delta)) / (a * log(1 + (1/a) * log(1 / delta))))
}


online_KR <- function(k, cvs, delta, pvalues, a=1, return_Vk = FALSE) {
  
  nb_rej_k <- sum(pvalues[1:k] <= cvs[1:k])
  sum_cvs <- sum(cvs[1:k])
  cte <- log_bar(delta, a)
  if(return_Vk) return(cte * (a + sum_cvs))
  else return(min(1, cte * (a + sum_cvs) / max(1, nb_rej_k)))
}

online_KR_optim <- function(a_max, k, cvs, delta, pvalues, return_Vk = FALSE) {
  
  # a <- round(a)
  # delta_ <- delta / ((pi**2 / 6) * a**2)
  
  kapa = (pi**2 / 6)
  a_seq <- 1:a_max
  delta_a_seq <- delta / (kapa * a_seq**2)
  cte_seq <- log_bar(delta_a_seq, a_seq)
  nb_rej_k <- sum(pvalues[1:k] <= cvs[1:k])
  sum_cvs <- sum(cvs[1:k])
  
  if(return_Vk) return(cte_seq * (a_seq + sum_cvs))
  else return(pmin(1, cte_seq * (a_seq + sum_cvs) / max(1, nb_rej_k)))
}

online_KR_U <- function(k, cvs, delta, pvalues, return_Vk = FALSE) {
  
  nb_rej_k <- sum(pvalues[1:k] <= cvs[1:k])
  
  if(nb_rej_k <= 1) {
    delta_ <- (6 * delta) / (pi**2)
    bound <- online_KR(k, cvs, delta_, pvalues, return_Vk = return_Vk)
    best_a <- 1
  }  
  else {
    bound_seq <- online_KR_optim(nb_rej_k, k, cvs, delta, pvalues, return_Vk)
    # optim <- optimize(online_KR_optim, c(1, nb_rej_k), k, cvs, delta, pvalues, return_Vk,
    #                   lower = 1, upper = nb_rej_k, maximum = FALSE)
    # 
    # bound <- optim$objective
    # best_a <- optim$minimum
    
    bound <- min(bound_seq)
    best_a <- which.min(bound_seq)
    
  }
  output <- list(bound = bound, a = round(best_a))
  return(output)
}


online_Freedman <- function(k, cvs, delta, pvalues, return_Vk = FALSE) {
  
  nb_rej_k <- sum(pvalues[1:k] <= cvs[1:k])
  sum_cvs <- sum(cvs[1:k])
  eps_ <- sapply(sum_cvs, eps, delta)
  cte = 2 * sqrt(eps_) * sqrt(sum_cvs) + eps_ 
  
  if(return_Vk) return(sum_cvs + cte)
  else return(min(1, (sum_cvs + cte) / nb_rej_k))
}


