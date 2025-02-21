generate_pvalues <- function(m, mu_1, pi0, seed=NULL) {
  #-----------------------------------------------------------------------------
  if (!is.null(seed)){
    set.seed(seed)
  }
  #-----------------------------------------------------------------------------
  
  true_indices <- rbinom(m, 1, 1-pi0)
  observations <- rnorm(m) + mu_1 * true_indices
  raw.pvalues <- pnorm(observations, lower.tail = FALSE)
  
  data <- list(true_indices = true_indices, observations = observations, raw = raw.pvalues)
  return(data)
}


generate_DU_pvalues <- function(m, pi0, seed = NULL) {
  
  return(c(rep(0, m * (1 - pi0)), runif(m * pi0)))
}


func_h <- function(x) {
  return(x * (log(x) - 1) + 1)
}

inverse <- function (f, lower = 1, upper = 9e10) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}

BH_k_hat <- function(alpha, pvalues) {
  
  k_hat=0
  m = length(pvalues)
  fdp_hat_topk <- (m * sort(pvalues)) / (1:m)
  set = which(fdp_hat_topk <= alpha) 
  if(length(set)>0) k_hat = max(set)
  
  return(k_hat)
}

m0_hat_vec_Simes <- function(t, pvalues, delta){
  # V_t <- sum(pvalues > t)
  # sorted_pvalues <- sort(pvalues)
  m = length(pvalues)
  V_t <- m - t
  return(V_t / (1- (pvalues[t]/delta)))
}

m0_hat_vec_DKW <- function(t, pvalues, delta){
  # V_t <- sum(pvalues > t)
  # sorted_pvalues <- sort(pvalues)
  m = length(pvalues)
  V_t <- m - t
  C <- log(1/delta) / 2
  return( (( sqrt(C) / (2 * (1 - pvalues[t]))) + sqrt( (C / (4*(1 - pvalues[t])**2)) + (V_t / (1 - pvalues[t])) ))**2 )
  
}

m0_hat_vec_KR <- function(t, pvalues, delta) {
  # V_t <- sum(pvalues > t)
  # sorted_pvalues <- sort(pvalues)
  m = length(pvalues)
  V_t <- m - t
  C <- log(1/delta) / log(1 + log(1/delta))
  
  return((C + V_t) / (1 - C*pvalues[t]))
}

m0_hat_vec_Well <- function(t, pvalues, delta) {
  # V_t <- sum(pvalues > t)
  # sorted_pvalues <- sort(pvalues)
  m = length(pvalues)
  V_t <- m - t
  C_t <- 2*log(pi**2 / (6*delta)) + 4*log(1 + log(1/pvalues[t], base = 2))
  
  if (pvalues[t] == 1){
    return(m)
  }
  else {
    return( (( sqrt(C_t * pvalues[t]) / (sqrt(2) * (1 - pvalues[t]))) + sqrt( (C_t / (2*(1 - pvalues[t])**2)) + (V_t / (1 - pvalues[t])) ) )**2 )
  }
}

#---------------------------------------------------------------------------------------------------

topk_Wellner <- function(k, delta, pvalues, adaptive = FALSE, hatm0 = NULL, return_Vk = FALSE) {
  
  if (adaptive) {
    nb_pvalues = hatm0
  }
  else {
    nb_pvalues = length(pvalues)
  }
  p_topk <- sort(pvalues)[k]
  c = pi**2 / (6 *delta)
  value = (2 * log(c) + 4 * log(1 + log(1/p_topk, base = 2))) / (nb_pvalues * p_topk)
  h_inv <- inverse(func_h, lower = 1, upper = value * 1e4)
  const = h_inv(value)$root
  bound = const * p_topk * nb_pvalues 
  
  if (return_Vk) return(bound)
  else return(min(1, bound / k))
}

topk_KR <- function(k, delta, pvalues, adaptive = FALSE, hatm0 = NULL, return_Vk = FALSE) {
  
  if (adaptive) {
    nb_pvalues = hatm0
  }
  else {
    nb_pvalues = length(pvalues)
  }
  
  p_topk <- sort(pvalues)[k]
  const = log(1 / delta) / log(1 + log(1/delta))
  
  if (return_Vk) return(const * (1 + nb_pvalues * p_topk))
  else return(min(1, (const * (1 + nb_pvalues * p_topk)) / k))
}

topk_Simes <- function(k, delta, pvalues, adaptive = FALSE, hatm0 = NULL, return_Vk = FALSE) {
  
  if (adaptive) {
    nb_pvalues = hatm0
  }
  else {
    nb_pvalues = length(pvalues)
  }
  
  p_topk <- sort(pvalues)[k]
  
  if (return_Vk) return((nb_pvalues * p_topk) / delta)
  else return(min(1, (nb_pvalues * p_topk) / (k * delta )))
}

topk_DKW <- function(k, delta, pvalues, adaptive = FALSE, hatm0 = NULL, return_Vk = FALSE) {
  
  if (adaptive) {
    nb_pvalues = hatm0
  }
  else {
    nb_pvalues = length(pvalues)
  }
  
  p_topk <- sort(pvalues)[k]
  
  if (return_Vk) return(p_topk * nb_pvalues + sqrt(nb_pvalues * 0.5 * log(1 / delta)))
  else return(min(1, (p_topk * nb_pvalues + sqrt(nb_pvalues * 0.5 * log(1 / delta))) / k))
}

topk_hybrid_KR_Well <- function(k, delta, pvalues, adaptive = FALSE, hatm0_KR= null, hatm0_Well = null, return_Vk = FALSE) {
  
  KR_bound <- topk_KR(k, delta / 2, pvalues, adaptive, hatm0_KR, return_Vk)
  Well_bound <- topk_Wellner(k, delta / 2, pvalues, adaptive, hatm0_Well, return_Vk)
  
  if (return_Vk) return(min(KR_bound, Well_bound))
  else return(min(1, min(KR_bound, Well_bound)))
  
}

#-----------------------------------------------------------------------------------------------------

topk_Simes_BH <- function(alpha, delta, adaptive = FALSE, hatpi0 = NULL,  k_hat = NULL, return_Vk = FALSE) {
  
  if(!adaptive) hatpi0 = 1
  
  if (return_Vk) return(k_hat * hatpi0 * (alpha/delta))
  else return(min(1, hatpi0 * (alpha/delta)))
}

topk_DKW_BH <- function(alpha, delta, m, k_hat, adaptive = FALSE, hatpi0 = NULL, return_Vk = FALSE) {
  
  if (!adaptive) hatpi0 = 1
  bound = hatpi0 * alpha + (sqrt(hatpi0 * m * 0.5 * log(1/delta)) / max(1, k_hat))
  
  if (return_Vk) return(k_hat * bound)
  else return(min(1, bound))
}

topk_KR_BH <- function(alpha, delta, k_hat, adaptive = FALSE, hatpi0 = NULL, return_Vk = FALSE) {
  
  if (!adaptive) hatpi0 = 1
  bound =  (log(1/delta) / log(1 + log(1/delta))) * (hatpi0 * alpha + (1 / max(1, k_hat)))
  
  if(return_Vk) return(k_hat * bound)
  else return(min(1, bound))
}

topk_Wellner_BH <- function(alpha, delta, m, k_hat, adaptive = FALSE, hatpi0 = NULL, return_Vk = FALSE) {
  
  if (!adaptive) hatpi0 = 1
  h_inv <- inverse(func_h)
  c = pi**2 / (6 * delta)
  bound = hatpi0 * alpha * h_inv((2 * log(c) + 4 * log(1 + log(m/(alpha * max(1, k_hat)), base = 2))) / (alpha * max(1, k_hat) * hatpi0))$root
  
  if(return_Vk) return(k_hat * bound)
  else return(min(1, bound))
}

topk_hybrid_KR_Well_BH <- function(alpha, delta, m, k_hat, adaptive = FALSE, hatpi0_KR = NULL, hatpi0_Well = NULL, return_Vk = FALSE) {
  
  KR_bound <- topk_KR_BH(alpha, delta / 2, k_hat, adaptive, hatpi0_KR, return_Vk)
  Well_bound <- topk_Wellner_BH(alpha, delta / 2, m, k_hat, adaptive, hatpi0_Well, return_Vk)
  
  if (return_Vk) return(min(KR_bound, Well_bound))
  else return(min(1, min(KR_bound, Well_bound)))
  
}

#-----------------------------------------------------------------------------------------------------


m0_hat_Simes <- function(pvalues, delta){
  m = length(pvalues)
  index_vec <- which(sort(pvalues) < delta)
  if (length(index_vec) == 0) return(m)
  else {
    m0_hat_vec <- sapply(index_vec, m0_hat_vec_Simes, sort(pvalues), delta)
    return(min(m, min(m0_hat_vec))) 
  }
}

m0_hat_DKW <- function(pvalues, delta) {
  m = length(pvalues)
  m0_hat_vec <- sapply(1:m, m0_hat_vec_DKW, sort(pvalues), delta)
  return(min(m, min(m0_hat_vec)))
}

m0_hat_KR <- function(pvalues, delta) {
  m = length(pvalues)
  C <- log(1/delta) / log(1 + log(1/delta))
  index_vec <- which(sort(pvalues) < (1 / C))
  if (length(index_vec) == 0) return(m)
  else {
    m0_hat_vec <- sapply(index_vec, m0_hat_vec_KR, sort(pvalues), delta)
    return(min(m, min(m0_hat_vec)))
  }

}

m0_hat_Well <- function(pvalues, delta) {
  m = length(pvalues)
  m0_hat_vec <- sapply(1:m, m0_hat_vec_Well, sort(pvalues), delta)
  return(min(m, min(m0_hat_vec)))
}