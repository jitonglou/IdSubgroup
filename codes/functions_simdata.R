#### Generate measurements ####
f_tpint = function(t_up, t_low, Lambda_k=NULL,
                   L_ik=NULL, Y_ik_mean=NULL, eta_k=NULL){
  val_L = ifelse(is.null(eta_k), 1, exp(c(L_ik(Y_ik_mean) %*% eta_k)))
  val_L*(Lambda_k(t_up) - Lambda_k(t_low))
}

f_tpsolve = function(t, t0, Lambda_k=NULL, val1=NULL,
                     L_ik=NULL, Y_ik_mean=NULL, eta_k=NULL){
  val_L = ifelse(is.null(eta_k), 1, exp(c(L_ik(Y_ik_mean) %*% eta_k)))
  Lambda_k(t) - Lambda_k(t0) - val1/val_L
}


eps_covmat_block = function(p, t_all_vec, t_next, b, Omega_next){
  mat_coef = exp(-((t_all_vec - t_next)/b)^2)
  
  (t_all_vec %>% lapply(function(t_exist){
    (Omega(t_exist) + Omega_next)/2
  }) %>% do.call(rbind, .)) * rep(mat_coef, each=p)
}

gen_nextmk = function(k_idx, t_next, X_i, Omega=NULL, 
                      g_inv_k=NULL, beta_k=NULL, f_k=NULL, phi_k=NULL,
                      t_all_vec=NULL, eps_all_vec=NULL, eps_covmat_joint=NULL, b){
  if (length(t_all_vec) == 0){
    # if t_next is the first observation of this subject
    eps_covmat_joint = Omega(t_next)
    
    eps_next_vec = mvrnorm(n=1, mu=c(0,p), Sigma=eps_covmat_joint)
    eps_next = eps_next_vec[k_idx]
    
    mu_next = g_inv_k(c(X_i %*% beta_k(t_next)) + eps_next)
    
    Y_next = f_k(mu_next, phi_k)
  } else {
    Omega_next = Omega(t_next)
    # the covariance matrix of the joint distribution of eps_all_vec and eps_next is
    # (eps_covmat_joint covmat_block )
    # (t(covmat_block)  Omega_next   )
    
    covmat_block = eps_covmat_block(p, t_all_vec, t_next, b, Omega_next)
    
    # sample eps_next_vec given eps_all_vec using the conditional distribution
    eps_next_vec = mvrnorm(n = 1, 
                           mu = t(covmat_block) %*% solve(eps_covmat_joint) %*% eps_all_vec %>% c, 
                           Sigma = Omega_next - t(covmat_block) %*% solve(eps_covmat_joint) %*% covmat_block
    )
    eps_next = eps_next_vec[k_idx]
    
    mu_next = g_inv_k(c(X_i %*% beta_k(t_next)) + eps_next)
    
    Y_next = f_k(mu_next, phi_k)
    
    eps_covmat_joint = eps_covmat_joint %>% 
      cbind(covmat_block) %>%
      rbind(cbind(t(covmat_block), Omega(t_next)))
  }
  return(list(eps_covmat_joint=eps_covmat_joint, eps_next_vec=eps_next_vec,
              eps_next=eps_next, mu_next=mu_next, Y_next=Y_next))
}


# tau_i=t_max = 12; i=1; p=2; X_i=X[i,]
# gen_1st_tp %>%
#   mapply(k = 1:p, tau_i=list(tau_i), X_i=list(X_i), 
#          gamma_k=lapply(1:p, function(k){gamma_p[,k]}), 
#          Lambda_k=Lambda_p)

gen_1st_tp = function(tau_i, X_i, gamma_k, Lambda_k){
  ## generate first tp for each marker
  
  t_ik = NA
  Z_i = X_i[-1]
  U = runif(1) # u smaller, root bigger
  # U = 0.5
  t_up = tau_i
  t_low = 0
  val1 = -log(U)*exp(-c(Z_i %*% gamma_k))
  val2 = f_tpint(t_up,t_low,Lambda_k)
  
  if (val1 <= val2){
    # must exists a root
    fit = brent(f_tpsolve, a=t_low, b=t_up, t0=t_low,
                Lambda_k=Lambda_k, val1=val1)
    t_ik = fit$root
  } # end of if val1 < val2 at beginning  
  
  return(t_ik)
}

gen_nexttp = function(tau_i, X_i, gamma_k, Lambda_k,
                      L_pre_ik=NULL, L_ik=NULL, eta_k=NULL,
                      t_now=NULL, t_ik=NULL, Y_ik=NULL){
  t_next = NA
  Z_i = X_i[-1]
  U = runif(1) # u smaller, root bigger
  # U = 0.5
  
  val1 = -log(U)*exp(-c(Z_i %*% gamma_k))
  
  t_pre = t_ik[t_ik > t_now - L_pre_ik] # last is t_now
  Y_pre = Y_ik[t_ik > t_now - L_pre_ik]
  n_pre = length(t_pre) # n_pre>=1
  
  if (n_pre==1){
    t_low_vec = t_now
    t_up_vec = t_now + L_pre_ik        
  } else {
    t_low_vec = c(t_now, t_pre[-n_pre]+L_pre_ik)
    t_up_vec = c(t_pre+L_pre_ik)        
  }
  
  ## case 1, t_next between t_now and min(t_now+L_pre_ik, tau_i)
  for (j in 1:n_pre){
    t_low = t_low_vec[j]
    # t_up = t_up_vec[j] # t_up_vec[j] might be > tau_i
    t_up = min(c(t_up_vec[j], tau_i))
    Y_ik_mean = mean(Y_pre[j:n_pre])
    
    val2 = f_tpint(t_up, t_low, Lambda_k, 
                   L_ik, Y_ik_mean, eta_k)
    
    if (val1 <= val2){
      fit = brent(f_tpsolve, a=t_low, b=t_up, t0=t_low,
                  Lambda_k=Lambda_k, val1=val1,
                  L_ik=L_ik, Y_ik_mean=Y_ik_mean, eta_k=eta_k)
      
      t_next = fit$root
      val1 = 0
      break # break j loop
    } else {
      val1 = val1-val2 # >0
      t_low = t_up
      next # iteration j+1
    } # end of ifelse(val1 <= val2)
  } # end of j loop
  
  ## case 2, t_next between min(t_now+L_pre_ik, tau_i) and tau_i
  if (val1 > 0 & val1 <= f_tpint(t_up=tau_i,t_low,Lambda_k)){
    # root > t_now + L_pre_ik 
    # t_low = t_now + L_pre_ik
    
    fit = brent(f_tpsolve, a=t_low, b=tau_i, t0=t_low,
                Lambda_k=Lambda_k, val1=val1,
                L_ik=NULL, Y_ik_mean=NULL, eta_k=NULL)
    
    t_next = fit$root
    # end of solving  t, which is > t_now + L_pre_ik and < tau_i
  }
  
  ## otherwise, no solution, return t_next=NA
  
  return(t_next)
}

gen_tp_v3 <- function(p, tau_i, X_i, gamma_p, Lambda_p,
                      Omega=NULL, g_inv_p=NULL, beta_p=NULL, f_p=NULL, phi_p=NULL,
                      L_pre_ip=NULL, L_ip=NULL, eta_p=NULL, b=0.01){
  ## 2020/7/17, generate eps for all markers of a subject together
  #{
  # 0. all results are NULL
  # 1. gen first tp vector (potential tp vector)
  # while not every element in the tp vector is NA, do
  # {
  # 2. select the smallest tp, generate the observations (given previous tp, eps, mk) at this tp
  # 3. generate next tp for this mk, and update potential tp vector
  # }
  # 4. return results
  #}
  
  ## step 1
  # may merge gen_1st_tp() to the general function gen_nexttp() later
  t_potential_vec = gen_1st_tp %>%
    mapply(tau_i=list(tau_i), X_i=list(X_i), 
           gamma_k=lapply(1:p, function(k){gamma_p[,k]}), 
           Lambda_k=Lambda_p) # a vec of len p (element can be NA)
  
  # intialize outputs
  t_ip = Y_ip = eps_ip = mu_ip = vector("list", length=p)
  t_all_vec = eps_all_vec = eps_covmat_joint = NULL
  
  while((!(is.na(t_potential_vec))) %>% any){
    # print(t_potential_vec)
    
    ## step 2
    k_idx = which.min(t_potential_vec)
    
    t_next = t_potential_vec[k_idx]
    t_ip[[k_idx]] = append(t_ip[[k_idx]], t_next)
    
    # result_next has 5 sublists: eps_covmat_joint, eps_next_vec, eps_next, mu_next, Y_next
    # smaller b -> smaller correlation between different tps -> eps_covmat_joint tend to be positive definite
    result_next = gen_nextmk(k_idx, t_next, X_i, Omega, 
      g_inv_k=g_inv_p[[k_idx]], beta_k=beta_p[[k_idx]], 
      f_k=f_p[[k_idx]], phi_k=phi_p[[k_idx]],
      t_all_vec, eps_all_vec, eps_covmat_joint, b)
    
    eps_ip[[k_idx]] = append(eps_ip[[k_idx]], result_next$eps_next)
    
    mu_ip[[k_idx]] = append(mu_ip[[k_idx]], result_next$mu_next)
    
    Y_ip[[k_idx]] = append(Y_ip[[k_idx]], result_next$Y_next)
    
    # update save all tps and eps to calculate the covariance matrix for all eps at tps
    t_all_vec = append(t_all_vec, t_next)
    eps_all_vec = append(eps_all_vec, result_next$eps_next_vec) # vector, length=p*length(t_all_vec)
    eps_covmat_joint = result_next$eps_covmat_joint # matrix, nrow=ncol=p*length(t_all_vec)
    
    ## step 3
    t_potential_vec[k_idx] = gen_nexttp(tau_i, X_i, 
      gamma_k=gamma_p[,k_idx], Lambda_k=Lambda_p[[k_idx]],
      L_pre_ik=L_pre_ip[[k_idx]], L_ik=L_ip[[k_idx]], eta_k=eta_p[,k_idx],
      t_now=t_next, t_ik=t_ip[[k_idx]], Y_ik=Y_ip[[k_idx]])
  } # end of condition while(t_potential_vec %>% is.na %>% not %>% any)
  return(list(t_ip, Y_ip, eps_ip, mu_ip))
}


