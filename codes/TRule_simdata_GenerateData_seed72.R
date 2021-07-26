rm(list=ls())

#### Libraries ####
library(statmod) # gauss.quad()
library(clusterGeneration) # genPositiveDefMat()
library(MASS) # mvrnorm(), select is masked from dplyr
library(purrr) # map_dfr()

library(dplyr) # easier data wrangling 
library(pracma) # brentDekker(), newtonRaphson(), ridders()
library(survival) # coxph()
library(rlist) # list.append

#### Set directories ####
jobname = "TRule_simdata_GenerateData"
rseed = 72; rsd = paste0("seed",rseed) # random seed used to generate one data set
casename = paste(jobname,rsd,sep="_")

maindir = "/pine/scr/j/i/jitong/thesis/TRule_model" # project directory
codesdir = file.path(maindir,"codes") # code directory
wsdir = file.path(maindir,"rdata") # R workspace directory

#### Load files ####
source(file=file.path(codesdir,"functions_simdata.R"))
source(file=file.path(codesdir,"functions_estparams.R"))


#### Settings of variables and parameters ####
## Following settings are described in the page 1 of the supplementary material
N = 10000 # number of subjects
p = 2 # number of health markers
m = 2 # number of pre-treatment covariates
t_max = 12 # maximum observation time for each subject
t.grids = 0:12 # time points at which estimate parameters

## R functions used to generate measurements of health markers
f_p = list(
  function(mu, phi){rnorm(1,mean=mu,sd=sqrt(phi))}, # normal
  function(mu, phi){rbinom(1,1,prob=mu)} # bernoulli
)

phi_p = c(0.1, NA) # dispersion parameters

## Inverse functions of canonical links and their derivatives for each health marker
g_inv_p = list(
  function(x){x}, # normal
  function(x){exp(x)/(1+exp(x))} # bernoulli
)

g_inv_deriv_p = list(
  function(x){rep(1,length(x))}, # normal
  function(x){exp(x)/(1+exp(x))^2} # bernoulli
)

## Baseline intensity functions
Lambda_p = list(
  function(u){1*u},
  function(u){1.2*u}
)

## A m*p matrix, each column is gamma parameters of the corresponding counting process
gamma_p = matrix(rep(c(0.5,0.25),p), ncol=p, byrow=F)

## A matrix of p columns, each column is eta parameters of the corresponding counting process
eta_p =  matrix(c(0.3, -0.1,
                  0.3, -0.1), ncol=p, byrow=F)

## functions used to calculate observed health history
L_ip = list(
  function(y){ifelse(is.null(y), return(c(0,0)), return(c(1,y)))},
  function(y){ifelse(is.null(y), return(c(0,0)), return(c(1,y)))}
)

L_pre_ip = rep(3,p) # length of time to calculate observed health history

## beta parameters (at time t)
beta_p = list(
  function(t){c(-1.36+t/10, sin(0.76+t), cos(-0.3+t))},
  function(t){c(cos(-0.25+t), 0.37+t/10, sin(-0.68+t))}
)

## Covariance structure of multivariate latent processes
Omega_0 = matrix(c(0.5,-0.25,-0.25,0.5),nrow=p)
Omega = function(t, scale_vec = sqrt(c(1/20,1/20))){
  Omega_0 +
    diag(scale_vec) %*% matrix(c(sin(t+2),cos(t-0.5),
                                 cos(t-0.5),sin(t+3)
    ), nrow=p) %*% diag(scale_vec)
}
b = 0.5 # a parameter used in generating Cov(epsilon(t),epsilon(s))

#### Generate data ####
## Generate pre-treatment covariates of subjects, X_i, a N*(m+1) matrix, each row for each subject
set.seed(rseed)
X = data.frame(int=1, age=rnorm(N, 0, 1/3), sex=rbinom(N, 1, 0.5)-0.5) %>% as.matrix

# ## A list of length N. Each element for each subject.
# In each element, there are four sublists of length p. 
# The first sublist records measurement time points for all health markers.
# The second sublist records measurements for all health markers.
# The third sublist records latent Gaussian processes for all health markers.
# The fourth sublist records the conditional means of measurements for all health markers.
obs_all = vector("list", length=N)
set.seed(rseed)
for (i in 1:N){
  results = try(
    gen_tp_v3(p, tau_i=t_max, X_i=X[i,], gamma_p, Lambda_p,
              Omega, g_inv_p, beta_p, f_p, phi_p,
              L_pre_ip, L_ip, eta_p, b),
    silent = TRUE)
  if (class(results) == "list"){
    obs_all[[i]] = results
  } else {
    obs_all[[i]] = vector("list", length=4)
  }
  if (i%%1000==0){print(paste(i, "is done"))}
}

## Restructure the above list by health marker
obs_p_all = 1:p %>% lapply(function(k){
  obs_all %>% lapply(function(obs_i){
    obs_i %>% lapply(function(obs_ip){
      obs_ip[[k]]
    }) # length 4
  }) # length N
}) # length p


#### Convert lists to inputs (data frames) ####
time_Xikj_all_list = 1:p %>%
  lapply(function(k){
    data.frame(i = rep(obs_p_all[[k]] %>% length %>% seq_len,
                       obs_p_all[[k]] %>% sapply(function(z){z[[1]] %>% length})),
               k = k,
               j = obs_p_all[[k]] %>% lapply(function(z){z[[1]] %>% length %>% seq_len}) %>% unlist,
               time = obs_p_all[[k]] %>% lapply(function(z){z[[1]]}) %>% unlist,
               value = obs_p_all[[k]] %>% lapply(function(z){z[[2]]}) %>% unlist,
               eps = obs_p_all[[k]] %>% lapply(function(z){z[[3]]}) %>% unlist,
               mu = obs_p_all[[k]] %>% lapply(function(z){z[[4]]}) %>% unlist
    )
  })

## Calculate health history information from inputs
time_Xikj_all_df_list = mapply(FUN=recent_summary, df=time_Xikj_all_list, t_win=L_pre_ip,
                               SIMPLIFY=FALSE)

#### Calculate sample variance for each health marker at each time point ####
var_vec_eps_lists = var_vec_lists = vector("list",length=length(t.grids))
for (mmm in 1:length(t.grids)){
  ttt = t.grids[mmm]
  var_vec_eps = var_vec = rep(0, p)
  for (kkk in 1:p){
    y_kkk = time_Xikj_all_list[[kkk]] %>%
      closest_Obs_eps(ttt, kkk)
    var_vec[kkk] = var(y_kkk[,2])
    var_vec_eps[kkk] = var(y_kkk[,3])
  }
  var_vec_lists[[mmm]] = var_vec
  var_vec_eps_lists[[mmm]] = var_vec_eps
}
names(var_vec_eps_lists) = names(var_vec_lists) = paste0("t=",t.grids)

#### Calculate sample correlation for each pair of health markers at each time point ####
cor_mat_eps_lists = cor_mat_lists = vector("list",length=length(t.grids))
for (mmm in 1:length(t.grids)){
  ttt = t.grids[mmm]
  cor_mat_eps = cor_mat = matrix(0, nrow=p, ncol=p)
  for (kkk in 1:(p-1)){
    for (lll in ((kkk+1):p)){
      y_kkk = time_Xikj_all_list[[kkk]] %>%
        semi_join(time_Xikj_all_list[[lll]], by="i") %>%
        closest_Obs_eps(ttt, kkk) # df(i, value_closest, eps_closest)
      y_lll = time_Xikj_all_list[[lll]] %>%
        semi_join(time_Xikj_all_list[[kkk]], by="i") %>%
        closest_Obs_eps(ttt, lll)
      cor_mat_eps[kkk,lll] = cor(y_kkk[,3],y_lll[,3]) # eps correlation
      cor_mat[kkk,lll] = cor(y_kkk[,2],y_lll[,2]) # obs correlation
    }
  }
  cor_mat_lists[[mmm]] = cor_mat
  cor_mat_eps_lists[[mmm]] = cor_mat_eps
}


#### Save above data ####
save(N, p, m, X, obs_p_all, t.grids, t_max,
     var_vec_eps_lists, var_vec_lists, cor_mat_eps_lists, cor_mat_lists,
     time_Xikj_all_df_list,
     beta_p, gamma_p, Lambda_p, Omega_0, Omega, g_inv_p, g_inv_deriv_p, f_p, phi_p,
     L_pre_ip, L_ip, eta_p,
     file=file.path(wsdir,paste0(casename,"_workspace.RData")))

#### Data for Andersen-Gill model ####
## Rename some objects
N0 = N
Obs_allbm = time_Xikj_all_df_list
X_idx = data.frame(id=1:N0, X)

## number of measurements for each subject and each health marker
nobs_ik = obs_p_all %>%
  lapply(function(obs_k_all){
    obs_k_all %>%
      sapply(function(obs_ki){
        length(obs_ki[[1]])
      })
  }) %>%
  do.call(cbind, .)
nobs_idx = which(rowSums(nobs_ik) == 0)

## convert the original data to the format of survival data
df_coxph_list = lapply(1:p, function(k){
  lapply(1:N0, function(i){
    my_SurvData_v2(y=obs_p_all[[k]][[i]][[2]], t=obs_p_all[[k]][[i]][[1]], 
                   t_win=L_pre_ip[k], t_max) %>%
      mutate(id = i, age = X[i,2], sex = X[i,3]) # will return a data frame even measurements(y) or time points (t) is NULL
  }) %>%
    do.call(rbind, .) %>%
    return()
})


# ## Save survival data for Andersen-Gill model
# save(df_coxph_list,file=file.path(wsdir,paste0(casename, "_intensity_df.RData")))

#### Fit Andersen-Gill model ####
gamma_est_summary = df_coxph_list %>%
  lapply(function(df_coxph_k){
    coxph(Surv(start, stop, event) ~ (age + sex + L1 + L2) + cluster(id),
          df_coxph_k %>% filter(stop-start>1e-6 & !(id %in% nobs_idx))
    ) # remove some observations taken at almost the same time
  })
gamma_est_mat = mapply(function(s){s$coefficient}, 
                       gamma_est_summary) # a matrix of column p, each column for each health marker

#### Save intensity parameter estimates ####
save(gamma_est_mat,
     file = file.path(wsdir,paste0(casename, "_intensity_est.RData")))

