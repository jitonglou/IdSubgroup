rm(list=ls())

#### Libraries ####
library(dplyr) # easier data wrangling 
library(magrittr) # %<>%
library(purrr) # reduce()
library(Rcpp) # sourceCpp()
library(RcppArmadillo)
library(dendextend) # set()

#### Set directories ####
jobname =  "TRule_simdata_DistMat"
rseed = 72; rsd = paste0("seed",rseed) # random seed used to generate one data set
cs = paste0("case", 0)
casename = paste(jobname,cs,rsd,sep="_")

maindir = "/pine/scr/j/i/jitong/thesis/TRule_model" # project directory
codesdir = file.path(maindir,"codes") # code directory
wsdir = file.path(maindir,"rdata") # R workspace directory

#### Load files ####
load(file.path(wsdir,paste0("TRule_simdata_GenerateData_",rsd,"_workspace.RData")))
load(file.path(wsdir,paste0("TRule_simdata_est_",cs,"_",rsd,".RData"))) # load Beta, Var, Cov estimates, Beta_tps, Var_tps, Cov_tps
source(file=file.path(codesdir,"functions_estparams.R"))
sourceCpp(file.path(codesdir,"functions_DistMat.cpp"), showOutput = TRUE)

#### Model settings ####
Obs_allbm = time_Xikj_all_df_list # measurements and time points
N = nrow(X) # number of subjects
m = ncol(X)-1 # number of baseline variables
p = length(Obs_allbm) # number of health markers
tps = t.grids[4:12] # interested time points

logy_list = list(
  logy_norm_C, logy_bern_C
) # a list of log likelihood functions of health markers, defined in functions_DistMat.cpp
no_phi = c(FALSE, TRUE) # whether the distribution of each health marker not have the dispersion parameter phi

#### Caculate the estimates of latent effects epsilon(t) ####
for (mmm in 1:length(tps)){
  ttt = tps[mmm]
  print(c(mmm, ttt))

  ## Calculate the covariance matrix of latent variables at time point ttt 
  omega_est = matrix(ncol=p,nrow=p)
  omega_est[lower.tri(omega_est)] = Cov_tps[[mmm]]
  omega_est[upper.tri(omega_est)] = t(omega_est)[upper.tri(t(omega_est))]
  diag(omega_est) = Var_tps[[mmm]]

  mgh = mgauss.hermite(n=5, mu=rep(0,p), sigma=omega_est) # multivariate Gauss-Hermite quadrature nodes and weights
  
  ## Extract the closest measurements near the time point ttt
  closest_Obs_list = lapply(Obs_allbm, function(x){
    closest_Obs_allbm(x, ttt)
  })
  
  closest_Obs_mat = data.frame(i = 1:N)
  for (i in 1:p){
    closest_Obs_mat %<>% full_join(closest_Obs_list[[i]], by="i")
  }
  
  flag = closest_Obs_mat[,-1] %>% 
    apply(2, function(col){
      !is.na(col)
    }) %>%
    apply(1, any) # exclude subjects who have no measurement for all health markers
  
  closest_Obs_mat_filter = closest_Obs_mat[flag,]
  
  closest_Y_mat = closest_Obs_mat_filter %>%
    dply::select_at(-1) %>%
    as.matrix() # a matrix of p columns, each column for each health marker
  
  ## Calculate estimates of measurements
  X_filter = X[closest_Obs_mat_filter$i,]
  lf_mat = X_filter %*% matrix(Beta_tps[[mmm]], ncol=p, byrow=FALSE) # linear effect matrix, n*p matrix, each column is x_i^T*beta_k 
  Y_hat_mat = sapply(1:p, function(k){g_inv_p[[k]](lf_mat[,k])})

  ## Calculate estimates of dispersion parameters
  var_Y_est = apply((closest_Y_mat - Y_hat_mat)^2, 2, function(x){
    mean(x, na.rm=T)
  })
  phi = var_Y_est; phi[no_phi] = 0
  
  ## Calculate the posterior mean of epsilon_ik given Y_ik 
  epsilon_est_mat = matrix(nrow=nrow(X_filter), ncol=p)

  for (i in 1:nrow(X_filter)){
    logfyeps = apply(mgh$points, 1, function(z){
      logpp = sapply(1:p, function(k){
        logy_list[[k]](y=closest_Y_mat[i,k],lf=lf_mat[i,k]+z[k],par2=phi[k])
      }) # vector of length p
      out = sum(logpp, na.rm=T)
    }) # vector of length n^p
    logfw = logfyeps + log(mgh$weights) - logSumExp_C(logfyeps + log(mgh$weights))      
    epsilon_est_mat[i,] = colSums(exp(logfw) * mgh$points)
  }
  save(epsilon_est_mat, file=file.path(maindir, paste0(casename,"_tp_",ttt,".RData")))
} # end of mmm loop

#### Calculate averaged distance(similarity) matrix between estimated latent effects ####
for (mmm in 1:length(tps)){
  ttt = tps[mmm]
  load(file.path(maindir, paste0(casename,"_tp_",ttt,".RData")))
  omega_hat_mat = cov(epsilon_est_mat)
  MDist_temp = MDist_C(nrow(X_filter), epsilon_est_mat, solve(omega_hat_mat)) # calculate the Mahalanobis distance between estimated latent effects
  if (mmm == 1) {MDist_avg = MDist_temp} else {MDist_avg = MDist_avg + MDist_temp}
  rm(MDist_temp)
}
MDist_avg = sqrt(MDist_avg/length(tps))
MDist_avg[lower.tri(MDist_avg)] = t(MDist_avg)[lower.tri(t(MDist_avg))]

#### Create hierarchical clusters of subjects based on the distance matrix ####
hc = hclust(as.dist(MDist_avg), method="ward.D2") # hc$order is the row number of the distance matrix displayed on the dendrogram, from left to right
hcd = as.dendrogram(hclust(as.dist(MDist_avg), method="ward.D2"))
save(hc, hcd, file=file.path(wsdir,paste0(casename,"_dendrogram.RData")))

## Generate dendrogram
num = 4 # number of clusters
clus = cutree(hc, num) # group index of each subject
hcd_dend = hcd %>% set("branches_k_color", k = num)
dendpath = file.path(maindir, paste0(jobname,"_rectdend_",length(hc$order),
                                     "_subjects_",num,"_groups.pdf"))
pdf(file = dendpath, width=14, height=8.5)
hcd_dend %>% plot(main = paste("Dendrogram for",length(hc$order),"subjects"),
                  leaflab = "none")
dev.off()

## Check questionable group
closest_Obs_mat_filter[which(clus==4),] # subgroup 4 only has one subject, ID = 1814
lapply(Obs_allbm, function(df){
  df %>% filter(i==1814)
}) # this person has one outlier for health marker 1 at t = 6.744, so we excluded this person 
