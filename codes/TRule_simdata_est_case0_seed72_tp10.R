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
jobname =  "TRule_simdata_est"
rseed = 72; rsd = paste0("seed",rseed) # random seed used to generate one data set
cs = paste0("case", 0)
casename = paste(jobname,cs,rsd,sep="_")

maindir = "/pine/scr/j/i/jitong/thesis/TRule_model" # project directory
codesdir = file.path(maindir,"codes") # code directory
wsdir = file.path(maindir,"rdata") # R workspace directory

#### Load files ####
wsname = "TRule_simdata_GenerateData"
load(file.path(wsdir,paste0(wsname,"_",rsd,"_workspace.RData")))
load(file.path(wsdir,paste0(wsname,"_",rsd,"_intensity_est.RData"))) # gamma_est_mat(d1,d2), split_idx
source(file=file.path(codesdir,"functions_simdata.R"))
source(file=file.path(codesdir,"functions_estparams.R"))

#### Settings for parameter estimation ####
Obs_allbm = time_Xikj_all_df_list # measurements and time points
N0 = nrow(X) # number of subjects
m = ncol(X)-1 # number of baseline variables
p = length(Obs_allbm) # number of health markers
idx = 1:N0 # # index of subjects

## Gauss-Hermite quadrature nodes and weights
n1 = 20 # number of nodes of Gauss-Hermite quadrature for normal distribution
n2 = 10 # number of nodes of Gauss-Hermite quadrature for bivariate normal distribution
gh = gauss.hermite(n1) # univariate Gauss-Hermite quadrature nodes and weights
bigh = bigauss.hermite(n2) # bivariate Gauss-Hermite quadrature nodes and weights

## Specify kernel functions
k_epan = function(lambda, x, x0){
  t = abs(x-x0)/lambda
  return(ifelse(t<=1, 3/4*(1-t^2)/lambda, 0))
}
kernel.list = list(k_epan, k_epan)

## Specify bandwidths in kernel functions
bandwidth.vec = c(bw_beta=0.3, bw_var=0.2, bw_cov=0.2)

## Estimated gamma and eta parameters for the whole data
gamma.hat.all = gamma_est_mat[1:m,] # m*p matrix, each column for each temporal process
eta.hat.all = gamma_est_mat[(m+1):nrow(gamma_est_mat),]

## Interested time point
mmm = 11
ttt = t.grids[mmm]

#### Beta and variance estimators ####
## Beta initial values
beta.init = beta_p[1:p] %>%
  lapply(function(f){f(ttt)}) %>%
  unlist() %>%
  matrix(ncol=p)
beta_lower = rep(-10, m+1); beta_upper = rep(10, m+1)

## Var initial values
sigma.init = diag(Omega(ttt))
sigma_lower = 0.1; sigma_upper = 1.5

## Parameters for convergence of algorithm
factr = 1e1; abstol = 1e-3; maxit_out = 100; maxit_in = 100

## Parameter estimation
par_est_bv = vector("list",p) # result list, length p, each element for each health marker
for (kkk in 1:p){
  fname = paste(casename,"tp",ttt,"marker",kkk,sep="_")
  tempfile = file.path(maindir, paste0("reportOptim_",fname)) # temporary outputs
  sink(tempfile)
  par.est.bv = my_iter_k_v2(kernel.func.list = kernel.list,
                             an.vec = bandwidth.vec,
                             t = ttt, k = kkk,
                             Z = X, N_sub = nrow(X),
                             gh = gh, bigh = bigh,
                             time_Xikj_k = Obs_allbm[[kkk]],
                             gfun.inv.k = g_inv_p[[kkk]],
                             gfun.inv.deriv.k = g_inv_deriv_p[[kkk]],
                             gamma.hat.k = gamma.hat.all[,kkk],
                             eta.hat.k = eta.hat.all[,kkk],
                             beta_init = beta.init[,kkk], beta_lower = beta_lower, beta_upper = beta_upper,
                             sigma_init = sigma.init[kkk], sigma_lower = sigma_lower, sigma_upper = sigma_upper,
                             factr = factr, abstol = abstol,
                             maxit_out = maxit_out, maxit_in = maxit_in, verbose=TRUE)
  sink()
  all_iterOptim = readLines(tempfile)
  unlink(tempfile)
  write(all_iterOptim, file.path(maindir, paste0(fname,"_alliters.txt"))) # algorithm outputs
  
  par_est_bv[[kkk]] = par.est.bv
}

beta_est = sapply(par_est_bv, function(x){x$par[1:(m+1)]}) # beta estimator matrix, dim (m+1)*p, each column for each health marker
var_est = sapply(par_est_bv, function(x){x$par[m+2]}) # variance estimator vector, length p, each element for each health marker

#### Covariance estimators ####
## Cov initial values
cov.init = Omega(ttt)

## Parameters for convergence of algorithm
factr = 1e1; maxit_in = 500

## Parameter estimation
par_est_cov = vector("list",p*(p-1)/2) # result list, length p*(p-1)/2, each element for each pair of health markers
an.kl = bandwidth.vec[3] # bandwidth for estimating covariance parameters

for (kkk in 1:(p-1)){
  for (lll in ((kkk+1):p)){
    idx_vec = kkk*p-(p-lll)-(1+kkk)*kkk/2  
    ## initial value
    param.uplmt = sqrt(var_est[kkk]*var_est[lll])-1e-3
    param.lowlmt = -param.uplmt
    par_iter = param.uplmt * cov.init[kkk,lll]/sqrt(cov.init[kkk,kkk]*cov.init[lll,lll])

    ## Only include subjects who have measurements in the time interval 
    time_Yikj_temp = Obs_allbm[[kkk]] %>%
      filter(time > ttt-an.kl & time < ttt+an.kl)
    time_Yilj_temp = Obs_allbm[[lll]] %>%
      filter(time > ttt-an.kl & time < ttt+an.kl)
    
    time_Yikj = time_Yikj_temp %>%
      semi_join(time_Yilj_temp, by="i")
    time_Yilj = time_Yilj_temp %>%
      semi_join(time_Yikj_temp, by="i")

    ## estimate parameters
    fname = paste(casename,"tp",ttt,"cov",kkk,lll,sep="_")
    tempfile = file.path(maindir, paste0("reportOptim_",fname)) # temporary outputs
    sink(tempfile)
    par.est.cov = optim(par=par_iter, fn=obj.func.kl.v2, gr=grad.func.kl.v2, hessian=TRUE,
                        k=kkk, l=lll,
                        kernel.func.kl=kernel.list[[2]], an.kl=an.kl,
                        t=ttt, Z=X, N_sub = nrow(X),
                        t_X_1=time_Yikj, t_X_2=time_Yilj,
                        beta.k=beta_est[,kkk],
                        gfun.inv.k=g_inv_p[[kkk]],
                        gfun.inv.k.deriv=g_inv_deriv_p[[kkk]],
                        sigma.sq.k=var_est[kkk],
                        gamma.hat.k = gamma.hat.all[,kkk],
                        eta.hat.k = eta.hat.all[,kkk],
                        var.deriv.k=FALSE,
                        beta.l=beta_est[,lll],
                        gfun.inv.l=g_inv_p[[lll]],
                        gfun.inv.l.deriv=g_inv_deriv_p[[lll]],
                        sigma.sq.l=var_est[lll],
                        gamma.hat.l=gamma.hat.all[,lll],
                        eta.hat.l = eta.hat.all[,lll],
                        var.deriv.l=FALSE,
                        gh=gh, bigh=bigh,
                        lower=param.lowlmt, upper=param.uplmt,
                        method="L-BFGS-B", control=list(trace=6, maxit=maxit_in, factr=factr))
    sink()
    all_iterOptim = readLines(tempfile)
    unlink(tempfile)
    write(all_iterOptim, file.path(maindir,paste0(fname,"_alliters.txt"))) # algorithm outputs
    
    par_est_cov[[idx_vec]] = par.est.cov
  }}

cov_est = sapply(par_est_cov, function(x){x$par}) # covariance estimator vector, length p*(p-1)/2

#### Save all estimators ####
par_est_all = list(bv=par_est_bv, cov=par_est_cov)
save(par_est_all, file=file.path(wsdir, paste0(casename,"_tp",ttt,".RData")))

