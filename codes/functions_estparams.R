#### Extract measurements and latent Gaussian process near a time point ####
closest_Obs_eps <- function(obsdf, ttt, kkk, bvdf=NULL){
  out = obsdf %>% 
    filter(k==kkk) %>% 
    group_by(i) %>% 
    summarize(value_closest=value[which.min(abs(time-ttt))],
              eps_closest=eps[which.min(abs(time-ttt))]) %>%
    as.data.frame()
  if(!is.null(bvdf)){
    out = out %>%
      left_join(bvdf, by=c("i"="ID"))
  }
  return(out)
}

closest_Obs_allbm = function(obsdf, ttt){
  out <- obsdf %>%
    group_by_at(1) %>% # group by first column
    summarize(value_at_closest_t=value[which.min(abs(time-ttt))])
  return(out)
}

## Transform data and fit Andersen-Gill model ####
recent_summary = function(df, t_win){
  nn = nrow(df)
  L1 = rep(0, nn)
  L2 = rep(0, nn)
  for (ii in 1:nn){
    y_pre = df %>%
      filter(i == df$i[ii] & time < df$time[ii] & time >= df$time[ii]-t_win)
    if (nrow(y_pre) == 0){next} else {
      L1[ii] = 1
      L2[ii] = mean(y_pre$value)
    }
  }
  return(cbind(df %>% select(-eps, -mu),L1,L2))
}

my_SurvData_v2 = function(y=NULL, t=NULL, t_win=0, t_max){
  # y=obs_p_all[[k]][[i]][[2]];t=obs_p_all[[k]][[i]][[1]];t_win=L_pre_ip[k]
  n_y = length(y[t<t_max]); n_t = length(y[t<t_max]);
  if (n_y != n_t){
    return(print("Y and T don't match."))
  } else if (n_t == 0){
    return(data.frame(y=NA, event=0, start=0, stop=t_max, L1=0, L2=0))
  } else {
    y = y[t<t_max]
    t = t[t<t_max]
    
    df = data.frame(y=c(y,NA), 
                    event=c(rep(1,n_y),0), 
                    start=c(0,t),
                    t_ext=c(0,t)+t_win,
                    stop=c(t,t_max))
    
    output =  data.frame(y = df$y[1],
                         event = 1,
                         start = df$start[1],
                         stop = df$stop[1],
                         L1 = 0, L2 = 0)
    
    for (j in 2:nrow(df)){
      t_idx = which(df$t_ext > df$start[j] & df$t_ext < df$stop[j])
      n_idx = length(t_idx)
      if (df$start[j] < df$stop[j] - t_win){
        out_tmp = data.frame(y = c(rep(NA, n_idx), df$y[j]),
                             event = c(rep(0, n_idx), 1),
                             start = c(df$start[j], df$t_ext[t_idx]),
                             stop = c(df$t_ext[t_idx], df$stop[j]),
                             L1 = c(rep(1, n_idx), 0), 
                             L2 = c(rep(1, n_idx), 0))
        for (jj in 1:n_idx){
          out_tmp$L2[jj] = mean(df$y[t_idx[jj:n_idx]-1], na.rm=TRUE)
        }
      } else if (n_idx == 0){
        out_tmp =  data.frame(y = df$y[j],
                              event = 1,
                              start = df$start[j],
                              stop = df$stop[j],
                              L1 = 1, 
                              L2 = mean(df$y[which(df$start >= df$stop[j] - t_win & df$start < df$stop[j])-1],
                                        na.rm=TRUE))
      } else {
        out_tmp = data.frame(y = c(rep(NA, n_idx), df$y[j]),
                             event = c(rep(0, n_idx), 1),
                             start = c(df$start[j], df$t_ext[t_idx]),
                             stop = c(df$t_ext[t_idx], df$stop[j]),
                             L1 = rep(1, n_idx+1), 
                             L2 = rep(NA, n_idx+1))
        for (jj in 1:n_idx){
          # same number in df$start < df$t_ext if use subtraction on RHS
          # out_tmp$L2[jj] = mean(df$y[which(df$start >= df$t_ext[t_idx[jj]] - t_win & df$start < df$t_ext[t_idx[jj]])-1],
          #                       na.rm=TRUE)
          out_tmp$L2[jj] = mean(df$y[which(df$stop + t_win >= df$t_ext[t_idx[jj]] & df$stop < df$t_ext[t_idx[jj]])],
                                na.rm=TRUE)
        }
        out_tmp$L2[n_idx+1] = mean(df$y[which(df$t_ext[1:j]>=df$stop[j]) - 1],
                                   na.rm=TRUE)
      }
      output = rbind(output, out_tmp)
    } # end of for (j in 2:nrow(df))
    output$event[output$stop == t_max] = 0
    return(output)
  }
}

#### Gauss-Hermite quadrature ####
hermite <- function (points, z) {
  p1 <- 1/pi^0.4
  p2 <- 0
  for (j in 1:points) {
    p3 <- p2
    p2 <- p1
    p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) * p3
  }
  pp <- sqrt(2 * points) * p2
  c(p1, pp)
}

gauss.hermite <- function (points, iterlim = 50) {
  x <- w <- rep(0, points)
  m <- (points + 1)/2
  for (i in 1:m) {
    z <- if (i == 1) 
      sqrt(2 * points + 1) - 2 * (2 * points + 1)^(-1/6)
    else if (i == 2) 
      z - sqrt(points)/z
    else if (i == 3 || i == 4) 
      1.9 * z - 0.9 * x[i - 2]
    else 2 * z - x[i - 2]
    for (j in 1:iterlim) {
      z1 <- z
      p <- hermite(points, z)
      z <- z1 - p[1]/p[2]
      if (abs(z - z1) <= 1e-15) 
        break
    }
    if (j == iterlim) 
      warning("iteration limit exceeded")
    x[points + 1 - i] <- -(x[i] <- z)
    w[i] <- w[points + 1 - i] <- 2/p[2]^2
  }
  # r <- cbind(x * sqrt(2), w/sum(w))
  # colnames(r) <- c("Points", "Weights")
  r <- data.frame(points=x * sqrt(2), weights=w/sum(w))
  r
}

bigauss.hermite <- function(n, mu=c(0,0), sigma=diag(2), prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  
  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  return(list(points=pts, weights=wts))
  
  # ## R matrix
  # a <- sigma[1,1]
  # b <- sigma[1,2]
  # d <- sigma[2,2]
  # tau <- a+d
  # delta <- sqrt(a*d-b^2)
  # R.kl.11 <- (a+delta)/sqrt(tau+2*delta)
  # R.kl.12 <- R.kl.21 <- b/sqrt(tau+2*delta)
  # R.kl.22 <- (d+delta)/sqrt(tau+2*delta)
  # 
  # 
  # y.pts <- matrix(c(mu[1]+R.kl.11*pts[,1] + R.kl.21*pts[,2], mu[2]+R.kl.12*pts[,1] + R.kl.22*pts[,2]),
  #                 ncol=dm, byrow=F)
  # 
  # return(list(points=y.pts, weights=wts))
}

mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))
    stop("mu and sigma have nonconformable dimensions")
  
  dm  <- length(mu)
  gh  <- gauss.hermite(n)
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  
  ## rotate, scale, translate points
  eig <- eigen(sigma)
  if(dm>1) {rot <- eig$vectors %*% diag(sqrt(eig$values))
  } else {rot <- eig$vectors %*% sqrt(eig$values)
  }
  # rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}

#### Square root of a covariance maxtrix ####
mat.sqroot.2dsym <- function(a,d,b){
  tau <- a+d
  delta <- sqrt(a*d-b^2)
  R.11 <- (a+delta)/sqrt(tau+2*delta)
  R.12 <- R.21 <- b/sqrt(tau+2*delta)
  R.22 <- (d+delta)/sqrt(tau+2*delta)
  
  return(matrix(c(R.11,R.12,R.21,R.22),nrow=2,byrow=T))
}

#### Partial R.kl partial sigma.sq.k, sigma.sq.l, sigma.kl ####
partial.R.kl <- function(a,d,b){
  tau <- a+d
  delta <- sqrt(a*d-b^2)

  p.sigma.sq.k.11 <- ((1+1/2*d/delta)*sqrt(tau+2*delta) - (a+delta)/2/sqrt(tau+2*delta)*(1+d/delta))/(tau+2*delta)
  p.sigma.sq.k.21 <- p.sigma.sq.k.12 <- b*(-1/2)*((tau+2*delta)^(-3/2))*(1+d/delta)
  p.sigma.sq.k.22 <- ((1/2*d/delta)*sqrt(tau+2*delta) - (d+delta)/2/sqrt(tau+2*delta)*(1+d/delta))/(tau+2*delta)

  p.sigma.sq.l.11 <- ((1/2*a/delta)*sqrt(tau+2*delta) - (a+delta)/2/sqrt(tau+2*delta)*(1+a/delta))/(tau+2*delta)
  p.sigma.sq.l.21 <- p.sigma.sq.l.12 <- b*(-1/2)*((tau+2*delta)^(-3/2))*(1+a/delta)
  p.sigma.sq.l.22 <- ((1+1/2*a/delta)*sqrt(tau+2*delta) - (d+delta)/2/sqrt(tau+2*delta)*(1+a/delta))/(tau+2*delta)

  p.sigma.kl.11 <- ((-b/delta)*sqrt(tau+2*delta) - (a+delta)/sqrt(tau+2*delta)*(-b/delta))/(tau+2*delta)
  p.sigma.kl.21 <- p.sigma.kl.12 <- (sqrt(tau+2*delta)-b/sqrt(tau+2*delta)*(-b/delta))/(tau+2*delta)
  p.sigma.kl.22 <- ((-b/delta)*sqrt(tau+2*delta) - (d+delta)/sqrt(tau+2*delta)*(-b/delta))/(tau+2*delta)

  return(list(p.sigma.sq.k=matrix(c(p.sigma.sq.k.11,p.sigma.sq.k.12,p.sigma.sq.k.21,p.sigma.sq.k.22),nrow=2,byrow=T),
              p.sigma.sq.l=matrix(c(p.sigma.sq.l.11,p.sigma.sq.l.12,p.sigma.sq.l.21,p.sigma.sq.l.22),nrow=2,byrow=T),
              p.sigma.kl=matrix(c(p.sigma.kl.11,p.sigma.kl.12,p.sigma.kl.21,p.sigma.kl.22),nrow=2,byrow=T)
  ))
}


#### E.ik, l.k and their partial derivatives for beta.k, sigma.sq.k ####
E.ik.and.derivs <- function(Z.i, beta.k, gfun.inv.k, gfun.inv.deriv.k,
                            n, mu.k=0, sigma.sq.k, prune.fraction=NULL, var.deriv=FALSE, gh){
  u.k <- c(t(Z.i) %*% beta.k) + gh$points*sqrt(sigma.sq.k) + mu.k
  
  if(var.deriv){
    return(list(val = sum(gfun.inv.k(u.k) * gh$weights)
				, beta.k.deriv = sum(gfun.inv.deriv.k(u.k) * gh$weights) * t(Z.i)
				, sigma.sq.k.deriv = sum(gfun.inv.deriv.k(u.k) * gh$points/(2*sqrt(sigma.sq.k)) * gh$weights)
    ))
  } else {
    return(list(val = sum(gfun.inv.k(u.k) * gh$weights)
                , beta.k.deriv = sum(gfun.inv.deriv.k(u.k) * gh$weights) * t(Z.i)
    ))
  }
}

# use L.k as a function of t.ikj and X.ik_t.ikj, or as an input?
l.k.and.derivs <- function(k, K.an=kernel.func, an, t, Z, t_X_k, 
                           beta.k, gfun.inv.k, gfun.inv.deriv.k,
                           n, mu.k=0, sigma.sq.k, prune.fraction=NULL,
                           gamma.hat.k, eta.hat.k=NULL, L.k=NULL, var.deriv=FALSE, gh){

	# # beta.k: (m+1)*1 vector
	# # gamma.hat.k: m*1 vector
	# colnames(t_X_k) are i,k,j,time,value,L_ik^T  

  if (var.deriv){
    l.k <- matrix(0, ncol(Z))
    l.k.beta.k.deriv <- matrix(0, ncol(Z), ncol(Z))
    l.k.sigma.sq.k.deriv <- matrix(0, ncol(Z))
        
    # Z: N*(m+1) matrix, the row vector is (1, age, sex) for the ith patient 
    for (ii in unique(t_X_k[,1])){
      # Z.i: (m+1)*1 vector
      Z.i <- as.matrix(Z[ii,])
      E.ik.results <- E.ik.and.derivs(Z.i, beta.k, gfun.inv.k, gfun.inv.deriv.k, 
                                      n, mu.k, sigma.sq.k, prune.fraction, var.deriv, gh)
      t_X_ik = t_X_k %>% filter(i==ii) %>% as.matrix # not all subjects in Z_i have X_ik, matrix gurantee L is a vector  
	  
  	  for (j in 1:nrow(t_X_ik)){
          t.ikj <- t_X_ik[j,4] # time
          X.ik_t.ikj <- t_X_ik[j,5] # value
  		    L.ik_t.ikj = t_X_ik[j,6:ncol(t_X_k)]
          dN.ik_t.ikj <- exp(-c(t(Z.i[-1,]) %*% gamma.hat.k + t(L.ik_t.ikj) %*% eta.hat.k))
          
          # return a (m+1)*1 matrix
          l.k <- l.k + K.an(an, t.ikj, t) * Z.i * (X.ik_t.ikj - E.ik.results$val) * dN.ik_t.ikj 

          # return a (m+1)*(m+1) matrix
          l.k.beta.k.deriv <- l.k.beta.k.deriv + K.an(an, t.ikj, t) * Z.i %*% (- E.ik.results$beta.k.deriv) * dN.ik_t.ikj 
          
          # return a (m+1)*1 matrix
          l.k.sigma.sq.k.deriv <- l.k.sigma.sq.k.deriv + K.an(an, t.ikj, t) * Z.i * (- E.ik.results$sigma.sq.k.deriv) * dN.ik_t.ikj
      } # for (j in 1:nrow(t_X_ik))
    } # for (i in unique(t_X_k[,1]))
    return(list(val=l.k, beta.k.deriv=l.k.beta.k.deriv
                , sigma.sq.k.deriv=l.k.sigma.sq.k.deriv
    ))
  } else {
    l.k <- matrix(0, ncol(Z)) # a (m+1)*1 matrix
    l.k.beta.k.deriv <- matrix(0, ncol(Z), ncol(Z))
        
    # Z: N*(m+1) matrix, the row vector is (1, age, sex) for the ith patient 
    for (ii in unique(t_X_k[,1])){
      # Z.i: (m+1)*1 vector
      Z.i <- as.matrix(Z[ii,])
      E.ik.results <- E.ik.and.derivs(Z.i, beta.k, gfun.inv.k, gfun.inv.deriv.k, 
                                      n, mu.k, sigma.sq.k, prune.fraction, var.deriv, gh)
      t_X_ik = t_X_k %>% filter(i==ii) %>% as.matrix
      
  	  for (j in 1:nrow(t_X_ik)){
          t.ikj <- t_X_ik[j,4] # time
          X.ik_t.ikj <- t_X_ik[j,5] # value
      		L.ik_t.ikj = t_X_ik[j,6:ncol(t_X_k)]
          dN.ik_t.ikj <- exp(-c(t(Z.i[-1,]) %*% gamma.hat.k + t(L.ik_t.ikj) %*% eta.hat.k))
          
          # return a (m+1)*1 matrix
          l.k <- l.k + K.an(an, t.ikj, t) * Z.i * (X.ik_t.ikj - E.ik.results$val) * dN.ik_t.ikj 
          
          # return a (m+1)*(m+1) matrix
          l.k.beta.k.deriv <- l.k.beta.k.deriv + K.an(an, t.ikj, t) * Z.i %*% (- E.ik.results$beta.k.deriv) * dN.ik_t.ikj 
      }
    }
    return(list(val=l.k, beta.k.deriv=l.k.beta.k.deriv))
  } # end of else (fixed variances)
}

## E.ikk, l.kk and their partial derivatives for beta.k, sigma.sq.k ####

E.ikk.and.derivs <- function(Z.i, beta.k, gfun.inv.k, gfun.inv.deriv.k,
                             n, mu.k=0, sigma.sq.k, prune.fraction=NULL, var.deriv=FALSE, gh){
  if (var.deriv){
    u.k <- c(t(Z.i) %*% beta.k) + gh$points*sqrt(sigma.sq.k) + mu.k
    return(list(val = sum((gfun.inv.k(u.k))^2 * gh$weights)
                , beta.k.deriv = sum(2*gfun.inv.k(u.k) * gfun.inv.deriv.k(u.k) * gh$weights) * t(Z.i)
                , sigma.sq.k.deriv = sum(2*gfun.inv.k(u.k) * gfun.inv.deriv.k(u.k) * gh$points/(2*sqrt(sigma.sq.k)) * gh$weights)
    ))
  } else {return(NULL)}
}

l.kk.and.derivs <- function(k, K.an=kernel.func, an, t, Z, t_X_kk, 
                            beta.k, gfun.inv.k, gfun.inv.deriv.k,
                            n, mu.k=0, sigma.sq.k, prune.fraction=NULL,
                            gamma.hat.k, eta.hat.k=NULL, L.k=NULL, var.deriv=FALSE, gh){
  # Z: N*(m+1) matrix, the row vector is (1, age, sex) for the ith patient 
  
  if (var.deriv){
    l.kk <- 0
    l.kk.beta.k.deriv <- matrix(0, 1, ncol(Z))
    l.kk.sigma.sq.k.deriv <- 0
    
    # # beta.k: (m+1)*1 vector
    # # gamma.hat.k: m*1 vector

    for (ii in unique(t_X_kk[,1])){
      # Z.i: (m+1)*1 vector
      Z.i <- as.matrix(Z[ii,])
      E.ikk.results <- E.ikk.and.derivs(Z.i, beta.k, gfun.inv.k, gfun.inv.deriv.k,
                                        n, mu.k, sigma.sq.k, prune.fraction, var.deriv, gh)
      t_X_ik = t_X_kk %>% filter(i==ii) %>% as.matrix
      
  	  for (j in 1:nrow(t_X_ik)){
        t.ikj <- t_X_ik[j,4] # time
        X.ik_t.ikj <- t_X_ik[j,5] # value
    		L.ik_t.ikj = t_X_ik[j,6:ncol(t_X_kk)]
        dN.ik_t.ikj <- exp(-c(t(Z.i[-1,]) %*% gamma.hat.k + t(L.ik_t.ikj) %*% eta.hat.k))

    		for (jj in 1:nrow(t_X_ik)){
    		  if (jj != j){
    		  # if (jj == j){
    		    t.ikjj <- t_X_ik[jj,4] # time
      			X.ik_t.ikjj <- t_X_ik[jj,5] # value
      			L.ik_t.ikjj = t_X_ik[jj,6:ncol(t_X_kk)]
      			dN.ik_t.ikjj <- exp(-c(t(Z.i[-1,]) %*% gamma.hat.k + t(L.ik_t.ikjj) %*% eta.hat.k))

      			# return a scalar
      			l.kk <- l.kk + K.an(an, t.ikj, t) * K.an(an, t.ikjj, t)  * (X.ik_t.ikj * X.ik_t.ikjj - E.ikk.results$val) * dN.ik_t.ikj * dN.ik_t.ikjj

      			# return a 1*(m+1) matrix
      			l.kk.beta.k.deriv <- l.kk.beta.k.deriv + K.an(an, t.ikj, t) * K.an(an, t.ikjj, t)  * (- E.ikk.results$beta.k.deriv) * dN.ik_t.ikj * dN.ik_t.ikjj

      			# return a scalar
      			l.kk.sigma.sq.k.deriv <- l.kk.sigma.sq.k.deriv + K.an(an, t.ikj, t) * K.an(an, t.ikjj, t)  * (- E.ikk.results$sigma.sq.k.deriv) * dN.ik_t.ikj * dN.ik_t.ikjj
  		    } # if (jj != j)
    		}	# for (jj in 1:nrow(t_X_ik))	  
      } # for (j in 1:nrow(t_X_ik))
    } # for (i in unique(t_X_kk[,1]))
    return(list(val=l.kk, beta.k.deriv=l.kk.beta.k.deriv
                , sigma.sq.k.deriv=l.kk.sigma.sq.k.deriv
    ))
  } else {return(NULL)}
}

#### E.ikl, l.kl and their partial derivatives for beta.k, beta.l, sigma.sq.k, sigma.sq.l, sigma.kl ####
E.ikl.and.derivs <- function(Z.i, beta.1, gfun.inv.1, gfun.inv.deriv.1, 
                             beta.2, gfun.inv.2, gfun.inv.deriv.2,
                             n=10, mu.12=c(0,0), sigma.sq.1, sigma.sq.2, sigma.12, prune.fraction=NULL,
                             var.deriv.1=FALSE, var.deriv.2=FALSE, bigh){

  R.kl <- mat.sqroot.2dsym(sigma.sq.1, sigma.sq.2, sigma.12) # a 2x2 matrix
  u.k <- c(t(Z.i) %*% beta.1) + mu.12[1] + R.kl[1,1]*bigh$points[,1] + R.kl[2,1]*bigh$points[,2]
  u.l <- c(t(Z.i) %*% beta.2) + mu.12[2] + R.kl[1,2]*bigh$points[,1] + R.kl[2,2]*bigh$points[,2]
  
  p.R.kl <- partial.R.kl(sigma.sq.1, sigma.sq.2, sigma.12)
  p.R.kl.sigma.sq.k <- p.R.kl$p.sigma.sq.k
  p.R.kl.sigma.sq.l <- p.R.kl$p.sigma.sq.l
  
  p.R.kl.sigma.kl <- p.R.kl$p.sigma.kl
  partial.u.k.sigma.kl <- p.R.kl.sigma.kl[1,1]*bigh$points[,1] + p.R.kl.sigma.kl[2,1]*bigh$points[,2]
  partial.u.l.sigma.kl <- p.R.kl.sigma.kl[1,2]*bigh$points[,1] + p.R.kl.sigma.kl[2,2]*bigh$points[,2]
  
  rrr <- list(val = sum(gfun.inv.1(u.k) * gfun.inv.2(u.l) * bigh$weights),
              
              beta.k.deriv = sum(gfun.inv.deriv.1(u.k) * gfun.inv.2(u.l) * bigh$weights) * t(Z.i),
              beta.l.deriv = sum(gfun.inv.1(u.k) * gfun.inv.deriv.2(u.l) * bigh$weights) * t(Z.i),
              
              sigma.kl.deriv = sum((gfun.inv.deriv.1(u.k) * partial.u.k.sigma.kl * gfun.inv.2(u.l) 
                                    + gfun.inv.1(u.k) * gfun.inv.deriv.2(u.l) * partial.u.l.sigma.kl) * bigh$weights)
  )
  
  if (var.deriv.1){
    partial.u.k.sigma.sq.k <- p.R.kl.sigma.sq.k[1,1]*bigh$points[,1] + p.R.kl.sigma.sq.k[2,1]*bigh$points[,2]
    partial.u.l.sigma.sq.k <- p.R.kl.sigma.sq.k[1,2]*bigh$points[,1] + p.R.kl.sigma.sq.k[2,2]*bigh$points[,2]
    rrr <- list.append(rrr, sigma.sq.k.deriv = sum((gfun.inv.deriv.1(u.k) * partial.u.k.sigma.sq.k * gfun.inv.2(u.l)
                                                    + gfun.inv.1(u.k) * gfun.inv.deriv.2(u.l) * partial.u.l.sigma.sq.k) * bigh$weights)
    )
  }
  
  if(var.deriv.2){
    partial.u.k.sigma.sq.l <- p.R.kl.sigma.sq.l[1,1]*bigh$points[,1] + p.R.kl.sigma.sq.l[2,1]*bigh$points[,2]
    partial.u.l.sigma.sq.l <- p.R.kl.sigma.sq.l[1,2]*bigh$points[,1] + p.R.kl.sigma.sq.l[2,2]*bigh$points[,2]
    rrr <- list.append(rrr, sigma.sq.l.deriv = sum((gfun.inv.deriv.1(u.k) * partial.u.k.sigma.sq.l * gfun.inv.2(u.l)
                                                    + gfun.inv.1(u.k) * gfun.inv.deriv.2(u.l) * partial.u.l.sigma.sq.l) * bigh$weights)
    )
  }
  
  return(rrr)
}

l.kl.and.derivs.v2 <- function(k1, k2, K.an=kernel.func, an, t, Z, 
                               t_X_1, t_X_2,
                            beta.1, gfun.inv.1, gfun.inv.deriv.1, beta.2, gfun.inv.2, gfun.inv.deriv.2,
                            n=10, mu.12=c(0,0), sigma.sq.1, sigma.sq.2, sigma.12, prune.fraction=NULL,
                            gamma.hat.1, gamma.hat.2, eta.hat.1, eta.hat.2,
                            var.deriv.1=FALSE, var.deriv.2=FALSE, bigh){
  l.kl <- 0
  l.kl.beta.k.deriv <- matrix(0, 1, ncol(Z))
  l.kl.beta.l.deriv <- matrix(0, 1, ncol(Z))
  l.kl.sigma.kl.deriv <- 0
  
  if (var.deriv.1 & var.deriv.2){
    l.kl.sigma.sq.k.deriv <- 0
    l.kl.sigma.sq.l.deriv <- 0
    
    # Z: N*(m+1) matrix, the row vector is (1, age, sex) for the ith patient 
    for (ii in unique(t_X_1[,1])){
      Z.i <- as.matrix(Z[ii,])
      E.ikl.results <- E.ikl.and.derivs(Z.i, beta.1, gfun.inv.1, gfun.inv.deriv.1, 
                                        beta.2, gfun.inv.2, gfun.inv.deriv.2,
                                        n, mu.12, sigma.sq.1, sigma.sq.2, sigma.12, prune.fraction=NULL,
                                        var.deriv.1, var.deriv.2, bigh)
      t_X_ik = t_X_1 %>% filter(i==ii) %>% as.matrix
      t_X_il = t_X_2 %>% filter(i==ii) %>% as.matrix
      
      # t.k: a list containing all measured timepoints of the kth temporal process for the all patients 
      for (j in 1:nrow(t_X_ik)){
        t.ikj <- t_X_ik[j,4] # time
        X.ik_t.ikj <- t_X_ik[j,5] # value
        L.ik_t.ikj = t_X_ik[j,6:ncol(t_X_1)]
        dN.ik_t.ikj <- exp(-c(t(Z.i[-1,]) %*% gamma.hat.1 + t(L.ik_t.ikj) %*% eta.hat.1))
        
        for (jj in 1:nrow(t_X_il)){
          t.iljj <- t_X_il[jj,4]
          X.il_t.iljj <- t_X_il[jj,5]
          L.il_t.iljj = t_X_il[jj,6:ncol(t_X_2)]
          dN.il_t.iljj <- exp(-c(t(Z.i[-1,]) %*% gamma.hat.2 + t(L.il_t.iljj) %*% eta.hat.2))
          
          # return a scalar
          l.kl <- l.kl + K.an(an, t.ikj, t) * K.an(an, t.iljj, t)  * (X.ik_t.ikj * X.il_t.iljj - E.ikl.results$val) * dN.ik_t.ikj * dN.il_t.iljj
          
          # return a 1*(m+1) vector
          l.kl.beta.k.deriv <- l.kl.beta.k.deriv + K.an(an, t.ikj, t) * K.an(an, t.iljj, t)  * (- E.ikl.results$beta.k.deriv) * dN.ik_t.ikj * dN.il_t.iljj
          
          # return a 1*(m+1) vector
          l.kl.beta.l.deriv <- l.kl.beta.l.deriv + K.an(an, t.ikj, t) * K.an(an, t.iljj, t)  * (- E.ikl.results$beta.l.deriv) * dN.ik_t.ikj * dN.il_t.iljj
          
          # return a scalar
          l.kl.sigma.sq.k.deriv <- l.kl.sigma.sq.k.deriv + K.an(an, t.ikj, t) * K.an(an, t.iljj, t)  * (- E.ikl.results$sigma.sq.k.deriv) * dN.ik_t.ikj * dN.il_t.iljj
          
          # return a scalar
          l.kl.sigma.sq.l.deriv <- l.kl.sigma.sq.l.deriv + K.an(an, t.ikj, t) * K.an(an, t.iljj, t)  * (- E.ikl.results$sigma.sq.l.deriv) * dN.ik_t.ikj * dN.il_t.iljj
          
          # return a scalar
          l.kl.sigma.kl.deriv <- l.kl.sigma.kl.deriv + K.an(an, t.ikj, t) * K.an(an, t.iljj, t)  * (- E.ikl.results$sigma.kl.deriv) * dN.ik_t.ikj * dN.il_t.iljj
        } # for (jj in 1:nrow(t_X_il))
      } # for (j in 1:nrow(t_X_ik))
    } # for (i in 1:nrow(Z_12))
    return(list(val=l.kl, 
                beta.k.deriv=l.kl.beta.k.deriv, beta.l.deriv=l.kl.beta.l.deriv,
                sigma.sq.k.deriv=l.kl.sigma.sq.k.deriv, 
                sigma.sq.l.deriv=l.kl.sigma.sq.l.deriv,
                sigma.kl.deriv=l.kl.sigma.kl.deriv
    ))
  } else {
    # Z: N*(m+1) matrix, the row vector is (1, age, sex) for the ith patient 
    for (ii in unique(t_X_1[,1])){
      Z.i <- as.matrix(Z[ii,])
      E.ikl.results <- E.ikl.and.derivs(Z.i, beta.1, gfun.inv.1, gfun.inv.deriv.1, 
                                        beta.2, gfun.inv.2, gfun.inv.deriv.2,
                                        n, mu.12, sigma.sq.1, sigma.sq.2, sigma.12, prune.fraction=NULL,
                                        var.deriv.1, var.deriv.2, bigh)
      t_X_ik = t_X_1 %>% filter(i==ii) %>% as.matrix
      t_X_il = t_X_2 %>% filter(i==ii) %>% as.matrix
      
      # t.k: a list containing all measured timepoints of the kth temporal process for the all patients 
      for (j in 1:nrow(t_X_ik)){
        t.ikj <- t_X_ik[j,4] # time
        X.ik_t.ikj <- t_X_ik[j,5] # value
        L.ik_t.ikj = t_X_ik[j,6:ncol(t_X_1)]
        dN.ik_t.ikj <- exp(-c(t(Z.i[-1,]) %*% gamma.hat.1 + t(L.ik_t.ikj) %*% eta.hat.1))
        
        for (jj in 1:nrow(t_X_il)){
          t.iljj <- t_X_il[jj,4]
          X.il_t.iljj <- t_X_il[jj,5]
          L.il_t.iljj = t_X_il[jj,6:ncol(t_X_2)]
          dN.il_t.iljj <- exp(-c(t(Z.i[-1,]) %*% gamma.hat.2 + t(L.il_t.iljj) %*% eta.hat.2))
          
          # return a scalar
          l.kl <- l.kl + K.an(an, t.ikj, t) * K.an(an, t.iljj, t)  * (X.ik_t.ikj * X.il_t.iljj - E.ikl.results$val) * dN.ik_t.ikj * dN.il_t.iljj
          
          # return a 1*(m+1) vector
          l.kl.beta.k.deriv <- l.kl.beta.k.deriv + K.an(an, t.ikj, t) * K.an(an, t.iljj, t)  * (- E.ikl.results$beta.k.deriv) * dN.ik_t.ikj * dN.il_t.iljj
          
          # return a 1*(m+1) vector
          l.kl.beta.l.deriv <- l.kl.beta.l.deriv + K.an(an, t.ikj, t) * K.an(an, t.iljj, t)  * (- E.ikl.results$beta.l.deriv) * dN.ik_t.ikj * dN.il_t.iljj
          
          # return a scalar
          l.kl.sigma.kl.deriv <- l.kl.sigma.kl.deriv + K.an(an, t.ikj, t) * K.an(an, t.iljj, t)  * (- E.ikl.results$sigma.kl.deriv) * dN.ik_t.ikj * dN.il_t.iljj
        }
      }
    }
    return(list(val=l.kl, 
                beta.k.deriv=l.kl.beta.k.deriv, beta.l.deriv=l.kl.beta.l.deriv,
                sigma.kl.deriv=l.kl.sigma.kl.deriv
    ))
  }
  
}


#### Obj.k; only k; beta, cov ####
obj.func.k.v3 <- function(parmt, kernel.func.k, kernel.func.kk, an.k, an.kk,
                          t, k, Z, N_sub, gh, bigh, t_X_k, t_X_kk,
                          gfun.inv.k, gfun.inv.deriv.k, 
                          gamma.hat.k, eta.hat.k, 
                          sigma.sq.fix.k, var.deriv.k){
  
  par.beta = parmt[1:(m+1)] # a (m+1) vector
  
  par.var = sigma.sq.fix.k
  if (var.deriv.k){
    par.var = parmt[m+2]
  }
  
  l_k <- l.k.and.derivs(k, kernel.func.k, an.k, t, Z, t_X_k,
                        beta.k=par.beta, gfun.inv.k, gfun.inv.deriv.k,
                        n, mu.k=0, sigma.sq.k=par.var, prune.fraction=NULL,
                        gamma.hat.k, eta.hat.k, L.k=NULL, var.deriv.k, gh)$val
  if(var.deriv.k){
    l_kk <- l.kk.and.derivs(k, kernel.func.kk, an.kk, t, Z, t_X_kk,
                            beta.k=par.beta, gfun.inv.k, gfun.inv.deriv.k,
                            n, mu.k=0, sigma.sq.k=par.var, prune.fraction=NULL,
                            gamma.hat.k, eta.hat.k, L.k=NULL, var.deriv.k, gh)$val
  } else {l_kk <- 0}
  
  W = (c(t(l_k) %*% l_k) + l_kk^2)/N_sub^2  # sum up equations for beta.k and variance sigma.

  return(W)
} # obj.func.k.v3

grad.func.k.v3 <- function(parmt, kernel.func.k, kernel.func.kk, an.k, an.kk,
                           t, k, Z, N_sub, gh, bigh, t_X_k, t_X_kk,
                           gfun.inv.k, gfun.inv.deriv.k, 
                           gamma.hat.k, eta.hat.k, 
                           sigma.sq.fix.k, var.deriv.k){
  par.beta = parmt[1:(m+1)] # a (m+1) vector
  
  par.var = sigma.sq.fix.k
  if (var.deriv.k){
    par.var = parmt[m+2]
  }
  
  l.k.results <- l.k.and.derivs(k, kernel.func.k, an.k, t, Z, t_X_k, 
                                beta.k=par.beta, gfun.inv.k, gfun.inv.deriv.k,
                                n, mu.k=0, sigma.sq.k=par.var, prune.fraction=NULL,
                                gamma.hat.k, eta.hat.k, L.k=NULL, var.deriv.k, gh)
  
  par.beta.grad = c(2*(t(l.k.results$val) %*% l.k.results$beta.k.deriv))
  
  if (var.deriv.k){
    par.var.grad = 2*c(t(l.k.results$val) %*% l.k.results$sigma.sq.k.deriv)
    
    l.kk.results <- l.kk.and.derivs(k, kernel.func.kk, an.kk, t, Z, t_X_kk,
                                    beta.k=par.beta, gfun.inv.k, gfun.inv.deriv.k,
                                    n, mu.k=0, sigma.sq.k=par.var, prune.fraction=NULL,
                                    gamma.hat.k, eta.hat.k, L.k=NULL, var.deriv.k, gh)
    
    par.beta.grad = c(2*l.kk.results$val * l.kk.results$beta.k.deriv)
    
    par.var.grad = 2*l.kk.results$val * l.kk.results$sigma.sq.k.deriv
  } else {par.var.grad = 0}
  
  par.grad = c(par.beta.grad, par.var.grad[var.deriv.k])/N_sub^2
  return(par.grad)
} # grad.func.k.v3

#### Obj kk; only k; cov ####
obj.func.kk.v2 <- function(parmt, kernel.func.kk, an.kk,
                          t, k, Z, N_sub, gh, bigh, t_X_kk,
                          gfun.inv.k, gfun.inv.deriv.k, 
                          gamma.hat.k, eta.hat.k, 
                          beta.fix.k, var.deriv.k=TRUE){

  l_kk <- l.kk.and.derivs(k, kernel.func.kk, an.kk, t, Z, t_X_kk,
                          beta.k=beta.fix.k, gfun.inv.k, gfun.inv.deriv.k,
                          n, mu.k=0, sigma.sq.k=parmt, prune.fraction=NULL,
                          gamma.hat.k, eta.hat.k, L.k=NULL, var.deriv.k, gh)$val
  
  W = l_kk^2
  return(W)
} # obj.func.kk.v2

grad.func.kk.v2 <- function(parmt, kernel.func.kk, an.kk,
                           t, k, Z, N_sub, gh, bigh, t_X_kk,
                           gfun.inv.k, gfun.inv.deriv.k, 
                           gamma.hat.k, eta.hat.k, 
                           beta.fix.k, var.deriv.k=TRUE){
  l.kk.results <- l.kk.and.derivs(k, kernel.func.kk, an.kk, t, Z, t_X_kk,
                                  beta.k=beta.fix.k, gfun.inv.k, gfun.inv.deriv.k,
                                  n, mu.k=0, sigma.sq.k=parmt, prune.fraction=NULL,
                                  gamma.hat.k, eta.hat.k, L.k=NULL, var.deriv.k, gh)
  
  par.var.grad = 2*l.kk.results$val * l.kk.results$sigma.sq.k.deriv
  return(par.var.grad)
} # grad.func.kk.v2

#### Iteration; only k; beta, cov ####
my_iter_k_v2 = function(kernel.func.list, an.vec,
                     t, k, Z, N_sub, gh, bigh, time_Xikj_k,
                     gfun.inv.k, gfun.inv.deriv.k, 
                     gamma.hat.k, eta.hat.k,
                     beta_init, beta_lower, beta_upper,
                     sigma_init, sigma_lower, sigma_upper,
                     factr = 1e7, abstol = 1e-8, 
                     maxit_out = 10, maxit_in = 100, verbose=FALSE){
  
  n_iter = 0
  beta_iter = beta_init; sigma_iter = sigma_init;
  params_iter = c(beta_iter, sigma_iter)
  # params_diff_2norm  = params_iter^2 %>% sum %>% sqrt
  params_diff_2norm  = 100
  
  kernel.func.k <- kernel.func.list[[1]] # for equations l.k
  an.k <- an.vec[1]
  kernel.func.kk <- kernel.func.list[[2]] # for equations l.kk and l.kl
  an.kk <- an.vec[2]

  t_X_k = time_Xikj_k %>%
    filter(time > t-an.k & time < t+an.k)
  # if (nrow(t_X_k)==0){next}
  t_X_kk = time_Xikj_k %>%
    filter(time > t-an.kk & time < t+an.kk)
  
  while(n_iter < maxit_out & params_diff_2norm > abstol){
    n_iter = n_iter + 1
    print(paste(n_iter, "is started."))
    
    beta_fit = optim(fn=obj.func.k.v3,
                     gr=grad.func.k.v3,
                     hessian=TRUE,
                     kernel.func.k = kernel.func.k, kernel.func.kk = kernel.func.kk,
                     an.k = an.k, an.kk = an.kk,
                     t = t, k = k,
                     Z = Z, N_sub = N_sub,
                     gh = gh, bigh = bigh,
                     t_X_k = t_X_k, t_X_kk = t_X_kk,
                     gfun.inv.k = gfun.inv.k,
                     gfun.inv.deriv.k = gfun.inv.deriv.k,
                     gamma.hat.k = gamma.hat.k,
                     eta.hat.k = eta.hat.k,
                     par = beta_iter,
                     sigma.sq.fix.k = sigma_iter,
                     var.deriv.k = FALSE,
                     lower = beta_lower, upper = beta_upper,
                     method="L-BFGS-B", control=list(trace=6, maxit=maxit_in, factr=factr))
    beta_iter = beta_fit$par
    
    sigma_fit = optim(fn=obj.func.kk.v2,
                      gr=grad.func.kk.v2,
                      hessian=TRUE,
                      kernel.func.kk = kernel.func.kk,
                      an.kk = an.kk,
                      t = t, k = k,
                      Z = Z, N_sub = N_sub, 
                      gh = gh, bigh = bigh,
                      t_X_kk = t_X_kk,
                      gfun.inv.k = gfun.inv.k,
                      gfun.inv.deriv.k = gfun.inv.deriv.k,
                      gamma.hat.k = gamma.hat.k,
                      eta.hat.k = eta.hat.k,
                      par = sigma_iter,
                      beta.fix.k = beta_iter,
                      var.deriv.k = TRUE,
                      lower = sigma_lower, upper = sigma_upper,
                      method="L-BFGS-B", control=list(trace=6, maxit=maxit_in, factr=factr))
    sigma_iter = sigma_fit$par
    
    params_diff_2norm = (c(beta_iter,sigma_iter)-params_iter)^2 %>% sum %>% sqrt
    params_iter = c(beta_iter, sigma_iter)
    
    if (verbose){
      print(c("beta:", round(beta_iter,4)) %>% paste(collapse=" "))
      print(c("sigma:", round(sigma_iter,4)) %>% paste(collapse=" "))
      print(c("param diff:", round(params_diff_2norm,8)) %>% paste(collapse=" "))
    }
  } # while(n_iter < maxit & params_diff_2norm > abstol)
  
  convergence = ifelse(n_iter == maxit_out & params_diff_2norm > abstol, 1, 0)
  
  return(list(par = params_iter, 
              convergence = convergence, 
              iterations = n_iter, 
              params_diff = params_diff_2norm,
              beta_fit = beta_fit,
              sigma_fit = sigma_fit))
} # my_iter_k_v2

#### Obj kl; only k,l ####
obj.func.kl.v2 <- function(parmt, k, l, kernel.func.kl, an.kl, t, Z, N_sub,
                        t_X_1=time_Xikj, t_X_2=time_Xilj,
                        beta.k, gfun.inv.k, gfun.inv.k.deriv, 
                        beta.l, gfun.inv.l, gfun.inv.l.deriv,
                        n=10, mu.12=c(0,0), sigma.sq.k, sigma.sq.l, prune.fraction=NULL,
                        gamma.hat.k, gamma.hat.l, eta.hat.k, eta.hat.l,
                        var.deriv.k, var.deriv.l, gh, bigh){
  
  l_kl <- l.kl.and.derivs.v2(k, l, kernel.func.kl, an.kl, t, Z, t_X_1, t_X_2,
                              beta.k, gfun.inv.k, gfun.inv.k.deriv, beta.l, gfun.inv.l, gfun.inv.l.deriv,
                              n=10, mu.12=c(0,0), sigma.sq.k, sigma.sq.l, parmt, prune.fraction=NULL,
                              gamma.hat.k, gamma.hat.l, eta.hat.k, eta.hat.l,
                              var.deriv.k, var.deriv.l, bigh)$val
  
  return(l_kl^2)
}

grad.func.kl.v2 <- function(parmt, k, l, kernel.func.kl, an.kl, t, Z, N_sub,
                         t_X_1=time_Xikj, t_X_2=time_Xilj,
                         beta.k, gfun.inv.k, gfun.inv.k.deriv, 
                         beta.l, gfun.inv.l, gfun.inv.l.deriv,
                         n=10, mu.12=c(0,0), sigma.sq.k, sigma.sq.l, prune.fraction=NULL,
                         gamma.hat.k, gamma.hat.l, eta.hat.k, eta.hat.l,
                         var.deriv.k, var.deriv.l, gh, bigh){
  
  l.kl.results <- l.kl.and.derivs.v2(k, l, kernel.func.kl, an.kl, t, Z, t_X_1, t_X_2,
                                      beta.k, gfun.inv.k, gfun.inv.k.deriv, beta.l, gfun.inv.l, gfun.inv.l.deriv, 
                                      n=10, mu.12=c(0,0), sigma.sq.k, sigma.sq.l, parmt, prune.fraction=NULL,
                                      gamma.hat.k, gamma.hat.l, eta.hat.k, eta.hat.l,
                                      var.deriv.k, var.deriv.l, bigh)
  
  return(2*l.kl.results$val * l.kl.results$sigma.kl.deriv)
}
