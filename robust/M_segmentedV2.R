# File: M_segmentedV2.R
# Version: 2.1-2
# Date: 2023-04-27
# Author: Alesandro Sarnaglia <alessandro.sarnaglia@ufes.br>
# Depends: R (>= 1.8.0), lmtest, MASS
# Description: Functions to fit robust segmented regression.
# License: GPL (>= 2)
# URL: https://github.com/ffajardo64/segmented-regression


Plus <- function(t) return(t*Ind0(t))
Ind0 <- function(t) return(as.numeric(t>0))

##-------------------------------------------------------------------##
##                              seg_reg_fixed_psi function
##-------------------------------------------------------------------##

seg_reg_fixed_psi <- function(Y, X=NULL, Z, psi, method='classic',
                              last_slope_zero = FALSE, psi_func=NULL){
  n <- length(Y); q <- length(psi)
  U <- matrix(nrow = n,ncol = q)
  for(k in 1:q){
    U[,k]<- (Z - psi[k])*Ind0(Z - psi[k])
  }
  if(last_slope_zero){
    Z_ <- as.matrix(Z - U[,q]); U_ <- as.matrix(U[,-q] - U[,q])
  }else{
    Z_ <- as.matrix(Z); U_ <- as.matrix(U)
  }
  
  if(method=='classic') fit<-lm(Y ~ Z_ + U_) else fit<-rlm(Y ~ Z_ + U_, psi = psi_func)
  if(method == 'classic') return(sd(fit$residuals)) else return(mad(fit$residuals))
}

##-------------------------------------------------------------------##
##                              grid_search_seg function
##-------------------------------------------------------------------##

grid_search_seg <- function(Y, X=NULL, Z, nintervals, npsi, method='classic',
                            last_slope_zero = FALSE, psi_func=NULL){
  suppressPackageStartupMessages(if(!require(MASS)) install.packages("MASS", repos = "http://cran.us.r-project.org"))
  qs <- quantile(x = Z, probs = ((1:(nintervals-1))/nintervals))
  lst <- list()
  for (i in 1:npsi) {
    lst[[i]] <- qs
  }
  names(lst) <- paste('psi_',1:npsi,sep='')
  grd_F <- as.matrix(expand.grid(lst)); grd <- NULL
  for(i in 1:nrow(grd_F)){
    if(all(diff(grd_F[i,]) > 0)) grd <- rbind(grd, grd_F[i,])
  }
  var_errors_v <- NULL
  for(i in 1:nrow(grd)){
    var_errors_v[i] <- seg_reg_fixed_psi(Y = Y,X = X, Z = Z, psi = grd[i,],
                                         method = method, last_slope_zero = last_slope_zero,
                                         psi_func = psi_func)
  }
  return(cbind(grd, var_res = var_errors_v))
}


##-------------------------------------------------------------------##
##                              M_segmented_aux function
##-------------------------------------------------------------------##

M_segmented_aux <-function(Y, X=NULL, Z, psi0, method='classic',
                      last_slope_zero = FALSE, psi_func=NULL, boot_resample=FALSE){
  
  suppressPackageStartupMessages(if(!require(lmtest)) install.packages("lmtest", repos = "http://cran.us.r-project.org"))
  suppressPackageStartupMessages(if(!require(MASS)) install.packages("MASS", repos = "http://cran.us.r-project.org"))
  
  if(is.null(X)){
    g1<-10; g2<-.01; mZ <- min(Z); MZ <- max(Z)
    q<-length(psi0)
    n<-length(Y)
    psi<- psi0
    s<-0
    if(method!='classic' & is.null(psi_func)){psi_func = psi.huber; warning('Null psi function. Using Huber!')}
    while(g1 > .1 & g2 < .1 & s < 50){  ### If the convergence criterion is the psi update
      s<-s+1
      U<-V<-matrix(nrow = n,ncol = q)
      for(k in 1:q){
        U[,k]<- (Z - psi[k])*Ind0(Z - psi[k])
        V[,k]<- -Ind0(Z - psi[k])
      }
      if(last_slope_zero){
        Z_ <- as.matrix(Z - U[,q]); U_ <- as.matrix(U[,-q] - U[,q]); V_ <- as.matrix(V)
      }else{
        Z_ <- as.matrix(Z); U_ <- as.matrix(U); V_ <- as.matrix(V)
      }

      if(method=='classic') fit<-lm(Y ~ Z_ + U_ + V_) else fit<-rlm(Y ~ Z_ + U_ + V_, psi = psi_func)
      coeficientes<-fit$coefficients
      
      alpha<-coeficientes[grep('Z_', names(coeficientes))]
      beta<-coeficientes[grep('U_', names(coeficientes))]
      
      gama<-coeficientes[grep('V_', names(coeficientes))]
      
      if(length(beta) < length(gama)) beta <- c(beta, -alpha-sum(beta))
      zeta0<- coeficientes[1]
      psinovo<- (gama/beta)/(1.2^s) + psi
      
      g1 <- max(abs(gama/beta)/sqrt(s))
      g2 <- min(coeftest(fit)[grep('V_', rownames(coeftest(fit))),4])
      
      
      if(boot_resample){
        psinovo[psinovo > quantile(x = Z, probs=.95)] <- quantile(x = Z, probs=.95)
        psinovo[psinovo < quantile(x = Z, probs=.05)] <- quantile(x = Z, probs=.05)
      }
      
      if(!boot_resample){
        if(any(psinovo >= MZ | psinovo <= mZ)){
          s <- 51
          psinovo[psinovo >= MZ] <- MZ
          psinovo[psinovo <= mZ] <- mZ
        }
      }
      psi <- psinovo
    } # 
    w<-0
    if(s==50) {warning('Did not converge!'); w<-1}
    if(s==51) {warning('Psi outside bounds of Z!'); w<-1}
    
    reta <- zeta0 + alpha*Z #+ beta*(Z-psi)*(Z>psi)
    for(i in 1:q){reta <- reta + beta[i]*(Z - psi[i])*(Z > psi[i])}
    res <- Y - reta

    if(method == 'classic') sig <- sd(res) else sig <- mad(res)
    
    result<-c(zeta0,alpha,beta,psi,sig^2)
    
    names(result)<-c("zeta_0","alpha",kronecker('beta_',1:q,paste,sep=''),kronecker('psi_',1:q,paste,sep=''), 'sig2')
    
    return(list(coef=result, w=w, residuals=res, predict=reta))
  }
  else{
    stop('Not implemented yet!')
  }
}

##-------------------------------------------------------------------##
##                              M_segmented function
##-------------------------------------------------------------------##

M_segmented <- function(Y, X=NULL, Z, psi0, method='classic',
                        last_slope_zero = TRUE, n.boot.se=0,
                        psi_func=NULL){
  if(method!='classic' & is.null(psi_func)){psi_func = psi.huber; warning('Null psi function. Using Huber!')}
  n <- length(Y); psi0v <- NULL
  bb <- 1
  while (bb <= 15) {
    ii <- sample(x = 1:n, size = n, replace = T)
    fit_b <- tryCatch(M_segmented_aux(Y = Y[ii], X = X[ii], Z = Z[ii], psi0 = psi0, method = method,
                                      last_slope_zero = last_slope_zero,
                                      psi_func = psi_func, boot_resample=TRUE), error = function(e) NULL)
    if(!is.null(fit_b)){
      psi0v <- rbind(psi0v, fit_b$coef[grep('psi', names(fit_b$coef))])
      bb <- bb+1
    }
  }
  psi00 <- colMeans(psi0v)
  
  if(n.boot.se == 0){
    return(
      list(
        fit = M_segmented_aux(Y = Y, X = X, Z = Z, psi0 = psi00, method = method,
                              last_slope_zero = last_slope_zero,
                              psi_func = psi_func),
        entries = list(Y=Y, X=X, Z=Z, psi0=psi0, method=method,
                       last_slope_zero = last_slope_zero, n.boot.se=n.boot.se,
                       psi_func=psi_func))
      )
  }else{
    fit0 <- M_segmented_aux(Y = Y, X = X, Z = Z, psi0 = psi00, method = method,
                            last_slope_zero = last_slope_zero,
                            psi_func = psi_func)
    n <- length(Y)
    coef_star <- NULL
    for(b in 1:n.boot.se){
      r_star <- sample(x = fit0$residuals, size = n, replace = TRUE)
      Y_star <- fit0$predict + r_star
      fit_star <- M_segmented_aux(Y = Y_star, X = X, Z = Z,
                                  psi0 = psi00,
                                  method = method, last_slope_zero = last_slope_zero,
                                  psi_func = psi_func)
      coef_star <- rbind(coef_star, fit_star$coef)
    }
    res <- cbind(estimate = fit0$coef,
                 mean = colMeans(coef_star),
                 se = sqrt(colMeans(coef_star^2)-colMeans(coef_star)^2))
    rownames(res) <- names(fit0$coef)
    return(list(coef = res,
                fit = fit0,
                entries = list(Y=Y, X=X, Z=Z, psi0=psi0, psi00=psi00, method=method,
                               last_slope_zero = last_slope_zero, n.boot.se=n.boot.se,
                               psi_func=psi_func)))
  }
}

##-------------------------------------------------------------------##
##                      M_seg_boot_plateau function
##-------------------------------------------------------------------##

M_seg_boot_plateau <- function(fit, n.boot=200, psi0=NULL){
  cfs <- fit$fit$coef
  ff <- fit$fit
  n <-length(ff$residuals)
  zeta0 <- cfs['zeta_0']
  alpha <- cfs['alpha']
  beta <- cfs[grep(pattern = 'beta', names(cfs))]
  psi <- cfs[grep(pattern = 'psi', names(cfs))]
  
  if(is.null(psi0)) psi0 <- fit$entries$psi00
  q <- length(psi)
  reta <- zeta0 + alpha*fit$entries$Z
  for(i in 1:(q-1)) reta <- reta + beta[i]*(fit$entries$Z - psi[i])*(fit$entries$Z > psi[i])
  psi_q <- psi[q]
  psi_q_boot <- NULL
  for (b in 1:n.boot) {
    if(b%%50 == 0) cat('rep:',b,'\n')
    ii <- sample(1:n, n, TRUE)
    r_star <- ff$residuals[ii]
    Y_star <- reta + r_star
    
    fit_star <- M_segmented_aux(Y = Y_star, X = fit$entries$X, Z = fit$entries$Z,
                                psi0 = psi0,
                                method = fit$entries$method,
                                last_slope_zero = fit$entries$last_slope_zero,
                                psi_func = fit$entries$psi_func, boot_resample = FALSE)
    psi_q_boot[b] <- (fit_star$coef[grep(pattern = 'psi', names(fit_star$coef))])[q]
  }
  return((sum(psi_q_boot < psi_q)+1)/(n.boot+1))
}
