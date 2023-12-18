# File: simulation.R
# Version: 1.0-2
# Date: 2023-05-07
# Author: Alesandro Sarnaglia <alessandro.sarnaglia@ufes.br>
# Depends: R (>= 1.8.0), lmtest, MASS
# Description: Functions to fit robust segmented regression.
# License: GPL (>= 2)
# URL: https://github.com/ffajardo64/segmented-regression

rm(list=ls(all=TRUE))

Ind<-function(cond=TRUE){
  return(as.numeric(cond))
}

source('./M_segmentedV2.R')
#source('/home/ffajardo/Dropbox/felix_sarnaglia_fajardo/ALGORITMOS/new_simul/simulacao_mod_V_Alessandro.R')
#source('~/Dropbox/alessandro/research/my_papers/felix_sarnaglia_fajardo/ALGORITMOS/new_simul/m_segmented_V_Alessandro_slope3_zero.R')

zeta0<-11
alfa<-25
beta<-c(-21, -4) 
p <- 0.8
psi<-c(1.0, p*12) # 2.5, 5.0, 7.5
n<-600 # 300, 600
z<-sort(runif(n=n,min = 0,max = 12))
reta<- zeta0 + alfa*z +
  beta[1]*(z-psi[1])*(z>psi[1]) +
  beta[2]*(z-psi[2])*(z>psi[2])

wClas<-wRob<-eClas<-eRob<-0
ResClas<-ResRob.bsqr <- ResRob.huber <- ResRob.hampel<-NULL

MC<-1000

##-------------------------------------------------------------------##
##                              Contaminated scenario
##-------------------------------------------------------------------##

mc<-1
n_out <- 2 # Number of outliers

sample.tmp <- table(ifelse(z<psi[1],0,ifelse(psi[1]<=z & z<=psi[2],1,2)))['0']:table(ifelse(z<psi[1],0,ifelse(psi[1]<=z & z<=psi[2],1,2)))['1']

s1 <- sample(sample.tmp,n_out,replace = F) # Second segment

while (mc<=MC){
  if(mc%%100==0)
    cat('End experiment:',mc,'\n')
  ##-------------------------------------------------------------------##
  ##                    contaminated distributions
  ##-------------------------------------------------------------------##
  
  # e <- 0.9*rnorm(n*0.5,0,1) + 0.1*rt(n*0.5,df=4)
  e <- 0.9*rnorm(n*0.5,0,1) + 0.1*rcauchy(n*0.5)
  
  y <- reta + e
  
  res1<-tryCatch(M_segmented(Y = y, Z = z, psi0 = psi, last_slope_zero = T), error=function(e){NULL})
  res2<-tryCatch(M_segmented(Y = y, Z = z, psi0 = psi, last_slope_zero = T, method = 'robust',psi_func = psi.bisquare), error=function(e){NULL})
  res3<-tryCatch(M_segmented(Y = y, Z = z, psi0 = psi, last_slope_zero = T, method = 'robust',psi_func = psi.huber), error=function(e){NULL})
  res4<-tryCatch(M_segmented(Y = y, Z = z, psi0 = psi, last_slope_zero = T, method = 'robust',psi_func = psi.hampel), error=function(e){NULL})
  
  
  if((!is.null(res1))&(!is.null(res2))){
    ResClas <-rbind(ResClas, res1$fit$coef)
    ResRob.bsqr <-rbind(ResRob.bsqr, res2$fit$coef)
    ResRob.huber <-rbind(ResRob.huber, res3$fit$coef)
    ResRob.hampel <-rbind(ResRob.hampel, res4$fit$coef)
    wClas<-wClas+res1$w;wRob<-wRob+res2$w
    mc<-mc+1
  }
  else{
    cat('Fim Réplica:',mc,'Deu ruim\n')
    eClas<-eClas+as.numeric(is.null(res1))
    eRob<-eRob+as.numeric(is.null(res2))
  }
}

##-------------------------------------------------------------------##
##                    calculating averages
##-------------------------------------------------------------------##

classico <- rbind(round((colMeans(ResClas)-c(zeta0, alfa, beta, psi)), digits = 4),
                  round(sqrt(colMeans(ResClas^2)-colMeans(ResClas)^2+(c(zeta0,alfa,beta,psi)-colMeans(ResClas))^2), digits=4))
rownames(classico)<-c('Bias', 'RMSE')

robusto.bsqr <- rbind(round((colMeans(ResRob.bsqr)-c(zeta0, alfa, beta, psi)), digits = 4),round(sqrt(colMeans(ResRob.bsqr^2)-colMeans(ResRob.bsqr)^2+(c(zeta0,alfa,beta,psi)-colMeans(ResRob.bsqr))^2), digits=4))
rownames(robusto.bsqr)<-c('Bias', 'RMSE')

robusto.huber <- rbind(round((colMeans(ResRob.huber)-c(zeta0, alfa, beta, psi)), digits = 4),round(sqrt(colMeans(ResRob.huber^2)-colMeans(ResRob.huber)^2+(c(zeta0,alfa,beta,psi)-colMeans(ResRob.huber))^2), digits=4))
rownames(robusto.huber)<-c('Bias', 'RMSE')

robusto.hampel <- rbind(round((colMeans(ResRob.hampel)-c(zeta0, alfa, beta, psi)), digits = 4),round(sqrt(colMeans(ResRob.hampel^2)-colMeans(ResRob.hampel)^2+(c(zeta0,alfa,beta,psi)-colMeans(ResRob.hampel))^2), digits=4))
rownames(robusto.hampel)<-c('Bias', 'RMSE')

cat("Classical method\n");t(classico);cat("\nRobust method (Bisquare)\n");t(robusto.bsqr);cat("\nRobust method (Huber)\n");t(robusto.huber);cat("\nRobust method (Hampel)\n");t(robusto.hampel)

cbind(t(classico), t(robusto.bsqr), t(robusto.huber), t(robusto.hampel))