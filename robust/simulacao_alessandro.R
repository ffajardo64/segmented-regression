Ind<-function(cond=TRUE){
  return(as.numeric(cond))
}

source('m_segmented.R')

zeta0<-2
alfa<-0
beta<-1# 1, 2
psi<-2.5 # 2.5, 5.0, 7.5
n<-300 # 300, 600
z<-runif(n=n,min = 0,max = 10)
reta<- zeta0 + alfa*z + beta*(z-psi)*(z>psi)
MC<-1000

wClas<-wRob<-eClas<-eRob<-0
ResClas<-ResRob<-NULL
mc<-1
while (mc<=MC){
  cat('Fim R?plica:',mc)
  # y<- reta + rnorm(n=n,mean = 0,sd = 1) # Sem valores at?picos
  y<- reta + 1*rt(n = n,df = 3)/sqrt(3) # Com valores at?picos
  # fit1<-lm(y~z)
  # require(segmented)
  # segmentado<-segmented(obj = fit1,seg.Z = ~z,psi = psi)
  # summary(segmentado)
  # c(zeta0, alfa, beta, psi)
  res1<-tryCatch(M_segmented(Y = y, Z = z, psi0 = psi), error=function(e){NULL})
  res2<-tryCatch(M_segmented(Y = y, Z = z, psi0 = psi, method = 'robust'), error=function(e){NULL})
  if((!is.null(res1))&(!is.null(res2))){
    ResClas <-rbind(ResClas, res1$result)
    ResRob <-rbind(ResRob, res2$result)
    wClas<-wClas+res1$w;wRob<-wRob+res2$w
    mc<-mc+1
    cat('\n')
  }
  else{
    cat('Deu ruim\n')
    eClas<-eClas+as.numeric(is.null(res1))
    eRob<-eRob+as.numeric(is.null(res2))
  }
}

bias<-rbind(round((colMeans(ResClas)-c(zeta0, alfa, beta, psi)), digits = 4),
            round((colMeans(ResRob)-c(zeta0, alfa, beta, psi)), digits = 4))
rownames(bias)<-c('Classic', 'Robust')

rmse<-rbind(round(sqrt(colMeans(ResClas^2)-colMeans(ResClas)^2+(c(zeta0,alfa,beta,psi)-colMeans(ResClas))^2), digits=4),
            round(sqrt(colMeans(ResRob^2)-colMeans(ResRob)^2+(c(zeta0,alfa,beta,psi)-colMeans(ResRob))^2), digits=4))
rownames(rmse)<-c('Classic', 'Robust')
erros<-rbind(c(wClas, eClas), c(wRob, eRob))
rownames(erros)<-c('Classic', 'Robust'); colnames(erros)<-c('Warnings', 'Errors')

bias
rmse
erros
