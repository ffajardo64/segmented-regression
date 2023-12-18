Ind<-function(cond=TRUE){
  return(as.numeric(cond))
}

source('./m_segmented.R')

zeta0<-2
alfa<-0
beta<-1# 1, 2
psi<-5.0 # 2.5, 5.0, 7.5
n<-100 # 300, 600
z<-sort(runif(n=n,min = 0,max = 10))
reta<- zeta0 + alfa*z + beta*(z-psi)*(z>psi)
MC<-1000

n_out <- 5

sample.tmp <- table(ifelse(reta==zeta0,0,1))['0']:(table(ifelse(reta==zeta0,0,1))['0']+table(ifelse(reta==zeta0,0,1))['1'])

#s1 <- sample((1:n)[!1:n %in% sample.tmp],n_out,replace = F) # Primeiro segmento
s1 <- sample(sample.tmp,n_out,replace = F)# Segundo segmento




wClas<-wRob<-eClas<-eRob<-0
ResClas<-ResRob<-NULL
mc<-1;df <- 3
while (mc<=MC){
  if(mc%%100==0)
    cat('Fim Réplica:',mc,'\n')
  e <- rnorm(n=n, mean = 0, sd = 1) # Sem valores atípicos
  for(j in s1){
    e[j]<-e[j] + rbinom(1,20,0.5)
  }
  
  y <- reta + e # Sem valores at?picos
  # y<- reta + 1*rt(n = n,df = df)/sqrt(df/(df-2)) # Com valores at?picos
  # y<- reta + (rchisq(n = n,df = df) - df)/sqrt(2*df)# Com valores at?picos
  # fit1<-lm(y~z)
  # require(segmented)
  # segmentado<-segmented(obj = fit1,seg.Z = ~z,psi = psi)
  # summary(segmentado)
  # c(zeta0, alfa, beta, psi)

  res1<-tryCatch(M_segmented(Y = y, Z = z, psi0 = psi), error=function(e){NULL})
  res2<-tryCatch(M_segmented(Y = y, Z = z, psi0 = psi, method = 'robust'), error=function(e){NULL})
  if((!is.null(res1))&(!is.null(res2))){
    ResClas <-rbind(ResClas, res1$coef)
    ResRob <-rbind(ResRob, res2$coef)
    wClas<-wClas+res1$w;wRob<-wRob+res2$w
    mc<-mc+1
  }
  else{
    cat('Fim Réplica:',mc,'Deu ruim\n')
    eClas<-eClas+as.numeric(is.null(res1))
    eRob<-eRob+as.numeric(is.null(res2))
  }
}


erros<-rbind(c(wClas, eClas), c(wRob, eRob))
rownames(erros)<-c('Classic', 'Robust'); colnames(erros)<-c('Warnings', 'Errors')

classico <- rbind(round((colMeans(ResClas)-c(zeta0, alfa, beta, psi)), digits = 4),
      round(sqrt(colMeans(ResClas^2)-colMeans(ResClas)^2+(c(zeta0,alfa,beta,psi)-colMeans(ResClas))^2), digits=4))
rownames(classico)<-c('Bias', 'RMSE')


robusto <- rbind(round((colMeans(ResRob)-c(zeta0, alfa, beta, psi)), digits = 4),round(sqrt(colMeans(ResRob^2)-colMeans(ResRob)^2+(c(zeta0,alfa,beta,psi)-colMeans(ResRob))^2), digits=4))
rownames(robusto)<-c('Bias', 'RMSE')

# ------ * Figuras * ---------
plot(z,y,pch=19, cex=ifelse(1:n %in% s1, 1, 0.5),
     col=ifelse(1:n %in% s1, "red", "black"))
text(z+0.3,y, labels=ifelse(1:n %in% s1, s1, ""), cex=0.8)
lines(z,reta, lwd=2)
lines(z, res1$predict, col='red', lwd=2) # Classic
lines(z, res2$predict, col='blue', lwd=2) # Robust
legend(x = "bottomright",          # Position
       legend = c("Original","Classic", "Robust"),  # Legend texts
       #lty = c(1, 2),           # Line types
       col = c("black", "red", "blue"),           # Line colors
       lwd = 2)

# ------ * Resultados * ---------
print("Classical method");classico;print("Robust method");robusto
print("Erros"); erros
