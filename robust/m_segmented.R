M_segmented<-function(Y, X=NULL, Z, psi0, method='classic'){
  # require(lmtest)
  if(method!='classic') require(MASS)
  if(is.null(X)){
    g<-.01
    q<-length(psi0)
    n<-length(Y)
    psi<- psi0
    s<-0
    while((g<.1)&(s<50)){
      # cat(g,'\n')
      s<-s+1
      U<-V<-matrix(nrow = n,ncol = q)
      for(i in 1:n){
        for(k in 1:q){
          U[i,k]<- (Z[i] - psi[k])*Ind(Z[i] > psi[k])
          V[i,k]<- -Ind(Z[i] > psi[k])
        }
      }
      
      if(method=='classic') fit<-lm(Y ~ Z + U + V) else fit<-rlm(Y ~ Z + U + V, psi = psi.huber)
      coeficientes<-fit$coefficients
      # stderror<- sqrt(diag(vcov(fit)))
      beta<-coeficientes[3:(q+2)] #de inicio beta ate o final de beta
      gama<-coeficientes[(q+3):(2*q+2)] #de o final do beta (um depois, no caso come?o de gama) ate o final do vetor
      zeta0<- coeficientes[1]
      alfa<-coeficientes[2]
      psinovo<- (gama/beta) + psi 
      # g<- max(abs(gama/beta))
      g<-min(coeftest(fit)[(q+3):(2*q+2),4])
      psi <- psinovo
    }
    w<-0
    if(s==50) {warning('Did not converge!'); w<-1}
    result<-c(zeta0,alfa,beta,psi)
    
    reta <- zeta0 + alfa*Z + beta*(Z-psi)*(Z>psi)
    res <- Y - reta
    
    # se<-c(stderror)
    names(result)<-c("zeta_0","alfa",kronecker('beta_',1:q,paste,sep=''),kronecker('psi_',1:q,paste,sep=''))
    #names(se)<-c("zeta_0","alfa",kronecker('beta_',1:q,paste,sep=''),kronecker('psi_',1:q,paste,sep=''))
    return(list(coef=result, w=w, residuals=res, predict=reta))
  }
  
  else{
    stop('Not implemented yet!')
  }
}