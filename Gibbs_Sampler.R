library(mvtnorm)
library(MCMCpack)
library(Matrix)

funv <- function(v0,a,b,samples){
  
  vs <- c();
  vs[1] <- v0;
  
  fz <- function(z){ 
    
    (( z/(z+3))^(1/2))*(( trigamma(z/2) - trigamma((z+1)/2) -(2*(z+3))/((z*(z+1)^2))  )^(1/2))
  }
  

  for (i in 2:(samples+1)){
    
    v_prop <- rgamma(1,2,1);
    
    bvs <- exp(-(b*vs[i-1])/2)
    if (bvs < 1e-20) {bvs<-1e-20} else {bvs<-bvs}

    anum<-a^(v_prop/2-1)
    if (is.infinite(anum)) {anum = 10^100} else {anum = anum}
    
    adenom<-a^(vs[i-1]/2-1)
    if (is.infinite(adenom)) {adenom = 10^100} else {adenom = adenom} 
    
    rr <- ((anum)*(exp(-(b*v_prop)/2))*fz(v_prop)*dgamma(vs[i-1],2,1) )/ ((adenom)*(bvs)*fz(vs[i-1])*dgamma(v_prop, 2,1)  );
  
    r <- min(rr,1);
    
    u <- runif(1,0,1);
    
    
    if (u<r) { vs[i] <- v_prop;} else vs[i] <- vs[i-1];
    
  }
  
  sample.v <- vs[samples];
  
  return(sample.v);
  
}


obmvrstd = function(x,y,max.steps=1000){
  
  n <- length(x[,1])
  p <- length(x[1,])
  d <- length(y[1,])
  
  v <- 2
  
  betaSamples <- array(0,dim=c(max.steps,d,p))
  vsamples <- array(0,dim=c(max.steps))
  vsamples[1] <- v

  invSigma <- diag(1,d,d)
  w <- rep(1,p)
  ee <- 1
  lambda2 <- 1

  xtx <- t(x)%*%x
  beta <- t(ginv(xtx)%*%t(x)%*%y)
  
  k <- 0
  
  while (k < max.steps) {
    k <- k + 1
    
    if (k %% 10 == 0) {
      cat('Iteration:', k, "\r")
    }
  

    ###### Q and V ######
    q <- rep(NA,n)
    
    for (i in seq(n)){
      q[i] <- rgamma(1,(v+d)/2,rate=((t(as.matrix(t(y)[,i])-beta%*%x[i,])%*%invSigma%*%(as.matrix(t(y)[,i])-beta%*%x[i,]))+v)/2)
    }
    
    invQ <- diag(q)
    
    if (is.infinite(prod(q))) {pq = 10^50} else {pq = prod(q)}

#v <- funv(vsamples[k],pq,sum(q),2)
    v <- 2
    vsamples[k+1] <- v
   
    ####### TAU and LAMBDA ###########
    
    tau2 <- rep(NA,p)
    invtau2 <- rep(NA,p)
    for (j in seq(p)){
      
      tau2[j] <- rinvgamma(1,(d+1)/2,(1/w[j])+((t(beta[,j])%*%invSigma%*%beta[,j])/(2*lambda2)))
      invtau2[j] <- 1/tau2[j]
    }
    
    lb <- rep(NA,p)
    for (j in seq(p)){
      
      lb[j] <- (invtau2[j]*(t(beta[,j])%*%invSigma%*%beta[,j]))
      
    }
    
    lambda2 <- rinvgamma(1,((p*d)+1)/2,(1/ee)+(0.5*sum(lb)))
    invlambda2 <- 1/lambda2
    
    
    ######### Wj and ee ###########
    
    for (j in seq(p)){
      w[j] <- rinvgamma(1,1,1+(1/tau2[j]))
    }
    
    ee <- rinvgamma(1,1,1+(1/lambda2))
    
    ###### OMEGA ########
    
    Dtau <- diag(invtau2*invlambda2,p)
    
    invOmega <- t(x)%*%invQ%*%as.matrix(x) + Dtau
    eig<-eigen(invOmega)
    Omega<-eig$vectors%*%diag(1/eig$values)%*%t(eig$vectors)
    
    ####### MU #########
    
    mu <- Omega%*%t(x)%*%invQ%*%as.matrix(y)
    
    ####### PHI #########
    
    phi <- t(y)%*%invQ%*%as.matrix(y) - t(mu)%*%invOmega%*%as.matrix(mu)
  
    ######## Sigma #########
    
    Si <- tryCatch(riwish(n,phi), error=function(e) {Si=1})
    if (!is.matrix(Si)) {k=k-1;}
    else{
    eig2<-eigen(Si)
    invSigma<-eig2$vectors%*%diag(1/eig2$values)%*%t(eig2$vectors)
    }
    ######## BETA #########
    
    #Kr <- kronecker(Si,Omega)
    #eig4<-eigen(Kr)
    #sqrKr<-diag(sqrt(eig4$values))%*%t(eig4$vectors)
    #vecmu <- matrix(as.vector(mu),1,p*d)
    #nor<-matrix(rnorm(p*d,0,1),1,p*d)
    #vecbeta <- t(sqrKr)%*%t(nor) + t(vecmu)
    
    #trbeta <- matrix(c(t(vecbeta)),p,d)

    beta <- tryCatch(rmatrixnorm(t(mu),as.positive.definite(Si),as.positive.definite(Omega)), error=function(e) {beta=1})
    if (!is.matrix(beta)) {k=k-1; beta=betaSamples[k,,];} 
    else {betaSamples[k,,] <- beta}
    
  }
  list(betaSamples,vsamples)
}
