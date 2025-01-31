LW1 = function(y,alphas,betas,tau2s,xs,delta){
  n  = length(y)
  N  = length(xs)
  quants = array(0,c(n,4,3))
  parss<-array(0,c(N,3,n))
  h2 = 1-((3*delta-1)/(2*delta))^2
  a  = sqrt(1-h2)
  pars = cbind(alphas,betas,log(tau2s))
  xss<-NULL
  ws<-NULL
  ESS<-NULL
  #like = rep(0,n)
  par(mfrow=c(1,1))
  for (t in 1:n){
    #like[t] = sum(dnorm(y[t],0.0,exp(xs/2)))
    # Resampling
    mus     = pars[,1]+pars[,2]*xs
    mpar    = apply(pars,2,mean)
    vpar    = var(pars)
    ms      = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    weight  = dnorm(y[t],0.0,exp(mus/2),log=TRUE)
    weight1 = exp(weight-max(weight))
    k       = sample(1:N,size=N,replace=T,prob=weight1)
    # Propagating
    ms1 = ms[k,]+matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    xt   = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    w    = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=T,prob=w)
    xs   = xt[ind]
    pars = ms1[ind,]
    xss<-rbind(xss,xs)
    parss[,,t]<-pars
    ws<-rbind(ws,w)
    cv2<-var(w)/(mean(w)^2)
    ESS<-c(ESS,N/(1+cv2))
    ts.plot(ESS,xlim=c(1,n))
    # Storing quantiles
    quants[t,1,] = quantile(pars[,1],c(0.5,0.025,0.975))
    quants[t,2,] = quantile(pars[,2],c(0.5,0.025,0.975))
    quants[t,3,] = quantile(exp(pars[,3]),c(0.5,0.025,0.975))
    quants[t,4,] = quantile(exp(xs/2),c(0.5,0.025,0.975))  
  }
  return(list(quants=quants,parss=parss,xss=xss))
}
LWW = function(y,alphas,betas,tau2s,xs,delta,lev,M=75,method=1){
  n  = length(y)
  N  = length(xs)
  quants = array(0,c(n,4,3))
  parss<-array(0,c(N,3,n))
  h2 = 1-((3*delta-1)/(2*delta))^2
  a  = sqrt(1-h2)
  pars = cbind(alphas,betas,log(tau2s))
  xss<-NULL
  ws<-NULL
  ESS<-NULL
  #like = rep(0,n)
  par(mfrow=c(1,1))
  for (t in 1:n){
    #like[t] = sum(dnorm(y[t],0.0,exp(xs/2)))
    # Resampling
    mus     = pars[,1]+pars[,2]*xs
    mpar    = apply(pars,2,mean)
    vpar    = var(pars)
    ms      = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    weight  = dnorm(y[t],0.0,exp(mus/2),log=TRUE)
    #ww  <- plWav1j(weight,filter.number = 4, family = 'DaubLeAsymm', a = 1, b = 1, bet = 1,
    #               nu0 = 5,lamb0 = 1, M=75, Ne = 200, j0=lev, plot.EMPL = FALSE)
    #weight  <- ww$bayesrec5
    weight1 = exp(weight-max(weight))
    k       = sample(1:N,size=N,replace=T,prob=weight1)
    # Propagating
    ms1 = ms[k,]+matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    xt   = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    if(method==1){
      xst  <- plWav1j(xt,filter.number = 4, family = 'DaubLeAsymm', a = 1, b = 1, bet = 1,
                    nu0 = 5,lamb0 = 1, M=M, Ne = 200, j0=lev, plot.EMPL = FALSE)}
    else{
      xst  <- BAYES.THR(xt,alpha=1,beta=1,filter.number = 4, family = 'DaubLeAsymm',plotfn=FALSE,j0=lev,dev=var)
    }
    xt <- xst
    w    = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    #ww  <- plWav1j(w,filter.number = 4, family = 'DaubLeAsymm', a = 1, b = 1, bet = 1,
    #                    nu0 = 5,lamb0 = 1, M=75, Ne = 200, j0=lev, plot.EMPL = FALSE)
    #w   <- ww$bayesrec5
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=T,prob=w)
    xs   = xt[ind]
    #xsw  <- plWav1j(xs,filter.number = 4, family = 'DaubLeAsymm', a = 1, b = 1, bet = 1,
    #                nu0 = 5,lamb0 = 1, M=75, Ne = 200, j0=lev, plot.EMPL = FALSE)
    #xs <- xsw$bayesrec5
    pars = ms1[ind,]
    xss<-rbind(xss,xs)
    parss[,,t]<-pars
    ws<-rbind(ws,w)
    cv2<-var(w)/(mean(w)^2)
    ESS<-c(ESS,N/(1+cv2))
    ts.plot(ESS,xlim=c(1,n))
    # Storing quantiles
    quants[t,1,] = quantile(pars[,1],c(0.5,0.025,0.975))
    quants[t,2,] = quantile(pars[,2],c(0.5,0.025,0.975))
    quants[t,3,] = quantile(exp(pars[,3]),c(0.5,0.025,0.975))
    quants[t,4,] = quantile(exp(xs/2),c(0.5,0.025,0.975))  
  }
  return(list(quants=quants,parss=parss,xss=xss))
}
LWW2 = function(y,alphas,betas,tau2s,xs,delta,lev,M=2,method=1){
  n  = length(y)
  N  = length(xs)
  quants = array(0,c(n,4,3))
  parss<-array(0,c(N,3,n))
  h2 = 1-((3*delta-1)/(2*delta))^2
  a  = sqrt(1-h2)
  pars = cbind(alphas,betas,log(tau2s))
  xss<-NULL
  ws<-NULL
  ESS<-NULL
  #like = rep(0,n)
  par(mfrow=c(1,1))
  for (t in 1:n){
    #like[t] = sum(dnorm(y[t],0.0,exp(xs/2)))
    # Resampling
    mus     = pars[,1]+pars[,2]*xs
    mpar    = apply(pars,2,mean)
    vpar    = var(pars)
    ms      = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    weight  = dnorm(y[t],0.0,exp(mus/2),log=TRUE)
    weight1 = exp(weight-max(weight))
    k       = sample(1:N,size=N,replace=T,prob=weight1)
    # Propagating
    ms1 = ms[k,]+matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    xt   = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    if(method==1){
      xst  <- plWav1j(xt,filter.number = 4, family = 'DaubLeAsymm', a = 1, b = 1, bet = 1,
                    nu0 = 5,lamb0 = 1, M=M, Ne = 2, j0=lev, plot.EMPL = FALSE)}
    else{
      xst  <- BAYES.THR(xt,alpha=1,beta=1,filter.number = 4, family = 'DaubLeAsymm',plotfn=FALSE,j0=lev,dev=var)
    }
    xt <- xst
    w    = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=T,prob=w)
    xs   = xt[ind]
    pars = ms1[ind,]
    xss<-rbind(xss,xs)
    parss[,,t]<-pars
    ws<-rbind(ws,w)
    cv2<-var(w)/(mean(w)^2)
    ESS<-c(ESS,N/(1+cv2))
    ts.plot(ESS,xlim=c(1,n))
    # Storing quantiles
    quants[t,1,] = quantile(pars[,1],c(0.5,0.025,0.975))
    quants[t,2,] = quantile(pars[,2],c(0.5,0.025,0.975))
    quants[t,3,] = quantile(exp(pars[,3]),c(0.5,0.025,0.975))
    quants[t,4,] = quantile(exp(xs/2),c(0.5,0.025,0.975))  
  }
  return(list(quants=quants,parss=parss,xss=xss))
}

LW5 = function(y,alphas,betas,tau2s,xs,delta,lev,M=75){
  n  = length(y)
  N  = length(xs)
  quants = array(0,c(n,4,3))
  parss<-array(0,c(N,3,n))
  h2 = 1-((3*delta-1)/(2*delta))^2
  a  = sqrt(1-h2)
  pars = cbind(alphas,betas,log(tau2s))
  xss<-NULL
  ws<-NULL
  ESS<-NULL
  #like = rep(0,n)
  par(mfrow=c(1,1))
  for (t in 1:n){
    #like[t] = sum(dnorm(y[t],0.0,exp(xs/2)))
    # Resampling
    mus     = pars[,1]+pars[,2]*xs
    mpar    = apply(pars,2,mean)
    vpar    = var(pars)
    ms      = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    weight  = dnorm(y[t],0.0,exp(mus/2),log=TRUE)
    #ww  <- plWav1j(weight,filter.number = 4, family = 'DaubLeAsymm', a = 1, b = 1, bet = 1,
    #               nu0 = 5,lamb0 = 1, M=75, Ne = 200, j0=lev, plot.EMPL = FALSE)
    #weight  <- ww$bayesrec5
    weight1 = exp(weight-max(weight))
    k       = sample(1:N,size=N,replace=T,prob=weight1)
    # Propagating
    ms1 = ms[k,]+matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    xt   = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    xst  <- plWav1j(xt,filter.number = 4, family = 'DaubLeAsymm', a = 1, b = 1, bet = 1,
                    nu0 = 5,lamb0 = 1, M=M, Ne = 200, j0=lev, plot.EMPL = FALSE)
    xt <- xst
    w    = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    #ww  <- plWav1j(w,filter.number = 4, family = 'DaubLeAsymm', a = 1, b = 1, bet = 1,
    #                    nu0 = 5,lamb0 = 1, M=75, Ne = 200, j0=lev, plot.EMPL = FALSE)
    #w   <- ww$bayesrec5
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=T,prob=w)
    xs   = xt[ind]
    #xsw  <- plWav1j(xs,filter.number = 4, family = 'DaubLeAsymm', a = 1, b = 1, bet = 1,
    #                nu0 = 5,lamb0 = 1, M=75, Ne = 200, j0=lev, plot.EMPL = FALSE)
    #xs <- xsw$bayesrec5
    pars = ms1[ind,]
    xss<-rbind(xss,xs)
    parss[,,t]<-pars
    ws<-rbind(ws,w)
    cv2<-var(w)/(mean(w)^2)
    ESS<-c(ESS,N/(1+cv2))
    ts.plot(ESS,xlim=c(1,n))
    # Storing quantiles
    quants[t,1,] = quantile(pars[,1],c(0.5,0.025,0.975))
    quants[t,2,] = quantile(pars[,2],c(0.5,0.025,0.975))
    quants[t,3,] = quantile(exp(pars[,3]),c(0.5,0.025,0.975))
    quants[t,4,] = quantile(exp(xs/2),c(0.5,0.025,0.975))  
  }
  return(list(quants=quants,parss=parss,xss=xss))
}
LW6 = function(y,alphas,betas,tau2s,xs,delta,lev){
  n  = length(y)
  N  = length(xs)
  quants = array(0,c(n,4,3))
  parss<-array(0,c(N,3,n))
  h2 = 1-((3*delta-1)/(2*delta))^2
  a  = sqrt(1-h2)
  pars = cbind(alphas,betas,log(tau2s))
  xss<-NULL
  ws<-NULL
  ESS<-NULL
  #like = rep(0,n)
  par(mfrow=c(1,1))
  for (t in 1:n){
    #like[t] = sum(dnorm(y[t],0.0,exp(xs/2)))
    # Resampling
    mus     = pars[,1]+pars[,2]*xs
    mpar    = apply(pars,2,mean)
    vpar    = var(pars)
    ms      = a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    weight  = dnorm(y[t],0.0,exp(mus/2),log=TRUE)
    weight1 = exp(weight-max(weight))
    k       = sample(1:N,size=N,replace=T,prob=weight1)
    # Propagating
    ms1 = ms[k,]+matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    xt   = rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    xst     <- BAYES.THR(xt,alpha=1,beta=1,filter.number = 4, family = 'DaubLeAsymm',plotfn=FALSE,j0=lev,dev=var)
    xt <- xst
    w    = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=T,prob=w)
    xs   = xt[ind]
    pars = ms1[ind,]
    xss<-rbind(xss,xs)
    parss[,,t]<-pars
    ws<-rbind(ws,w)
    cv2<-var(w)/(mean(w)^2)
    ESS<-c(ESS,N/(1+cv2))
    ts.plot(ESS,xlim=c(1,n))
    # Storing quantiles
    quants[t,1,] = quantile(pars[,1],c(0.5,0.025,0.975))
    quants[t,2,] = quantile(pars[,2],c(0.5,0.025,0.975))
    quants[t,3,] = quantile(exp(pars[,3]),c(0.5,0.025,0.975))
    quants[t,4,] = quantile(exp(xs/2),c(0.5,0.025,0.975))  
  }
  return(list(quants=quants,parss=parss,xss=xss))
}
errorF<-function(y,vec){
  SMem<-matrix(NA,nrow=ncol(vec),ncol=1)
  rownames(SMem)<-names(vec)
  colnames(SMem)<-"RMSE"
  Serr<-NULL
  for(j in 1:ncol(vec)){
    for(i in 1:length(y)){
      Serr[i]<-(vec[i,j]-y[i])^2
    }
    SMem[j,]<-sqrt(sum(Serr)/length(y))
  }
  return(SMem)
}
plWav1j<-function(dat, filter.number = 4, family = "DaubLeAsymm", a = 10, b = 10, bet = 1,
                  nu0 = 5,lamb0 = 10, M = 10, Ne, j0 = nlevelsWT(vw), plot.EMPL = FALSE){
  require(wavethresh)
  set.seed(1321)
  vw<-wd(dat,filter.number=filter.number,family=family)
  sigma<-mad(accessD(vw,level=(nlevelsWT(vw)-1)))/0.6745
  pr1a<-NULL
  for (j in 0:(j0-1)){
    coefthr<-accessD(vw,level=j)
    pr1a[(2^j:(2^(j+1)-1))]<-rev(coefthr)
  }
  a<-a
  b<-b*sigma^2
  Oj<-NULL
  tau12<-1/rgamma(1,a+3,b+10)
  Psi<-NULL
  pr<-NULL
  for (j in 0:(j0-1)){
    if(j==0){
      pr[1]<-0.999999
    }else{
      pr[j+1]<-2^(-j*bet)
    }
    coefthr<-accessD(vw,level=j)
    Oj[(2^j:(2^(j+1)-1))]<-(sqrt(sigma^2+tau12)/sigma)*((1-pr[j+1])/pr[j+1])*exp(-0.5*(rev(coefthr)^2/sigma^2)*(tau12/(tau12+sigma^2)))
    Psi<-1/(1+Oj)
  }
  pi0<-c0<-NULL
  for(j in 0:(j0-1)){
    pi0[j+1]<-sum(Psi[(2^j:(2^(j+1)-1))])/length(Psi[(2^j:(2^(j+1)-1))])
    c0[j+1]<-max(0,sum(Psi[(2^j:(2^(j+1)-1))]*pr1a[(2^j:(2^(j+1)-1))]^2)/(sigma^2*sum(Psi[(2^j:(2^(j+1)-1))]))-1)
  }
  d<-dx<-matrix(0,ncol=length(vw$D),nrow=Ne)
  nu0<-nu0
  lamb0<-lamb0
  S<-cbind(rep(nu0/2,Ne),rep(lamb0*nu0/2,Ne))
  sig1<-1/rgamma(Ne,nu0/2,lamb0*nu0/2)
  #prF<-pi0
  #stF<-c0
  Psi1F<-matrix(0,ncol=length(vw$D),nrow = Ne)
  Oj1F<-matrix(0,ncol=length(vw$D),nrow = Ne)
  pr1b<-NULL
  M<-M
  prF<-matrix(0,nrow=M,ncol=j0)
  stF<-matrix(0,nrow=M,ncol=j0)
  prF[1,]<-pi0
  stF[1,]<-c0
  Psi1p<-NULL
  sig1p<-NULL
  for(m in 2:M){
    for(j in 0:(j0-1)){
      coefthr<-accessD(vw,level=j)
      pr1b[(2^j:(2^(j+1)-1))]<-rev(coefthr)
      for(i in (2^j:(2^(j+1)-1))){
        #PL-step
        Oj1F[,i]<-(sqrt(1+stF[m-1,j+1])/1)*((1-prF[m-1,j+1])/prF[m-1,j+1])*exp(-0.5*(pr1b[i]^2/sig1)*(1/(1+stF[m-1,j+1]^(-1))))
        Psi1F[,i]<-1/(1+Oj1F[,i])
        d[,i]<-rnorm(Ne,Psi1F[,i]*stF[m-1,j+1]/(1+stF[m-1,j+1])*pr1b[i],sd=sqrt(Psi1F[,i]*sig1*stF[m-1,j+1]/(1+stF[m-1,j+1])))
        w<-dnorm(pr1b[i],d[,i],sd=sqrt(sig1))
        k<-sample(1:Ne,size=Ne,replace=TRUE,prob=w)
        S[,1]<-S[k,1]+1/2
        S[,2]<-S[k,2]+(pr1b[i]-d[,i])^2/2  # pag 112 UPCcourse-handouts.pdf
        sig1<-1/rgamma(Ne,S[,1],S[,2])
        d[,i]<-rnorm(Ne,Psi1F[,i]*stF[m-1,j+1]/(1+stF[m-1,j+1])*pr1b[i],sd=sqrt(Psi1F[,i]*sig1*stF[m-1,j+1]/(1+stF[m-1,j+1])))
        #dx[,i]<-sample(d[,i],size=Ne,replace=TRUE,prob=w)
        Psi1p[i]<-sample(Psi1F[,i],size=1,replace=FALSE,prob=rep(1/Ne,Ne)) 
        #Psi1p[i]<-sum(Psi1F[,i]/Ne)
        sig1p<-mean(sig1)
      }
      #E-step
      prF[m,j+1]<-sum(Psi1p[(2^j:(2^(j+1)-1))])/length(Psi1p[(2^j:(2^(j+1)-1))])
      stF[m,j+1]<-max(0,sum(Psi1p[(2^j:(2^(j+1)-1))]*pr1a[(2^j:(2^(j+1)-1))]^2)/(sig1p*sum(Psi1p[(2^j:(2^(j+1)-1))]))-1)
    }
  }
  bayes3<-vw
  Ex5F1<-NULL
  for(i in 1:length(pr1b)){
    Ex5F1[i]<-mean(d[,i])}
  for (j in 0:(j0-1)){
    bayes3<-putD(bayes3,level=j,v=rev(Ex5F1[(2^j:(2^(j+1)-1))]))
  }
  bayesrec5<-wr(bayes3)
  D<-accessC(vw,level=0)
  if (plot.EMPL == TRUE) {
    x <- seq(1, length(dat))/length(dat)
    par(mfrow = c(1, 2))
    plot(x, dat, type = "l", ylab = "(a) Datos")
    plot(x, bayesrec5, type = "l", ylab = "Bayes_Shrink", ylim = c(min(dat), max(dat)))
  }
  return(bayesrec5=bayesrec5)
  #lisc<-list(bayesrec5=bayesrec5,pr1a=pr1a,Ex5F1=Ex5F1,Psi1p=Psi1p,prF=prF,stF=stF,D=D)
  return(lisc)
}
