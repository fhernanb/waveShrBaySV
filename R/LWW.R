#' LWW
#'
#' Esta funcion calcula el modelo de volatilidad estocástica basada en el filtro de partículas de Liu & West
#' el algoritmo incorpora los pasos de empuje bayesianos basados en la transformación wavelet.
#' plWav1j, metodología propuesta para la eliminación de ruido aditivo basado en la transformación wavelet.
#' BAYES.THR, metodología de  Abramovich et al. (1998) para la eliminación de ruido aditivo basado en la transformación wavelet.
#'
#' @param y representa la serie de observaciones reales.
#' @param alphas representa la matriz de observaciones de contraste.
#' @param betas bla bla bla, Omar debe completar esto.
#' @param tau2s bla bla bla, Omar debe completar esto.
#' @param xs bla bla bla, Omar debe completar esto.
#' @param delta bla bla bla, Omar debe completar esto.
#' @param lev  bla bla bla, Omar debe completar esto.
#' @param M  bla bla bla, Omar debe completar esto.
#' @param Ne  bla bla bla, Omar debe completar esto.
#' @param method  bla bla bla, Omar debe completar esto.
#'
#' @return Esta funcion retorna bla bla bla, Omar debe completar esto..
#'
#' @example examples/examples_LWW.R
#'
#' @author Omar Rios Saavedra, \email{orioss@unal.edu.co}
#'
#' @export
#' @importFrom wavethresh BAYES.THR
LWW = function(y,alphas,betas,tau2s,xs,delta,lev,M=75,Ne=20,method=1){
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
  like = rep(0,n)
  like2 = rep(0,n)
  #par(mfrow=c(1,1))
  for (t in 1:n){
    like[t] = sum(dnorm(y[t],0.0,exp(xs/2)))
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
                    nu0 = 50,lamb0 = 1, M = M, Ne = Ne, j0=lev, plot.EMPL = FALSE)}
    else{
      xst  <- BAYES.THR(xt,alpha=1,beta=1,filter.number = 4, family = 'DaubLeAsymm',plotfn=FALSE,j0=lev,dev=var)
    }
    like2[t] = sum(dnorm(y[t],0.0,exp(xst/2)))
    if(like2[t]>like[t]){
      xt <- xst}else{
      xt <- xt
      }
    #xt <- xst
    w    = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=T,prob=w)
    xs   = xt[ind]
    pars = ms1[ind,]
    xss<-rbind(xss,xs)
    parss[,,t]<-pars
    #ws<-rbind(ws,w)
    #cv2<-var(w)/(mean(w)^2)
    #ESS<-c(ESS,N/(1+cv2))
    #ts.plot(ESS,xlim=c(1,n))
    # Storing quantiles
    quants[t,1,] = quantile(pars[,1],c(0.5,0.025,0.975))
    quants[t,2,] = quantile(pars[,2],c(0.5,0.025,0.975))
    quants[t,3,] = quantile(exp(pars[,3]),c(0.5,0.025,0.975))
    quants[t,4,] = quantile(exp(xs/2),c(0.5,0.025,0.975))
  }
  return(list(quants=quants,parss=parss,xss=xss))
}

