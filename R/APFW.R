#' APFW
#'
#' Esta funcion calcula el modelo de volatilidad estocástica basada en el filtro auxiliar de partículas (APF).
#' El algoritmo incorpora los pasos de empuje bayesianos basados en la transformación wavelet.
#' El algoritmo asume que los parámetros del modelo son desconocidos.
#'
#' @param y  representa la serie de observaciones reales.
#' @param ms1  vector de parámetros (\eqn{\alpha}, \eqn{\beta} y \eqn{\tau^2}).
#' @param xs  partículas iniciales de la variable latente a partir de la distribución a priori.
#' @param levj  nivel de resolución en la transformación wavelet.
#' @param M  parámetro de la función bayeShrinkPL.
#' @param Ne  parámetro de la función bayeShrinkPL.
#' @param method  1 o 2, método de eliminación de ruido a partir de la transformación wavelet method = 1 (bayeShrinkPL), method = 2 (BAYES.THR).
#'
#' @return Esta funcion retorna los cuantiles (2.5%, 50% y 97.5%) de las estimaciones de la volatilidad estocástica
#' con el método de filtrado de partículas basado en wavelets. Además también retorna la función de verosimilitud
#'
#' @example examples/examples_APFW.R
#'
#' @author Omar Rios Saavedra, \email{orioss@unal.edu.co}
#'
#' @export
#' @importFrom stats dnorm plot.ts quantile var
APFW<-function(y,ms1,xs,levj,M=5,Ne=5,method=1){
  n = length(y)
  N = length(xs)
  like = rep(0,n)
  like2 = rep(0,n)
  quants = array(0,c(n,1,3))
  xss<-NULL
  ws<-NULL
  #ESS<-NULL
  #par(mfrow=c(1,1))
  for (t in 1:n){
    like[t] = sum(dnorm(y[t],0.0,exp(xs/2)))
    weight  = dnorm(y[t],0.0,exp(xs/2),log=TRUE)
    weight1 = exp(weight-max(weight))
    k       = sample(1:N,size=N,replace=T,prob=weight1)
    xt      = rnorm(N,ms1[1]+ms1[2]*xs[k],exp(ms1[3]/2))
    if(method==1){
      xst  <- bayeShrinkPL(xt,filter.number = 4, family = 'DaubLeAsymm', M = M, Ne = Ne, j0 = levj, plot.bayeShrinkPL = FALSE)}
    else{
      xst  <- BAYES.THR(xt,alpha=1,beta=1,filter.number = 4, family = 'DaubLeAsymm', plotfn = FALSE, j0 = levj, dev = var)
    }
    like2[t] = sum(dnorm(y[t],0.0,exp(xst/2)))
    if(like2[t]>like[t]){
    xt <- xst#}else{
    xt <- xt}
    xt <- xst
    w    = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=TRUE,prob=w)
    xs   = xt[ind]
    xss<-rbind(xss,xs)
    ws<-rbind(ws,w)
    #cv2<-var(w)/(mean(w)^2)
    #ESS<-c(ESS,N/(1+cv2))
    #ts.plot(ESS,xlim=c(1,n))
    quants[t,1,] = quantile(exp(xs/2),c(0.5,0.025,0.975))
  }
  return(list(quants=quants,xss=xss))
}
