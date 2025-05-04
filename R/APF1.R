#' APF1
#'
#' Esta funcion calcula la volatilidad latente de un conjunto de observaciones utilizando
#' el filtro auxiliar de partículas (APF). El algoritmo asume que los parámetros del modelo son conocidos.
#'
#' @param y representa la serie de observaciones reales.
#' @param ms1 vector de parámetros (\eqn{\alpha}, \eqn{\beta} y \eqn{\tau^2}).
#' @param xs partículas iniciales de la variable latente a partir de la distribución a priori.
#'
#' @return Esta funcion retorna los cuantiles (2.5%, 50% y 97.5%) de las estimaciones de la volatilidad estocástica.
#' Además también retorna la función de verosimilitud
#'
#' @example examples/examples_APF1.R
#'
#' @author Omar Rios Saavedra, \email{orioss@unal.edu.co}
#'
#' @export
#' @importFrom stats dnorm quantile ts.plot var
APF1<-function(y,ms1,xs){
  n = length(y)
  N = length(xs)
  quants = array(0,c(n,1,3))
  xss<-NULL
  ws<-NULL
  #ESS<-NULL
  like = rep(0,n)
  #par(mfrow=c(1,1))
  for (t in 1:n){
    like[t] = sum(dnorm(y[t],0.0,exp(xs/2)))
    weight  = dnorm(y[t],0.0,exp(xs/2),log=TRUE)
    weight1 = exp(weight-max(weight))
    k       = sample(1:N,size=N,replace=T,prob=weight1)
    xt      = rnorm(N,ms1[1]+ms1[2]*xs[k],exp(ms1[3]/2))
    w    = dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    w    = exp(w-max(w))
    ind  = sample(1:N,size=N,replace=TRUE,prob=w)
    xs   = xt[ind]
    xss<-rbind(xss,xs)
    #ws<-rbind(ws,w)
    #cv2<-var(w)/(mean(w)^2)
    #ESS<-c(ESS,N/(1+cv2))
    #ts.plot(ESS,xlim=c(1,n))
    quants[t,1,] = quantile(exp(xs/2),c(0.5,0.025,0.975))
  }
  return(list(quants=quants,xss=xss,like=like))
}
