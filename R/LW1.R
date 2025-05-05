#' LW1
#'
#' Función de filtro de partículas con aprendizaje de parámetros enfoque Liu & West
#' Basada en el algoritmo de Lopes & Tsay (2011)
#'
#' @param y representa la serie de observaciones reales.
#' @param alphas representa los valores iniciales para el parámetro de reversión de la media en el proceso de volatilidad estocástica.
#' @param betas representa los valores iniciales para el parámetros de persistencia de volatilidad.
#' @param tau2s representa los valorese iniciales para la varianza de la variable latente (volatilidad estocástica).
#' @param xs partículas iniciales de la variable latente a partir de la distribución a priori.
#' @param delta constante de ponderación para el aprendizaje de parámetros en el algoritmo Liu & West (2001).
#'
#' @return Esta funcion retorna los cuantiles (2.5%, 50% y 97.5%) de
#' las estimaciones de la volatilidad estocástica,
#' y sus parámetros (\eqn{\alpha}, \eqn{\beta} y \eqn{\tau^2}).
#'
#' @example examples/examples_LW1.R
#'
#' @author Omar Rios Saavedra, \email{orioss@unal.edu.co}
#'
#' @export
#' @importFrom stats dnorm quantile ts.plot var
LW1 <- function(y, alphas, betas, tau2s, xs, delta) {
  n <- length(y)
  N <- length(xs)
  quants <- array(0,c(n,4,3))
  parss<-array(0,c(N,3,n))
  h2 <- 1-((3*delta-1)/(2*delta))^2
  a  <- sqrt(1-h2)
  pars <- cbind(alphas,betas,log(tau2s))
  xss<-NULL
  ws<-NULL
  #ESS<-NULL
  #like <- rep(0,n)
  #par(mfrow=c(1,1))
  for (t in 1:n){
    #like[t] <- sum(dnorm(y[t],0.0,exp(xs/2)))
    # Resampling
    mus     <- pars[,1]+pars[,2]*xs
    mpar    <- apply(pars,2,mean)
    vpar    <- var(pars)
    ms      <- a*pars+(1-a)*matrix(mpar,N,3,byrow=T)
    weight  <- dnorm(y[t],0.0,exp(mus/2),log=TRUE)
    weight1 <- exp(weight-max(weight))
    k       <- sample(1:N,size=N,replace=T,prob=weight1)
    # Propagating
    ms1 <- ms[k,]+matrix(rnorm(3*N),N,3)%*%chol(h2*vpar)
    xt   <- rnorm(N,ms1[,1]+ms1[,2]*xs[k],exp(ms1[,3]/2))
    w    <- dnorm(y[t],0.0,exp(xt/2),log=TRUE)-weight[k]
    w    <- exp(w-max(w))
    ind  <- sample(1:N,size=N,replace=T,prob=w)
    xs   <- xt[ind]
    pars <- ms1[ind,]
    xss<-rbind(xss,xs)
    parss[,,t]<-pars
    #ws<-rbind(ws,w)
    #cv2<-var(w)/(mean(w)^2)
    #ESS<-c(ESS,N/(1+cv2))
    #ts.plot(ESS,xlim=c(1,n))
    # Storing quantiles
    quants[t,1,] <- quantile(pars[,1],c(0.5,0.025,0.975))
    quants[t,2,] <- quantile(pars[,2],c(0.5,0.025,0.975))
    quants[t,3,] <- quantile(exp(pars[,3]),c(0.5,0.025,0.975))
    quants[t,4,] <- quantile(exp(xs/2),c(0.5,0.025,0.975))
  }
  return(list(quants=quants,parss=parss,xss=xss))
}
