#' bayeShrinkPL
#'
#' Esta funcion elimina el ruido aditivo de una serie de observaciones por medio de un algoritmo bayesiano basado en el aprendizaje de partículas.
#' La función utiliza un método shrinkage bayesiano basado en particle learning para suavizar el valor de los
#' coeficientes de la transformación wavelet de las observaciones.
#'
#' @param dat Serie de observaciones a ingresar a las que se le debe eliminar el ruido aditivo.
#' @param filter.number Parámetro de la transformación wavelet que indica en número de momentos de desvanecimiento.
#' @param family Familia de la transformación wavelet c('DaubExPhase','DaubLeAsymm','Coiflets',...).
#' @param M  Número de repeticiones en el proceso de maximización del algoritmo.
#' @param Ne  Número de partículas el proceso secuencial Monte Carlo.
#' @param j0  Nivel de resolución de la transformación wavelet.
#' @param plot.bayeShrinkPL Gráfica de comparación entre la serie de observaciones y serie libre de ruido.
#'
#' @return Esta funcion retorna la serie de observaciones libre de ruido al nivel de resolución especificado.
#'
#' @example examples/examples_bayeShrinkPL.R
#'
#' @author Omar Rios Saavedra, \email{orioss@unal.edu.co}
#'
#' @export
#' @importFrom stats mad rgamma
#' @importFrom graphics par
#' @importFrom wavethresh putD accessD wr accessC nlevelsWT
bayeShrinkPL <- function(dat, filter.number = 4, family = "DaubLeAsymm", M = 10, Ne = 100, j0 = nlevelsWT(vw), plot.bayeShrinkPL = FALSE){
  #a = 10; b = 10; bet = 1; nu0 = 5; lamb0 = 10
  #set.seed(1321)
  vw <- wavethresh::wd(dat,filter.number=filter.number,family=family)
  sigma<-mad(wavethresh::accessD(vw,level=(nlevelsWT(vw)-1)))/0.6745
  pr1a<-NULL
  for (j in 0:(j0-1)){
    coefthr<-accessD(vw,level=j)
    pr1a[(2^j:(2^(j+1)-1))]<-rev(coefthr)
  }
  a<-10
  b<-10*sigma^2
  Oj<-NULL
  tau12<-1/rgamma(1,a+3,b+10)
  Psi<-NULL
  pr<-NULL
  for (j in 0:(j0-1)){
    if(j==0){
      pr[1]<-0.999999
    }else{
      pr[j+1]<-2^(-j*1)
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
  nu0<-5
  lamb0<-10
  S<-cbind(rep(nu0/2,Ne),rep(lamb0*nu0/2,Ne))
  sig1<-1/rgamma(Ne,nu0/2,lamb0*nu0/2)
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
        Psi1p[i]<-sample(Psi1F[,i],size=1,replace=FALSE,prob=rep(1/Ne,Ne))
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
    bayes3<-wavethresh::putD(bayes3,level=j,v=rev(Ex5F1[(2^j:(2^(j+1)-1))]))
  }
  bayesrec5<-wavethresh::wr(bayes3)
  D<-wavethresh::accessC(vw,level=0)
  if (plot.bayeShrinkPL == TRUE) {
    x <- seq(1, length(dat))/length(dat)
    par(mfrow = c(2, 1))
    plot(x, dat, type = "l", ylab = "Datos con ruido")
    plot(x, bayesrec5, type = "l", ylab = "Bayes_Shrink", ylim = c(min(dat), max(dat)))
  }
  return(bayesrec5=bayesrec5)
}
