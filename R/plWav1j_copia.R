#' plWav1j_copia
#'
#' Esta funcion calcula bla bla bla.
#'
#' @param dat bla bla bla, Omar debe completar esto.
#' @param filter.number bla bla bla, Omar debe completar esto.
#' @param family bla bla bla, Omar debe completar esto.
#' @param a bla bla bla, Omar debe completar esto.
#' @param b bla bla bla, Omar debe completar esto.
#' @param bet bla bla bla, Omar debe completar esto.
#' @param nu0  bla bla bla, Omar debe completar esto.
#' @param lamb0  bla bla bla, Omar debe completar esto.
#' @param M  bla bla bla, Omar debe completar esto.
#' @param Ne  bla bla bla, Omar debe completar esto.
#' @param j0  bla bla bla, Omar debe completar esto.
#' @param plot.EMPL  bla bla bla, Omar debe completar esto.
#'
#' @return Esta funcion retorna bla bla bla, Omar debe completar esto..
#'
#' @example examples/examples_plWav1j_copia.R
#'
#' @author Omar Rios Saavedra, \email{orioss@unal.edu.co}
#'
#' @export
#' @importFrom stats mad rgamma
#' @importFrom graphics par
plWav1j_copia <- function(dat, filter.number = 4, family = "DaubLeAsymm", a = 10, b = 10, bet = 1,
                  nu0 = 5,lamb0 = 10, M = 10, Ne, j0 = nlevelsWT(vw), plot.EMPL = FALSE){
  #require(wavethresh) Omar, no use esta notacion
  set.seed(1321)
  vw <- wavethresh::wd(dat,filter.number=filter.number,family=family)
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
  lisc<-list(bayesrec5=bayesrec5,pr1a=pr1a,Ex5F1=Ex5F1,Psi1p=Psi1p,prF=prF,stF=stF,D=D)
  return(lisc)
}
