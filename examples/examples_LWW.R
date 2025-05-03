# Ejemplo de aplicación del filtro de Liu West con empuje basado en wavelet
\dontrun{
library(wavethresh)

rlike <- function(x){rnorm(1,0,exp(x/2))}
n     =  2^10
alpha =  -0.004
beta  =  0.985
tau2  =  0.1
tau   = sqrt(tau2)
y1     = rep(0,n)
x1     = rep(0,n)
x1[1]  = alpha/(1-beta)
y1[1]  = rlike(x1[1])
set.seed(116)
for (t in 2:n){
  x1[t] = rnorm(1,alpha+beta*x1[t-1],tau)
  y1[t] = rlike(x1[t])
}

alpha.true<-alpha
beta.true<-beta
tau2.true<-tau2

N=2^11
delta<-0.975
lev<-3:6

m0<-0.0;C0<-0.1;sC0<-sqrt(C0)
ealpha<-alpha;valpha<-0.01
ephi<-beta;vphi<-0.01
nu<-3;lambda<-tau2

xs<-rnorm(N,m0,sC0)

alphas<-rnorm(N,ealpha,sqrt(valpha))
betas<-rnorm(N,ephi,sqrt(vphi))
tau2s<-1/rgamma(N,nu/2,nu*lambda/2)
rest2<-LWW(y1,alphas,betas,tau2s,xs,delta,lev[1],M=2,method=1)
rest2
}
