#Ejemplo de estimación de volatilidad estocástica con parámetros conocidos

rlike<-function(x){rnorm(1,0,exp(x/2))}
n     =  2^11
alpha = -0.00645
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

#dev.off()
alpha.true<-alpha
beta.true<-beta
tau2.true<-tau2

theta<-c(alpha,beta,tau2)
m0<-0.0;C0<-0.1;sC0<-sqrt(C0)
N=2^10
set.seed(123);xs<-rnorm(N,m0,sC0)

vol1<-APF1(y1,theta,xs)
mvol1<-vol1$quants[,1,1]
plot.ts(y1,ylab='Retornos')
plot.ts(exp(x1/2),ylab='Volatilidad latente')
plot.ts(mvol1,ylab='Volatilidad estimada')
lines(exp(x1/2),ylim=c(-0.01,6),col='red')
