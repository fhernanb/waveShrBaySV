rlike<-function(x){rnorm(1,0,exp(x/2))}
n     =  2^10
alpha = -0.00645
beta  =  0.985
tau2  =  0.1
tau   = sqrt(tau2)
y     = rep(0,n)
x     = rep(0,n)
x[1]  = alpha/(1-beta)
y[1]  = rlike(x[1])
set.seed(116)
for (t in 2:n){
  x[t] = rnorm(1,alpha+beta*x[t-1],tau)
  y[t] = rlike(x[t])
}
par(mfrow=c(2,1));plot.ts(y);plot.ts(exp(x/2))
