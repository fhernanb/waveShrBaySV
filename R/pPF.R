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

#sigma2 <- 0.1

#m0<-0.0;C0<-0.1;sC0<-sqrt(C0)
#ealpha<-0;valpha<-1#10
#ephi<-0;vphi<-1
#nu<-3;lambda<-sigma2 #tau2 # 

xs<-rnorm(N,m0,sC0)
#xst  <- plWav1j(xs2,filter.number = 4, family = 'DaubLeAsymm', a = 1, b = 1, bet = 1,
#                nu0 = 5,lamb0 = 1, M=75, Ne = 200, j0=lev[4], plot.EMPL = FALSE)
#xs <- xst
alphas<-rnorm(N,ealpha,sqrt(valpha))
betas<-rnorm(N,ephi,sqrt(vphi))
tau2s<-1/rgamma(N,nu/2,nu*lambda/2)

set.seed(8642);rest1<-LW1(y,alphas,betas,tau2s,xs,delta)
set.seed(8642);rest2<-LWW(y,alphas,betas,tau2s,xs,delta,lev[1],M=15,method=1)
set.seed(8642);rest3<-LWW(y,alphas,betas,tau2s,xs,delta,lev[4],M=15,method=1) #lev[4]=6
set.seed(8642);rest4<-LWW(y,alphas,betas,tau2s,xs,delta,lev[3],M=40,method=1)
set.seed(8642);rest5<-LWW(y,alphas,betas,tau2s,xs,delta,lev[4],M=40,method=1)
set.seed(8642);rest6<-LWW(y,alphas,betas,tau2s,xs,delta,10,method = 2)


mvolp1<-rest1$quants[,4,1]
mvolp2<-rest2$quants[,4,1]
mvolp3<-rest3$quants[,4,1]
mvolp4<-rest4$quants[,4,1]
mvolp5<-rest5$quants[,4,1]
mvolp6<-rest6$quants[,4,1]

vec42<-data.frame(mvolp1,mvolp2,mvolp3,mvolp4,mvolp5,mvolp6)
errorF(exp(x/2),vec42)

x11()
par(mfrow=c(3,2))
for(i in 1:6){
  plot.ts(vec42[,i],xlab="time",ylab="",main="",lwd=2)
  lines(exp(x/2),ylim=c(-0.01,6),col='red')
}

malpha1<-rest1$quants[,1,1]
lalpha1<-rest1$quants[,1,2]
ualpha1<-rest1$quants[,1,3]
malpha2<-rest2$quants[,1,1]
lalpha2<-rest2$quants[,1,2]
ualpha2<-rest2$quants[,1,3]
malpha3<-rest3$quants[,1,1]
lalpha3<-rest3$quants[,1,2]
ualpha3<-rest3$quants[,1,3]
malpha4<-rest4$quants[,1,1]
lalpha4<-rest4$quants[,1,2]
ualpha4<-rest4$quants[,1,3]
malpha5<-rest5$quants[,1,1]
lalpha5<-rest5$quants[,1,2]
ualpha5<-rest5$quants[,1,3]
malpha6<-rest6$quants[,1,1]
lalpha6<-rest6$quants[,1,2]
ualpha6<-rest6$quants[,1,3]

vecAl<-data.frame(malpha1,malpha2,malpha3,malpha4,malpha5,malpha6)
errorF(alpha.true,vecAl)

par(mfrow=c(1,6))
ts.plot(malpha1,ylim=range(lalpha1,ualpha1),ylab="",main=expression(alpha))
lines(lalpha1,lwd=1.25)
lines(ualpha1,lwd=1.25)
abline(h=alpha.true,col=2,lwd=2)
ts.plot(malpha2,ylim=range(lalpha2,ualpha2),ylab="",main=expression(alpha))
lines(lalpha2,lwd=1.25)
lines(ualpha2,lwd=1.25)
abline(h=alpha.true,col=2,lwd=2)
ts.plot(malpha3,ylim=range(lalpha3,ualpha3),ylab="",main=expression(alpha))
lines(lalpha3,lwd=1.25)
lines(ualpha3,lwd=1.25)
abline(h=alpha.true,col=2,lwd=2)
ts.plot(malpha4,ylim=range(lalpha4,ualpha4),ylab="",main=expression(alpha))
lines(lalpha4,lwd=1.25)
lines(ualpha4,lwd=1.25)
abline(h=alpha.true,col=2,lwd=2)
ts.plot(malpha5,ylim=range(lalpha5,ualpha5),ylab="",main=expression(alpha))
lines(lalpha5,lwd=1.25)
lines(ualpha5,lwd=1.25)
abline(h=alpha.true,col=2,lwd=2)
ts.plot(malpha6,ylim=range(lalpha6,ualpha6),ylab="",main=expression(alpha))
lines(lalpha6,lwd=1.25)
lines(ualpha6,lwd=1.25)
abline(h=alpha.true,col=2,lwd=2)

mbeta1<-rest1$quants[,2,1];mean(mbeta1)
lbeta1<-rest1$quants[,2,2]
ubeta1<-rest1$quants[,2,3]
mbeta2<-rest2$quants[,2,1]
lbeta2<-rest2$quants[,2,2]
ubeta2<-rest2$quants[,2,3]
mbeta3<-rest3$quants[,2,1]
lbeta3<-rest3$quants[,2,2]
ubeta3<-rest3$quants[,2,3]
mbeta4<-rest4$quants[,2,1]
lbeta4<-rest4$quants[,2,2]
ubeta4<-rest4$quants[,2,3]
mbeta5<-rest5$quants[,2,1]
lbeta5<-rest5$quants[,2,2]
ubeta5<-rest5$quants[,2,3]
mbeta6<-rest6$quants[,2,1]
lbeta6<-rest6$quants[,2,2]
ubeta6<-rest6$quants[,2,3]

vecBt<-data.frame(mbeta1,mbeta2,mbeta3,mbeta4,mbeta5,mbeta6)
errorF(beta.true,vecBt)

par(mfrow=c(1,6))
ts.plot(mbeta1,ylim=range(lbeta1,ubeta1),ylab="",main=expression(beta))
lines(lbeta1,lwd=1.25,col='red')
lines(ubeta1,lwd=1.25,col='red')
abline(h=beta.true,col=2,lwd=2)
ts.plot(mbeta2,ylim=range(lbeta2,ubeta2),ylab="",main=expression(beta))
lines(lbeta2,lwd=1.25,col='red')
lines(ubeta2,lwd=1.25,col='red')
abline(h=beta.true,col=2,lwd=2)
ts.plot(mbeta3,ylim=range(lbeta3,ubeta3),ylab="",main=expression(beta))
lines(lbeta3,lwd=1.25,col='red')
lines(ubeta3,lwd=1.25,col='red')
abline(h=beta.true,col=2,lwd=2)
ts.plot(mbeta4,ylim=range(lbeta4,ubeta4),ylab="",main=expression(beta))
lines(lbeta4,lwd=1.25,col='red')
lines(ubeta4,lwd=1.25,col='red')
abline(h=beta.true,col=2,lwd=2)
ts.plot(mbeta5,ylim=range(lbeta5,ubeta5),ylab="",main=expression(beta))
lines(lbeta5,lwd=1.25,col='red')
lines(ubeta5,lwd=1.25,col='red')
abline(h=beta.true,col=2,lwd=2)
ts.plot(mbeta6,ylim=range(lbeta6,ubeta6),ylab="",main=expression(beta))
lines(lbeta6,lwd=1.25,col='red')
lines(ubeta6,lwd=1.25,col='red')
abline(h=beta.true,col=2,lwd=2)

mtau21<-rest1$quants[,3,1] # apply(exp(rest1$parss[,3,]),2,median) #
ltau21<-rest1$quants[,3,2] # apply(exp(rest1$parss[,3,]),2,quant025) # 
utau21<-rest1$quants[,3,3] # apply(exp(rest1$parss[,3,]),2,quant975) #
mtau22<-rest2$quants[,3,1]
ltau22<-rest2$quants[,3,2]
utau22<-rest2$quants[,3,3]
mtau23<-rest3$quants[,3,1]
ltau23<-rest3$quants[,3,2]
utau23<-rest3$quants[,3,3] 
mtau24<-rest4$quants[,3,1]
ltau24<-rest4$quants[,3,2]
utau24<-rest4$quants[,3,3] 
mtau25<-rest5$quants[,3,1]
ltau25<-rest5$quants[,3,2]
utau25<-rest5$quants[,3,3]
mtau26<-rest6$quants[,3,1]
ltau26<-rest6$quants[,3,2]
utau26<-rest6$quants[,3,3] 

vecTa<-data.frame(mtau21,mtau22,mtau23,mtau24,mtau25,mtau26)
errorF(tau2.true,vecTa)

par(mfrow=c(1,6))
ts.plot(mtau21,ylim=range(ltau21,utau21),ylab="",main=expression(tau^2))
lines(ltau21,lwd=1.25)
lines(utau21,lwd=1.25)
abline(h=tau2.true,col=2,lwd=2)
ts.plot(mtau22,ylim=range(ltau22,utau22),ylab="",main=expression(tau^2))
lines(ltau22,lwd=1.25)
lines(utau22,lwd=1.25)
abline(h=tau2.true,col=2,lwd=2)
ts.plot(mtau23,ylim=range(ltau23,utau23),ylab="",main=expression(tau^2))
lines(ltau23,lwd=1.25)
lines(utau23,lwd=1.25)
abline(h=tau2.true,col=2,lwd=2)
ts.plot(mtau24,ylim=range(ltau24,utau24),ylab="",main=expression(tau^2))
lines(ltau24,lwd=1.25)
lines(utau24,lwd=1.25)
abline(h=tau2.true,col=2,lwd=2)
ts.plot(mtau25,ylim=range(ltau25,utau25),ylab="",main=expression(tau^2))
lines(ltau25,lwd=1.25)
lines(utau25,lwd=1.25)
abline(h=tau2.true,col=2,lwd=2)
ts.plot(mtau26,ylim=range(ltau26,utau26),ylab="",main=expression(tau^2))
lines(ltau26,lwd=1.25)
lines(utau26,lwd=1.25)
abline(h=tau2.true,col=2,lwd=2)

saveRDS(rest1,file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest1.RData")
rest1<-readRDS(file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest1.RData")
saveRDS(rest2,file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest2.RData")
rest2<-readRDS(file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest2.RData")
saveRDS(rest3,file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest3.RData")
rest3<-readRDS(file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest3.RData")
saveRDS(rest4,file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest4.RData")
rest4<-readRDS(file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest4.RData")
saveRDS(rest5,file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest5.RData")
rest5<-readRDS(file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest5.RData")
saveRDS(rest6,file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest6.RData")
rest6<-readRDS(file="C:/Users/Omar/OneDrive/Favorites/Documents/vector/NuPF_Par/rest6.RData")



