# Ejemplo filtro de partículas con parámetros desconocidos (Liu & West (2001))

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

alpha.true <- alpha
beta.true <- beta
tau2.true <- tau2

# Número de partículas y ponderación de la mixtura
N = 2^11
delta <- 0.975

# Valores iniciales y distribuciones a priori
m0 <- 0.0; C0 <- 0.1; sC0 <- sqrt(C0)
ealpha <- alpha; valpha <- 0.1
ephi <- beta; vphi <- 0.1
nu <- 4; lambda <- tau2

xs <- rnorm(N, m0, sC0)
alphas <- rnorm(N, ealpha, sqrt(valpha))
betas <- rnorm(N, ephi, sqrt(vphi))
tau2s <- 1/rgamma(N, nu/2, nu*lambda/2)

# Estimación de la volatilidad latente y parámetros
rest1 <- LW1(y1, alphas, betas, tau2s, xs, delta)
mvolp1 <- rest1$quants[,4,1]

# Volatilidad latente y estimación
plot.ts(y1, ylab = 'Retornos')
plot.ts(exp(x1/2), ylim = c(-0.01, 4), ylab = 'Volatilidad latente', col = 'red')
plot.ts(mvolp1, ylim = c(-0.01, 4), ylab = 'Volatilidad estimada')
lines(exp(x1/2), col = 'red')
legend('topleft', legend = c('Volatilidad estimada', 'Volatilidad latente'),
       lty = 1, col = c('black', 'red'), bty = 'n')

# Parámetro de la media
malpha1 <- rest1$quants[,1,1]
lalpha1 <- rest1$quants[,1,2]
ualpha1 <- rest1$quants[,1,3]
ts.plot(malpha1, ylim = range(lalpha1, ualpha1), main = 'LW', ylab = expression(mu))
lines(lalpha1, lwd = 1.25, col = 'gray')
lines(ualpha1, lwd = 1.25, col = 'gray')
abline(h = alpha.true, col = 2, lwd = 1)

# Parámetro de persistencia de la volatilidad
mbeta1 <- rest1$quants[,2,1]
lbeta1 <- rest1$quants[,2,2]
ubeta1 <- rest1$quants[,2,3]
ts.plot(mbeta1, ylim = range(lbeta1,ubeta1), main = 'LW', ylab = expression(phi))
lines(lbeta1, lwd = 1.25, col = 'gray')
lines(ubeta1, lwd = 1.25, col = 'gray')
abline(h=beta.true, col = 2, lwd = 1)

# Parámetro de la varianza de la volatilidad latente
mtau21 <- rest1$quants[,3,1]
ltau21 <- rest1$quants[,3,2]
utau21 <- rest1$quants[,3,3]
ts.plot(mtau21, ylim = range(ltau21, utau21), main = 'LW', ylab = expression(tau^2))
lines(ltau21, lwd = 1.25, col = 'gray')
lines(utau21, lwd = 1.25, col = 'gray')
abline(h=tau2.true, col = 2, lwd = 1)

