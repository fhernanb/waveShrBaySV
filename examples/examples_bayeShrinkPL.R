# Ejemplo de eliminación de ruido

library(wavethresh)

y <- DJ.EX()$doppler

s <- sd(y)
SNR <- 7
e <- rnorm(length(y),mean = 0, sd = s/SNR)
YNoi <- y + e

YDNoi<-bayeShrinkPL(YNoi)

plot.ts(y,ylab='Datos originales')
plot.ts(YNoi,ylab='Datos con ruido aditivo')
plot.ts(YDNoi,ylab='BayeShrink')
