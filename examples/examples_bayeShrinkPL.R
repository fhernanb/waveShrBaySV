# Ejemplo de eliminación de ruido

y <- wavethresh::DJ.EX()$blocks
s <- sd(y)
SNR <- 7
eN <- rnorm(length(y),mean = 0, sd = s/SNR)

YNoi<-y + eN

bayeShrinkPL(YNoi, plot.bayeShrinkPL = TRUE)
