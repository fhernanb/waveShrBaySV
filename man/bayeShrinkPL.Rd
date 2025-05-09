% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayeShrinkPL.R
\name{bayeShrinkPL}
\alias{bayeShrinkPL}
\title{bayeShrinkPL}
\usage{
bayeShrinkPL(
  dat,
  filter.number = 4,
  family = "DaubLeAsymm",
  M = 10,
  Ne = 100,
  j0 = nlevelsWT(vw),
  plot.bayeShrinkPL = FALSE
)
}
\arguments{
\item{dat}{Serie de observaciones a ingresar a las que se le debe eliminar el ruido aditivo.}

\item{filter.number}{Parámetro de la transformación wavelet que indica en número de momentos de desvanecimiento.}

\item{family}{Familia de la transformación wavelet c('DaubExPhase','DaubLeAsymm','Coiflets',...).}

\item{M}{Número de repeticiones en el proceso de maximización del algoritmo.}

\item{Ne}{Número de partículas el proceso secuencial Monte Carlo.}

\item{j0}{Nivel de resolución de la transformación wavelet.}

\item{plot.bayeShrinkPL}{Gráfica de comparación entre la serie de observaciones y serie libre de ruido.}
}
\value{
Esta funcion retorna la serie de observaciones libre de ruido al nivel de resolución especificado.
}
\description{
Esta funcion elimina el ruido aditivo de una serie de observaciones por medio de un algoritmo bayesiano basado en el aprendizaje de partículas.
La función utiliza un método shrinkage bayesiano basado en particle learning para suavizar el valor de los
coeficientes de la transformación wavelet de las observaciones.
}
\examples{
# Ejemplo de eliminación de ruido

library(wavethresh)

y <- DJ.EX()$doppler

s <- sd(y)
SNR <- 9
e <- rnorm(length(y),mean = 0, sd = s/SNR)
YNoi <- y + e

YDNoi<-bayeShrinkPL(YNoi)

plot.ts(y,ylab='Datos originales')
plot.ts(YNoi,ylab='Datos con ruido aditivo')
plot.ts(YDNoi,ylab='BayeShrinkPL')
}
\author{
Omar Rios Saavedra, \email{orioss@unal.edu.co}
}
