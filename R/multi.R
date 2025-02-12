#' multi
#'
#' Esta funcion sirve para obtener la multiplicacion de dos numeros reales.
#'
#' @param x A number.
#' @param y A number.
#'
#' @return Esta funcion retorna un numero que corresponde a la sum of \code{x} and \code{y}.
#'
#' @example examples/examples_multi.R
#'
#' @details
#' Esta funcion sirva para bla bla bla
#' en este parrafo Omar debe dar todos los detalles tecnicos de la funcion
#' aqui es donde explica lo que considere que se necesita.
#'
#' @author Valentina Hurtado Sepúlveda, \email{vhurtados@unal.edu.co}
#'
#' @export
multi <- function(x, y) {
  res <- myaux(x, y)
  res
}
#' @importFrom stats rnorm
myaux <- function(x, y) {
  muestra <- rnorm(n=10)
  x * y + mean(muestra)
}
