#' errorF
#'
#' Esta funcion se encarga del cálculo de la raíz cuadrada del error cuadrático medio.
#' La función facilita la comparación entre varios errores de un mismo procedimiento.
#'
#' @param y representa la serie de observaciones reales.
#' @param vec representa la matriz de observaciones de contraste. Cada columna representa una variable.
#'
#' @return Esta funcion retorna una matriz de la raiz del error cuadrático medio para vectores de contraste ingresados en vec.
#'
#' @example examples/examples_errorF.R
#'
#' @author Omar Rios Saavedra, \email{orioss@unal.edu.co}
#'
#' @export
errorF <- function(y, vec) {
  SMem <- matrix(NA, nrow = ncol(vec), ncol = 1)
  rownames(SMem) <- names(vec)
  colnames(SMem) <- "RMSE"
  Serr <- NULL
  for (j in 1:ncol(vec)) {
    for (i in 1:length(y)) {
      Serr[i] <- (vec[i, j] - y[i])^2
    }
    SMem[j, ] <- sqrt(sum(Serr)/length(y))
  }
  return(SMem)
}
