
# waveShrBaySV <img src="docs/logo.png" align="right" alt="" width="120" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/fhernanb/waveShrBaySV/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fhernanb/waveShrBaySV/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

El objetivo de waveShrBaySV es eliminar el ruido aditivo en series de tiempo con error estocástico 
endógeno a través de la técnica de contracción (shrinkage) de los coeficientes de la transformación 
wavelet. Así mismo, se implementa en el modelo de volatilidad estocástica la metodología de 
eliminación de ruido a partir de un proceso de filtro de partículas con empuje basado en wavelets.

## Installation

You can install the development version of waveShrBaySV like so:

``` r
# install.packages("pak")
pak::pak("fhernanb/waveShrBaySV")
```

You can visit the [package
website](https://fhernanb.github.io/waveShrBaySV/) to explore the vignettes
(articles) and function reference.
