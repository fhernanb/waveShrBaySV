# Example for two-dimensions
x <- c(1, 6)
v <- matrix(c(5, 3,
              1, 4), ncol=2, byrow=TRUE)
errorF(y=x, vec=v)

# Example for three-dimensions
x <- c(1, 6, 5)
v <- matrix(c(5, 3, 1,
              1, 4, 3,
              0, 2, 2), ncol=3, byrow=TRUE)
errorF(y=x, vec=v)

