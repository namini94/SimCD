#' Symmetric Kullback-Leibler divergence between samples drawn from random variables
#'
#' This function calculates the symmetric Kullback-Leibler divergence between two
#' empirical distriburions.
#'
#' @param x1 a matrix, where each row corresponds to draws from a random distribution
#' @param x2 a matrix, where each row corresponds to draws from a random distribution
#' and should have the same number of rows as \code{x1}
#' @return A vector with the length equal to the number of rows of input matrices,
#' @examples
#' x1 <- matrix(0, 2, 100)
#' x1[1,] <- rnorm(100, 0, 1)
#' x1[2,] <- rnorm(100, -1, 2)
#' x2 <- matrix(0, 2, 100)
#' x2[1,] <- rnorm(100, 1, 0)
#' x2[1,] <- rnorm(100, 1, 3)
#' kl <- KLsym(x1, x2)
#' @export
KLsym <- function(x1, x2)
  # Calculates KL divergence between the rows of two matrices x1 and x2
  # Output is a vector of distances with the same length as the number of rows of x1 and x2
{
  n <- nrow(x1)
  if (nrow(x2)!=n) stop("Number of rows should be the same.")
  nbins <- 100
  c <- 1.5
  KLdist <- NULL
  for (i in 1:n)
  {
    #find outliers
    q <- quantile(c(x1[i, ], x2[i, ]))
    Q1 <- q[2]
    Q3 <- q[4]
    IQD <- Q3 - Q1
    xmin <- max(Q1 - c*IQD,0)
    xmax <- Q3 + c*IQD
    idx <- which(x1[i,]>=xmin & x1[i,]<=xmax)
    f1 <- (hist(x1[i, idx], breaks = seq(xmin, xmax, length.out = nbins), plot = FALSE)$counts)/length(idx) + 1e-10
    idx <- which(x2[i,]>=xmin & x2[i,]<=xmax)
    f2 <- (hist(x2[i, idx], breaks = seq(xmin, xmax, length.out = nbins), plot = FALSE)$counts)/length(idx) + 1e-10
    KLdist <- c(KLdist, sum((f1-f2)*log(f1/f2)))
  }
  return(KLdist)
}