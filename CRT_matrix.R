CRT_matrix <- function(x,r)
    ## Chinese Restaurant Table distribution with matrix input and matrix output
    ## Siamak Zamani
    ## Created Jul 2016    
{
    dyn.load("CRT_matrix.so")
    if (sum(dim(x)!=dim(r))>0) stop('Inputs of CRT_matrix should have the same dimensions')
    dx <- dim(x)
    M <- dx[1]
    N <- dx[2]
    x <- c(x)
    r <- c(r)
    L <- rep(0, length(x))
    
    out <- .C("CRT_matrix", x=as.double(x), 
              r=as.double(r), M=as.integer(M), N=as.integer(N), L=as.integer(L))
    return(matrix(out$L, M, N))
}
