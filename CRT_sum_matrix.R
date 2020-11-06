CRT_sum_matrix <- function(x,r)
    ## Chinese Restaurant Table distribution with matrix input and vector output
    ## Siamak Zamani
    ## Created Dec 2017    
{
    dyn.load("CRT_sum_matrix")
    K <- dim(x)[1]
    J <- dim(x)[2]
    Lsum <- rep(0L,J)
    out <- .C("CRT_sum_matrix", x=as.double(c(x)), 
              r=as.double(r), K=as.integer(K), J=as.integer(J), Lsum=as.integer(Lsum))
    return(out$Lsum)
}