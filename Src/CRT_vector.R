CRT_vector <- function(x, r)
## Chinese Restaurant Table distribution with vector input and output arguments
## Siamak Zamani
## Created Jan 2016
{   
    dyn.load("CRT_vector")
    L <- rep(0L, length(x))
    out <- .C("CRT_vector", x=as.double(x), 
              r=as.double(r), Lenx=as.integer(length(x)), L=as.integer(L))
    return(out$L)
}