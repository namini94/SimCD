CRT_sum <- function(x,r)
## Chinese Restaurant Table distribution with vector input and scalar output
## Siamak Zamani
## Created Jan 2016    
{
    dyn.load("CRT_sum")
    Lsum <- 0L
    out <- .C("CRT_sum", x=as.double(x), 
              r=as.double(r), Lenx=as.integer(length(x)), Lsum=as.integer(Lsum))
    return(out$Lsum)
}