logcosh <- function(x)
{
    return(abs(x)-log(2)+log1p(exp(-2*abs(x))))
}
