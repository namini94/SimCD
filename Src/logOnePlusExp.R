logOnePlusExp <- function(x)
# Mingyuan Zhou
{
    dex <- x<0
    y <- matrix(0,dim(x)[1],dim(x)[2])
    y[dex] <- log1p(exp(x[dex]))
    y[!dex] <- x[!dex] + log1p(exp(-x[!dex]))
    return(y)
}