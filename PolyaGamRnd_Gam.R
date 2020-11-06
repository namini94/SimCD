PolyaGamRnd_Gam <- function(a, c, Truncation=6L)
# Coverted to R from MATLAB code by Mingyuan Zhou
# Siamak Zamani, June 2016
{
    realmin <- 2.2251e-308
    c <- abs(c)
    idx <- a!=0       # 'a' should not be zero
    a <- a[idx]
    c <- c[idx]
    
    if (!is.null(a))
    {
        idx.c <- c!=0
        c1 <- c[idx.c]
        a1 <- a[idx.c]
        
        xmeanfull <- a/4
        xmeanfull[idx.c] <- a1*tanh(c1/2)/(2*c1)
        xvarfull <- a/24
        
        idx.c1 <- c>=1e-3
        c1 <- c[idx.c1]
        a1 <- a[idx.c1]
        xvarfull[idx.c1] <- 0.5*exp(log(a1)-3*log(c1)+log(-expm1(-2*c1)-2*c1*exp(-c1))-log(1+exp(-2*c1)+2*exp(-c1)))
        
        idx.c1 <- c<1e-3
        c1 <- c[idx.c1]
        a1 <- a[idx.c1]
        
        xvarfull[idx.c1] <- 0.5*exp(log(a1)+ pmax(-3*log(c1)+log(-expm1(-2*c1)-2*c1*exp(-c1))-log(1+exp(-2*c1)+2*exp(-c1)), -log(12)-2*logcosh(c1/2), na.rm = T))
        
        if (Truncation>1)
        {
            temp <- matrix(((1:(Truncation-1))-0.5)^2, length(c), Truncation-1, byrow = T) + matrix(c^2/(4*pi^2), length(c), Truncation-1, byrow = F)
            xmeantruncate <- 1/(2*pi^2)*a*rowSums(1/temp)
            xmean <- pmax(xmeanfull - xmeantruncate,0)
            xvartruncate <- 1/4/pi^4*a*rowSums(1/(temp)^2)
            
            x <- 1/2/pi^2*rowSums(matrix(rgamma(length(a)*(Truncation-1), rep(a,Truncation-1)), length(a), Truncation-1)/temp)
            
            xvar <- pmax(xvarfull - xvartruncate,0)
            
            dex1 <- xvarfull>=(xvartruncate+realmin)
            
            x[!dex1] <- x[!dex1] + xmean[!dex1]
            if (sum(dex1)>0)
            {
                x[dex1] <- x[dex1] + rgamma(sum(dex1), xmean[dex1]^2/xvar[dex1],rate = 1) * (xvar[dex1]/xmean[dex1])
            }
            
        }else{
                cc <- xmeanfull/pmax(xvarfull,realmin)
                aa <- xmeanfull*cc
                x <- randg(aa)/pmax(cc,realmin)
        }
        
        temp <- x
        x <- rep(0,length(idx))
        x[idx] <- temp
            
        
    }else 
        x <- NULL
    return(x)
}