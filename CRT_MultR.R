CRT_MultR <- function(y,Phi,Theta)
      
{
    dyn.load("CRT_MultR")
    V <- dim(Phi)[1]
    K <- dim(Phi)[2]
    J <- dim(y)[2]
    lkj <- rep(0L,K*J)
    lvk <- rep(0L,K*V)
    out <- .C("CRT_MultR", y=as.double(c(y)), 
              Phi=as.double(c(Phi)), Theta=as.double(c(Theta)), V=as.integer(V),
              K=as.integer(K), J=as.integer(J), lkj=as.integer(lkj), lvk=as.integer(lvk))
    lkj <- matrix(out$lkj,K,J)
    lvk <- matrix(out$lvk,V,K)
    return(list(lkj=lkj,lvk=lvk))
}