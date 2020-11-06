dirrnd <- function(A)
{
    G <- rgamma(length(A),c(A),rate = 1)
    if (!is.null(dim(A))){
        G <- matrix(G,dim(A)[1],dim(A)[2])
        G <- t(t(G)/colSums(G))
    }else{
        G <- G/sum(G)
    }
    return(G)
}