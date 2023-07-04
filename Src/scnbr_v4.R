scnbr_v4 <- function(counts, X, Z, K, ncores, Burnin = 1000L, Collections = 1000L, PGTruncation = 10L, randtry = 2017)
{
    set.seed(randtry)
    y <- as.matrix(counts)
    idx.nz <- rowSums(y) > -1
    y <- y[idx.nz,]
    V <- dim(y)[1]       # V
    J <- dim(y)[2]     # J
    # hyperparameters
    a0 <- b0 <- 1
    e0 <- f0 <- 0.01
    
    P <- dim(X)[1]      # X is P by J
    Q <- dim(Z)[1]      # Z is Q by V
    # Initialization
    Beta <- matrix(0,P,V)
    Phi <- matrix(1/V,V,K)
    Theta <- matrix(0,K,J)
    Delta <- matrix(0,Q,J)
    rj <- rep(1,J)
    h <- 1
    gammat <- rep(1,K)
    alpha <- rep(1,P)
    eta <- rep(1,Q)
    omegat <- matrix(0,V,J)
    s <- matrix(0,V,2*Collections)
    u <- matrix(0,V,1)
    Psi <- t(Beta) %*% X + t(Z) %*% Delta + Phi %*% Theta
    
    realmin <- .Machine$double.xmin
    maxlogLike <- -Inf
    out <- list(Phi=0*Phi,Beta=Beta,Delta=Delta,rj=0*rj,Theta=Theta, nama=s , kl=u)
    iterMax <- Burnin+Collections
    
    # set up the cluster
    #cl <- parallel::makeCluster(ncores)           # Windows
     cl <- parallel::makeForkCluster(ncores)     # Forking
    doParallel::registerDoParallel(cl)
    
    for (iter in 1:iterMax)
    {
        cat(iter, '\n')
        
        # sample r_j
        ell <- CRT_matrix(y,matrix(rj,V,J,byrow = T))
        temp <- logOnePlusExp(Psi)
        rj <- rgamma(J, shape = e0+colSums(ell), rate = h+colSums(temp))
        
        # sample h
        h <- rgamma(1, shape = e0+J*e0, rate = f0+sum(rj))
        
        for (v in 1:V){
            # sample omega[v,]
            omegat[v,] <- PolyaGamRnd_Gam(y[v,]+rj, Psi[v,], Truncation = PGTruncation)
        }
        
        # sample Beta
        Beta <- foreach(v=1:V,.combine = cbind) %dopar% {
            omega <- omegat[v,]
            sigmak <- X %*% diag(omega) %*% t(X)
            diag(sigmak) <- diag(sigmak) + pmax(alpha,1e-3)
            errMSG <- sum(eigen(sigmak)$values<=0)
            if (errMSG == 0){
                invchol <- chol(sigmak)
            }else{
                count_0 <- 0
                while (errMSG != 0)
                {
                    diag(sigmak) <- diag(sigmak) + 10^(count_0-6)
                    errMSG <- sum(eigen(sigmak)$values<=0)
                    count_0 <- count_0 + 1
                }
                invchol <- chol(sigmak)
            }
            invchol <- solve(invchol)
            muk1 <- X %*% t(0.5*(y[v,]-rj)-omega*(t(Phi[v,]) %*% Theta + t(Z[,v]) %*% Delta))
            if (iter>Burnin){
             # out$mean <- (invchol %*% t(invchol) %*% muk1) +  out$mean 
            }
             invchol %*% (rnorm(P)+t(invchol) %*% muk1)
            
          
            #print(sigmak)
            #print(muk)
        }
        
        
        # sample Beta_mean
       # s <- foreach(v=1:V,.combine = cbind) %dopar% {
        #  omega <- omegat[v,]
         # sigmak <- X %*% diag(omega) %*% t(X)
        #  diag(sigmak) <- diag(sigmak) + pmax(alpha,1e-3)
         # errMSG <- sum(eigen(sigmak)$values<=0)
        #  if (errMSG == 0){
         #   invchol <- chol(sigmak)
        #  }else{
         #   count_0 <- 0
          #  while (errMSG != 0)
           # {
            #  diag(sigmak) <- diag(sigmak) + 10^(count_0-6)
            #  errMSG <- sum(eigen(sigmak)$values<=0)
            #  count_0 <- count_0 + 1
            #}
            #invchol <- chol(sigmak)
          #}
          #invchol <- solve(invchol)
          #muk1 <- X %*% t(0.5*(y[v,]-rj)-omega*(t(Phi[v,]) %*% Theta + t(Z[,v]) %*% Delta))
          
          #invchol %*% t(invchol) %*% muk1 

        #}
        
        # sample Beta_var
        #u <- foreach(v=1:V,.combine = cbind) %dopar% {
          #omega <- omegat[v,]
          #sigmak <- X %*% diag(omega) %*% t(X)
          #diag(sigmak) <- diag(sigmak) + pmax(alpha,1e-3)
          #errMSG <- sum(eigen(sigmak)$values<=0)
          #if (errMSG == 0){
           # invchol <- chol(sigmak)
          #}else{
           # count_0 <- 0
            #while (errMSG != 0)
            #{
             # diag(sigmak) <- diag(sigmak) + 10^(count_0-6)
              #errMSG <- sum(eigen(sigmak)$values<=0)
              #count_0 <- count_0 + 1
            #}
            #invchol <- chol(sigmak)
          #}
          #invchol <- solve(invchol)
          #muk1 <- X %*% t(0.5*(y[v,]-rj)-omega*(t(Phi[v,]) %*% Theta + t(Z[,v]) %*% Delta))
          
          #invchol %*% t(invchol) 
          
        #}
        
        # sample Phi    
        Phi <- foreach(v=1:V,.combine = rbind) %dopar%  {   
            omega <- omegat[v,]
            sigmak <- Theta %*% diag(omega) %*% t(Theta)
            diag(sigmak) <- diag(sigmak) + 1
            errMSG <- sum(eigen(sigmak)$values<=0)
            if (errMSG == 0){
                invchol <- chol(sigmak)
            }else{
                count_0 <- 0
                while (errMSG != 0)
                {
                    diag(sigmak) <- diag(sigmak) + 10^(count_0-6)
                    errMSG <- sum(eigen(sigmak)$values<=0)
                    count_0 <- count_0 + 1
                }
                invchol <- chol(sigmak)
            }
            invchol <- solve(invchol)
            muk <- Theta %*% t(0.5*(y[v,]-rj)-omega*(t(Beta[,v]) %*% X + t(Z[,v]) %*% Delta))
            t(invchol %*% (rnorm(K)+t(invchol) %*% muk))
        }
        
        # sample Theta
        Theta <- foreach(j=1:J,.combine = cbind) %dopar% {
            omega <- omegat[,j]
            sigmak <- t(Phi) %*% diag(omega) %*% Phi
            diag(sigmak) <- diag(sigmak) + gammat
            errMSG <- sum(eigen(sigmak)$values<=0)
            if (errMSG == 0){
                invchol <- chol(sigmak)
            }else{
                count_0 <- 0
                while (errMSG != 0)
                {
                    diag(sigmak) <- diag(sigmak) + 10^(count_0-6)
                    errMSG <- sum(eigen(sigmak)$values<=0)
                    count_0 <- count_0 + 1
                }
                invchol <- chol(sigmak)
            }
            invchol <- solve(invchol)
            muk <- t(Phi) %*% t(0.5*(y[,j]-rj[j])-omega*(t(X[,j]) %*% Beta + t(Delta[,j]) %*% Z))
            invchol %*% (rnorm(K)+t(invchol) %*% muk)
        }
        
        
        # sample Delta
        Delta <- foreach(j=1:J,.combine = cbind) %dopar% {
            omega <- omegat[,j]
            sigmak <- Z %*% diag(omega) %*% t(Z)
            diag(sigmak) <- diag(sigmak) + pmax(eta,1e-3)
            errMSG <- sum(eigen(sigmak)$values<=0)
            if (errMSG == 0){
                invchol <- chol(sigmak)
            }else{
                count_0 <- 0
                while (errMSG != 0)
                {
                    diag(sigmak) <- diag(sigmak) + 10^(count_0-6)
                    errMSG <- sum(eigen(sigmak)$values<=0)
                    count_0 <- count_0 + 1
                }
                invchol <- chol(sigmak)
            }
            invchol <- solve(invchol)
            muk <- Z %*% t(0.5*(y[,j]-rj[j])-omega*(t(X[,j]) %*% Beta + t(Theta[,j]) %*% t(Phi)))
            invchol %*% (rnorm(Q)+t(invchol) %*% muk)
        }
        
        Psi <- t(Beta) %*% X + t(Z) %*% Delta + Phi %*% Theta
        
        # sample gammat
        gammat <- rgamma(K, e0+J/2, f0+rowSums(Theta^2)/2)
        
        # Sample alpha, eta
        alpha <- rgamma(P, 1 + V/2, rate = 1 + 0.5*rowSums(Beta^2))
        eta <- rgamma(Q, 1 + J/2, 1 + 0.5*rowSums(Delta^2))
        
        gc()
        
        if (iter>Burnin)
        {
            out$Beta <- Beta + out$Beta
            out$nama[,iter-Burnin] <- exp(Beta[1,])
            out$nama[,iter-Burnin+Collections] <- exp(Beta[1,]+Beta[2,])
            #out$me <- s + out$me
            #out$va <- u + out$va
            out$Delta <- Delta + out$Delta
            out$rj <- rj + out$rj
            
            logLike <- sum(lgamma(rj+t(y))) - V*sum(lgamma(rj)) - 
                sum(y * logOnePlusExp(-Psi)) -
                sum(matrix(rj,V,J,byrow = T) * logOnePlusExp(Psi))
            
            if(logLike > maxlogLike){
                maxlogLike <- logLike
                out$Phi <- Phi
                out$Theta <- Theta
            }
        }
        
    }
    
    parallel::stopCluster(cl)
    #out$me <- out$me/Collections
    #out$va <- out$va/Collections
    out$kl <- KLsym(out$nama[,1:Collections],out$nama[,(Collections+1):(2*Collections)])
    out$Delta <- out$Delta/Collections
    out$Beta <- out$Beta/Collections
    out$rj <- out$rj/Collections
    return(out)
}
