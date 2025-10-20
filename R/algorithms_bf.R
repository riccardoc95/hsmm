
backward.forward <- function(fit, Gamma, Pi, K, M){
  
  ld <- rep(M, K)
  n <- nrow(fit)
  
  psi <- alpha <- Beta <- array(dim = c(n, sum(ld)))
  C	<- array(dim = n)
  
  post.pi.aug <- array(dim = c(n, sum(ld)))
  post.bi.pi.aug <- array(dim = c(n-1, sum(ld), sum(ld)))
  
  psi[1,] <- Pi
  psi_per_fit <- psi[1,]*fit[1,]
  C[1] <- sum(psi_per_fit)
  alpha[1,] <- (psi_per_fit)/C[1] 
  
  for(i in 2:n){
    psi[i, ] <- alpha[i-1,]%*%Gamma[i-1,,]
    psi_per_fit <- psi[i,]*fit[i,]
    C[i] <- sum(psi_per_fit)
    alpha[i,] <- (psi_per_fit)/C[i] 
    
  }
  
  l <- sum(log(C))
  
  Beta[n,] <- 1/C[n]
  
  for(i in (n-1):1){
    Beta[i, ] <- c(Gamma[i,,]%*%(fit[i+1,]*Beta[i+1,]))/C[i]
  }
  
  alphaBeta <- alpha*Beta
  post.pi.aug = alphaBeta/rowSums(alphaBeta)
  
  for(h in 1:sum(ld)){
    for(k in 1:sum(ld)){
      post.bi.pi.aug[,h,k] <- alpha[1:(n-1),h]*Gamma[,h,k]*fit[2:n,k]*Beta[2:n,k] 
    }
  }
  return(list(l = l, post.pi.aug = post.pi.aug, post.bi.pi.aug = post.bi.pi.aug))
}
