
vdp_svd <- function(object){
  # This function uses singular values decomposition to compute variance-decomposition
  # proportions for a design matrix X
  # Input: object is a design matrix X
  # Output: vdp
  
  s <- svd(object)
  # matrix based on singular values
  ms <- s$v  %*% diag((1/s$d)^2) %*% t(s$v)
  mu <- s$d
  
  msdim <- dim(ms)
  pom <- matrix(NA, nrow = msdim[1], ncol = msdim[2])
  
  for (i in 1:msdim[1]){
    for (j in 1:msdim[2]){
      
      pom[i,j]<-(s$v[i,j]^2/mu[j]^2)/diag(ms)[i]
      
    }
  }
  ci<- max(mu) / mu
  res <- as.matrix(cbind(mu, ci, t(pom)))
  return(res)
}