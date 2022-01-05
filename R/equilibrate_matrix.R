equilibrate_matrix <- function(object){
  # Input is a design matrix X
  s_i <- 1 / sqrt(t(object) %*% object)
  
  eX <- object %*% diag(1 / sqrt(diag(t(object) %*% object)))
  return(eX)
}