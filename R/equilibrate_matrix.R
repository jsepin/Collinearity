equilibrate_matrix <- function(object){
  # Input is a design matrix X
  eX <- object %*% diag(1 / sqrt(diag(t(object) %*% object)))
  colnames(eX) <- colnames(object)
  return(eX)
}