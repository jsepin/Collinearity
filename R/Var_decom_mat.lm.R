Var_decom_mat.lm <- function(object){
  ## Input: a lm object 
  ## Output: the variance decompossition matrix along with the equilibrated condition indeces
  
  mat_X <- model.matrix(object)
  eq_mat_X <- round(vdp_svd(equilibrate_matrix(mat_X)), 3)
  
  return(eq_mat_X)
}