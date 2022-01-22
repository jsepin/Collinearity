Var_decom_mat.lm <- function(object, equilibration = TRUE){
  ## Input: a lm object 
  ## Output: the variance decompossition matrix along with the equilibrated condition indeces
  
  mat_X <- model.matrix(object)
  
  if(equilibration){ 
    
    eq_mat_X <- round( vdp_svd( equilibrate_matrix( mat_X )), 3)
    
  } else{
    
    eq_mat_X <- round( vdp_svd( mat_X ), 3)
  }
  colnames(eq_mat_X) <- c("mu" ,"cond_ind", colnames(mat_X))
  return(eq_mat_X)
}