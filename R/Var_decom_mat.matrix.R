Var_decom_mat.matrix <- function(object, equilibration = TRUE){
  ## Input: a matrix object 
  ## Output: the variance decompossition matrix along with condition indices
  ## If equilibration is TRUE, the equilibrated condition indices are provided
  
  mat_X <- object
  
  if(equilibration){ 
    
    eq_mat_X <- round( vdp_svd( equilibrate_matrix( mat_X )), 3)
    
  } else{
    
    eq_mat_X <- round( vdp_svd( mat_X ), 3)
  }
  colnames(eq_mat_X) <- c("mu" ,"cond_ind", colnames(mat_X))
  
  return(eq_mat_X)
}