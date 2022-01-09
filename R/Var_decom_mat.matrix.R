Var_decom_mat.matrix <- function(object, single_val_decomp = TRUE){
  ## Input: a matrix object 
  ## Output: the variance decompossition matrix along with the equilibrated condition indeces
  
  mat_X <- object
  
  if(single_val_decomp){ 
    
    eq_mat_X <- round( vdp_svd( equilibrate_matrix( mat_X )), 3)
  } else{
    
    cond_ind <- round( sqrt( cond_num( eq_mat_X, values = T)$`condition index`), 2)
    eq_mat_X <- cbind(cond_ind, eq_mat_X[, -1])
    colnames(eq_mat_X) <- c("cond_ind", colnames( mat_X ))
    
    }
  
  
  return(eq_mat_X)
}