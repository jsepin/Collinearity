cond_num.lm <- function(object, values = FALSE){
  regressors <- model.matrix( object )
  
  regressors_Tr <- t( regressors )
  reg_eigen_val <- eigen( regressors_Tr %*% regressors )$values
  
  cond_number <- max( reg_eigen_val) / min( reg_eigen_val )
  
  vec_mat <- eigen( regressors_Tr %*% regressors )$vectors
  rownames(vec_mat) <- colnames(model.matrix(object))
  
  reg_eigen_val_mat <- matrix(reg_eigen_val, ncol = 1)
  
  cond_ind <- max( reg_eigen_val ) / reg_eigen_val
  
  out_l <- list("condition number" = cond_number,
                "Eigen values" = reg_eigen_val_mat,
                "vectors" = vec_mat,
                "condition index" = cond_ind)
  
  if(values){
    
    return(out_l)}else{
      
      return(cond_number)                    
      
    }
}