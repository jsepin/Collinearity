# R: Function to determine sample size or power
#' @param alpha Significance level
#' @param n Sample size. Default is NULL.
#' @param power Power. Default is NULL.
#' @param Delta Relevant difference
#' @param sigma Experimental error
#' @param p Number of parameters in whole model. Default is 3.
#' @param voilen Variable Of Interest length. Default is 0.5.
#' @param trouble Diagonal entry of the inverted squared equilibrated design matrix. Default is 1.
#' @param m lm model. Default is NULL.
#' @param voi Variable of interest. Default is NULL.
#'
#' @return Returns the power or sample size

copowerlm <- function(power=NULL,n = NULL ,alpha=0.05, Delta, sigma, p=3 , voilen=0.5, trouble = 1,m = NULL, voi = NULL){
  # if lm model is given
  if(!is.null(m)){
    X <- model.matrix(m)
    # search voi
    if(!(voi %in% colnames(X)) || is.null(voi)){
      stop(paste0("voi not found! \n Try: ", paste0(gsub("[`]","",colnames(X)),collapse = ", ") ) )
    }
    E <- equilibrate_matrix(X)
    trouble <- solve(t(E) %*% E)[voi,voi]
    voilen  <- var(X[,voi]) + mean(X[,voi])^2
    sigma <- stats::sigma(m)
    p <- ncol(X)
    if(is.null(power)){n <- nrow(X)}
  }
  
  # power function
  myFpower <- function(alpha, n, Delta, sigma, p , voilen, trouble){
    ncp <- ((n *voilen*  Delta^2)/(sigma^2))  * (1/trouble)
    df1 <- 1
    df2 <- n-p
    Fcrit <-  qf(1-alpha,df1,df2)
    power <- 1-pf(Fcrit,df1,df2,ncp)
    return(power)
  }
  # Power is wanted:
  if(is.null(power)){
    conditions <- expand.grid("n"=n,"alpha"=alpha,
                              "Delta"=Delta, "sigma"=sigma,
                              "p"=p, "voilen"=voilen, "trouble"=trouble)
    re <-apply(conditions, 1, function(x) {
      n <- x[1]; alpha <- x[2]; Delta <- x[3]; sigma <- x[4]; p <- x[5];
      voilen <- x[6]; trouble <- x[7]
      # power calculation
      tryCatch({
        result <- myFpower(alpha, n,Delta, sigma, p, voilen, trouble) },
        warning = function(w){result <<- NA; print(warnings()) },
        error  = function(e){result <<- NA; print("Error! Return NA") }
      )
      return(result)
      })
    # give back
    out <- list("conditions" = conditions, "power" = re)
  }
  
  # n is wanted:
  else if(is.null(n)){
  conditions <- expand.grid("power"=power,"alpha"=alpha,
                            "Delta"=Delta, "sigma"=sigma,
                            "p"=p, "voilen"=voilen, "trouble"=trouble)
  
  re <-apply(conditions, 1, function(x) {
    x <- unlist(x) 
    power <- x[1]; alpha <- x[2]; Delta <- x[3]; sigma <- x[4]; p <- x[5];
    voilen <- x[6]; trouble <- x[7]
    
    #loss function
    f <- function(power, alpha, n, Delta, sigma, p, voilen, trouble){
      return(power - myFpower(alpha, n, Delta, sigma, p, voilen, trouble))
    }
    # uniroot
    tryCatch({
      result <- uniroot(f, interval = c(p+1,1000000), extendInt="yes", tol = 0.00001,
                        power = power, alpha = alpha, Delta=Delta,p=p,
                        sigma = sigma, voilen=voilen, trouble = trouble
      )$root },
      warning = function(w){result <<- NA; print(summary(warnings()) )},
      error  = function(e){result <<- NA; print("Error! Return NA")}
    )
    return(result)
  })
  
  # give back
  out <- list("conditions" = conditions, "n" = re)
  }
  
  return(out)
}


