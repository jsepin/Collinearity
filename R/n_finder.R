# R: Function to determine sample size or power
#' @param alpha Significance level
#' @param n Sample size
#' @param power Power
#' @param Delta Relevant difference
#' @param sigma Experimental error
#' @param p Number of parameters in whole model
#' @param voilen Variable Of Interest length
#' @param trouble Diagonal entry of the inverted squared equilibrated design matrix
#'
#' @return Returns the power or sample size

n_finder <- function(power=NULL,n = NULL ,alpha=0.05, Delta, sigma, p=3 , voilen=0.5, trouble = 1){
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
        warning = function(w){result <<- NA},
        error  = function(e){result <<- NA}
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
      warning = function(w){result <<- NA},
      error  = function(e){result <<- NA}
    )
    return(result)
  })
  
  # give back
  out <- list("conditions" = conditions, "n" = re)
  }
  
  return(out)
}

