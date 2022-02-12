\name{Var_decom_mat.matrix}
\alias{Var_decom_mat.matrix}

\title{ Variance decomp
}
\description{Variance decomp}
\usage{
Var_decom_mat.matrix(object, equilibration)
}
\arguments{
  \item{object}{a matrix}
  \item{equilibration}{Logical. If TRUE, the design matrix is equilibrated first. Default value is TRUE}
  
}
\author{Georgios Kazantzidis, Malgorzata Roos}
\examples{
n <- 20  
set.seed(1237) 
cnst <- rep(1,n)
pred_1 <- runif(n, -30, 40)  
pred_2 <- runif(n, 30, 100)  
pred_3 <- runif(n, -80, -10) 
pred_4 <- runif(n, -1, 5) 
pred_5 <- runif(n, -100, -80) 
pred_6 <- runif(n, -20, 20) 
pred_7 <- runif(n, 0, 35) 
pred_8 <- runif(n, 0, 35) 
# Predictor 9 is exactly collinear
pred_9 <-  pred_1 + pred_2 + pred_3 + pred_4 + pred_5 + pred_6 + pred_7 + pred_8 + rnorm(n, mean = 0, sd = .005)

eps_y <- rnorm(n, mean = 0, sd = .5)  # noise for response 
y_fin <- pred_1 * 2 + pred_2 * 3 + pred_3 * 1 - pred_4 * 7 + pred_5 * .5 + pred_6 * 2 + pred_7 * 2 + pred_8 * 2 + pred_9 + eps_y
f_1 <- y_fin ~ pred_1 + pred_2 + pred_3 + pred_4 + pred_5 + pred_6 + pred_7 + pred_8 + pred_9 + cnst - 1

d_dat <- data.frame(y_fin, pred_1, pred_2, pred_3, pred_4, 
                    pred_5, pred_6, pred_7, pred_8, pred_9)
                    
mod_9 <- lm(f_1, data = d_dat)
mod_9_XX <- model.matrix(mod_9)
Var_decom_mat.matrix(mod_9_XX)
}

