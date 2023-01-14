\name{copowerlm}
\alias{copowerlm}
\title{Power calculations for adjustment of collinearity}
\description{Compute the power or determine the sampe size to obtain a target power.}
\usage{copowerlm(power=NULL, n = NULL, alpha=0.05, Delta, sigma, p=3 , voilen=0.5, trouble = 1)}
\arguments{
  \item{power}{power of test (1 minus Type II error probability)}
  \item{n}{number of observations in total}
  \item{alpha}{significance level (Type I error probability)}
  \item{Delta}{relevant effect of variable of interest}
  \item{sigma}{experimental error}
  \item{p}{number of parameters in whole model}
  \item{voilen}{variable of interest length}
  \item{trouble}{diagonal entry of the inverted squared equilibrated design matrix}
}
\details{}
\value{}
\references{}
\author{Jerome Sepin, Malgorzata Roos}
\examples{
library(Collinearity)
# when a lm model is available the relevant parts are extracted
m <- lm(data = iris, Sepal.Length ~.)
voi <- "Petal.Length"
copowerlm(power = 0.8,  n = NULL, alpha = 0.05,
          Delta = c(0.1, 0.2, 0.3), m = m, voi = voi)$n
copowerlm(power = NULL, n = NULL, alpha = 0.05,
          Delta = c(0.1, 0.2, 0.3), m = m, voi = voi)$power

# all components can also be given individually
sigma <- stats::sigma(m) # 0.3068261
p <- ncol(model.matrix(m)) # 6
voilen <- var(model.matrix(m)[,voi]) + mean(model.matrix(m)[,voi])^2 # 17.23884
E <- equilibrate_matrix(model.matrix(m))
trouble <- solve(t(E) %*% E)[voi,voi] # 128.8317
copowerlm(power = 0.8, n = NULL, alpha = 0.05, Delta = c(0.1, 0.2, 0.3),
          sigma = sigma, p = p , voilen = voilen, trouble = trouble )
}

