\name{copowerlm}
\alias{copowerlm}
\title{Power calculations for adjustment of collinearity}
\description{Compute the power or determine the sampe size to obtain a target power.}
\usage{copowerlm(power=NULL, n=NULL, alpha=0.05, Delta, sigma, p=3 , voilen=0.5, trouble=1, m=NULL, voi=NULL)}
\arguments{
  \item{power}{Power of test (1 minus Type II error probability)}
  \item{n}{Number of observations in total}
  \item{alpha}{Significance level (Type I error probability)}
  \item{Delta}{Relevant effect of variable of interest}
  \item{sigma}{Experimental error}
  \item{p}{Number of parameters in whole model}
  \item{voilen}{Variable of interest length}
  \item{trouble}{Diagonal entry of the inverted squared equilibrated design matrix}
  \item{m}{lm model. Default is NULL.}
  \item{voi}{Variable of interest. Default is NULL.}
}
\details{
The underlying core of the function is the following equation:
\deqn{
\text{power}=1-\texttt{pf}\Bigl( 
\texttt{q=}\texttt{qf}\left(\texttt{p=}1-\alpha, \\
\texttt{df1=}1,\texttt{ df2=}n-p\right), \\
\texttt{df1=}1,\texttt{ df2=}n-p, \\
\texttt{ncp=}\left[\sum_{i=1}^{n}\boldsymbol{X}[i,j]^2\right]\cdot\frac{\Delta^2}{\hat\sigma^2 }\cdot\frac{1}{\left[\left(\boldsymbol{E^\top E}\right)^{-1}\right][j,j]} 
\Bigr)
}
or
\deqn{
\text{power}=424 \\
\textcolor{white}{\text{power}}=424\\
\quad=424
}
with \eqn{\Delta=\texttt{Delta}}, \eqn{\hat\sigma=\texttt{sigma}}, \eqn{\left[\left(\boldsymbol{E^\top E}\right)^{-1}\right][j,j]=\texttt{trouble}} where \eqn{j} represents the position of the variable of interest (\eqn{\texttt{voi}}).
Furthermore, \eqn{\boldsymbol{E}} represents the equilibrated design matrix \eqn{\boldsymbol{X}}.
In systems without collinearity as it is for example in univariate regression models (also without intercept), \eqn{\texttt{trouble}} reduces to 1.

Moreover, it can be shown that the following holds:
\deqn{
\text{E}\left(\sum_{i=1}^{n}X_{ij}^2\right)=n\left(\text{Var}(X_{ij})+\text{E}(X_{ij})^2\right)
}
Thus, \eqn{\texttt{voilen}} is taken to be \eqn{\text{Var}(\boldsymbol{X}[i,j])+\text{E}(\boldsymbol{X}[i,j])^2}. Note, if \eqn{\texttt{voi}} is binary, \eqn{X_{ij}} is Bernoulli distributed with probability to success \eqn{p} which is then also equals \eqn{\texttt{voilen}}.
}
\value{}
\references{}
\author{Jerome Sepin, Malgorzata Roos}
\examples{
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
trouble <- solve(t(E) \%*\% E)[voi,voi] # 128.8317
copowerlm(power = 0.8, n = NULL, alpha = 0.05, Delta = c(0.1, 0.2, 0.3),
          sigma = sigma, p = p , voilen = voilen, trouble = trouble )

# the power with the data currently at hand can also be determined
copowerlm(power = NULL, n = dim(iris)[1], alpha = 0.05, Delta = c(0.1, 0.2, 0.3),
          sigma = sigma, p = p , voilen = voilen, trouble = trouble )
          
}

