\name{cond_num.matrix}
\alias{cond_num.matrix}
\title{Condition number}
\description{Calculates the condition number of a design matrix from a linear model.}
\usage{
cond_num(object)
}
\arguments{object
  \item{object}{A linear model or a matrix}
   \item{values}{Logical If TRUE, all the eigen value of the matrix are provided} 
}
\value{numeric. The condition number of the object provided. If values = TRUE, a list with the condition number, the eigen values, the eigen vectors and the condition indexes.}
\references{Montgomery, Douglas C., Elizabeth A. Peck, and G. Geoffrey Vining. Introduction to linear regression analysis. John Wiley and Sons, 2021.}
\author{Georgios kazantzidis \email{george_kazantzidis@outlook.com}}

