\name{lmcopri}
\alias{lmcopri}
\title{Collinearity fingerprint}
\description{Compute the Wald statistics and the condition numbers of bootstrap samples to check the robustness of the results.}
\usage{lmcopri(m, B = 100, returndata = FALSE, ncon = NULL, main = "Collinearity Fingerprint", xdi = 0, ydi = 0, alpha = 0.05, cex.vd = 1, cex.main = 2, cex.sub = 1, cex.prop  = 1)}
\arguments{
  \item{m} {Model fitted}
  \item{B} {Number of bootstrap samples to perform. Default is 100.}
  \item{returndata} {Logical or data frame. If TRUE, the bootstrap data frame used for plotting is returned. This data frame can then be given as input so that bootstrapping, which may be time consuming does not have to be executed again. Default is FALSE.}
  \item{ncon} {Number of condition indices plotted}
  \item{main} {Title of the plot}
  \item{xdi} {Horizontal position of the variance decomposition proportions. Default is 0.}
  \item{ydi} {Amount to increase ylim. Default is 0.}
  \item{alpha} {Significance level (Type I error probability). Default is 0.05.}
  \item{cex.vd} {The magnification to be used for the variance decomposition proportions relative to the current setting of cex. Default is 1.}
  \item{cex.main} {The magnification to be used for main titles relative to the current setting of cex. Default is 2.}
  \item{cex.sub} {The magnification to be used for sub-titles relative to the current setting of cex. Default is 1.}
  \item{cex.prop} {The magnification to be used for the proportions of significant results relative to the current setting of cex. Default is 1.}
}
\details{}
\value{}
\references{}
\author{Jerome Sepin, Malgorzata Roos}

