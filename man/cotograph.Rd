\name{cotograph}
\alias{cotograph}
\title{Collinearity graph}
\description{Compute the Wald statistics and the condition numbers of bootstrap samples to check the robustness of the results.}
\usage{cotograph(m, voi, equilibrate = FALSE, main = "Graph: Relation of explanatory variables\n Multivariable fitted Model", cex.node = 1, cex.tovoi = 1, cex.main = 1, col.edge.line = "blue", col.edge.text = "blue", col.node.voi = "green", col.node.nonvoi = "lightblue", radius_circle = 0.2, subR2 = FALSE, lwdtovoi = 1, mar = c(2.1, 2.1,2.1, 2.1))}
\arguments{
  \item{m}{Model fitted}
  \item{voi}{Variable of interest specified as character}
  \item{equilibrate}{If TRUE, the design matrix is equilibrated first. Default value is FALSE}
  \item{main}{Title for the plot.}
  \item{cex.node}{Text size of nodes. Default is 1.}
  \item{cex.tovoi}{Text size of edge to voi-node. Default is 1.}
  \item{col.edge.line}{Color of edges to voi-node. Default is blue.}
  \item{col.edge.text}{Color of edges text to voi-node. Default is blue.}
  \item{col.node.voi}{Color of the voi node. Default is green.}
  \item{col.node.nonvoi}{Color of non-voi node(s). Default is lightblue.}
  \item{radius_circle}{Size of node. Default is 0.2.}
  \item{subR2}{Logical. If TRUE, the R2 values of the non-voi nodes are plotted as well. Default value is FALSE.}
  \item{lwdtovoi}{Maximum linewidth of edges. Default is 5.}
  \item{mar}{Margin of the plot. c(2.1, 2.1,2.1, 2.1) }
}
\details{}
\value{}
\references{}
\author{Jerome Sepin, Malgorzata Roos}
\examples{
dev.off(); par(mfrow = c(2,1))
m <- lm(data = iris, Sepal.Length ~.)
cotograph(m=m, voi = "Petal.Length")

# equilibrated design matrix
cotograph(m=m, voi = "Petal.Length", equilibrate = TRUE)

# with sub-R2 values added
dev.off()
m <- glm(data = iris, I(iris$Species=="versicolor") ~. ,family = binomial)
cotograph(m=m, voi = "Petal.Length", subR2 = TRUE)
}