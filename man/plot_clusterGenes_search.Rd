% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shinyModules.R
\name{plot_clusterGenes_search}
\alias{plot_clusterGenes_search}
\title{scClustViz plot: Plot within-cluster gene expression highlighting selected genes}
\usage{
plot_clusterGenes_search(sCVd, selClust, GOI)
}
\arguments{
\item{sCVd}{The sCVdata object.}

\item{selClust}{A named character vector representing the cluster to be
displayed. If \code{\link{labelCellTypes}} has been run, pass the
appropriate element of \code{attr(Clusters(sCV),"ClusterNames")} to this
argument to show both cluster number and cell type label in the legend.}

\item{GOI}{A character vector of gene names to highlight.}
}
\description{
This function makes a scatterplot of gene detection rate vs. mean detected
gene abundance, highlighting specified genes.
}
\examples{
\dontrun{
pdf("filepath.pdf",width=12,height=7)
plot_clusterGenes_search(sCVd=sCVdata,
                         selClust="1",
                         GOI=c("Actb","Sox2"))
dev.off()
}

}
