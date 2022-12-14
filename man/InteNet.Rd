% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InteNet.R
\name{InteNet}
\alias{InteNet}
\title{Integrate network data from single-cell RNA-seq and ATAC-seq}
\usage{
InteNet(RNA_net, ATAC_net, parallel.cores = 2, verbose = TRUE)
}
\arguments{
\item{RNA_net}{Network data for RNA datasets. Produced by the \code{ConNetGNN} function.}

\item{ATAC_net}{Network data for ATAC datasets. Produced by the \code{ConNetGNN} function.}

\item{parallel.cores}{Number of processors to use when doing the calculations in parallel (default: \code{2}). If \code{parallel.cores=0}, then it will use all available core processors unless we set this argument with a smaller number.}

\item{verbose}{Gives information about each step. Default: \code{TRUE}.}
}
\value{
A list.
}
\description{
For the SNARE-seq dataset, a droplet-based method to simultaneously profile gene expression and chromatin accessibility in each of thousands of single nuclei,
the \code{InteNet} function can integrate network data of scRNA-seq data and scATAC-seq data (results of the \code{ConNetGNN} function) to into a gene-cell network.
}
\details{
InteNet


The scATAC-seq dataset needs to be converted into a gene activity matrix according to the process of \code{Signac}(\code{https://satijalab.org/signac/articles/snareseq.html}).
The subsequent process is consistent with the scRNA-seq dataset. The \code{InteNet} function then integrates the network data of RNA-seq data and ATAC-seq data into a gene-cell network.
With integrated network data as input, \code{scPathway} and \code{cpGModule} functions will infer pathway activity score matrix and gene modules supported by single-cell multi-omics.
}
\examples{
require(ActivePathways)
require(parallel)
data(RNA_net)
data(ATAC_net)
\dontrun{
RNA_ATAC_IntNet<-InteNet(RNA_net,ATAC_net,parallel.cores=1)
}

# View data
data(RNA_ATAC_IntNet)
summary(RNA_ATAC_IntNet)
}
