##' InteNet
##'
##' @title Integrate network data from single-cell RNA-seq and ATAC-seq
##'
##' @description
##' For the SNARE-seq dataset, a droplet-based method to simultaneously profile gene expression and chromatin accessibility in each of thousands of single nuclei,
##' the \code{InteNet} function can integrate network data of scRNA-seq data and scATAC-seq data (results of the \code{ConNetGNN} function) to into a gene-cell network.
##'
##' @param RNA_net Network data for RNA datasets. Produced by the \code{ConNetGNN} function.
##' @param ATAC_net Network data for ATAC datasets. Produced by the \code{ConNetGNN} function.
##' @param parallel.cores Number of processors to use when doing the calculations in parallel (default: \code{2}). If \code{parallel.cores=0}, then it will use all available core processors unless we set this argument with a smaller number.
##' @param verbose Gives information about each step. Default: \code{TRUE}.
##'
##' @details
##' The scATAC-seq dataset needs to be converted into a gene activity matrix according to the process of \code{Signac}(\code{https://satijalab.org/signac/articles/snareseq.html}).
##' The subsequent process is consistent with the scRNA-seq dataset. The \code{InteNet} function then integrates the network data of RNA-seq data and ATAC-seq data into a gene-cell network.
##' With integrated network data as input, \code{scPathway} and \code{cpGModule} functions will infer pathway activity score matrix and gene modules supported by single-cell multi-omics.
##'
##'
##' @return A list.
##'
##' @importFrom ActivePathways merge_p_values
##' @importFrom parallel makeCluster
##' @importFrom parallel clusterEvalQ
##' @importFrom parallel parLapply
##' @importFrom parallel stopCluster
##'
##' @export
##'
##' @examples
##' require(ActivePathways)
##' require(parallel)
##' data(RNA_net)
##' data(ATAC_net)
##' \dontrun{
##' RNA_ATAC_IntNet<-InteNet(RNA_net,ATAC_net,parallel.cores=1)
##' }
##'
##' # View data
##' data(RNA_ATAC_IntNet)
##' summary(RNA_ATAC_IntNet)

InteNet<-function(RNA_net,ATAC_net,parallel.cores=2,verbose=TRUE){
  if(!isLoaded("ActivePathways")){
    stop("The package ActivePathways is not available!")
  }

  if(!isLoaded("parallel")){
    stop("The package parallel is not available!")
  }

  if(all(colnames(ATAC_net[["gene_cell_network"]])==colnames(RNA_net[["gene_cell_network"]]))==FALSE){
    stop("The cells of scRAN-seq data and scATAC-seq data do not match!")
  }

  Integrate_m<-list()
  ver<-c("Integrate cell network \n","Integrate gene network \n","Integrate cell-gene network \n")

  #cell
  if(verbose==TRUE){
    cat("Merge cell network \n")
  }
  cl <- makeCluster(parallel.cores)
  clusterEvalQ(cl,library(ActivePathways))
  px<-match(colnames(RNA_net[["cell_network"]]),colnames(ATAC_net[["cell_network"]]))
  inte_cell<-parLapply(cl,1:nrow(RNA_net[["cell_network"]]),function(j,inter_ATAC,inter_RNA){
      ret<-merge_p_values(as.matrix(cbind(as.numeric(inter_RNA[j,]), as.numeric(inter_ATAC[j,]))), method='Brown')
      return(ret)
  },ATAC_net[["cell_network"]][px,px],RNA_net[["cell_network"]])
  stopCluster(cl)
  inte_cell<-do.call("rbind",inte_cell)
  colnames(inte_cell)<-colnames(RNA_net[["cell_network"]])
  row.names(inte_cell)<-row.names(RNA_net[["cell_network"]])
  Integrate_m[["cell_network"]]<-inte_cell

  #gene
  if(verbose==TRUE){
    cat("Merge gene network \n")
  }
  allgenes<-unique(c(colnames(ATAC_net[["gene_network"]]),colnames(RNA_net[["gene_network"]])))
  rg<-intersect(colnames(ATAC_net[["gene_network"]]),colnames(RNA_net[["gene_network"]]))
  pp<-match(rg,colnames(ATAC_net[["gene_network"]]))
  ATAC_jj<-ATAC_net[["gene_network"]][pp,pp]
  pp<-match(rg,colnames(RNA_net[["gene_network"]]))
  RNA_jj<-RNA_net[["gene_network"]][pp,pp]
  cl <- makeCluster(parallel.cores)
  clusterEvalQ(cl,library(ActivePathways))
  inte_cell<-parLapply(cl,1:nrow(RNA_jj),function(j,inter_ATAC,inter_RNA){
    ret<-merge_p_values(as.matrix(cbind(as.numeric(inter_RNA[j,]), as.numeric(inter_ATAC[j,]))), method='Brown')
    return(ret)
  },ATAC_jj,RNA_jj)
  stopCluster(cl)
  inte_cell<-do.call("rbind",inte_cell)
  colnames(inte_cell)<-colnames(RNA_jj)
  row.names(inte_cell)<-row.names(RNA_jj)

  gene_network<-matrix(0,nrow = length(allgenes),ncol = length(allgenes))
  colnames(gene_network)<-allgenes
  row.names(gene_network)<-allgenes

  pp<-match(colnames(inte_cell),colnames(gene_network))
  gene_network[pp,pp]<-inte_cell

  yjg<-setdiff(colnames(ATAC_net[["gene_network"]]),rg)
  pp<-match(yjg,colnames(ATAC_net[["gene_network"]]))
  pp1<-match(yjg,colnames(gene_network))
  gene_network[pp1,pp1]<-ATAC_net[["gene_network"]][pp,pp]

  yjg<-setdiff(colnames(RNA_net[["gene_network"]]),rg)
  pp<-match(yjg,colnames(RNA_net[["gene_network"]]))
  pp1<-match(yjg,colnames(gene_network))
  gene_network[pp1,pp1]<-RNA_net[["gene_network"]][pp,pp]
  Integrate_m[["gene_network"]]<-gene_network

  #gene-cell
  if(verbose==TRUE){
    cat("Merge gene-cell network \n")
  }
  px<-match(colnames(RNA_net[["gene_cell_network"]]),colnames(ATAC_net[["gene_cell_network"]]))
  ATAC_net[["gene_cell_network"]]<-ATAC_net[["gene_cell_network"]][,px]
  rna_m<-RNA_net[["gene_cell_network"]]
  atac_m<-ATAC_net[["gene_cell_network"]]

  mi<-min(rna_m)
  ma<-max(rna_m)
  rna_m<-(rna_m-mi)/(ma-mi)

  mi<-min(atac_m)
  ma<-max(atac_m)
  atac_m<-(atac_m-mi)/(ma-mi)

  gene_cell_network<-matrix(0,nrow = length(allgenes),ncol = ncol(rna_m))
  colnames(gene_cell_network)<-colnames(rna_m)
  row.names(gene_cell_network)<-allgenes

  pp<-match(rg,row.names(atac_m))
  ATAC_jj<-atac_m[pp,]

  pp<-match(rg,row.names(rna_m))
  RNA_jj<-rna_m[pp,]

  cl <- makeCluster(parallel.cores)
  clusterEvalQ(cl,library(ActivePathways))
  inte_cell<-parLapply(cl,1:nrow(RNA_jj),function(j,inter_ATAC,inter_RNA){
    ret<-merge_p_values(as.matrix(cbind(as.numeric(inter_RNA[j,]), as.numeric(inter_ATAC[j,]))), method='Brown')
    return(ret)
  },ATAC_jj,RNA_jj)
  stopCluster(cl)
  inte_cell<-do.call("rbind",inte_cell)
  colnames(inte_cell)<-colnames(RNA_jj)
  row.names(inte_cell)<-row.names(RNA_jj)

  pp<-match(row.names(inte_cell),row.names(gene_cell_network))
  gene_cell_network[pp,]<-inte_cell

  yjg<-setdiff(row.names(rna_m),rg)
  pp<-match(yjg,row.names(rna_m))
  pp1<-match(yjg,row.names(gene_cell_network))
  gene_cell_network[pp1,]<-rna_m[pp,]

  yjg<-setdiff(row.names(atac_m),rg)
  pp<-match(yjg,row.names(atac_m))
  pp1<-match(yjg,row.names(gene_cell_network))
  gene_cell_network[pp1,]<-atac_m[pp,]
  Integrate_m[["gene_cell_network"]]<-gene_cell_network

  x<-Integrate_m[[1]]
  x[upper.tri(x)] <- 0
  s<-Integrate_m[[1]]
  s[lower.tri(s)] <- 0
  s<-t(s)
  cell_m<-(x+s)/2
  b<-cell_m+t(cell_m)
  diag(b)<-diag(cell_m)
  Integrate_m[[1]]<-b

  x<-Integrate_m[[2]]
  x[upper.tri(x)] <- 0
  s<-Integrate_m[[2]]
  s[lower.tri(s)] <- 0
  s<-t(s)
  cell_m<-(x+s)/2
  b<-cell_m+t(cell_m)
  diag(b)<-diag(cell_m)
  Integrate_m[[2]]<-b

  return(Integrate_m)
}
