##' cpGModule
##'
##' @title Identify cell phenotype activated gene module
##'
##' @description Mining activated gene modules in specific cell phenotype.
##'
##' @param network.data Network data constructed by the \code{ConNetGNN} function.
##' @param cellset A vector of cell id. The specified cell set, which will be used as the restart set.
##' @param nperm Number of random permutations. Default: \code{100}.
##' @param cut.pvalue The threshold of P-value, and genes below this threshold are regarded as gene modules activated by the cell set. Default: \code{0.01}.
##' @param cut.fdr The threshold of false discovery rate (FDR), and genes below this threshold are regarded as gene modules activated by the cell set. Default: \code{0.05}.
##' @param parallel.cores Number of processors to use when doing the calculations in parallel (default: \code{2}). If \code{parallel.cores=0}, then it will use all available core processors unless we set this argument with a smaller number.
##' @param rwr.gamma Restart parameter. Default: \code{0.7}.
##' @param normal_dist Whether to use pnorm to calculate P values. Default: \code{TRUE}.Note that if normal_dist is FALSE, we need to increase nperm (we recommend 100).
##' @param verbose Gives information about each step. Default: \code{TRUE}.
##'
##' @details
##' The \code{cpGModule} function takes a user-defined cell set as a restart set to automatically
##' identify activated gene modules. A perturbation analysis was used to calculate a significant P-value for each gene.
##' The \code{Benjamini & Hochberg (BH)} method was used to adjust the P-value to obtain the FDR.
##' Genes with a significance level less than the set threshold are considered as cell phenotype activated gene modules.
##'
##' @return
##' A data frame contains four columns:
##' \describe{
##'   \item{Genes}{Gene ID.}
##'   \item{AS}{Activity score.}
##'   \item{Pvalue}{Significant P-value.}
##'   \item{FDR}{False discovery rate.}
##' }
##'
##' @importFrom parallel makeCluster
##' @importFrom parallel clusterEvalQ
##' @importFrom parallel parLapply
##' @importFrom parallel stopCluster
##' @importFrom stats p.adjust
##' @importFrom stats na.omit
##'
##' @export
##'
##' @examples
##' require(parallel)
##' require(stats)
##'
##' # Load the result of the ConNetGNN function.
##' data(ConNetGNN_data)
##' data(Hv_exp)
##'
##' # Construct the cell set corresponding to 0h.
##' index<-grep("0h",colnames(Hv_exp))
##' cellset<-colnames(Hv_exp)[index]
##' cpGModule_data<-cpGModule(ConNetGNN_data,cellset,nperm=10,parallel.cores=1)


cpGModule<-function(network.data,cellset,nperm=100,cut.pvalue=0.01,cut.fdr=0.05,parallel.cores=2,rwr.gamma=0.7,normal_dist=TRUE,verbose=TRUE){
  if(!isLoaded("parallel")){
    stop("The package parallel is not available!")
  }

  if(!isLoaded("stats")){
    stop("The package stats is not available!")
  }


  cg_net<-network.data[[3]]
  cg_net<-apply(cg_net,2,function(x){
    if(all(x==0)){
      return(x)
    }else{
      return(x/max(x))
    }
  })

  c_net<-network.data[[1]]
  diag(c_net)<-0

  g_net<-network.data[[2]]
  diag(g_net)<-0

  merge1<-cbind(c_net,t(cg_net))
  merge2<-cbind(cg_net,g_net)
  cell_gene_network<-rbind(merge1,merge2)

  cell_n<-dim(c_net)[1]
  cellindex<-c(1:cell_n)
  geneindex<-(cell_n+1):nrow(cell_gene_network)

  if (verbose) {
    cat("Start RWR  \n")
  }

  pp<-match(cellset,row.names(cell_gene_network))
  names(pp)<-cellset
  pp<-na.omit(pp)

  diag.D <- apply(cell_gene_network,1,sum);
  diag.D[diag.D==0] <- Inf;
  inv.diag.D <- 1/diag.D;
  nadjM <-cell_gene_network*inv.diag.D
  rm(diag.D)
  rm(cell_gene_network)

  resW <- RWR(nadjM, pp, gamma=rwr.gamma)
  Actscores<-resW[geneindex]

  if (verbose) {
    cat("Start perturbation analysis \n")
  }

  cl <- makeCluster(parallel.cores)
  rdmatrix<-parLapply(cl,1:nperm,function(r,cellindex,pp,nadjM,rwr.gamma,geneindex){
    samplei<-sample(cellindex,size=length(pp))
    names(samplei)<-row.names(nadjM)[samplei]
    resW_rd <- RWR(nadjM, samplei, gamma=rwr.gamma)
    return(resW_rd[geneindex])
  },cellindex,pp,nadjM,rwr.gamma,geneindex)
  stopCluster(cl)
  rdmatrix<-do.call("cbind",rdmatrix)
  
  gc()
  
  pvalue<-NULL
  for(j in 1:length(Actscores)){
    ifelse(normal_dist==TRUE, pvalue[j]<-pnorm(Actscores[j],mean = mean(rdmatrix), sd = sd(rdmatrix),lower.tail = F),
           pvalue[j]<-sum(rdmatrix>=Actscores[j])/length(rdmatrix))
  }

  fdr<-p.adjust(pvalue,"BH",length(pvalue))

  GeneModule<-data.frame(Genes=names(Actscores),AS=Actscores,Pvalue=pvalue,FDR=fdr)
  GeneModule<-GeneModule[which(GeneModule$FDR<=cut.fdr&GeneModule$Pvalue<=cut.pvalue),]

  return(GeneModule)
}

