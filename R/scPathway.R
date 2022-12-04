##' scPathway
##'
##' @title Infer pathway activation score matrix at single-cell resolution
##'
##' @description Calculate pathway activity score of single-cell by random walk with restart (RWR).
##'
##' @param network.data The input network data is the result from the \code{ConNetGNN} function.
##' @param gmt.path Pathway database in \code{GMT} format.
##' @param pathway.min Minimum size (in genes) for pathway to be considered. Default: \code{10}.
##' @param pathway.max Maximum size (in genes) for database gene sets to be considered. Default: \code{500}.
##' @param nperm Number of random permutations. Default: \code{10}. We recommend setting it to 100 times.
##' @param parallel.cores Number of processors to use when doing the calculations in parallel (default: \code{2}). If \code{parallel.cores=0}, then it will use all available core processors unless we set this argument with a smaller number.
##' @param rwr.gamma Restart parameter. Default: \code{0.7}.
##' @param verbose Gives information about each step. Default: \code{TRUE}.
##'
##' @details
##' The \code{scPathway} function integrates the results of ConNetGNN into a gene-cell association network.
##' The genes included in each pathway are used as a restart set in the gene-cell association network to calculate the strength of its association with each cell through \code{RWR}.
##' Perturbation analysis was performed to remove noise effects in the network and to obtain the final single-cell pathway activity score matrix.
##'
##' @return A matrix of single-cell pathway activity score.
##'
##' @importFrom utils txtProgressBar
##' @importFrom utils setTxtProgressBar
##' @importFrom parallel makeCluster
##' @importFrom parallel clusterExport
##' @importFrom parallel parLapply
##' @importFrom parallel stopCluster
##' @importFrom stats na.omit
##'
##' @export
##'
##' @examples
##' require(parallel)
##' require(utils)
##' # Load the result of the ConNetGNN function.
##' data(ConNetGNN_data)
##' kegg.path<-system.file("extdata", "KEGG_human.gmt", package = "scapGNN")
##' # We recommend the use of a compiler.
##' # The compiler package can be used to speed up the operation.
##' # library(compiler)
##' # scPathway<- cmpfun(scPathway)
##' scPathway_data<-scPathway(ConNetGNN_data,gmt.path=kegg.path,
##'                           pathway.min=25,nperm=6,parallel.cores=1)
##'


scPathway<-function(network.data,gmt.path=NULL,pathway.min=10,pathway.max=500,nperm=10,parallel.cores=2,rwr.gamma=0.7,verbose=TRUE){
  if(!isLoaded("utils")){
    stop("The package utils is not available!")
  }

  if(!isLoaded("parallel")){
    stop("The package parallel is not available!")
  }

  cg_net<-network.data[[3]]
  cg_net<-t(apply(cg_net,1,function(x){
    if(all(x==0)){
      return(x)
    }else{
      return(x/max(x))
    }
  }))

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

  pathway_list<-load_path_data(gmt.path)

  del<-NULL
  for(i in 1:length(pathway_list)){
    inte<-length(intersect(pathway_list[[i]],row.names(cell_gene_network)))
    if(inte<pathway.min|inte>pathway.max){
      del<-c(del,i)
    }
  }
  pathway_list<-pathway_list[-del]

  if(verbose){
    cat(paste(length(pathway_list),"pathways are used for RWR","\n",sep = " "))
  }

  cl <- makeCluster(parallel.cores)
  #clusterEvalQ(cl,library(scapGNN))
  clusterExport(cl,'RWR')


  saPW_matrix<-NULL

  if(verbose){
    cat("start perturbation \n")
  }

  pb <- txtProgressBar(style=3)
  for(i in 1:length(pathway_list)){
    setTxtProgressBar(pb, i/length(pathway_list))

    pp<-match(pathway_list[[i]],row.names(cell_gene_network))
    names(pp)<-pathway_list[[i]]
    pp<-na.omit(pp)

    resW <- RWR(cell_gene_network, pp, gamma=rwr.gamma)

    rdmatrix<-parLapply(cl,1:nperm,function(r,geneindex,pp,cell_gene_network,rwr.gamma,cellindex){
      samplei<-sample(geneindex,size=length(pp))
      names(samplei)<-row.names(cell_gene_network)[samplei]
      resW_rd <- RWR(cell_gene_network, samplei, gamma=rwr.gamma)
      return(resW_rd[cellindex])
    },geneindex,pp,cell_gene_network,rwr.gamma,cellindex)
    rdmatrix<-do.call("rbind",rdmatrix)

    pvalue<-NULL
    for(j in 1:cell_n){
      pvalue[j]<-sum(rdmatrix[,j]>=resW[j])/nperm
    }
    saPW_matrix<-rbind(saPW_matrix,pvalue)
  }
  stopCluster(cl)
  row.names(saPW_matrix)<-names(pathway_list)
  colnames(saPW_matrix)<-colnames(c_net)

  saPW_matrix <- 1-saPW_matrix
  return(saPW_matrix)
}
