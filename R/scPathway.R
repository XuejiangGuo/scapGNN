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
##' @param nperm Number of random permutations. Default: \code{50}. We recommend the setting of 100.
##' @param parallel.cores Number of processors to use when doing the calculations in parallel (default: \code{2}). If \code{parallel.cores=0}, then it will use all available core processors unless we set this argument with a smaller number.
##' @param rwr.gamma Restart parameter. Default: \code{0.7}.
##' @param normal_dist Whether to use pnorm to calculate P values. Default: \code{TRUE}.Note that if normal_dist is FALSE, we need to increase nperm (we recommend 100).
##' @param seed Random number generator seed.
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
##'                           pathway.min=25,nperm=2,parallel.cores=1)
##'


scPathway<-function(network.data,gmt.path=NULL,pathway.min=10,pathway.max=500,nperm=50,parallel.cores=2,rwr.gamma=0.7,normal_dist=TRUE,seed=1217,verbose=TRUE){
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

  rm(c_net)
  rm(g_net)
  rm(cg_net)
  rm(network.data)

  diag.D <- apply(cell_gene_network,1,sum);
  diag.D[diag.D==0] <- Inf;
  inv.diag.D <- 1/diag.D;
  nadjM <-cell_gene_network*inv.diag.D
  rm(diag.D)
  rm(cell_gene_network)

  pathway_list<-load_path_data(gmt.path)

  del<-NULL
  path_length <- NULL
  for(i in 1:length(pathway_list)){
    intgenes <- intersect(pathway_list[[i]],row.names(nadjM))
    inte<-length(intgenes)
    if(inte<pathway.min|inte>pathway.max){
      del<-c(del,i)
    }else{
      path_length <- c(path_length,inte)
      pathway_list[[i]] <- intgenes
    }
  }

  if(is.null(del)==FALSE){
    pathway_list<-pathway_list[-del]
  }

  path_length_uniq <- unique(path_length)

  if(verbose){
    cat(paste(length(pathway_list),"pathways are used for RWR","\n",sep = " "))
  }



  if(verbose){
    cat("Start perturbation \n")
  }

  cl <- makeCluster(parallel.cores)

  clusterExport(cl,'RWR')
  clusterExport(cl,'set.seed')
  set.seed(seed)
  seed_v <- sample(1:(nperm+500),nperm)

  pb <- txtProgressBar(style=3)
  rd_m <- list()
  for(i in 1:length(path_length_uniq)){
    setTxtProgressBar(pb, i/length(path_length_uniq))

    cd <- path_length_uniq[i]
    rdmatrix<-parLapply(cl,1:nperm,function(r,geneindex,nadjM,rwr.gamma,cellindex,cd,seed_v){
      set.seed(seed_v[r])
      samplei<-sample(geneindex,size=cd)
      names(samplei)<-row.names(nadjM)[samplei]

      resW_rd <- RWR(nadjM, samplei, gamma=rwr.gamma)
      return(resW_rd[cellindex])
    },geneindex,nadjM,rwr.gamma,cellindex,cd,seed_v)
    rd_m[[i]]<-do.call("rbind",rdmatrix)

    gc()
  }
  stopCluster(cl)
  rm(rdmatrix)

  if(verbose){
    cat("\nCalculate the scores \n")
  }

  pb <- txtProgressBar(style=3)
  saPW_matrix <- lapply(1:length(pathway_list),function(i){
    setTxtProgressBar(pb, i/length(pathway_list))

    index<-match(pathway_list[[i]],row.names(nadjM))
    names(index)<-pathway_list[[i]]
    resW <- RWR(nadjM, index, gamma=rwr.gamma)

    pp <- which(path_length_uniq==path_length[i])
    pvalue<-NULL
    for(j in 1:cell_n){
      ifelse(normal_dist==TRUE, pvalue[j]<-1-pnorm(resW[j],mean = mean(rd_m[[pp]][,j]), sd = sd(rd_m[[pp]][,j]),lower.tail = F),
             pvalue[j]<- 1-sum(rd_m[[pp]][,j]>=resW[j])/nperm)
    }
    return(pvalue)
  })
  saPW_matrix<-do.call("rbind",saPW_matrix)

  row.names(saPW_matrix)<-names(pathway_list)
  colnames(saPW_matrix)<-colnames(nadjM)[cellindex]
  saPW_matrix <- signif(saPW_matrix,2)
  return(saPW_matrix)
}
