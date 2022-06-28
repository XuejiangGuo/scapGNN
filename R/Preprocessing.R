##' Preprocessing
##'
##' @title Data preprocessing
##'
##' @description This function is to prepare data for the \code{ConNetGNN} function.
##'
##' @param data The input data should be a data frame or a matrix where the rows are genes and the columns are cells. The \code{seurat} object are also accepted.
##' @param verbose Gives information about each step. Default: \code{TRUE}.
##'
##' @details
##' The function is able to interface with the \code{seurat} framework. The process of \code{seurat} data processing refers to \code{Examples}.
##' The input data should be containing hypervariable genes and log-transformed. Left-truncated mixed Gaussian (LTMG) modeling to calculate gene
##' regulatory signal matrix. Positively correlated gene-gene and cell-cell are used as the initial gene correlation matrix and cell correlation matrix.
##'
##' @return A list:
##' \describe{
##'   \item{orig_dara}{User-submitted raw data, rows are highly variable genes and columns are cells.}
##'   \item{cell_features}{Cell feature matrix.}
##'   \item{gene_features}{Gene feature matrix.}
##'   \item{ltmg_matrix}{Gene regulatory signal matrix for LTMG.}
##'   \item{cell_adj}{The adjacency matrix of the cell correlation network.}
##'   \item{gene_adj}{The adjacency matrix of the gene correlation network.}
##' }
##'
##' @importFrom coop tpcor
##' @importFrom methods new
##'
##' @export
##'
##' @examples
##' \dontrun{
##' # Load dependent packages.
##' require(coop)
##'
##' # Seurat data processing.
##' require(Seurat)
##'
##' # Load the PBMC dataset (Case data for seurat)
##' pbmc.data <- Read10X(data.dir = "../data/pbmc3k/filtered_gene_bc_matrices/hg19/")
##'
##' # Our recommended data filtering is that only genes expressed as non-zero in more than
##' # 1% of cells, and cells expressed as non-zero in more than 1% of genes are kept.
##' # In addition, users can also filter mitochondrial genes according to their own needs.
##' pbmc <- CreateSeuratObject(counts = pbmc.data, project = "case",
##'                                     min.cells = 3, min.features = 200)
##' pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
##' pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
##'
##' # Normalizing the data.
##' pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")
##'
##' # Identification of highly variable features.
##' pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)
##'
##' # Run Preprocessing.
##' Prep_data <- Preprocessing(pbmc)
##'
##' }
##'
##' # Users can also directly input data
##' # in data frame or matrix format
##' # containing highly variable genes.
##' data("Hv_exp")
##' Prep_data <- Preprocessing(Hv_exp[1:300,])

Preprocessing<-function(data,verbose=TRUE){
  if(!isLoaded("coop")){
    stop("The package coop is not available!")
  }

  if(sum(grep("Seurat",class(data)))!=0){
    exp<-data@assays[["RNA"]]@data
    hvgi<-match(data@assays[["RNA"]]@var.features,row.names(exp))
    exp<-as.matrix(exp[hvgi,])
  }else{
    exp<-as.matrix(data)
  }


  cell_features<-t(apply(exp, 2, function(x){
    return(x/max(x))
  }))

  gene_features<-t(apply(exp, 1, function(x){
    return(x/max(x))
  }))

  if (verbose) {
    cat("Run LTMG \n")
  }
  LTMG_Object<- new(Class = 'LTMG', InputData =  exp)

  if(verbose){
    suppressMessages(LTMG <- RunLTMG(LTMG_Object,Gene_use = 'all', verbose=verbose)@OrdinalMatrix)
  }

  if (verbose) {
    cat("Calculate cell correlation matrix \n")
  }
  cellnet <- coop::tpcor(cell_features)
  cellnet[cellnet<0]<-0
  diag(cellnet)<-0

  cellnet_p<-t(apply(cellnet,1,function(x){
    len<-length(which(x>0))
    x1<-rep(0,length(x))

    for(j in 1:length(x)){
      if(x[j]>0){
        p<-length(which(x>x[j]))/len
        if(p==0){
          x1[j]<-0.00001
        }else{
          x1[j]<-p
        }
      }
    }
    x1[x1>0.05]<-0
    return(x1)
  }))

  x<-cellnet_p
  x[upper.tri(x)] <- 0
  s<-cellnet_p
  s[lower.tri(s)] <- 0
  s<-t(s)
  s[s!=0]<-1
  x[x!=0]<-1
  cellnet_p<-s+x
  cellnet_p[cellnet_p==2]<-1
  cellnet_p<-cellnet_p+t(cellnet_p)
  cellnet<-cellnet*cellnet_p

  if (verbose) {
    cat("Calculate gene correlation matrix \n")
  }
  genenet <- coop::tpcor(gene_features)
  genenet[genenet<0]<-0
  diag(genenet)<-0

  genenet_p<-t(apply(genenet,1,function(x){
    len<-length(which(x>0))
    x1<-rep(0,length(x))

    for(j in 1:length(x)){
      if(x[j]>0){
        p<-length(which(x>x[j]))/len
        if(p==0){
          x1[j]<-0.00001
        }else{
          x1[j]<-p
        }
      }
    }
    x1[x1>0.05]<-0
    return(x1)
  }))

  x<-genenet_p
  x[upper.tri(x)] <- 0
  s<-genenet_p
  s[lower.tri(s)] <- 0
  s<-t(s)
  s[s!=0]<-1
  x[x!=0]<-1
  genenet_p<-s+x
  genenet_p[genenet_p==2]<-1
  genenet_p<-genenet_p+t(genenet_p)
  genenet<-genenet*genenet_p

  if (verbose) {
    cat("Done  \n")
  }

  results<-list(orig_dara=exp,cell_features=cell_features,gene_features=gene_features,
                ltmg_matrix=LTMG,cell_adj=cellnet,gene_adj=genenet)

  return(results)
}
