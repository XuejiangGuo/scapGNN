##' ConNetGNN
##'
##' @title Construct association networks for gene-gene, cell-cell, and gene-cell based on graph neural network (GNN)
##'
##' @description This function implements a graph neural network with two autoencoders. 1. AutoEncoder (AE) based on deep neural network:
##' Infer latent associations between genes and cells. 2. Graph AutoEncoder (GAE) based on graph convolutional neural network: Construct
##' association networks for gene-gene, cell-cell.
##'
##' @param Prep_data The input data is the result from the \code{Preprocessing} function.
##' @param python.path The path to a Python binary. If python.path="default", the program will use the current system path to python.
##' @param miniconda.path The path in which miniconda will be installed. If the \code{python.path} is NULL and conda or miniconda is not installed in the system, the program will automatically install miniconda according to the path specified by \code{miniconda.path}.
##' @param AE.epochs The number of epoch for the deep neural network (AE). Default: \code{1000}.
##' @param AE.learning.rate Initial learning rate of AE. Default: \code{0.001}.
##' @param AE.reg.alpha The LTMG regularized intensity. Default: \code{0.5}.
##' @param use.VGAE Whether to use Variational Graph AutoEncoder (VGAE). Default: \code{TRUE}.
##' @param GAE.epochs The number of epoch for the GAE. Default: \code{300}.
##' @param GAE.learning.rate Initial learning rate of GAE. Default: \code{0.01}.
##' @param cell_val_ratio For GAE that construct cell-cell association networks, the proportion of edges that are extracted as the validation set. Default: \code{0.05}.
##' @param gene_val_ratio As with parameter \code{cell_val_ratio}, it is simply applied with the construction of gene-gene association networks.
##' @param parallel Whether to use multiple processors to run GAE. Default: \code{FALSE} When \code{parallel=TRUE} (default), tow processors will be used to run GAE.
##' @param seed Random number generator seed.
##' @param verbose Gives information about each step. Default: \code{TRUE}.
##' @param GPU.use Whether to use GPU for GNN modules. Default: \code{FALSE}. If GPU.use=TRUE, CUDA needs to be installed.
##'
##' @details
##' The \code{ConNetGNN} function establishes a graph neural network (GNN) framework to mine latent relationships between genes and cells and within themselves.
##' This framework mainly includes two capabilities: \itemize{
##' \item 1.Deep neural network-based AutoEncoder inferring associations between genes and cells and generating gene features and cell features for the GAE.
##' \item 2.The GAE takes the gene feature and cell feature as the node features of the initial gene correlation network and cell correlation network,
##' and constructs the gene association network and cell association network through the graph convolution process.
##' }
##'
##' The GNN is implemented based on \code{pytorch}, so an appropriate python environment is required:
##' \itemize{
##' \item python >=3.9.7
##' \item pytorch >=1.10.0
##' \item sklearn >=0.0
##' \item scipy >=1.7.3
##' \item numpy >=1.19.5
##' }
##'
##' If the user has already configured the python environment, the path of the python binary file can be directly entered into \code{python.path}.
##' If the parameter \code{python.path} is NULL, the program will build a miniconda environment called \code{scapGNN_env} and configure python.
##' We also provide environment files for conda: \code{/inst/extdata/scapGNN_env.yaml}. Users can install it with the command: \code{conda env create -f scapGNN_env.yaml}.
##'
##' @return A list:
##' \describe{
##'   \item{cell_network}{Constructed cell association network.}
##'   \item{gene_network}{Constructed gene association network.}
##'   \item{cell_gene_network}{Constructed gene-cell association network.}
##' }
##'
##' @importFrom reticulate conda_version
##' @importFrom reticulate miniconda_path
##' @importFrom reticulate install_miniconda
##' @importFrom reticulate use_python
##' @importFrom reticulate py_module_available
##' @importFrom reticulate py_config
##' @importFrom reticulate conda_list
##' @importFrom reticulate py_available
##' @importFrom reticulate source_python
##' @importFrom parallel makeCluster
##' @importFrom parallel clusterEvalQ
##' @importFrom parallel parLapply
##' @importFrom parallel stopCluster
##'
##' @export
##'
##' @examples
##' require(coop)
##' require(reticulate)
##' require(parallel)
##' # Data preprocessing
##' data("Hv_exp")
##' Prep_data <- Preprocessing(Hv_exp[1:300,])
##' \dontrun{
##' # Specify the python path
##' ConNetGNN_data <- ConNetGNN(Prep_data,python.path="../miniconda3/envs/scapGNN_env/python.exe")
##' }
##'


ConNetGNN<-function(Prep_data,python.path=NULL,miniconda.path = NULL,AE.epochs=1000,AE.learning.rate=0.001,AE.reg.alpha=0.5,use.VGAE=TRUE,
                   GAE.epochs = 300,GAE.learning.rate = 0.01,cell_val_ratio=0.05, gene_val_ratio=0.05,parallel=FALSE,seed=125,GPU.use=FALSE,verbose=TRUE){
  if(!isLoaded("reticulate")){
    stop("The package reticulate is not available!")
  }

  if(parallel){
    if(!isLoaded("parallel")){
      stop("The package parallel is not available!")
    }
  }

  GAE_function<-NULL
  AE_function<-NULL


  ####python
  if(is.null(python.path)){
    condav<-conda_version()
    if(is.character(condav)==TRUE){
      cat(paste(condav,"is available!  \n",sep = " "))

      conda_env<-conda_list()
      p<-which(conda_env[,1]=="scapGNN_env")
      if(length(p)>0){
        cat("The environment scapGNN_env already exists!  \n")
        python.path<-conda_env[p,2]
      }else{
        #creat scapGNN_env
        cat("The environment scapGNN_env not found, will be created!  \n")
        create_scapGNN_env()
        conda_env<-conda_list()
        p2<-which(conda_env[,1]=="scapGNN_env")
        python.path<-conda_env[p2,2]
      }
    }else{
      #install conda
      cat("No conda or miniconda detected, miniconda will be created through the reticulate R package!  \n")
      if (is.null(miniconda.path)) {
        miniconda.path <- reticulate::miniconda_path()
      }
      status <- tryCatch(
        reticulate::install_miniconda(path = miniconda.path),
        error = function(e) {
          return(TRUE)
        }
      )
      if (isTRUE(status)) {
        stop(
          "Error during the installation of miniconda. Please see the website of the ",
          "miniconda for more details",
          call. = FALSE
        )
      }
      cat("Create the environment scapGNN_env!  \n")
      create_scapGNN_env()
      conda_env<-conda_list()
      p2<-which(conda_env[,1]=="scapGNN_env")
      python.path<-conda_env[p2,2]
    }
  }

  if(!py_available()){
    if(python.path=="default"){
      python.path <- Sys.which("python")
    }
    use_python(python.path, required = T)
  }


  if(!py_module_available("torch")){
    instPyModule("pytorch")
  }
  if(!py_module_available("sklearn")){
    instPyModule("sklearn")
  }
  if(!py_module_available("scipy")){
    instPyModule("scipy")
  }


  HVexp<-as.matrix(Prep_data[[1]])
  row.names(HVexp)<-NULL
  colnames(HVexp)<-NULL

  cell_features<-as.matrix(Prep_data[[2]])
  row.names(cell_features)<-NULL
  colnames(cell_features)<-NULL

  gene_features<-as.matrix(Prep_data[[3]])
  row.names(gene_features)<-NULL
  colnames(gene_features)<-NULL

  LTMG<-as.matrix(Prep_data[[4]])
  row.names(LTMG)<-NULL
  colnames(LTMG)<-NULL


  cell_adj<-as.matrix(Prep_data[[5]])
  row.names(cell_adj)<-NULL
  colnames(cell_adj)<-NULL

  gene_adj<-as.matrix(Prep_data[[6]])
  row.names(gene_adj)<-NULL
  colnames(gene_adj)<-NULL

  orig_adj<-list(cell_adj,gene_adj)


  if (verbose) {
    cat("Run AutoEncoder  \n")
  }

  if(GPU.use==TRUE){
	source_python(system.file("python", "AutoEncoder_GPU.py", package = "scapGNN"))
  }else{
	source_python(system.file("python", "AutoEncoder.py", package = "scapGNN"))
  }

  AE_data<-AE_function(cell_features=cell_features,
                 gene_features=gene_features,
                 exp=HVexp,ltmg_m=LTMG,DNN_epochs=AE.epochs,
                 DNN_learning_rate=AE.learning.rate,
                 reg_alpha=AE.reg.alpha,seed=seed)

  if (verbose) {
    cat(paste("loss:",AE_data[[4]],"\n",sep = " "))
  }

  if (verbose) {
    cat("Run Graph AutoEncoder  \n")
  }

  if(parallel==TRUE){
    ncores<-2
    cl <- makeCluster(ncores)
    clusterEvalQ(cl,library(reticulate))
    clusterEvalQ(cl,library(scapGNN))
    rt<-c(cell_val_ratio,gene_val_ratio)
    GAE_data<-parLapply(cl,1:2,function(i,AE_data,orig_adj,use.VGAE,GAE.epochs,GAE.learning.rate,seed,python.path,rt){
      use_python(python.path, required = T)

	  if(GPU.use==TRUE){
		source_python(system.file("python", "GraphAutoEncoder_GPU.py", package = "scapGNN"))
	  }else{
		source_python(system.file("python", "GraphAutoEncoder.py", package = "scapGNN"))
	  }

      res<-GAE_function(net_m=orig_adj[[i]],feature_m=AE_data[[i+1]],
                        use_model=use.VGAE,GAE_epochs=GAE.epochs,
                      GAE_learning_rate=GAE.learning.rate,seed=seed,ratio_val=rt[i])
      return(res)
    },AE_data,orig_adj,use.VGAE,GAE.epochs,GAE.learning.rate,seed,python.path,rt)
    stopCluster(cl)
  }else{
	if(GPU.use==TRUE){
		source_python(system.file("python", "GraphAutoEncoder_GPU.py", package = "scapGNN"))
	}else{
		source_python(system.file("python", "GraphAutoEncoder.py", package = "scapGNN"))
	}

    if (verbose) {
      cat("Construct cell-cell association network  \n")
    }
    cell_gae<-GAE_function(net_m=orig_adj[[1]],feature_m=AE_data[[2]],
                      use_model=use.VGAE,GAE_epochs=GAE.epochs,
                      GAE_learning_rate=GAE.learning.rate,seed=seed,ratio_val=cell_val_ratio)
    if (verbose) {
      cat("Construct gene-gene association network \n")
    }
    gene_gae<-GAE_function(net_m=orig_adj[[2]],feature_m=AE_data[[3]],
                      use_model=use.VGAE,GAE_epochs=GAE.epochs,
                      GAE_learning_rate=GAE.learning.rate,seed=seed,ratio_val=gene_val_ratio)
    GAE_data<-list(cell_gae,gene_gae)
  }

  if (verbose) {
    cat(paste("loss:",GAE_data[[1]][[2]],"\n",sep = " "))
    cat(paste("loss:",GAE_data[[2]][[2]],"\n",sep = " "))
  }

  cell_net<-GAE_data[[1]][[1]]
  gene_net<-GAE_data[[2]][[1]]
  cg_net<-AE_data[[1]]

  colnames(cg_net)<-colnames(Prep_data[[1]])
  row.names(cg_net)<-row.names(Prep_data[[1]])

  colnames(cell_net)<-colnames(Prep_data[[1]])
  row.names(cell_net)<-colnames(Prep_data[[1]])

  colnames(gene_net)<-row.names(Prep_data[[1]])
  row.names(gene_net)<-row.names(Prep_data[[1]])

  if (verbose) {
    cat("Done  \n")
  }

  results<-list(cell_network=cell_net,gene_network=gene_net,gene_cell_network=cg_net)
  return(results)
}
