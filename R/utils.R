##' isLoaded
##'
##' @title The internal functions of the \code{scapGNN} package
##' @description Determine if the package is loaded.
##' @param name Internal parameters.
isLoaded <- function(name) {
  (paste("package:", name, sep="") %in% search()) ||
    (name %in% loadedNamespaces())
}

##' create_scapGNN_env
##'
##' @title Create the create_scapGNN_env environment on miniconda
##' @description The internal functions of the \code{scapGNN} package.
create_scapGNN_env<-function(){
  status <- tryCatch(
    reticulate::conda_create(
      envname = "scapGNN_env",
      packages = "python==3.9.7"
    ),
    error = function(e) {
      return(TRUE)
    }
  )
  if (isTRUE(status)) {
    stop(
      "Error during the creation of conda environment. Please see the website of the ",
      "miniconda for more details",
      call. = FALSE
    )
  }
}

##' instPyModule
##'
##' @title Install the pyhton module through the reticulate R package
##' @description The internal functions of the \code{scapGNN} package.
##' @param module Internal parameters.
##' @importFrom reticulate py_install
instPyModule<-function(module){
  status<-tryCatch(py_install(module,method="conda"),error = function(e) {
    return(TRUE)
  })
  if (isTRUE(status)) {
    stop(
      paste("Error during the installation. Please see the website of the python module:",module,
            "for more details",sep = " "),
      call. = FALSE
    )
  }
}


##' load_path_data
##' @title load pathway or gene set's gmt file
##' @description
##' The internal functions of the \code{scapGNN} package.
##'
##' file format:
##' 1. first index: pathway's name or ID.
##' 2. second index: pathway's url or others, it dosen't matter.
##' 3. third to all: gene symbols in pathway.
##'
##' @param gmt_file_path Internal parameters.
##' @export
##' @return a list
load_path_data = function(gmt_file_path){
  tmp = readLines(gmt_file_path)
  gsets = list()
  for(i in 1:length(tmp)){
    t = strsplit(tmp[i],'\t')[[1]]
    genes = t[3:length(t)]
    genes = genes[which(genes != "")]
    gsets[[t[1]]] = genes
  }
  return (gsets)
}

##' Function that performs a random Walk with restart (RWR) on a given graph
##' @param W : adjacency matrix of the graph
##' @param ind.positives : indices of the "core" positive examples of the graph. They represent to the indices of W corresponding to the positive examples
##' @param gamma : restart parameter (def: 0.6)
##'
##' @return a list with three elements:
##' - p : the probability at the steady state
##' - ind.positives : indices of the "core" positive examples of the graph (it is equal to the same
##'                  input parameter
##' - n.iter : number of performed iterations
##' @export
##' @return a vector

RWR <- function(W, ind.positives, gamma=0.6) {
  n <- nrow(W);
  names.var <- rownames(W);
  diag.D <- apply(W,1,sum);
  diag.D[diag.D==0] <- Inf;
  inv.diag.D <- 1/diag.D;

  M <-W*inv.diag.D # M = D^-1 * W

  n <- nrow(M)
  p0 <- p <- numeric(n)
  names(p) <- names(p0) <- rownames(W)
  rm(W)
  n.positives <- length(ind.positives)
  p0[ind.positives] <- 1/n.positives

  p <- p0;
  for (t in 1:5000) {
    pold <- p
    p <- ((1-gamma) * as.vector(pold %*% M)) + gamma * p0
    if (sum(abs(p-pold)) < 1e-6){
      break()
    }
  }
  return(p)
}
