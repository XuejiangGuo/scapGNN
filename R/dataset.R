#' Results of ConNetGNN() for scATAC-seq data from SNARE-seq dataset
#'
#' @description   A list to store the gene association network of scATAC-seq data.
#' Case data from the SNARE-seq dataset.
#'
#' @format a list of three adjacency matrices.
#' @examples
#' data(ATAC_net)
"ATAC_net"


#' Results of ConNetGNN() for scRNA-seq data from SNARE-seq dataset
#'
#' @description   A list to store the gene association network of scRNA-seq data.
#' Case data from the SNARE-seq dataset.
#'
#' @format a list of three adjacency matrices.
#' @examples
#' data(RNA_net)
"RNA_net"


#' Results of InteNet() for integrating scRNA-seq and scATAC-seq data.
#'
#' @description  An integrated network of scRNA-seq and scATAC-seq data from SNARE-seq.
#'
#' @format a list of three adjacency matrices.
#' @examples
#' data(RNA_ATAC_IntNet)
"RNA_ATAC_IntNet"


#' Single-cell gene expression profiles
#'
#' @description   A log-transformed gene-cell matrix containing highly variable features.
#'
#' @format a matrix.
#' @examples
#' data(Hv_exp)
"Hv_exp"


#' The results of ConNetGNN() function
#'
#' @description   Results of ConNetGNN() function with Hv_exp as input.
#'
#' @format a list.
#' @examples
#' data(ConNetGNN_data)
"ConNetGNN_data"


#' Single cell pathway activity matrix
#'
#' @description   Results of scPathway() function.
#'
#' @format a matrix.
#' @examples
#' data(scPathway_data)
"scPathway_data"


#' Cell-activated gene modules under the 0-hour phenotype
#'
#' @description   Results of cpGModule() function.
#'
#' @format a list.
#' @examples
#' data(H9_0h_cpGM_data)
"H9_0h_cpGM_data"


#' Cell-activated gene modules under the 24-hour phenotype
#'
#' @description   Results of cpGModule() function.
#'
#' @format a list.
#' @examples
#' data(H9_24h_cpGM_data)
"H9_24h_cpGM_data"


#' Cell-activated gene modules under the 36-hour phenotype
#'
#' @description   Results of cpGModule() function.
#'
#' @format a list.
#' @examples
#' data(H9_36h_cpGM_data)
"H9_36h_cpGM_data"
