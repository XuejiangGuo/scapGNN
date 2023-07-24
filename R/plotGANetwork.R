##' plotGANetwork
##'
##' @title Visualize gene association network graph of a gene module or pathway at the specified cell phenotype
##'
##' @description
##' Based on the gene set input by the user, \code{plotGANetwork} functional draws the gene association network in the specified cell phenotype.
##' The node size in the network reflects the activation strength of the gene. The thickness of the edge indicates the strength of interaction between genes.
##'
##' @param network.data Network data constructed by the \code{ConNetGNN} function.
##' @param cellset A vector of cell id. A cell set corresponding to the specified cell phenotype.
##' @param geneset A vector of gene id. A gene module or pathway.
##' @param rwr.gamma Restart parameter. Default: \code{0.7}.
##' @param vertex.colors The fill color of the vertex. The number of colors should match the number of cell phenotypes. If \code{NULL (default)}, the system will automatically assign colors.
##' @param vertex.size The size of the vertex. Default: \code{10}.
##' @param vertex.label.cex The font size for vertex labels. Default: \code{0.8}.
##' @param vertex.label.dist The distance of the label from the center of the vertex. If it is 0 then the label is centered on the vertex. Default: \code{1}.
##' @param vertex.label.color The color of the labels. Default: \code{black}.
##' @param edge.width The width of the edge. This does not affect the relative size of the edge weights. Default: \code{5}.
##' @param margin The amount of empty space below, over, at the left and right of the plot, it is a numeric
##' vector of length four. Usually values between 0 and 0.5 are meaningful, but negative values
##' are also possible, that will make the plot zoom in to a part of the graph. If it is shorter than
##' four then it is recycled. Default: \code{0}.
##' @param layout Either a function or a numeric matrix. It specifies how the vertices will be placed on the plot. For details, please refer to the \code{igraph}Package. Default: \code{layout_as_star}.
##' @param main A main title for the plot.
##' @param plotgraph Whether to draw the picture. Default: \code{TRUE}. If \code{FALSE}, the image will not be displayed but the network data will be returned in the igraph data format.
##'
##' @return A graph or list.
##'
##' @importFrom igraph graph_from_adjacency_matrix
##' @importFrom igraph induced_subgraph
##' @importFrom igraph layout_as_star
##' @importFrom igraph E
##' @importFrom stats na.omit
##'
##' @export
##'
##' @examples
##' require(igraph)
##'
##' # Load the result of the ConNetGNN function.
##' data(ConNetGNN_data)
##'
##' data("Hv_exp")
##' index<-grep("0h",colnames(Hv_exp))
##' cellset<-colnames(Hv_exp)[index]
##' pathways<-load_path_data(system.file("extdata", "KEGG_human.gmt", package = "scapGNN"))
##' geneset<-pathways[[which(names(pathways)=="Tight junction [PATH:hsa04530]")]]
##' plotGANetwork(ConNetGNN_data,cellset,geneset,main = "Tight junction [PATH:hsa04530]")
##'

plotGANetwork<-function(network.data,cellset,geneset,rwr.gamma=0.7,vertex.colors=NULL,vertex.size=10,
                       vertex.label.cex=0.8,vertex.label.dist= 1,vertex.label.color="black",
                       edge.width=5,margin=0,layout=layout_as_star,main=NULL,plotgraph=TRUE){
  if(!isLoaded("igraph")){
    stop("The package igraph is not available!")
  }

  cg_net<-network.data[[3]]
  cg_net<-apply(cg_net,2,function(x){
    return(x/max(x))
  })

  c_net<-network.data[[1]]
  diag(c_net)<-0

  g_net<-network.data[[2]]
  diag(g_net)<-0

  merge1<-cbind(c_net,t(cg_net))
  merge2<-cbind(cg_net,g_net)
  cell_gene_network<-rbind(merge1,merge2)

  pp<-match(cellset,row.names(cell_gene_network))
  names(pp)<-cellset
  pp<-na.omit(pp)

  resW <- RWR(cell_gene_network, pp,gamma=rwr.gamma)
  pp<-match(geneset,row.names(cell_gene_network))
  pp<-na.omit(pp)
  geneW<-resW[pp]

  net<-graph_from_adjacency_matrix(cell_gene_network,mode="undirected",weighted=TRUE,diag=TRUE)

  subnet<-induced_subgraph(net,names(geneW))

  if(plotgraph){
    plot(subnet,vertex.size=geneW*2000*vertex.size,vertex.label.cex=vertex.label.cex,vertex.label.dist= vertex.label.dist,
         vertex.label.color=vertex.label.color,vertex.color=vertex.colors,
         edge.width=E(net)$weight*edge.width,margin=margin,layout=layout,main=main)
  }else{
    return(list(graph=subnet,Node_weight=geneW))
  }

}
