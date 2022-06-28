##' plotMulPhenGM
##'
##' @title Visualize gene association network graph for activated gene modules under multiple cell phenotypes
##'
##' @description For multiple cell phenotypes, the \code{plotMulPhenGM} function will display the activated gene modules for each phenotype and show the connection and status of genes in different cell phenotypes.
##'
##' @param data.list a list. Each element represents the \code{cpGModule} function result of a cell phenotype and the names of the lists are the corresponding cell phenotype.
##' @param network.data Network data constructed by the \code{ConNetGNN} function.
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
##' @param layout Either a function or a numeric matrix. It specifies how the vertices will be placed on the plot. For details, please refer to the \code{igraph} Package. Default: \code{layout_with_lgl}.
##' @param legend.position This places the legend on the inside of the plot frame at the given location. See the \code{legend()} function for details.
##' @param legend.cex The font size of legend. Default: \code{1.5}.
##' @param legend.pt.cex Expansion factor(s) for the points. Default: \code{3}.
##' @param plotgraph Whether to draw the picture. Default: \code{TRUE}. If \code{FALSE}, the image will not be displayed but the network data will be returned in the igraph data format.
##'
##' @details
##' If a gene is significantly activated in more than one cell phenotype, we call it a co-activated gene. These co-activated genes are shown on the sector diagram.
##' Each interval of the sector diagram represents the activation strength of the gene in this cell phenotype relative to other cell phenotypes.
##'
##' @return A graph or list.
##'
##' @importFrom igraph graph_from_adjacency_matrix
##' @importFrom igraph as_adjacency_matrix
##' @importFrom igraph layout_with_lgl
##' @importFrom igraph graph.union
##' @importFrom igraph E
##' @importFrom graphics legend
##' @importFrom grDevices rainbow
##' @importFrom stats na.omit
##'
##' @export
##'
##' @examples
##' require(igraph)
##' require(grDevices)
##' # Load the result of the ConNetGNN function.
##' data(ConNetGNN_data)
##' # Obtain cpGModule results for each cell phenotype.
##' data(H9_0h_cpGM_data)
##' data(H9_24h_cpGM_data)
##' data(H9_36h_cpGM_data)
##' data.list<-list(H9_0h=H9_0h_cpGM_data,H9_24h=H9_24h_cpGM_data,H9_36h=H9_36h_cpGM_data)
##' plotMulPhenGM(data.list,ConNetGNN_data)
##'

plotMulPhenGM<-function(data.list,network.data,vertex.colors=NULL,vertex.size=10,vertex.label.cex=0.8,
                        vertex.label.dist= 1,vertex.label.color="black",edge.width=5,
                        margin=0,layout=layout_with_lgl,legend.position="bottomright",legend.cex=1.5,legend.pt.cex = 3,plotgraph=TRUE){
  if(!isLoaded("igraph")){
    stop("The package igraph is not available!")
  }

  if(!isLoaded("graphics")){
    stop("The package graphics is not available!")
  }

  g_net<-network.data[[2]]
  diag(g_net)<-0

  qc<-NULL
  for(i in 1:length(data.list)){
    if(nrow(data.list[[i]])==0){
      qc<-c(qc,i)
    }
  }
  if(is.null(qc)==FALSE){
    data.list<-data.list[-qc]
  }


  if(is.null(vertex.colors)){
    vertex.colors<-rainbow(length(data.list))
  }


  ig<-NULL
  qc<-NULL
  ig_p<-NULL
  for(i in 1:length(data.list)){
    if(length(data.list[[i]]$Genes)==1){
      ig<-c(ig,data.list[[i]]$Genes)
      qc<-c(qc,i)
      ig_p<-names(data.list)[i]
    }
  }
  if(is.null(qc)==FALSE){
    data.list<-data.list[-qc]
  }

  allgenes<-NULL
  allsubg<-list()
  for(i in 1:length(data.list)){
    pp<-match(data.list[[i]]$Genes,row.names(g_net))
    subg<-g_net[pp,pp]

    allsubg[[i]]<-graph_from_adjacency_matrix(subg,mode="undirected",weighted=TRUE,diag=TRUE)
    allgenes<-c(allgenes,data.list[[i]]$Genes)

  }
  ggraph<-graph.union(allsubg)

  alladj<-as_adjacency_matrix(ggraph)
  for(i in 1:length(data.list)){
    alladj1<-as_adjacency_matrix(ggraph,attr=paste("weight",i,sep = "_"),sparse=FALSE)
    alladj1[is.na(alladj1)]<-1
    alladj<-alladj*alladj1
  }

  allgenes<-unique(allgenes)

  node.atlist<-list()
  for(i in 1:length(allgenes)){
    gene<-allgenes[i]
    values<-NULL
    for(j in 1:length(data.list)){
      pp<-which(data.list[[j]]$Genes==gene)
      if(length(pp)>0){
        values<-c(values,data.list[[j]]$AS[pp])
      }else{
        values<-c(values,0)
      }
    }
    values<-values/sum(values)
    names(values)<-names(data.list)
    node.atlist[[i]]<-values
  }
  names(node.atlist)<-allgenes

  if(length(ig)>=1){
    for(j in 1:length(ig)){
      pp<-which(names(node.atlist)==ig[j])
      if(pp>0){
        temp<-0
        names(temp)<-ig_p[j]

        for(i in 1:length(node.atlist)){
          if(i==pp){
            temp1<-1
            names(temp1)<-ig_p[j]
            node.atlist[[i]]<-c(node.atlist[[i]],temp1)
          }
          node.atlist[[i]]<-c(node.atlist[[i]],temp)
        }
      }else{
        temp<-0
        names(temp)<-ig_p[j]
        for(i in 1:length(node.atlist)){
          node.atlist[[i]]<-c(node.atlist[[i]],temp)
        }

        temp1<-rep(0,length(node.atlist[[1]]))
        names(temp1)<-names(node.atlist[[1]])
        temp1[length(temp1)]<-1
        node.atlist[[i+1]]<-temp1

        alladj<-cbind(alladj,rep(0,nrow(alladj)))
        alladj<-rbind(alladj,rep(0,ncol(alladj)))
        row.names(alladj)[nrow(alladj)]<-ig[j]
        colnames(alladj)[ncol(alladj)]<-ig[j]
      }

    }

  }

  net<-graph_from_adjacency_matrix(as.matrix(alladj),mode="undirected",weighted=TRUE,diag=TRUE)

  phen<-NULL
  if(length(names(vertex.colors))>0){
    for(i in 1:length(node.atlist)){
      pp<-match(names(vertex.colors),names(node.atlist[[i]]))
	  qc<-which(is.na(pp)==TRUE)
      pp<-pp[-qc]
	  vertex.colors<-vertex.colors[-qc]
      node.atlist[[i]]<-node.atlist[[i]][pp]
	  phen<-c(phen,names(node.atlist[[i]]))
    }
  }
  phen<-unique(phen)
  pp<-match(names(vertex.colors),phen)
  qc<-which(is.na(pp)==TRUE)
  if(length(qc)>0){
    vertex.colors<-vertex.colors[-qc]
  }

  if(plotgraph){
    plot(ggraph, vertex.shape="pie", vertex.pie=node.atlist, vertex.pie.color=list(vertex.colors),
         vertex.size=vertex.size,vertex.label.cex=vertex.label.cex,vertex.label.dist= vertex.label.dist,
         vertex.label.color=vertex.label.color,
         edge.width=E(net)$weight*edge.width,margin=margin,layout=layout)
    legend(legend.position, legend=factor(names(node.atlist[[1]]),levels = names(node.atlist[[1]])), bty = "n", cex =legend.cex,
           pt.cex = legend.pt.cex, pch=20, col = vertex.colors , horiz = FALSE)
  }else{
    return(list(graph=net,Node_weight=node.atlist))
  }
}
