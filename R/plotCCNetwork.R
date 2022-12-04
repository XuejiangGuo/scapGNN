##' plotCCNetwork
##'
##' @title Visualize cell cluster association network graph
##'
##' @description
##' The \code{plotCCNetwork} function takes cells belonging to the same phenotype as a cluster.
##' When cell phenotypes are not provided, the \code{plotCCNetwork} functions identify cell clusters based on edge betweenness.
##' Cell interactions between cell clusters are merged into one edge by mean.
##' The thickness of the edge indicates the strength of interaction between cell clusters.
##'
##' @param network.data The input network data is the result from the \code{ConNetGNN} function.
##' @param cell_id A vector of cell phenotype.
##' @param cell_cluster A binary value. Whether to automatically identify cell clusters based on edge betweenness. Default: \code{FALSE}.
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
##' @param layout Either a function or a numeric matrix. It specifies how the vertices will be placed on the plot. For details, please refer to the \code{igraph}Package. Default: \code{layout_with_lgl}.
##' @param legend.cex The font size of legend. Default: \code{1.5}.
##' @param legend.pt.cex Expansion factor(s) for the points. Default: \code{3}.
##' @param proportion This parameter specifies what percentage of edges to display (edges are sorted by their weight in descending order). Default: \code{1}, all edges are used.
##' @param plotgraph Whether to draw the picture. Default: \code{TRUE}. If \code{FALSE}, the image will not be displayed but the network data will be returned in the igraph data format.
##'
##' @return Graph or network data.
##'
##' @importFrom igraph cluster_edge_betweenness
##' @importFrom igraph graph_from_adjacency_matrix
##' @importFrom igraph layout_with_lgl
##' @importFrom igraph E
##' @importFrom igraph V
##' @importFrom graphics legend
##' @importFrom grDevices rainbow
##'
##' @export
##'
##' @examples
##' require(igraph)
##' require(graphics)
##'
##' data(ConNetGNN_data)
##'
##' # Construct the cell phenotype vector.
##' cell_id<-colnames(ConNetGNN_data[["cell_network"]])
##' temp<-unlist(strsplit(cell_id,"_"))
##' cell_phen<-temp[seq(2,length(temp)-1,by=3)]
##' names(cell_id)<-cell_phen
##' head(cell_id)
##' plotCCNetwork(ConNetGNN_data,cell_id,edge.width=10)


plotCCNetwork<-function(network.data,cell_id=NULL,cell_cluster=FALSE,vertex.colors=NULL,
                        vertex.size=10,vertex.label.cex=0.8,vertex.label.dist= 1,
                        vertex.label.color="black",edge.width=5,margin=0,layout=layout_with_lgl,
                        legend.cex=1.5,legend.pt.cex = 3,proportion=1,plotgraph=TRUE){
  if(!isLoaded("igraph")){
    stop("The package igraph is not available!")
  }
  
  c_net<-network.data[[1]]
  
  if(cell_cluster==TRUE){
    net<-graph_from_adjacency_matrix(c_net,mode="undirected",weighted=TRUE,diag=TRUE)
    net_c<-cluster_edge_betweenness(net,directed=FALSE)
    cluster_v<-net_c$names
    names(cluster_v)<-paste("cluster",net_c$membership,sep = "_")
    
    cluster_v<-cluster_v[match(cluster_v,colnames(c_net))]
    
    c_p<-names(table(names(cluster_v)))
    cc_amatrix<-NULL
    for(i in 1:length(c_p)){
      pp<-which(names(cluster_v)==c_p[i])
      if(length(pp)==1){
        cc_amatrix<-cbind(c_net[,pp],cc_amatrix)
      }else{
        cc_amatrix<-cbind(apply(c_net[,pp],1,mean),cc_amatrix)
      }
      
    }
    
    cc_amatrix1<-NULL
    for(i in 1:length(c_p)){
      pp<-which(names(cluster_v)==c_p[i])
      
      if(length(pp)==1){
        cc_amatrix1<-cbind(cc_amatrix[pp,],cc_amatrix1)
      }else{
        cc_amatrix1<-cbind(apply(cc_amatrix[pp,],2,mean),cc_amatrix1)
      }
      
      
    }
    
    diag(cc_amatrix1)<-0
    row.names(cc_amatrix1)<-c_p
    colnames(cc_amatrix1)<-c_p
    
    ccnp<-graph_from_adjacency_matrix(cc_amatrix1,mode="undirected",weighted=TRUE,diag=TRUE)
    
    if(is.null(cell_id)){
      if(plotgraph){
        plot(ccnp,vertex.size=vertex.size,vertex.label.cex=vertex.label.cex,vertex.label.dist= vertex.label.dist,
             vertex.label.color=vertex.label.color,edge.width=E(ccnp)$weight*edge.width,margin=margin,layout=layout)
      }else{
        return(ccnp)
      }
    }else{
      phen1<-names(table(names(cell_id)))
      if(is.null(vertex.colors)){
        vertex.colors<-rainbow(length(phen1))
      }
      names(vertex.colors)<-phen1
      
      
      n_v<-list()
      
      for(i in 1:length(c_p)){
        pp<-which(names(cluster_v)==c_p[i])
        cv<-cluster_v[pp]
        pp<-match(cv,cell_id)
        cp<-cell_id[pp]
        vg<-table(names(cp))
        v1<-rep(0,length(phen1))
        v1[match(names(vg),phen1)]<-vg
        n_v[[i]]<-v1/sum(v1)
      }
      names(n_v)<-c_p
      
      if(plotgraph){
        plot(ccnp, vertex.shape="pie", vertex.pie=n_v, vertex.pie.color=list(vertex.colors),
             vertex.size=vertex.size,vertex.label.cex=vertex.label.cex,vertex.label.dist= vertex.label.dist,
             vertex.label.color=vertex.label.color,edge.width=E(net)$weight*edge.width,margin=margin,layout=layout)
        legend("bottomright", legend=factor(phen1,levels = phen1), bty = "n", cex =legend.cex,
               pt.cex = legend.pt.cex, pch=20, col = vertex.colors , horiz = FALSE)
      }else{
        return(ccnp)
      }
    }
    
  }else{
    
    
    c_p<-names(table(names(cell_id)))
    
    if(is.null(vertex.colors)){
      vertex.colors<-rainbow(length(c_p))
      names(vertex.colors)<-c_p
    }
    
    vertex.colors<-vertex.colors[match(names(vertex.colors),c_p)]
    
    cc_amatrix<-NULL
    for(i in 1:length(c_p)){
      pp<-which(names(cell_id)==c_p[i])
      cc_amatrix<-cbind(apply(c_net[,pp],1,mean),cc_amatrix)
    }
    
    cc_amatrix1<-NULL
    for(i in 1:length(c_p)){
      pp<-which(names(cell_id)==c_p[i])
      cc_amatrix1<-cbind(apply(cc_amatrix[pp,],2,mean),cc_amatrix1)
    }
    
    diag(cc_amatrix1)<-0
    row.names(cc_amatrix1)<-c_p
    colnames(cc_amatrix1)<-c_p
    
    vg<-unique(as.vector(cc_amatrix1))
    vg<-sort(vg,decreasing = T)
    a<-vg[floor(length(vg)*proportion)]
    cc_amatrix1[cc_amatrix1<=a]<-0
    ccnp<-graph_from_adjacency_matrix(cc_amatrix1,mode="undirected",weighted=TRUE,diag=TRUE)
    
    
    if(plotgraph){
      plot(ccnp,vertex.size=vertex.size,vertex.label.cex=vertex.label.cex,vertex.label.dist= vertex.label.dist,
           vertex.label.color=vertex.label.color,edge.width=E(ccnp)$weight*edge.width,margin=margin,layout=layout,vertex.color=vertex.colors)
    }else{
      return(ccnp)
    }
  }
  
}