## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(scapGNN)

## ----echo = F, out.width = "100%"---------------------------------------------
knitr::include_graphics("../inst/extdata/flow_diagram1.png")

## ----echo = F, out.width = "100%"---------------------------------------------
knitr::include_graphics("../inst/extdata/flow_diagram2.png")

## -----------------------------------------------------------------------------
# Users can also directly input data in a data frame or matrix format which contains hypervariable genes and is log-transformed.
data("Hv_exp")
Prep_data <- Preprocessing(Hv_exp,verbose=FALSE)
summary(Prep_data)

## -----------------------------------------------------------------------------
# View the content of the ConNetGNN() results.
data(ConNetGNN_data)
summary(ConNetGNN_data)

## ----fig.width = 15,fig.height = 7--------------------------------------------
data(scPathway_data)
scPathway_data[1:3,1:3]


## ----echo = F, out.width = "100%"---------------------------------------------
knitr::include_graphics("../inst/extdata/heatmap.png")

## ----fig.width = 7,fig.height = 7---------------------------------------------
library(igraph)

# Load data.
data(ConNetGNN_data)
data("Hv_exp")

# Construct cell set.
index<-grep("0h",colnames(Hv_exp))
cellset<-colnames(Hv_exp)[index]

# Construct gene set.
pathways<-load_path_data(system.file("extdata", "KEGG_human.gmt", package = "scapGNN"))
geneset<-pathways[[which(names(pathways)=="Tight junction [PATH:hsa04530]")]]

plotGANetwork(ConNetGNN_data,cellset,geneset,vertex.label.dist=1.5,main = "Tight junction [PATH:hsa04530]")

## ----fig.width = 7,fig.height = 7---------------------------------------------
require(igraph)
require(graphics)

data(ConNetGNN_data)

# Construct the cell phenotype vector.
cell_id<-colnames(ConNetGNN_data[["cell_network"]])
temp<-unlist(strsplit(cell_id,"_"))
cell_phen<-temp[seq(2,length(temp)-1,by=3)]
names(cell_id)<-cell_phen
head(cell_id)
plotCCNetwork(ConNetGNN_data,cell_id,edge.width=10)

## -----------------------------------------------------------------------------
require(parallel)
require(stats)

# Load the result of the ConNetGNN function.
data(ConNetGNN_data)
data(Hv_exp)

# Construct the cell set corresponding to 0h.
index<-grep("0h",colnames(Hv_exp))
cellset<-colnames(Hv_exp)[index]
H9_0h_cpGM_data<-cpGModule(ConNetGNN_data,cellset,parallel.cores=1)
summary(H9_0h_cpGM_data)

## ----fig.width = 7,fig.height = 7---------------------------------------------
library(igraph)

# Load data.
data(ConNetGNN_data)
data("Hv_exp")
data("H9_0h_cpGM_data")

# Construct cell set.
index<-grep("0h",colnames(Hv_exp))
cellset<-colnames(Hv_exp)[index]

# Construct gene set.
geneset<-H9_0h_cpGM_data$Genes

plotGANetwork(ConNetGNN_data,cellset,geneset,vertex.label.dist=1.5,main = "Gene network of 0h cells activated gene module")

## ----fig.width = 7,fig.height = 7---------------------------------------------
require(igraph)
require(grDevices)
# Load the result of the ConNetGNN function.
data(ConNetGNN_data)
# Obtain cpGModule results for each cell phenotype.
data(H9_0h_cpGM_data)
data(H9_24h_cpGM_data)
data(H9_36h_cpGM_data)
data.list<-list(H9_0h=H9_0h_cpGM_data,H9_24h=H9_24h_cpGM_data,H9_36h=H9_36h_cpGM_data)
plotMulPhenGM(data.list,ConNetGNN_data,margin=-0.05)

