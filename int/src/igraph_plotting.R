# This script plots igraph networks with highlighted genes
# written by tiago lubiana 20/08/2018

# load data ----
load('./int/results/igraph_dengue.Rdata')
predicted_degs <- read.csv('./int/data/gene_candidates_top_dengue')

library(data.table)
library(igraph)
library(netresponse)
library(tibble)
library(dplyr)

vertices_g <- as_ids(V(g))

vertices_df <- as.data.frame(vertices_g)
colnames(vertices_df) <- 'nodes'
vertices_df$predicted <- ifelse(vertices_df$nodes %in%  predicted_degs$gene, 'yes', 'no' ) 

edgelist <- as.data.frame(as_edgelist(g, names = TRUE))

net <- graph_from_data_frame(d=edgelist, vertices=vertices_df, directed=F) 

V(net)[predicted=='yes']

net <- simplify(net, remove.multiple = F, remove.loops = T) 
plot(net)

colrs <- c("gray80","tomato")

E(net)$edge.color <- "black"
l <- layout_with_fr(g)
V(net)$size <- 4
V(net)$color <- colrs[as.factor(V(net)$predicted)]
plot(net, )

svg('./imgs/network_for_project.svg',)
plot(net,rescale=F, layout=l*0.098, edge.arrow.size=.4, edge.curved=.1, vertex.label=NA) 
dev.off()
