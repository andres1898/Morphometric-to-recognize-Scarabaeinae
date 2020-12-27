## Andres Ordoñez, Juliette Gualdrón & Daniel Rafael Miranda-Esquivel
## Escuela de Biología, Universidad Industrial de Santander
## Geometric and traditional morphometrics for differentiating Dung Beetles
## Part: Traditional morphometrics analyses
## R version 3.6.3
## Dic 2020

# 1. Load packeges ####
library(dendextend)
library(stringr)
#library(phytools)

# Read PCA scores
Head_MG <- read.csv("../Results/2020-12-27_Geometric_results/Head/head_PCAscores.csv", row.names = 1)
Protibia_MG <- read.csv("../Results/2020-12-27_Geometric_results/Protibia/protibia_PCAscores.csv", row.names = 1)
Head_MT <- read.csv("../Results/2020-12-27_Traditional_results/Head/PCA_scores.csv", row.names = 1)
Protibia_MT <- read.csv("../Results/2020-12-27_Traditional_results/Protibia/PCA_scores.csv", row.names = 1)

# make a list
Components <- list(Head_MG = Head_MG, Protibia_MG = Protibia_MG,
                   Head_MT = Head_MT, Protibia_MT = Protibia_MT)


# Read csv table of information
info_head <- read.csv("../Data/head_information_table.csv")
info_protibia <- read.csv("../Data/protibia_information_table.csv")

list_info <- list(protibia = info_protibia, head = info_head)

# Create distance matrices
Component_dist <- lapply(Components, dist)

# Make a cluster with UPGMA
Clusters <- lapply(Component_dist, hclust, method = "average")

# clusters as dendograms
Dendrograms <- lapply(Clusters, as.dendrogram)

# common nodes
j <- dendlist(HeadMG = Dendrograms$Head_MG, ProtibiaMG = Dendrograms$Protibia_MG, 
              HeadMT = Dendrograms$Head_MT, ProtibiaMT = Dendrograms$Protibia_MT)

h <- dist.dendlist(j)

h1 <- h - nnodes(Dendrograms$Head_MG)
h1 <- as.matrix(h1)
h1 <- as.data.frame(h1)

write_csv(x = h1, path = "../Results/Nodes_Comparison.csv")
