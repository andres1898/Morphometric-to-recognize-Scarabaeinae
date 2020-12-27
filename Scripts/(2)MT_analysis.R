## Andres Ordoñez, Juliette Gualdrón & Daniel Rafael Miranda-Esquivel
## Escuela de Biología, Universidad Industrial de Santander
## Geometric and traditional morphometrics for differentiating Dung Beetles
## Part: Traditional morphometrics analyses
## R version 3.6.3
## Dic 2020

# 1. Load packeges ####

library(geomorph)
library(psych)
library(tidyverse)
library(simstudy)
library(jtools)
library(nabor)
library(DMwR)
library(factoextra)
library(FactoMineR)
library(fpc)
library(dendextend)
library(NbClust)
library(cluster)
library(shapes)
library(ggplot2)
library(Morpho)
library(ellipse)
library(stringr)
library(fpc)
library(colorspace)
library(phangorn)

# Read measures ####

head_measures <- read.csv("../Data/head_measures.csv", header = T)
protibia_measures <- read.csv("../Data/protibia_measures.csv", header = T)
#
list_structures <- list(protibia = protibia_measures, head = head_measures)

#### Read info matrix
info_head <- read.csv("../Data/head_information_table.csv")
info_protibia <- read.csv("../Data/protibia_information_table.csv")
#
list_info <- list(protibia = info_protibia, head = info_head)

# list with the names of results directories
dir.create(paste0("../Results/", Sys.Date(), "_Traditional_results"), showWarnings = FALSE)
dir.create(paste0("../Results/", Sys.Date(), "_Traditional_results/Protibia"), showWarnings = FALSE)
dir.create(paste0("../Results/", Sys.Date(), "_Traditional_results/Head"), showWarnings = FALSE)

results <- c(paste0("../Results/", Sys.Date(), "_Traditional_results/Protibia/"), paste0("../Results/", Sys.Date(), "_Traditional_results/Head/"))

# Correlation ####

for (i in 1:2) {
pdf(file = paste0(results[i], "correlation.pdf"))
pairs.panels(list_structures[i][[1]][,-1], stars = T, ellipses = F)
dev.off()
}

#BoxViolinplot ####

for (f in 1:2) {
  
  pdf(paste0(results[f], "Boxplot_measures.pdf"))  
  
  for (i in 2:ncol(list_structures[f][[1]])) {
    
  data <- list_structures[f][[1]]
    p <- ggplot(data , aes(x = list_info[f][[1]][,2], y = data[, i])) + geom_violin() + 
      geom_boxplot(width = 0.2) +
      labs(title = paste0("medida de : ", colnames(list_structures[f][[1]])[i])) +
      xlab("Morfotipo") +
      ylab("longitud en mm") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 0.9))
  
    plot(p)
  }
  dev.off()
}


#PCA ####
colores <- c("#F68656", "#F0BA70", "#F46E01", "#F0DF66", "#DA2B27", "#ad6989", "#303960", "#a8df65")

for (i in 1:2){
  pca <- prcomp(list_structures[i][[1]][,-1])
  pca_scores <- as_data_frame(pca$x)
  rownames(pca$x) <- list_info[i][[1]]$Morfotipo
  pca_imp <- summary(pca)
  
  pdf(paste0(results[i], "PCA_ellipse_measures.pdf"), paper = "a4r")
    p <- ggplot(pca_scores, aes(PC1, PC2, color = list_info[i][[1]]$Morfotipo)) +
      geom_point() +
      stat_ellipse(type = "t") +
      xlab(paste0("PC1 = ", round(pca_imp$importance[2]*100, digits = 1), "%")) +
      ylab(paste0("PC2 = ", round(pca_imp$importance[5]*100, digits = 1), "%")) +
      # scale_shape_manual(values = FIGURAS) +
      scale_color_manual(values = colores) +
      theme(legend.position = "right",
            axis.text=element_text(size=14),
            axis.title=element_text(size=18),
            panel.background = element_rect("white"),
            panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                            colour = "gray90"), 
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "gray90"))
    plot(p)
  dev.off()
  
  write.csv(data.frame(pca_imp$importance), file = paste0(results[i], "PCA_importance_measures.csv"))
  
  write.csv(data.frame(pca$x), file = paste0(results[i], "PCA_scores.csv"))
}


# 8. Canonical Variable Analyses -CVA- #####

#CVA for test morphotype and locality recognition
for (i in 1:length(list_structures)) {
  
  for (j in 2) {
  
    #make the CVA, elimnate the group with one entry (individuo 19 from Santo Domingo, Cantagallo)
    cva_structure <- CVA(list_structures[i][[1]][,-1], list_info[i][[1]][j][[1]], 
                         weighting = T, plot = TRUE, 
                         rounds = 0, cv = TRUE, p.adjust.method = "none")
    
    #plot the CVA space and save
    p <- ggplot(as.data.frame(cva_structure$CVscores), 
                aes(x = cva_structure$CVscores[,1], 
                    y = cva_structure$CVscores[,2],
                    color = as.factor(row.names(cva_structure$CVscores)))) + 
      geom_point() +
      stat_ellipse(type = "t") +
      xlab(paste0("CV1 = ", 
                  round(cva_structure$Var[1,2], digits = 1), "%")) +
      ylab(paste0("CV2 = ", 
                  round(cva_structure$Var[2,2], digits = 1), "%")) +
      labs(color = names(list_info[i][[1]][j])) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14)) +
      scale_color_manual(values = colores)
    
    pdf(file = paste0(results[i], names(list_structures)[i], "_CVA_", names(list_info[i][[1]][j]), ".pdf"), paper = "a4r")
    plot(p)
    dev.off()
    
    #save the table of frecuencies hits
    tmp_cva <- print(cva_structure)
    
    write.csv(file = paste0(results[i], names(list_structures)[i], "_CVA_", names(list_info[i][[1]][j]), ".csv"), x = tmp_cva$table )
  } 
}

# 9. Cluster UPGMA ####
# colLab function of Joris Meys (http://stackoverflow.com/) modified by Ambrosio Torres to coloring edgePar and labels of a dendrogram ####

colLab <- function(n) {
  if(is.leaf(n)) {
    a <- attributes(n)
    # clusMember - a vector designating leaf grouping
    # labelColors - a vector of colors for the above grouping
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "edgePar") <- list(col = labCol)
    attr(n, "label") <- NULL ;
  }
  n
}

# 9.1 by morphotypes PCA-scores####
for (i in 1:length(list_structures)){
  pca <- prcomp(list_structures[i][[1]][,-1])
  pca_scores <- as_data_frame(pca$x)
  rownames(pca$x) <- list_info[i][[1]]$Morfotipo
  datavalores <- pca_scores
  row.names(datavalores) <- list_info[i][[1]]$ID
  
  distvalores <- dist(datavalores) 
  clust_valores <- hclust(distvalores, method = "complete") 
  #labelColors <- c("#F8766D", "#BB9D00", "#00B81F", "#00C0B8", "#00A5FF", "#E76BF3", "#FF6C90", "brown")
  labelColors <- colores
  
  clusMember <- rep(NA,length(rownames(datavalores)))
  
  for (j in 1:length(list_info[i][[1]]$Morfotipo)) {
    
    if (list_info[i][[1]]$Morfotipo[j] == "C. cyanellus") {
      clusMember[j] <- 1
    }
    if (list_info[i][[1]]$Morfotipo[j] == "C. juvencus") {
      clusMember[j] <- 2
    }
    if (list_info[i][[1]]$Morfotipo[j] == "C. septemmaculatus") {
      clusMember[j] <- 3
    }
    if (list_info[i][[1]]$Morfotipo[j] == "C. subhyalinus") {
      clusMember[j] <- 4
    }
    if (list_info[i][[1]]$Morfotipo[j] == "C. Glaphyrocanthon") {
      clusMember[j] <- 5
    }
    if (list_info[i][[1]]$Morfotipo[j] == "Dichotomius") {
      clusMember[j] <-6
    }
    if (list_info[i][[1]]$Morfotipo[j] == "Eurysternus") {
      clusMember[j] <- 7
    }
    if (list_info[i][[1]]$Morfotipo[j] == "Ontherus") {
      clusMember[j] <- 8
    }
  }
  
  names(clusMember) <- rownames(datavalores)
  
  clusDendro <- as.dendrogram(as.hclust(clust_valores))  
  
  clusDendro <- dendrapply(clusDendro, colLab)
  
  #make graph and save it
  pdf(file = paste0(results[i], names(list_structures)[i], "_Dendrogram_morphotypes.pdf"), paper = "a4r", width = 595, height = 842)
  #png(filename = paste0(results[i], names(list_structures)[i], "_Dendrogram_morphotypes.png"), width = 420*4, height = 290*3, res = 200)
  par(mar=c(0.5, 4.5, 1, 1))
  par(oma= c(0,0,0,0))
  par(lwd=3)
  plot(clusDendro,horiz=F,axes=T, ylab= "Euclidean distance", cex.axis=1.3,cex.lab=1.7)
  par(lwd=1)
  legend("topright", pch= 21, pt.bg=labelColors, 
         legend=expression(italic("C. cyanellus"), italic("C. juvencus"),italic("C. septemmaculatus"), italic("C. subhyalinus"), italic("C. Glaphyrocanthon"), italic("Dichotomius"), italic("Eurysternus"), italic("Ontherus")), pt.cex = .9, cex=.6)
  #dev.off()
  dev.off()
}

# 10. K-means analyses ####
############# optimal cluster by Gap statistic
for (i in 1:length(list_structures)) {
  pdf(file = paste0(results[i], names(list_structures)[i], "_GAP_kmeans.pdf"), paper = "a4r")
  plot(fviz_nbclust(list_structures[i][[1]][,-1], kmeans, method = "gap_stat", nboot = 500))
  dev.off()
}
# Head has 2 k-groups 
# Protibia has 6 k-groups