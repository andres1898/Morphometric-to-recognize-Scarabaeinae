## Andres Ordoñez, Juliette Gualdrón & Daniel Rafael Miranda-Esquivel
## Escuela de Biología, Universidad Industrial de Santander
## Geometric and traditional morphometrics for differentiating Dung Beetles
## Part: Geometric morphometrics  analyses
## R version 3.6.3
## Dic 2020

# 1. Load packeges ####

library(geomorph)
library(shapes)
library(ggplot2)
library(Morpho)
library(ellipse)
library(stringr)
library(factoextra)
library(FactoMineR)
library(fpc)
library(dendextend)
library(NbClust)
library(cluster)
library(colorspace)
library(phangorn)


# 2. Create directories for Geometrical Morphometric ######

dir.create("../Results")
dir.create(paste0("../Results/", Sys.Date(), "_Geometric_results"), showWarnings = FALSE)


# 3. Read data ####

######## read as tps and align
head <- readland.tps("../Data/head_aligned.tps", specID = "ID"); head_alg <- gpagen(head)
protibia <- readland.tps("../Data/protibia_aligned.tps", specID = "ID"); protibia_alg <- gpagen(protibia)

# create a list that contains tps
list_structures <- list(protibia = protibia_alg, head = head_alg)

######## read information matrix
info_head <- read.csv("../Data/head_information_table.csv", header = TRUE)
info_protibia <- read.csv("../Data/protibia_information_table.csv", header = TRUE)

# create a list that contains info matrices
list_information <- list(protibia = info_protibia, head = info_head)

# list with the names of results directories
dir.create(paste0("../Results/", Sys.Date(), "_Geometric_results/Protibia"), showWarnings = FALSE)
dir.create(paste0("../Results/", Sys.Date(), "_Geometric_results/Head"), showWarnings = FALSE)

resultados <- c(paste0("../Results/", Sys.Date(), "_Geometric_results/Protibia/"), paste0("../Results/", Sys.Date(), "_Geometric_results/Head/"))

# 4. Generalized Procrustes Analyses ####
for (i in 1:length(list_structures)) { #for each strucutre into the list. Saves the results, doesn't plot it
  pdf(file = paste0(resultados[i], names(list_structures)[i], "_GPA.pdf"), paper = "a4r")
    plotAllSpecimens(list_structures[i][[1]]$coords) #ver el gpa y todos align
  dev.off()
}

# 5.1 Principal Component Analyses with deformation plates ####

#create an empty list to fill with the result behind
pca_structures <- list(protibia = NA, head = NA)

for (i in 1:length(list_structures)) { #for each strucutre into the list. Saves the results, doesn't plot it
  pdf(file = paste0(resultados[i], names(list_structures)[i], "_PCA_deformation_Grid.pdf"), paper = "a4r")
  pca_structures[i][[1]] <- plotTangentSpace(list_structures[i][[1]]$coords, warpgrids = TRUE, legend = TRUE, groups = list_information[i][[1]][,"Morfotipo"])
  pca_structures[i][[1]] <- plotTangentSpace(list_structures[i][[1]]$coords, warpgrids = TRUE, legend = TRUE, groups = list_information[i][[1]][,"Morfotipo"])
  dev.off()
}

# 5.2 Principal Component Analyses with ellipse ####

colores <- c("#F68656", "#F0BA70", "#F46E01", "#F0DF66", "#DA2B27", "#ad6989", "#303960", "#a8df65")

for (i in 1:length(list_structures)) { #for each strucutre into the list. Saves the results, doesn't plot it
  
  pca_scores <- as.data.frame(pca_structures[i][[1]]$pc.scores) #individual pca scores as dataframe

  pdf(file = paste0(resultados[i], names(list_structures)[i], "_PCA_ellipses.pdf"), paper = "a4r")
  
  for (j in c(2,3)) {  #corresponde to category in info matrices
    p <- ggplot(pca_scores, aes(PC1, PC2, color = list_information[i][[1]][,j])) +
      geom_point() +
      stat_ellipse(type = "t") +
      xlab(paste0("PC1 = ", 
                  round(pca_structures[i][[1]]$pc.summary$importance[2]*100, digits = 1), "%")) +
      ylab(paste0("PC2 = ", 
                  round(pca_structures[i][[1]]$pc.summary$importance[5]*100, digits = 1), "%")) +
      labs(color = colnames(list_information[i][[1]])[j]) +
      theme(axis.text=element_text(size=14),
            axis.title=element_text(size=18),
            panel.background = element_rect("white"),
            panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                            colour = "gray90"), 
            panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                            colour = "gray90"))+ 
      scale_color_manual(values = colores)
      
    plot(p)
  }
  dev.off()
  
  # save components deformation
  ref <- mshape(list_structures[i][[1]]$coords)
  pdf(paste0(resultados[i], "PCA_deformation", names(list_structures)[i], ".pdf") )
  plotRefToTarget(ref, pca_structures[i][[1]]$pc.shapes$PC1min)
  title("PC1min")
  plotRefToTarget(ref, pca_structures[i][[1]]$pc.shapes$PC1max)
  title("PC1max")
  plotRefToTarget(ref, pca_structures[i][[1]]$pc.shapes$PC2min)
  title("PC2min")
  plotRefToTarget(ref, pca_structures[i][[1]]$pc.shapes$PC2max)
  title("PC2max")
  
  dev.off()

  # save components
  write.csv(file = paste0(resultados[i], names(list_structures)[i], "_PCAscores.csv"), x = pca_scores)

}


# 6. Deformation grills ####

# separate and analyses morphotypes
list_grills <- list(protibia = NA, head = NA) #list to fill

for (i in 1:length(list_structures)) { 
  structures <- list_structures[i][[1]]
  reference <- mshape(structures$coords) #for each structure
  #separate morphotypes
  C_juvencus <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. juvencus"]) 
  C_septemmaculatus <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. septemmaculatus"])
  C_subhyalinus <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. subhyalinus"])
  C_cyanellus <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. cyanellus"])
  C_Glaphyrocanthon <- mshape(structures$coords[,,list_information[i][[1]][,2]=="C. Glaphyrocanthon"])
  Dichotomius <- mshape(structures$coords[,,list_information[i][[1]][,2]=="Dichotomius"])
  Eurysternus <- mshape(structures$coords[,,list_information[i][[1]][,2]=="Eurysternus"])
  Onthophagus <- mshape(structures$coords[,,list_information[i][[1]][,2]=="Onthophagus"])
  
  list_grills[i][[1]] <- list(Reference = reference, C_juvencus = C_juvencus, 
                         C_septemmaculatus = C_septemmaculatus, C_subhyalinus = C_subhyalinus,
                         C_cyanellus = C_cyanellus, C_Glaphyrocanthon = C_Glaphyrocanthon, 
                         Dichotomius = Dichotomius, Eurysternus = Eurysternus, Onthophagus = Onthophagus)
}

# plot deformations for each structures and morphotypes respect to the reference
for (i in 1:length(list_structures)) { #for each structure
  
  pdf(file = paste0(resultados[i], names(list_structures)[i], "_deformation_grills.pdf"))
  
  for (j in 2:8) { #for each morphotype
    
    TPSgrid(TT = list_grills[i][[1]][1][[1]], YY = list_grills[i][[1]][j][[1]], 
            xbegin=-999, ybegin = -999, xwidth = -999,
            opt=1, ext=0.1, ngrid = 20, cex=1, pch=20, 
            col=1, mag=2.5, axes3=FALSE) #deformation for reference to each morphotype
    title(names(list_grills[i][[1]][j]))
    
  }
  
  dev.off()
  
}

# 7. Allometry ####

#make a list with all information
list_dataFrames_size_shape <- list(protibia = NA, head = NA) #list to fill

for (i in 1:length(list_structures)) { #for each strucutre
    
  tmp_size <- list_structures[i][[1]]$Csize #extract the centroide size
  tmp_shape <- list_structures[i][[1]]$coords #extract the shape coordinates
  list_dataFrames_size_shape[i][[1]] <- geomorph.data.frame(size = tmp_size, shape = tmp_shape,
                                                            groups = list_information[i][[1]]$Morfotipo)
}

#regresion between shape and other variables
for (i in 1:length(list_structures)) {
  cab_anova_size_spp <- procD.lm(shape ~ size*groups, data = list_dataFrames_size_shape[i][[1]]) #ANOVA-procrustes size vs shape, groups and shape*groups
  tmp <- summary(cab_anova_size_spp) #the result
  
  write.csv(file = paste0(resultados[i], names(list_structures)[i], "_allometry.csv"), as.data.frame(tmp$table)) #save the statistic values
  #
  pdf(file = paste0(resultados[i], names(list_structures)[i], "_allometry.pdf"))
      
  plot(cab_anova_size_spp, type = "regression", predictor = as.numeric(list_dataFrames_size_shape[i][[1]]$groups), reg.type = "RegScore")
  title("Shape vs Morphotypes")
  
  plot(cab_anova_size_spp, type = "regression", predictor = list_dataFrames_size_shape[i][[1]]$size, reg.type = "RegScore")
  title("Shape vs Size")
  
  plot(cab_anova_size_spp, type = "regression", predictor = as.numeric(list_dataFrames_size_shape[i][[1]]$groups)*list_dataFrames_size_shape[i][[1]]$size, reg.type = "RegScore")
  title("Shape vs Size*Morphotypes")
  
  dev.off()
  #
  
  plotlm <- plot(cab_anova_size_spp, type = "regression", predictor = list_dataFrames_size_shape[i][[1]]$size, reg.type = "RegScore")
  
  DF_ggplot <- data.frame(Reg_shape = plotlm$PredLine,
                          spp = list_dataFrames_size_shape[i][[1]]$groups,
                          CS = list_dataFrames_size_shape[i][[1]]$size)
  
  p <- ggplot(DF_ggplot, aes(x = CS, y = Reg_shape, colour = spp)) + 
    geom_point()
  
  pdf(file = paste0(resultados[i], names(list_structures)[i], "_Reg_shapeVSsize.pdf"))
      plot(p)
  dev.off()
  
}

# 8. Canonical Variable Analyses -CVA- #####

#CVA for test morphotype and locality recognition
for (i in 1:length(list_structures)) {
  
  for (j in c(2)) {
  
    #make the CVA, elimnate the group with one entry (individuo 19 from Santo Domingo, Cantagallo)
    cva_structure <- CVA(list_structures[i][[1]]$coords, list_information[i][[1]][j][[1]], 
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
      labs(color = names(list_information[i][[1]][j])) +
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14)) +
      scale_color_manual(values = c("black", "green", "blue", "cyan", "red", "pink", "yellow", "gray40"))
    
    pdf(file = paste0(resultados[i], names(list_structures)[i], "_CVA_", names(list_information[i][[1]][j]), ".pdf"), paper = "a4r")
      plot(p)
    dev.off()
    
    #save the table of frecuencies hits
    tmp_cva <- print(cva_structure)
    
    write.csv(file = paste0(resultados[i], names(list_structures)[i], "_CVA_", names(list_information[i][[1]][j]), ".csv"), x = tmp_cva$table )
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

# 9.1 by morphotypes Procrustes Distances####
for (i in 1:length(list_structures)){
  
  datavalores <- two.d.array(protibia_alg$coords)
  row.names(datavalores) <- list_information[i][[1]]$ID
  
  distvalores <- head_alg$procD
  clust_valores <- hclust(distvalores, method = "complete") 
  labelColors <- colores
  
  clusMember <- rep(NA, length(rownames(datavalores)))
  
  for (j in 1:length(list_information[i][[1]]$Morfotipo)) {
    
    if (list_information[i][[1]]$Morfotipo[j] == "C. cyanellus") {
      clusMember[j] <- 1
    }
    if (list_information[i][[1]]$Morfotipo[j] == "C. juvencus") {
      clusMember[j] <- 2
    }
    if (list_information[i][[1]]$Morfotipo[j] == "C. septemmaculatus") {
      clusMember[j] <- 3
    }
    if (list_information[i][[1]]$Morfotipo[j] == "C. subhyalinus") {
      clusMember[j] <- 4
    }
    if (list_information[i][[1]]$Morfotipo[j] == "C. Glaphyrocanthon") {
      clusMember[j] <- 5
    }
    if (list_information[i][[1]]$Morfotipo[j] == "Dichotomius") {
      clusMember[j] <-6
    }
    if (list_information[i][[1]]$Morfotipo[j] == "Eurysternus") {
      clusMember[j] <- 7
    }
    if (list_information[i][[1]]$Morfotipo[j] == "Ontherus") {
      clusMember[j] <- 8
    }
  }
  
  names(clusMember) <- rownames(datavalores)
  
  clusDendro <- as.dendrogram(as.hclust(clust_valores))  
  
  clusDendro <- dendrapply(clusDendro, colLab)
  
  #make graph and save it
  
  pdf(file = paste0(resultados[i], names(list_structures)[i], "_Dendrogram_ProcDist_morphotypes.pdf"), paper = "a4r")
  par(mar=c(0.5, 4.5, 1, 1))
  par(oma= c(0,0,0,0))
  par(lwd=2)
  plot(clusDendro,horiz=F,axes=T, ylab= "Euclidean distance", cex.axis=1.3,cex.lab=1.7)
  par(lwd=1)
  legend("topright", pch= 21, pt.bg=labelColors, 
         legend=expression(italic("C. cyanellus"), italic("C. juvencus"),italic("C. septemmaculatus"), italic("C. subhyalinus"), italic("C. Glaphyrocanthon"), italic("Dichotomius"), italic("Eurysternus"), italic("Ontherus")), pt.cex = .9, cex=.6)
  dev.off()
}

# 10. K-means analyses ####
############# optimal cluster by Gap statistic
for (i in 1:length(list_structures)) {
  tmp_data <- two.d.array(list_structures[i][[1]]$coords)
  pdf(file = paste0(resultados[i], names(list_structures)[i], "_GAP_kmeans.pdf"), paper = "a4r")
  plot(fviz_nbclust(tmp_data, kmeans, method = "gap_stat", nboot = 500))
  dev.off()
}
# head has 8 k-groups 
# Protibia has 5 k-groups