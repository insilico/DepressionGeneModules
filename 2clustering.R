# clustering code

library("dendextend")
library(WGCNA)
# library(dynamicTreeCut)
# source("http://bioconductor.org/biocLite.R")
# biocLite("GSVA")
library(GSVA)
library(reshape2)
library(grid)
library(broom)
library(ggplot2)
library(ade4)
library(leaps) 
library(dplyr)
# install.packages("flashClust")
library(flashClust)

rm(list = ls())
load("0.8genes.filtered.corrected.Rdata")

pathnames <- list.files(pattern="[.]R$", path=paste(getwd(),"/funcs", sep =""), full.names=TRUE);
sapply(pathnames, FUN=source);
sab.covs <- covs.short[,c("sex", "age", "batch")]
sab.covs$batch <- as.factor(sab.covs$batch)
genes <- colnames(rnaSeq)
rnaseq_counts <- rnaSeq
expr <- as.matrix(t(rnaseq_counts))

thres = 0.2
inputmat = rnaSeq
regs = genes
corrThresh = thres
inputmat <- cor(inputmat)
#   corrThresh <- 0.4
powCorr <- inputmat^1
adjacency <- abs(powCorr) > corrThresh
TOM <- TOMsimilarity(adjacency+0, TOMType = "unsigned");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, cutHeight = 0.95, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

mod_sizes <- table(dynamicMods)
length(mod_sizes)
mod_sizes

# # gets all cluster members except first (-1) because first is "unassigned genes"
numMods <- sort(unique(names(dynamicMods)))[-1]
names(dynamicMods) <- colnames(adjacency)
modNumbers <- unique(dynamicMods)
names(dynamicMods) <- colnames(adjacency)
gene.mod.list <- lapply(modNumbers, function(x) {
  cluster_genes <- dynamicMods[dynamicMods == x]
  names(cluster_genes)
})

gene.mod.df <- melt(gene.mod.list)
rownames(gene.mod.df) <- gene.mod.df[,1]

write.csv(gene.mod.df, paste("wgcna.", thres, ".group.new.csv", sep = ""))
