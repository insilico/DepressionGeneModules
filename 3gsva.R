rm(list = ls())
load("0.8genes.filtered.corrected.Rdata")

# load all the functions written in the folder funcs
pathnames <- list.files(pattern="[.]R$", path=paste(getwd(),"/funcs", sep =""), full.names=TRUE);
sapply(pathnames, FUN=source);
sab.covs <- covs.short[,c("sex", "age", "batch")]
sab.covs$batch <- as.factor(sab.covs$batch)
genes <- colnames(rnaSeq)

# read data set
rnaseq_counts <- rnaSeq
expr <- as.matrix(t(rnaseq_counts))

# read in modules as gene sets
group.genes <- read.csv("wgcna.0.3.group.new.csv", stringsAsFactors = F, row.names = 1)
colnames(group.genes) <- c("gene", "dynamicMods")
group.genes$dynamicMods <- as.factor(group.genes$dynamicMods)
c7_genes_sets <- list()
c7_genes_sets$genesets <- split(group.genes[,"gene"], group.genes$dynamicMods)
c7_genes_sets$geneset.names <- names(c7_genes_sets$genesets)

num.genes <- sapply(c7_genes_sets$genesets, length)
names(c7_genes_sets$genesets) <- paste("Module", names(c7_genes_sets$genesets), sep = "")
# run GSVA to get a gene set by subject matrix of enrichment scores
res <- gsva(expr, c7_genes_sets$genesets, method="ssgsea", verbose=T, rnaseq=T)
rownames(res) <- c7_genes_sets$geneset.names
num.modules <- nrow(res)
rownames(res) <- paste("Module", rownames(res), sep = "")
modules <- rownames(res)
# read phenotypes
phenotypes <- phenos
levels(phenotypes) <- c(0,1)
phenos.df <- data.frame(phenos)


# ------------------------------------------------------------------------------------
# Create star/spider plot
# ------------------------------------------------------------------------------------
collapse.genes <- t(res)
cases <- names(phenos)[phenos == "MDD"]
ctrls <- names(phenos)[phenos == "HC"]
collapse.genes.cases <- collapse.genes[cases,]
collapse.genes.ctrls <- collapse.genes[ctrls,]
save(collapse.genes, num.genes, cases, ctrls, collapse.genes.cases, collapse.genes.ctrls, 
     group.genes, num.modules, modules, phenotypes, phenos.df, sab.covs, file = "enrichmentProfile.Rdata")




