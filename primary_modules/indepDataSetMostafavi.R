library(GSVA)
library(broom) # function tidy
rm(list = ls())
# system.time(mostafaviMDD <- read.csv("mostafaviMDD.csv", row.names = 1)) # takes a minute or two
rowSums(is.na(mostafaviMDD))
# ignore last row because 8638 missing values
mostafaviMDD <- mostafaviMDD[-nrow(mostafaviMDD), ]


covariates <- read.csv("Biological_and_hidden_factors.csv", stringsAsFactors = F, row.names = 1)
clinic.variables <- read.table("Clinical_variables.txt")
mdd.status <- read.table("Dx_Case_status.txt", header = T)
mdd.status$MDDstatus <- as.factor(mdd.status$MDDstatus)
levels(mdd.status$MDDstatus) <- c("HC", "MDD") 
# check this with the author; on readme: 
# This file cotain disease status (MDD=1, Control=0) 
# for each of the 922individuals analyzed in this study.
group.genes <- read.csv("wgcna.0.3.group.new.csv", stringsAsFactors = F, row.names = 1)
colnames(group.genes) <- c("gene", "dynamicMods")
group.genes$dynamicMods <- as.factor(group.genes$dynamicMods)


mostafaviGenes <- intersect(colnames(mostafaviMDD), rownames(group.genes))
mostafaviMDDcut <- mostafaviMDD[, mostafaviGenes]
mostafaviCheck <- mostafaviMDDcut[, 1:20]
# mostafaviCheck <- mostafaviMDDcut[, 1:20]
boxplot(t(mostafaviCheck))
for (i in 1:20){
  print(sd(mostafaviCheck[,i]))
}
sum(is.na(mostafaviMDDcut))

sds <- apply(mostafaviMDDcut, 2, sd)
hist(sds)
group.genes <- group.genes[mostafaviGenes, ]
# group.genes <- data.frame(value =  names(dynamicMods), dynamicMods, stringsAsFactors = F)
c7_genes_sets <- list()
c7_genes_sets$genesets <- split(group.genes[,"gene"], group.genes$dynamicMods)
c7_genes_sets$geneset.names <- names(c7_genes_sets$genesets)

sapply(c7_genes_sets$genesets, length)
names(c7_genes_sets$genesets) <- paste("Module", names(c7_genes_sets$genesets), sep = "")
# run GSVA to get a gene set by subject matrix of enrichment scores
res <- gsva(as.matrix(t(mostafaviMDDcut)), c7_genes_sets$genesets, method="ssgsea", verbose=T, rnaseq=T)
rownames(res) <- c7_genes_sets$geneset.names
num.modules <- nrow(res)
rownames(res) <- paste("Module", rownames(res), sep = "")
modules <- rownames(res)
collapse.genes <- t(res)

phenos.df <- mdd.status[, 2, drop = F]
rownames(phenos.df) <- mdd.status[, 1]

# covariates$SmokeBefore <- as.factor(covariates$SmokeBefore)

linear.reg <- function(my.formula, data, y, subj.id = rownames(data), fam = gaussian){
  linear_result <- matrix(NaN, ncol(data), 3)
  colnames(linear_result) <- c("Intercept", "Effect", "pvalue")
  rownames(linear_result) <- colnames(data)
  ptm <- proc.time()
  for (i in 1:(ncol(data))){
    tryCatch({
      model.df <- data.frame(genei = data[subj.id,i], y = y[subj.id,], covariates[subj.id, c("Sex", "BMI_CURRENT", "SmokeBefore")])
      my.linear <- tidy(glm(as.formula(my.formula), data =  model.df, family = fam))
      my.linear.genei <- my.linear[1:2,]
      linear_result[i,1:2]  <- my.linear.genei$estimate
      linear_result[i,3] <- my.linear.genei$p.value[2]
    }, 
    warning = function(w){
      print(paste('fitted probabilities numerically 0 or 1 occurred at i =', i))
    }
    )
  }
  duration <- proc.time() - ptm
  print(duration)
  
  linear_result <- data.frame(na.omit(linear_result))
  linear_result$p.value.corrected.FDR <- p.adjust(linear_result$pvalue, "BH")
  linear.batch <- data.frame(linear_result[order(linear_result[,3],decreasing=F),])
  print(head(linear.batch, 23))
  linear.batch
}

phq.nobatch.lin <- linear.reg(my.formula = "y ~ . ", data = collapse.genes[, c("Module5", "Module17")], 
                              y = clinic.variables[, "Phq_Tot", drop = F], 
                                subj.id = rownames(collapse.genes))
pheno.nobatch.lin <- linear.reg(my.formula = "y ~ . - Sex - BMI_CURRENT", data = collapse.genes, y = phenos.df, 
                                subj.id = rownames(collapse.genes), 
                                fam = "binomial")
pheno.nobatch.lin <- linear.reg(my.formula = "y ~ . ", data = collapse.genes[, c("Module5", "Module17")], y = phenos.df, 
                                subj.id = rownames(collapse.genes), fam = "binomial")

pheno.nobatch.lin <- linear.reg(my.formula = "y ~ . ", data = collapse.genes[, ], y = phenos.df, 
                                subj.id = rownames(collapse.genes), fam = "binomial")


# subj.id = rownames(collapse.genes)
# tidy(glm(MDDstatus ~ ., 
#     data = cbind(MDDstatus = phenos.df[subj.id,], collapse.genes, 
#                  covariates[subj.id, c("Sex", "BMI_CURRENT")]), 
#     family = "binomial"))






