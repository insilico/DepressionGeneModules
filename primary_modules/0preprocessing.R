# source("http://bioconductor.org/biocLite.R")
# biocLite("edgeR")
# install.packages("abodOutlier")
# install.packages("statmod")
# install.packages("pca3d")
# source("https://bioconductor.org/biocLite.R")
# biocLite("sva")
library(abodOutlier) 
library(edgeR)
library(dplyr)
library(reshape2)
library(pca3d)
library(ade4)
library(grid)
library(ggplot2)
library(broom)
# library(sva)

# rm(list=ls())
source("projections.R")
load('combined_cohorts_data.Rdata') # raw counts
subj_ids_combined <- colnames(combined)
rna_ids_combined <- rownames(combined)
rnaSeq <- t(combined)
subject.attrs <- read.csv("Demographic_symptom.csv", stringsAsFactors=FALSE)
attrs.df <- tbl_df(subject.attrs)
phenos.df <- attrs.df %>% 
  filter(X %in% colnames(combined)) %>%
  dplyr::select(X, Diag)
phenos_cc <- phenos.df[match(colnames(combined), phenos.df$X), ]
phenos <- as.factor(ifelse(phenos.df$Diag == "HC", 0, 1))
colnames(phenos_cc) <- c("IID", "ccstatus")
pheno_ids <- phenos_cc$IID
pheno <- as.factor(phenos_cc$ccstatus)
names(pheno) <- pheno_ids
ppp <- ggplot() + coord_fixed() + labs(x="Comp1, Axis1", y="Comp2, Axis2") +
  geom_hline(yintercept=0, col="darkgrey") + 
  geom_vline(xintercept=0, col="darkgrey")



# -----------------------------------------------------------------------------
# Create DGEList object y and remove outlier genes
# -----------------------------------------------------------------------------
y <- DGEList(counts = t(rnaSeq), group = pheno)
keep <- rowSums(cpm(y) > 10) >= 16
summary(keep)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
head(y$samples,10)
head(t(rnaSeq))[,1:10]
plotMDS(y)
plot(y$samples$norm.factors)
text(1:159,y$samples$norm.factors,labels = rownames(rnaSeq))
boxplot(y$samples$norm.factors)


# -----------------------------------------------------------------------------
# Remove the outlier subjects
# -----------------------------------------------------------------------------
rnaSeqFiltered <- rnaSeq[,keep] # removed low counts genes
mycpm <- cpm(y, log = T)
# # # Angle-based outlier detection:
# # mycpm <- cpm(t(rnaSeqFiltered), log = T)
abof <- abod(t(mycpm), method = "knn", k = 10)
names(abof) <- rownames(rnaSeq)
head(sort(abof),5)
abof <- abof[order(abof)]
x <- names(abof)

pdf("logCPM.pdf", height = 5, width = 6)
plot(1:159, abof, main = "LogCPM", col = ifelse(x%in%c("AE396","AN736"), "red", "black"), xlab = "Subjects")
# # text(1:nrow(rnaSeq), abof, color = c("red", rep("black", nrow(rnaSeq))))
# save(abof, file = "abof.Rdata")
dev.off()

sum(rnaSeqFiltered["AE396",])
sum(rnaSeqFiltered["AN736",])

plot(rowSums(rnaSeqFiltered))
plot(rowSums(rnaSeq))
text(1:159,rowSums(rnaSeqFiltered),labels = rownames(rnaSeqFiltered))
plot(rnaSeqFiltered["AE396",])
plot(rnaSeqFiltered["AN736",])
sum(rnaSeq["AN736",])
text(1:159,rowSums(rnaSeq),labels = rownames(rnaSeq))
plot(rnaSeq["AN736",], type = "l")
# condensed.rnaSeqFiltered <- rnaSeqFiltered[rownames(rnaSeqFiltered)!="AN736",]
condensed.rnaSeq <- rnaSeqFiltered[!rownames(rnaSeqFiltered)%in%c("AN736","AE396"),]
condensed.pheno <- pheno[!names(pheno)%in%c("AN736","AE396")]
condensed.cohort <- c(rep("1", ncol(cohort1)), 
                      rep("2", ncol(cohort2)), 
                      rep("3", ncol(cohort3)),
                      rep("4", ncol(cohort4)-2))

y <- DGEList(counts = t(condensed.rnaSeq), group = condensed.pheno)

# -----------------------------------------------------------------------------
# Create new DGEList object y with 2 fewer subjects
# -----------------------------------------------------------------------------
# y <- DGEList(counts = t(condensed.rnaSeq), group = condensed.pheno)
# keep <- rowSums(cpm(y) > 4) >= 10
# summary(keep)
# y <- y[keep, , keep.lib.sizes = FALSE]
# y <- calcNormFactors(y) # TMM normalization
# plotMDS(y)
# plotMDS(y, pch = condensed.cohort)
# 
# plot(y$samples$norm.factors)
# text(1:157,y$samples$norm.factors,labels = rownames(condensed.rnaSeq))
# boxplot(y$samples$norm.factors)



# # estimate dispersion (square root of dispersion is the coeffcient of biological variation BCV)
# design <- model.matrix(~condensed.cohort + condensed.pheno)
# rownames(design) <- colnames(y)
# y <- estimateDisp(y, design, robust=TRUE)
# y$common.dispersion
# # coefficient of biological variation
# sqrt(y$common.dispersion)
# plotBCV(y)
# fit <- glmQLFit(y, design, robust=TRUE)
# plotQLDisp(fit)
# qlf <- glmQLFTest(fit, coef=2:3)
# topTags(qlf)
# FDR <- p.adjust(qlf$table$PValue, method="BH")
# sum(FDR < 0.05)
# qlf <- glmQLFTest(fit)
# topTags(qlf)
# top <- rownames(topTags(qlf))
# cpm(y)[top,1:10]
# 
# 
# # total number of genes significantly up-regulated or down-regulated at 5% FDR
# summary(dt <- decideTestsDGE(qlf))
# isDE <- as.logical(dt)
# DEnames <- rownames(y)[isDE]
# plotSmear(qlf, de.tags=DEnames)
# abline(h=c(-1,1), col="blue")
# 
# imp.genes <- cbind(DEnames)
# write(imp.genes, file = "DE.genes.csv")
# sessionInfo()



# -----------------------------------------------------------------------------
# LogCPM, plot PCA and save data
# -----------------------------------------------------------------------------
# mycpm <- cpm(t(rnaSeq), log = T)

mycpm <- cpm(y, log = T)
# text(1:159, abof, labels = names(abof))
my.pca <- dudi.pca(t(mycpm), nf = 10, scale = F, scannf = F, center = T)
varExplained <- as.integer(my.pca$eig/sum(my.pca$eig)*100)
ppp + geom_point(data=my.pca$li, aes(x=Axis1, y=Axis2, color = as.factor(condensed.cohort)), size = 3.5) + 
  geom_text(data=my.pca$li, aes(x=Axis1, y=Axis2, label=rownames(my.pca$li)))

ppp + geom_point(data=my.pca$li, 
                 aes(x=Axis1, y=Axis2,shape  = as.factor(condensed.cohort), color = condensed.pheno), 
                 size = 3.5) 
cpm.corrected <- removeBatchEffect(cpm(t(condensed.rnaSeq), log = T), condensed.cohort)
my.pca <- dudi.pca(t(cpm.corrected), nf = 4, scale = F, scannf = F, center = T)
# ppp + geom_point(data=my.pca$li, aes(x=Axis1, y=Axis2, color = as.factor(condensed.cohort)), size = 3.5) + 
#   geom_text(data=my.pca$li, aes(x=Axis1, y=Axis2, label=rownames(my.pca$li)))
ppp + geom_point(data=my.pca$li, 
                 aes(x=Axis1, y=Axis2,shape  = as.factor(condensed.cohort), color = condensed.pheno), 
                 size = 3.5) 
boxplot(cpm.corrected, las=2, cex.axis=0.75, main="Log CPM After Outliers and Batch Effect Removal",
        col=condensed.pheno)
boxplot(mycpm, las=2, cex.axis=0.75, main="Log CPM Before Outliers and Batch Effect Removal",
        col=pheno)

filtered.logcpm.uncorrected <- mycpm
# save(filtered.logcpm.uncorrected, file = "filtered.logcpm.uncorrected.Rdata")
filtered.logcpm.corrected <- cpm.corrected
save(filtered.logcpm.corrected, file = "filtered.logcpm.corrected.Rdata")



