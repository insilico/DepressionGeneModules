rm(list = ls())
load('filtered.logcpm.corrected.Rdata')
library(psych)
library(stats)
# -----------------------------------------------------------------------------
#
# filtering
#
# -----------------------------------------------------------------------------
thr <- 0.8
geneLowVarianceFilter <- function(dataMatrix, percentile=0.1) {
  variances <- apply(as.matrix(dataMatrix), 1, var)
  threshold <- quantile(variances, c(percentile))
  mask <- apply(dataMatrix, 1, function(x) var(x) > threshold)
  fdata <- dataMatrix[mask, ]
  # return the row mask and filtered data
  list(mask=mask, fdata=fdata)
}
geneLowValueFilter <- function(dataMatrix, percentile=0.1) {
  # Remove gene profiles with low absolute values in dataMatrix. Returns:
  # 1) a logical vector mask identifying gene expression profiles in dataMatrix
  #    that have absolute expression levels in the lowest 10% of the data set.
  # 2) a data matrix containing filtered expression profiles.
  threshold <- quantile(as.matrix(dataMatrix), c(percentile))
  mask <- apply(dataMatrix, 1, function(x) all(abs(x) < threshold))
  fdata <- dataMatrix[!mask, ]  
  # return the row mask and filtered data
  list(mask=!mask, fdata=fdata)
}
geneHighCoefOfVarFilter <- function(dataMatrix, coefOfVarThreshold=0.1) {
  coefOfVars <- apply(dataMatrix, 1, function(x) { sd(x) / abs(mean(x)) })
  # the smaller the threshold, the higher the experimental effect relative 
  # to the measurement precision
  # filter the data matrix
  fdata <- dataMatrix[coefOfVars < coefOfVarThreshold, ]
  # return the filtered data
  list(coefOfVars = coefOfVars, fdata=fdata)
}

filtered <- geneHighCoefOfVarFilter(filtered.logcpm.corrected, thr)
filtered.combined <- data.frame(filtered$fdata)

num.genes <- nrow(filtered.combined)
genes <- row.names(filtered.combined)
rnaSeq <- t(filtered.combined)
subjects <- row.names(rnaSeq)
my_subjs <- subjects


# -----------------------------------------------------------------------------
#
# import phenotypes:
#
# -----------------------------------------------------------------------------
covs <- read.table("Demographic_symptom.csv", sep = ",", header = T, row.names = 1)

covs.short <- covs[my_subjs,]
phenos <- factor(covs.short$Diag)
names(phenos) <- rownames(covs.short)
phenos <- phenos[my_subjs]
rnaSeq <- rnaSeq[my_subjs,]
num.genes <- ncol(rnaSeq)
save(num.genes, rnaSeq, my_subjs, phenos, covs.short,
     file = paste(thr, "genes.filtered.corrected.Rdata", sep = ""))

