#!/usr/bin/Rscript
#
# 4 cohorts.

library(edgeR)
library(dplyr)
library(reshape2)
library(limma)

#rm(list = ls())

# get raw counts - positive integers - from HTSeq
all.cohorts <- read.table("libr_mdd_rnaseq.tab", 
                          sep = "\t", header = TRUE, stringsAsFactors = FALSE, 
                          row.names = 1)
all.genes <- rownames(all.cohorts)

subj.total.counts <- colSums(all.cohorts)
mean(subj.total.counts) 

dim(all.cohorts) # 25369 x159

# remove low count genes 
y <- DGEList(counts = all.cohorts)
keep <- rowSums(cpm(y) > 10) >= 16
combined <- all.cohorts[keep,] # removed low counts genes
#combined <- all.cohorts[rowSums(all.cohorts) > 0, ]
all.subjs <- colnames(combined)
dim(combined)  # 10142 x 159

subj.total.counts <- colSums(all.cohorts)
mean(subj.total.counts) # 22,165,139 average counts per sample

# batch ids
cohort.ids <- read.table("cohort_subjects.tab", header = TRUE, sep = " ",
                         stringsAsFactors = FALSE)

#batch_strings <- subject.attrs$batch
#batch_strings[batch_strings=="1st"] <- 1 
#subject.attrs$X

# remove outliers
#outliers <- c("AN736", "AE396", "AI270")
combined <- combined[,-which(colnames(combined)=="AN736")]
combined <- combined[,-which(colnames(combined)=="AE396")]
#combined <- combined[,-which(colnames(combined)=="AI270")]

cohort.ids <- cohort.ids[-which(cohort.ids[,2]=="AN736"),]
cohort.ids <- cohort.ids[-which(cohort.ids[,2]=="AE396"),]
#cohort.ids <- cohort.ids[-which(cohort.ids[,2]=="AI270"),]

subj.total.counts <- colSums(combined)
hist(subj.total.counts)
mean(subj.total.counts) # 22M average counts per sample

# get matching phenotypes
subject.attrs <- read.csv("Demographic_symptom.csv", stringsAsFactors = FALSE)
attrs.df <- tbl_df(subject.attrs)
phenos.df <- attrs.df %>% 
  filter(X %in% colnames(combined)) %>%
  dplyr::select(X, Diag)
phenos.df <- phenos.df[match(colnames(combined), phenos.df$X), ]
phenos <- ifelse(phenos.df$Diag == "HC", 0, 1)

# transform raw counts
# CPM + log2 transform
combined.cpm <- cpm(combined, log = F)
combined.cpm.corrected <- removeBatchEffect(x = as.matrix(combined.cpm), 
                                            batch = cohort.ids$COHORT)
combined.cpm  <- combined.cpm.corrected
