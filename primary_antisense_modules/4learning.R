# install.packages("corrplot")
library(corrplot)

rm(list = ls())

load("0.8genes.filtered.corrected.Rdata")
load('enrichmentProfile.Rdata')
corrplot(cor(collapse.genes), order = "hclust")

psqi.df <- covs.short[,"psqi_score", drop = F]
madrs.df <- read.csv("R01_madrs_data.csv", row.names = 1, header = T)
madrs.df <- madrs.df[,-(1:3)]
madrs <- madrs.df[,"Total", drop = F]
madrs <- na.omit(madrs)
cases <- names(phenos)[phenos == "MDD"]
ctrls <- names(phenos)[phenos == "HC"]
psqi <- psqi.df
psqitype <- gsub("psqi_score_","",colnames(psqi))
psqi <- na.omit(psqi)
shaps.df <- covs.short[,"shaps_score", drop = F]
shaps <- shaps.df
shapstype <- gsub("shaps_score_","",colnames(shaps))
shaps <- na.omit(shaps)
subj.shaps.genes <- intersect(row.names(collapse.genes), row.names(shaps))
subj.psqi.genes <- intersect(row.names(collapse.genes), row.names(psqi.df))
subj.madrs.genes <- intersect(row.names(collapse.genes), row.names(madrs))
subj.madrs.cases <- intersect(cases, row.names(madrs))
subj.shaps.cases <- intersect(subj.shaps.genes, cases)
subj.psqi.cases <- intersect(subj.psqi.genes, cases)

diag.split <- ggplot(covs.short, aes(x = batch, fill = Diag)) +
  geom_bar(stat='Count', position='dodge') + scale_fill_brewer(palette="Dark2") #+ 
pdf(file = "BatchDiag.pdf", height = 4, width = 3.5)
diag.split 
dev.off()
# ------------------------------------------------------------------------------------
# Logistic regression with Phenotypes, PSQI, SHAPS
# ------------------------------------------------------------------------------------

linear.reg <- function(my.formula, data, y, subj.id, fam = gaussian){
  linear_result <- matrix(NaN, num.modules, 4)
  colnames(linear_result) <- c("Intercept", "Effect", "sdError", "pvalue")
  rownames(linear_result) <- modules
  ptm <- proc.time()
  for (i in 1:(num.modules)){
    tryCatch({
      model.df <- data.frame(genei = data[subj.id,i], y = y[subj.id,], sab.covs[subj.id,])
      my.linear <- tidy(glm(as.formula(my.formula), data =  model.df, family = fam))
      my.linear.genei <- my.linear[1:2,]
      linear_result[i,1:2]  <- my.linear.genei$estimate
      linear_result[i,3]  <- my.linear.genei$std.error[2]
      linear_result[i,4] <- my.linear.genei$p.value[2]
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
  linear.batch <- data.frame(linear_result[order(linear_result$pvalue,decreasing=F),])
  print(head(linear.batch, 10))
  linear.batch
}

pheno.nobatch.lin <- linear.reg(my.formula = "y ~ . - batch", 
                                data = collapse.genes, y = phenos.df, 
                                subj.id = rownames(collapse.genes), fam = "binomial")

psqi.batch.lin <- linear.reg(my.formula = "y ~ .", data = collapse.genes, y = psqi, subj.id = subj.psqi.cases)
shaps.batch.lin <- linear.reg(my.formula = "y ~ . ", data = collapse.genes, y = shaps, subj.id = subj.shaps.genes)
pheno.batch.lin <- linear.reg(my.formula = "y ~ . ", data = collapse.genes, y = phenos.df, subj.id = rownames(collapse.genes), fam = "binomial")
madrs.batch.lin <- linear.reg(my.formula = "y ~ . ", data = collapse.genes, y = madrs, subj.id = subj.madrs.genes)

scaled.madrs.batch.lin <- linear.reg(my.formula = "y ~ . ", data = collapse.genes, y = scale(madrs), subj.id = subj.madrs.genes)

my.modules <- rownames(madrs.batch.lin)
my.modules <- gsub("odule", "", my.modules)
plot(-log(madrs.batch.lin$p.value.corrected.FDR), xaxt = "n", 
     xlab='Module', ylab = "-log(p.adjust)", main = "Modules’ Associations with MADRAS")
abline(h = -log(0.05))
abline(h = -log(0.1), type = 2)
text(21.5,3.2, "5% FDR Threshold")
text(21.5,-log(0.1)+0.2, "10% FDR Threshold")
axis(1, at=1:nrow(madrs.batch.lin), labels=my.modules)

madrs.plot <- madrs.batch.lin
madrs.plot$signi <- -log(madrs.plot$p.value.corrected.FDR)
madrs.plot$Module <- my.modules
madrs.plot$signilevel <- as.factor((madrs.plot$p.value.corrected.FDR < 0.05) + 
                                     (madrs.plot$p.value.corrected.FDR < 0.1))
#install.packages("ggthemes")
library("ggthemes")

cairo_pdf(filename = "moduleAsso.pdf", height = 5, width = 10)
ggplot(madrs.plot, aes(x = Module, y = signi, alpha = signilevel)) + geom_point(size = 2) +
     xlab('Module') + ylab("-log(adjusted.p)") + ggtitle("Modules’ Associations with MADRS") + 
  geom_hline(yintercept = -log(0.05)) + geom_hline(yintercept = -log(0.1), linetype = 2) + 
  annotate("text", 1.9, -log(0.05) + 0.15, label = "0.05 FDR", size = 4)+
  annotate("text", 1.78, -log(0.1) + 0.15, label = "0.1 FDR", size = 4) + 
  scale_x_discrete(limits=paste("M", 1:23, sep = "")) + 
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 14), title = element_text(size = 16),
        panel.background = element_blank(), #panel.grid.minor.y = element_line(size=3),
        panel.grid.major = element_line(colour = "grey90"),
        # plot.margin = margin(half_line, half_line, half_line, half_line))
        plot.margin=unit(c(1, 1, 0.5, 0.5), "lines")) 
dev.off()

madrs.plot.ordered <- madrs.plot[order(madrs.plot$p.value.corrected.FDR, decreasing = T), ]
madrs.plot.ordered$Module <- paste0("Module ", gsub("M", "", madrs.plot.ordered$Module))
madrs.plot.ordered$Module <- factor(madrs.plot.ordered$Module, 
                                    levels = madrs.plot.ordered$Module)

pImp <- ggplot(madrs.plot.ordered, aes(y = Module, x = signi, alpha = signilevel)) + geom_point(size = 2) +
  xlab('') + ylab("") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5), 
        axis.title = element_text(size = 14), title = element_text(size = 16),
        panel.background = element_blank(), #panel.grid.minor.y = element_line(size=3),
        panel.grid.major = element_line(colour = "grey90"),
        # plot.margin = margin(half_line, half_line, half_line, half_line))
        plot.margin=unit(c(1, 1, 0.5, 0.5), "lines")) 

pdf("moduleImportances.pdf", height = 4, width = 2.5)
pImp
dev.off()

# important modules:
important.mods <- rownames(madrs.batch.lin)[madrs.batch.lin$p.value.corrected.FDR < 0.05]
# sapply(c7_genes_sets$genesets[important.mods], length)

madrs.collapse.lin <- glm(y ~ ., data = data.frame(genes = collapse.genes[subj.madrs.genes,], y = madrs[subj.madrs.genes,]))
tidy(madrs.collapse.lin)
collapse.lin <- glm(phenos ~ ., data = data.frame(collapse.genes, phenos), family = "binomial")
tidy(collapse.lin)    

write.csv(madrs.batch.lin, 'madrs.WGCNA.0.3.csv')


# Linear regression on each gene:
# madrs.batch.genes <- linear.reg(my.formula = "y ~ . ", data = rnaSeq, y = madrs, subj.id = subj.madrs.genes)
linear_result <- matrix(NaN, ncol(rnaSeq), 3)
colnames(linear_result) <- c("Intercept", "Effect", "pvalue")
rownames(linear_result) <- colnames(rnaSeq)
ptm <- proc.time()
data = rnaSeq
y = phenos.df
subj.id = rownames(rnaSeq)
fam = "binomial"
# y = madrs
# subj.id = subj.madrs.genes
# fam = "gaussian"

for (i in 1:(ncol(rnaSeq))){
# for (i in 1:100){
  # print(i)
  tryCatch({
    model.df <- data.frame(genei = data[subj.id,i], y = y[subj.id,], sab.covs[subj.id,])
    my.linear <- tidy(glm(y ~ ., data =  model.df, family = fam))
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
print(head(linear.batch, 10))
linear.batch$Significance <- -log(linear.batch$pvalue)
linear.batch.module <- merge(group.genes, linear.batch, by = "row.names")
linear.batch.module <- linear.batch.module[order(linear.batch.module$Significance, decreasing = T),]

write.csv(linear.batch.module, 'individual.logistic.csv')

# i = 17
imp.vec <- vector(mode = "numeric", length = num.modules)
size.vec <- vector(mode = "numeric", length = num.modules)
names(imp.vec) <- paste("Module", 1:num.modules, sep = "")
names(size.vec) <- paste("Module", 1:num.modules, sep = "")

for (i in 1:num.modules){
  genes.in.modi <- group.genes[group.genes$dynamicMods==i,1]
  imp.vec[i] <- sum(linear.batch[genes.in.modi, "Significance"])  
  size.vec[i] <- length(genes.in.modi)
}

imp.size.df <- data.frame(TotalSignificance = imp.vec, ClusterSize = size.vec, SingleSignificance = imp.vec/size.vec)
summary.df <- merge(imp.size.df, madrs.batch.lin, by = "row.names")
summary.df <- summary.df[order(summary.df$pvalue),]
# madrs.batch.lin$TotalSignificanceMADRS <- imp.vec
# madrs.batch.lin$ClusterSize <- size.vec
# madrs.batch.lin$ActualSignificance <- madrs.batch.lin$TotalSignificanceMADRS/madrs.batch.lin$ClusterSize
summary.df$pvalueSignificance <- -log(summary.df$pvalue)
tidy(cor.test(summary.df$ClusterSize, summary.df$pvalueSignificance))

cor(summary.df$pvalueSignificance, summary.df$SingleSignificance, method = "spearman")

k <- 100
hyper.df <- data.frame(mod = 1:num.modules, pval = NA)

for (i in 1:num.modules){
  bestgenes <- table(linear.batch.module[1:k, "dynamicMods"])
  bestgenes <- data.frame(t(bestgenes))
  bestgenes <- bestgenes[,-1]
  bestgenes[,1] <- paste("Module", 1:num.modules, sep = "")
  m <- summary.df[summary.df$Row.names == paste("Module", i, sep = ""), "ClusterSize"] # number of genes from M5; 
  n <- nrow(group.genes) - m; # number of genes from other modules
  x <- bestgenes[bestgenes$Var2 == paste("Module", i, sep = ""), "Freq"] -1
  x[x<0] <- 0
  hyper.df[i, "pval"] <- phyper(x, m, n, k, lower.tail = F)
  
}

hyper.df[order(hyper.df$pval),]
# write.csv(summary.df, file = "summaryGLM.csv")

cv.elastic.regr <- cv.glmnet(collapse.genes[subj.madrs.genes, ], madrs[subj.madrs.genes, ], alpha = 0.7)
plot(cv.elastic.regr)
best.lamb <- cv.elastic.regr$lambda.1se
elastic.coef <- predict(cv.elastic.regr, type = "coefficients", s = best.lamb)
t(elastic.coef)

##### Calculate centralities for each gene:
library(igraph)
inputmat = rnaSeq
corrThresh = 0.2
inputmat <- cor(inputmat)
powCorr <- inputmat^1
adjacency <- abs(powCorr) > corrThresh
adj.plot <- adjacency
diag(adj.plot) <- 0
graph <- graph.adjacency(adj.plot, weighted=TRUE, mode="lower")
gene.centralities <- eigen_centrality(graph)$vector
linear.batch.module <- read.csv('individual.logistic.csv')
imp.modules <- as.numeric(gsub("Module", "", rownames(madrs.batch.lin)))
# i = 5
par(mfrow = c(2,2))
idx <- 0
important.plots <- list()

idx <- idx + 1
linear.i <- linear.batch.module[linear.batch.module$dynamicMods %in% imp.modules[1:2],]
rownames(linear.i) <- linear.i$gene
# top 6 modules
linear.i <- linear.batch.module
centralities.i <- sort(gene.centralities[rownames(linear.i)], decreasing = T)
head(centralities.i)
genesi <- names(centralities.i)
genesi <- genesi[!startsWith(genesi, "LOC")]
imp.df <- data.frame(centrality = centralities.i[genesi], 
                     Module = factor(group.genes[genesi, "dynamicMods"]),
                     imp.score = -log(linear.i[genesi, "p.value.corrected.FDR"]),
                     rownames = genesi)
# important.plots[[idx]] <-
imp.df$modf = factor(imp.df$Module, levels=imp.modules)
imp.df$gene <- rownames(imp.df)
library(grDevices)
# pdf(file = "moduleCentralities.pdf", width = 7, height = 7)
# cairo_pdf(filename = "moduleCentralities175.pdf", width = 7, height = 4.5)
# cairo_pdf(filename = "moduleCentralitiesall.pdf", width = 7, height = 9)

ggplot(imp.df, aes(x = centrality, y = imp.score, 
                     color = Module)) +
  geom_text(aes(label = gene), size = 1.5) + geom_hline(aes(yintercept = -log(0.05)), linetype = 3) +
  # geom_text(aes(0.05, -log(0.05), label = "0.05 FDR Threshold", vjust = -1))+
  annotate("text", 0.1, -log(0.05) + 0.1, label = "0.05 FDR", size = 2) +
  geom_hline(aes(yintercept = -log(0.1)), linetype = 2) +
  annotate("text", 0.1, -log(0.1) + 0.1, label = "0.1 FDR", size = 2) +
  theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
  geom_smooth(method='lm', se = F) + facet_wrap( ~ modf, ncol=4) + theme_bw() + 
  labs(x = "Centrality", y = "-Log(adjusted P-value)",
       title = paste("Importance Score vs. Centrality of Individual Genes per Module", sep = ""))


  dev.off()

  
linear.i <- linear.batch.module

rownames(linear.i) <- linear.i$gene
# top 6 modules
centralities.i <- sort(gene.centralities[rownames(linear.i)], decreasing = T)

genesi <- names(centralities.i)
genesi <- genesi[!startsWith(genesi, "LOC")]
imp.df <- data.frame(centrality = centralities.i[genesi], 
                     Module = factor(group.genes[genesi, "dynamicMods"]),
                     imp.score = -log(linear.i[genesi, "p.value.corrected.FDR"]),
                     rownames = genesi)
  
  # ggsave(g.cent, file = "moduleCentralities.svg")  
pvals.cent <- vector(mode = "numeric")  
coef.cent <- vector(mode = "numeric")  
rsquared <- vector(mode = "numeric")  
  for (i in imp.modules){
    mydf <- imp.df[imp.df$Module == i,]
    mylm <- lm(imp.score ~ centrality, data = mydf)
    pvals.cent <- c(pvals.cent, tidy(mylm)$p.value[2])
    coef.cent <- c(coef.cent, mylm$coefficients[2])
    rsquared <- c(rsquared, summary(mylm)$r.squared)
  }
names(pvals.cent) <- imp.modules
names(coef.cent) <- imp.modules
coef.cent
(mod.cent.df <- data.frame(module = imp.modules, coefficients = coef.cent, 
                          pvalue = pvals.cent, rsquared = rsquared, signi = seq(1,num.modules)))
write.csv(mod.cent.df, file = "moduleCentralitySign.csv")

summary(lm(rsquared ~ signi, data = mod.cent.df))
summary(lm(pvalue ~ signi, data = mod.cent.df))
summary(lm(coefficients ~ signi, data = mod.cent.df))



