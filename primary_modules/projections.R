# projections.R - Bill White - 2/27/16

library(pca3d)

# -----------------------------------------------------------------------------
plotMds <- function(inData, mainTitle="", subjectColors=rep("black", nrow(inData))) {
  d <- dist(inData) # euclidean distances between the rows
  fit <- cmdscale(d, eig=TRUE, k=2) # k is the number of dim
  # plot solution 
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main=mainTitle, type="n")
  text(x, y, labels=row.names(inData), cex=0.7, col=subjectColors)
}

# -----------------------------------------------------------------------------
plotPca <- function(inData, mainTitle="", subjectColors=rep("black", nrow(inData)), 
                    scaleFlag=FALSE, showLabelsFlag=FALSE, showLabelsText=NULL,
                    showLabelsColor="black") {
  pcaRes <- prcomp(inData, scale=scaleFlag)
  varExplained <- as.integer((pcaRes$sdev) ^ 2 / sum(pcaRes$sdev ^ 2) * 100)
  # pca2d(pcaRes, group=subjectColors, show.labels=showLabelsFlag, title=mainTitle, 
  #       show.ellipses=FALSE)
  pcaScores <- pcaRes$x[, 1:2]
  plot(pcaScores, asp=1, pch=16, col=subjectColors, main=mainTitle, 
       xlab=paste("PCA 1 - ", varExplained[1], "% variance explained", sep=""), 
       ylab=paste("PCA 2 - ", varExplained[2], "% variance explained", sep=""))
  if(showLabelsFlag) {
    if(is.null(showLabelsText)) {
      text(pcaScores[, 1], pcaScores[, 2], colnames(inData), col=showLabelsColor)
    } else {
      text(pcaScores[, 1] + 2, pcaScores[, 2], showLabelsText, col=showLabelsColor)
    }
  }
  #biplot(pcaRes, 1:2)
}

# -----------------------------------------------------------------------------
geneCaseCtrlBoxplot <- function(expr, case.control, gene.name) {
  boxplot.data <- list(controls=as.numeric(expr[gene.name, case.control == 0]), 
                       cases=as.numeric(expr[gene.name, case.control == 1]))
  boxplot(boxplot.data, main=gene.name)
}

# -----------------------------------------------------------------------------
removeSubjects <- function(expr, subject.phenos, subject.ids) {
  if(length(subject.ids) > 0) {
    expr <- expr[, -which(colnames(expr) %in% subject.ids)]
  }
  subj.ids.new <- colnames(expr)
  gene.ids.new <- rownames(expr)
  phenos.new <- subject.phenos[which(subject.phenos[, 1] %in% subj.ids.new), ]
  phenos.new <- phenos.new[order(phenos.new[, 1]), ]
  pheno.ids.new <- as.character(phenos.new[, 1])
  phenos.new <- as.numeric(phenos.new[, 3])
  subj.colors <- ifelse(phenos.new == 1, "red", "black")
  condition <- factor(ifelse(phenos.new == 1, "case", "control"))
  list(expr=expr, subj.ids=subject.ids, gene.ids=gene.ids.new, 
       phenos=phenos.new, phenos.factor=condition, subj.colors=subj.colors)
}
