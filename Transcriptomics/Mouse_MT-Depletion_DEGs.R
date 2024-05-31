# Clear all variables in the environment
rm(list=ls(all=T))

# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(AnnotationHub)
library(org.Mm.eg.db)
library(clusterProfiler)
library(plyr)

# Set working directory
setwd("~/Mouse_Mito-depletion/Rresult")

# Define sample numbers
controlnum <- 2
casenum <- 2
num <- casenum + controlnum

# Load data
datad4 <- as.matrix(read.csv("../RNAdata/genematrix/gene_count_matrix.csv", row.names="gene_id"))
colnames(datad4) <- c("WT_iPS1", "WT_iPS2", "MD_D4_iPS1", "MD_D4_iPS2")
datad4 <- datad4[, c(3, 4, 1, 2)]
# Compare Day4 with WT
sampletitle <- "Mouse MT_D4 vs WT_iPS"
countData <- datad4
condition <- factor(c(rep("treat", casenum), rep("untreat", controlnum)))
colData <- data.frame(row.names=colnames(countData), condition)

# Differential expression analysis
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalize=TRUE)), by="row.names", sort=FALSE)
resdata <- resdata[!is.na(resdata$padj),]
resdata <- resdata[order(resdata$padj),]

# Gene annotation
PluripotencyGene <- read.table("~/PosdocProject/Denial/Human_Mito-depletion/12.5update/Pluripotency_Genes.txt", header = F)
MT_coding <- read.table("~/PosdocProject/Denial/Human_Mito-depletion/Rresult/Human mitochondria genome coding genes.csv", sep = ",", header = T)
MT_asso <- read.table("../Mouse_chromosome_related_mito_genelist.txt", sep = "\t", header = F)
resdata$genename <- sapply(resdata$Row.names, function(x) strsplit(x, "\\|")[[1]][2])
resdata$genename <- toupper(resdata$genename)
resdata$genetype <- NA
resdata$genetype[resdata$genename %in% toupper(c(MT_coding[, 2], "MT-CYTB"))] <- "MT_coding"
resdata$genetype[resdata$genename %in% toupper(MT_asso[, 1])] <- "MT_associated"
resdata$genetype[resdata$genename %in% toupper(PluripotencyGene[, 1])] <- "PluripotencyGenes"

# PCA plot for RNA-seq (Day4 only)
library(ggbiplot)
pcamat <- as.data.frame(t(resdata[, 8:11]))
rownames(pcamat) <- c("MT_Depletion_Day4_S1", "MT_Depletion_Day4_S2", "WT_Day0_S1", "WT_Day0_S2")
pcamatgroup <- c(rep("MT_Depletion", 2), rep("WT", 2))
pca_result <- prcomp(pcamat, scale = TRUE, center = FALSE)
pca_result$x <- scale(pca_result$x)

# Plot PCA
pdf(file = "Mouse_MT_Depletion_pcaplot.pdf", width = 5, height = 5)
ggbiplot(pca_result,
         var.axes = FALSE,
         obs.scale = 1,
         groups = pcamatgroup,
         ellipse = FALSE,
         circle = FALSE) +
  geom_text(aes(label = rownames(pcamat)), vjust = 1.5, size = 2) +
  scale_x_continuous(limits = c(-1.1, 1.1)) +
  scale_y_continuous(limits = c(-1.1, 1.1)) +
  theme_bw()
dev.off()

# Volcano plot and heatmap for differentially expressed genes
diff_gene_p <- resdata[resdata$padj < 0.05, ]
diff_gene_p$stat[diff_gene_p$log2FoldChange < 0] <- "Up"
diff_gene_p$stat[diff_gene_p$log2FoldChange > 0] <- "Down"
diff_gene_p$log2FoldChange <- -diff_gene_p$log2FoldChange
diff_gene_p <- diff_gene_p[abs(diff_gene_p$log2FoldChange) > 1.5, ]
write.table(diff_gene_p, sep = "\t", row.names = F, quote = F, "Mouse_MD_DEGs.txt")

# Volcano plot
resdata$log2FoldChange <- -resdata$log2FoldChange
data <- resdata
data$threshold <- as.factor(ifelse(data$padj < 0.05 & abs(data$log2FoldChange) >= 1.5, ifelse(data$log2FoldChange > 1.5, 'Up', 'Down'), 'Not'))
data$pos <- rowMeans(data[, 8:(8 + casenum - 1)])
data$neg <- rowMeans(data[, (8 + casenum):(8 + num - 1)])
data$genetype[!data$Row.names %in% diff_gene_p$Row.names] <- NA
data$genename[is.na(data$genetype)] <- NA

p <- ggplot(data, aes(x = log2(pos), y = log2(neg), colour = threshold)) + geom_point()
p1 <- p + geom_point(aes(color = data$threshold))
p2 <- p1 + scale_color_manual(values = c("green", "grey", "#f8766d")) + geom_point(alpha = 0.4, size = 1)
p3 <- p2 + theme_bw() + labs(x = "Case_genes", y = "Control_genes", title = sampletitle) + theme(plot.title = element_text(hjust = 0.5))
p4 <- p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p5 <- p4 + geom_abline(intercept = 0, slope = 1, lty = 4, col = "black", lwd = 0.6)
a <- cor.test(data$pos, data$neg, method = "pearson")
p6 <- p5 + annotate("text", x = 14, y = -0.5, label = paste("Pearson's correlation ", round(a$estimate, 3)), size = 3)
ggsave(p6, file = paste(sampletitle, "_correlation2.pdf", sep = ""), width = 8, height = 8)

# Heatmap for differentially expressed genes
datamat <- as.matrix(diff_gene_p[, c(8:11)])
rownames(datamat) <- diff_gene_p$genename
anno_gene2 <- diff_gene_p[!is.na(diff_gene_p$genetype), c(12, 13)]
anno_gene <- as.factor(anno_gene2$genetype)
anno_gene <- as.matrix(anno_gene)
anno_gene <- as.data.frame(anno_gene)
rownames(anno_gene) <- anno_gene2$genename
rownames(datamat)[rownames(datamat) == "NA"] <- ""

pdf(file = "Mouse_DEGs3.pdf", width = 8, height = 8)
pheatmap(mat = datamat, scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100), border_color = "white")
dev.off()
